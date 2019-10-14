#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>
#include <memory>

#include "cif++/Config.h"
#include "svm++.h"

using namespace std;

namespace cif {
extern int VERBOSE;
}

namespace svm
{

// --------------------------------------------------------------------

vector<pair<value_type,value_type>> Matrix::Scale()
{
	vector<pair<value_type,value_type>> result;
	
	// collect min/max pairs for each column
	for (auto& v: m_data)
	{
		if (result.size() < v.size())
			result.insert(result.end(), v.size() - result.size(),
				make_pair(numeric_limits<float>::max(), numeric_limits<float>::lowest()));
		
		for (size_t i = 0; i < v.size(); ++i)
		{
			if (result[i].first > v[i])
				result[i].first = v[i];
			if (result[i].second < v[i])
				result[i].second = v[i];
		}
	}
	
	// adjust the stored values
	for (auto& v: m_data)
	{
		for (size_t i = 0; i < v.size(); ++i)
			v[i] = (v[i] - result[i].first) / (result[i].second - result[i].first);
	}
	
	return result;
}

// --------------------------------------------------------------------
//	View is used to access just a subset of a Matrix, potentially reordered

class View : public Data
{
  public:
	View& operator=(const View&) = delete;

	View(const View& rhs)
		: m_m(rhs.m_m), m_p(rhs.m_p) {}
	View(const Data& m, const vector<size_t>& perm)
		: m_m(const_cast<Data&>(m)), m_p(const_cast<vector<size_t>&>(perm)) {}
	
	virtual Vector& operator[](size_t row)				{ return m_m[m_p[row]]; }
	virtual const Vector& operator[](size_t row) const	{ return m_m[m_p[row]]; }
	virtual size_t size() const							{ return m_p.size(); }

	void reserve(size_t size)							{ m_p.reserve(size); }
	void push_back(const Vector& v)
	{
		for (size_t i = 0; i < m_m.size(); ++i)
		{
			if (&m_m[i] == &v)
			{
				m_p.push_back(i);
				break;
			}
		}
	}

  private:
	Data&			m_m;
	vector<size_t>	m_p;
};

// --------------------------------------------------------------------
//	QMatrix 'contains' Q, which is y[i] * y[j] * K(i, j) for C_SVC e.g.
//	Actually, it is simply a wrapper and can be optimized eventually as is done in libsvm

template<SVMType S, KernelType K>
class QMatrix;

template<KernelType K>
class QMatrix<C_SVC, K>
{
  public:
	typedef Kernel<K>		kernel_type;

	QMatrix(const vector<int8_t>& y, const Data& x, const typename kernel_type::param_type& params)
		: m_y(y), m_K(x, params) {}
	
	double operator()(size_t i, size_t j) const
	{
		return m_y[i] * m_y[j] * m_K(i, j);
	}
	
  private:
	const vector<int8_t>&	m_y;
	kernel_type 		m_K;
};

// --------------------------------------------------------------------
//	Solver is used to solve the quadratic optimization problem.
//	Code here is based on the paper:
//
//	Journal of Machine Learning Research 6 (2005) 1889–1918
//	Submitted 04/05; Revised 10/05; Published 12/05
//	Working Set Selection Using Second Order Information
//	for Training Support Vector Machines
//	Rong-En Fan
//	Pai-Hsuen Chen
//	Chih-Jen Lin
//	
//	This is somewhat different from the code in libsvm (and more readable).

template<typename QM>
class Solver
{
  public:
	Solver(const QM& Q, double eps, double tau, double C)
		: m_Q(Q), m_eps(eps), m_tau(tau), m_C(C) {}
	
	tuple<vector<double>,double,double> Solve(const vector<label_type>& y)
	{
		size_t len = y.size();
		
		vector<double> A(len, 0);
		vector<double> G(len, -1);
		
		size_t iteration = 0;
		size_t max_iteration = 10000000;

		while (++iteration < max_iteration)
		{
			int32_t i, j;
			tie(i, j) = selectB(y, A, G);

			if (j == -1)
				break;
			
			auto a = m_Q(i, i) + m_Q(j, j) - 2 * y[i] * y[j] * m_Q(i, j);
			if (a <= 0)
				a = m_tau;
			
			auto b = -y[i] * G[i] + y[j] * G[j];
			auto oldAi = A[i], oldAj = A[j];
			A[i] += y[i] * b / a;
			A[j] -= y[j] * b / a;
			
			auto sum = y[i] * oldAi + y[j] * oldAj;
			if (A[i] > m_C)
				A[i] = m_C;
			if (A[i] < 0)
				A[i] = 0;
			A[j] = y[j] * (sum - y[i] * A[i]);
			if (A[j] > m_C)
				A[j] = m_C;
			if (A[j] < 0)
				A[j] = 0;
			A[i] = y[i] * (sum - y[j] * A[j]);
		
			auto deltaAi = A[i] - oldAi, deltaAj = A[j] - oldAj;
			for (size_t t = 0; t < len; ++t)
				G[t] += m_Q(t, i) * deltaAi + m_Q(t, j) * deltaAj;
		}
		
		if (cif::VERBOSE > 1)
			cout << "optimization finished in " << iteration << " iterations" << endl;
		
		double sum = 0, v = 0;
		size_t n = 0;
		for (size_t i = 0; i < len; ++i)
		{
			if (A[i] > 0 and A[i] < m_C)
			{
				sum += y[i] * G[i];
				++n;
			}
			v += A[i] * (G[i] - 1);
		}
		double rho = sum / n;
		double obj = v / 2;
		
		return make_tuple(A, rho, obj);
	}
	
	tuple<int32_t,int32_t> selectB(const vector<int8_t>& y, const vector<double>& A,
		const vector<double>& G)
	{
		// select i
		int32_t i = -1;
		double Gmax = numeric_limits<double>::lowest();
		for (size_t t = 0; t < y.size(); ++t)
		{
			if ((y[t] == +1 and A[t] < m_C) or
				(y[t] == -1 and A[t] > 0))
			{
				if (-y[t] * G[t] >= Gmax)
				{
					i = t;
					Gmax = -y[t] * G[t];
				}
			}
		}
		
		// select j
		int32_t j = -1;
		double Gmin = numeric_limits<double>::max();
		double objMin = numeric_limits<double>::max();
		for (size_t t = 0; t < y.size(); ++t)
		{
			if ((y[t] == +1 and A[t] > 0) or
				(y[t] == -1 and A[t] < m_C))
			{
				auto b = Gmax + y[t] * G[t];
				if (Gmin > -y[t] * G[t])
					Gmin = -y[t] * G[t];
				if (b > 0)
				{
					auto a = m_Q(i, i) + m_Q(t, t) - 2 * y[i] * y[t] * m_Q(i, t);
					if (a <= 0)
						a = m_tau;
					if (objMin > -(b * b) / a)
					{
						objMin = -(b * b) / a;
						j = t;
					}
				}
			}
		}
		
		if (Gmax - Gmin < m_eps)
			i = j = -1;
		
		return make_tuple(i, j);
	}
	
	
  private:
	const QM&	m_Q;
	double		m_eps, m_tau, m_C;
};

// --------------------------------------------------------------------
//	RBF Kernel implementation

double Kernel<RBF>::f(const Vector& x, const Vector& y, const param_type& params)
{
	double sum = 0;

	size_t l = x.size();
	if (l < y.size())
		l = y.size();

	for (size_t i = 0; i < l; ++i)
	{
		if (x.contains(i) and y.contains(i))
		{
			double d = x[i] - y[i];
			sum += d * d;
		}
		else if (x.contains(i))
		{
			sum += x[i] * x[i];
		}
		else if (y.contains(i))
		{
			sum += y[i] * y[i];
		}
	}
	
	return exp(-params.gamma * sum);
}

// --------------------------------------------------------------------

double Kernel<RBF>::operator()(std::size_t i, std::size_t j) const
{
	auto& xi = m_x[i];
	auto& xj = m_x[j];
	
	return exp(-m_gamma * (dot(xi, xi) + dot(xj, xj) - 2 * dot(xi, xj)));
}

// --------------------------------------------------------------------

ModelBase::ModelBase(zeep::xml::element& config)
{
	for (auto rho: config.find("rho"))
		m_rho.push_back(stod(rho->str()));
	
	for (auto c: config.find("class"))
	{
		ClassInfo ci = {
			static_cast<label_type>(stoi(c->get_attribute("label"))),
			stoull(c->get_attribute("nr_sv"))
		};
		
		m_class.push_back(ci);
	}

	for (auto sv: config.find("svm-node"))
	{
		auto sv_coef = sv->find("sv_coef");
		vector<double> sv_coef_v;
		for (auto svc: sv_coef)
			sv_coef_v.push_back(stod(svc->str()));
		m_sv_coef.push_back(sv_coef_v);

		Vector row;

		for (auto svv: sv->find("sv"))
		{
			size_t index = stoull(svv->get_attribute("index"));
			double v = stod(svv->str());
			row[index] = v;
		}
		
		m_sv.push_back(move(row));
	} 
}

ModelBase* ModelBase::Create(zeep::xml::element& n)
{
	auto config = &n;
	if (config->name() != "svm-config")
		config = n.find_first("svm-config");

	string svmType = config->get_attribute("svm-type");
	string kernelType = config->get_attribute("kernel-type");

	ModelBase* result = nullptr;

	if (svmType == SVMTypeTraits<C_SVC>::GetSVMTypeString())
	{
		if (kernelType == KernelTypeTraits<RBF>::GetKernelTypeString())
			result = new Model<RBF,C_SVC>(*config);
	}
	
	return result;
}

vector<double> ModelBase::PredictWithProbability(const Vector& v) const
{
	throw runtime_error("This SVM does not implement probabilities");
	return vector<double>();
}

zeep::xml::element* ModelBase::GetConfig()
{
	using zeep::xml::element;

	auto config = new element("svm-config");
	
	for (auto rho: m_rho)
	{
		auto e = new element("rho");
		e->content(to_string(rho));
		config->append(e);
	}
	
	for (auto ci: m_class)
	{
		auto e = new element("class");
		e->set_attribute("label", to_string(ci.label));
		e->set_attribute("nr_sv", to_string(ci.nr_sv));
		config->append(e);
	}
	
	for (size_t i = 0; i < m_sv.size(); ++i)
	{
		auto e = new element("svm-node");
		for (auto sv_coef: m_sv_coef[i])
		{
			auto svc = new element("sv_coef");
			svc->content(to_string(sv_coef));
			e->append(svc);
		}
		
		auto& row = m_sv[i];
		for (size_t j = 0; j < row.size(); ++j)
		{
			if (not row.contains(j))
				continue;
			
			auto sv = new element("sv");
			sv->set_attribute("index", to_string(j));
			sv->content(to_string(row[j]));
			e->append(sv);
		}
		
		config->append(e);
	}

	return config;
}

// --------------------------------------------------------------------

template<KernelType K>
zeep::xml::element* Model<K,C_SVC>::GetConfig()
{
	typedef SVMTypeTraits<C_SVC>	svm_type_traits;
	typedef KernelTypeTraits<K>		kernel_type_traits;
	
	auto config = ModelBase::GetConfig();
	config->set_attribute("svm-type", svm_type_traits::GetSVMTypeString());
	config->set_attribute("kernel-type", kernel_type_traits::GetKernelTypeString());
	
	zeep::xml::serializer sr(config);
	sr.serialize_element("params", m_params);
	
	if (m_params.probability)
	{
		sr.serialize_element("probA", m_probA);
		sr.serialize_element("probB", m_probB);
	}
	
	return config;
}

template<KernelType K>
void Model<K, C_SVC>::Predict(const Vector& x, vector<double>& sums) const
{
	size_t nr_class = this->m_class.size();
	
	vector<double> kvalue;
	for (auto& sv: this->m_sv)
		kvalue.push_back(kernel_type::f(x, sv, m_params));
				
	vector<size_t> start(nr_class);
	start[0] = 0;
	for (size_t i = 1; i < nr_class; ++i)
		start[i] = start[i - 1] + this->m_class[i - 1].nr_sv;

	sums.clear();
	sums.reserve(nr_class * (nr_class - 1));
	
	int p = 0;
	for (size_t i = 0; i < nr_class; ++i)
		for (size_t j = i + 1; j < nr_class; ++j)
		{
			double sum = 0;
			auto si = start[i];
			auto sj = start[j];
			auto ci = this->m_class[i].nr_sv;
			auto cj = this->m_class[j].nr_sv;
			
			for (size_t k = 0; k < ci; ++k)
				sum += this->m_sv_coef[si + k][j - 1] * kvalue[si + k];
			for (size_t k = 0; k < cj; ++k)
				sum += this->m_sv_coef[sj + k][i] * kvalue[sj + k];
			
			sum -= this->m_rho[p];
			
			sums.push_back(sum);

			++p;
		}
}

template<KernelType K>
label_type Model<K, C_SVC>::Predict(const Vector& x) const
{
	size_t nr_class = this->m_class.size();

	vector<double> sums;
	Predict(x, sums);

	vector<int> vote(nr_class, 0);
	
	int p = 0;
	for (size_t i = 0; i < nr_class; ++i)
		for (size_t j = i + 1; j < nr_class; ++j)
		{
			if (sums[p] > 0)
				++vote[i];
			else
				++vote[j];

			++p;
		}

	auto voted = distance(vote.begin(),
		max_element(vote.begin(), vote.end()));
	
	return this->m_class[voted].label;
}

// Platt's binary SVM Probablistic Output: an improvement from Lin et al.
tuple<double,double> SigmoidTrain(const vector<double>& dec_value, const vector<label_type>& labels)
{
	size_t l = dec_value.size();
	
	double prior1 = 0, prior0 = 0;

	for (auto l: labels)
		if (l > 0) ++prior1; else ++prior0;
	
	int max_iter = 100;			// Maximal number of iterations
	double min_step = 1e-10;	// Minimal step taken in line search
	double sigma = 1e-12;		// For numerically strict PD of Hessian
	double eps = 1e-5;
	double hiTarget = (prior1 + 1.0) / (prior1 + 2.0);
	double loTarget = 1 / (prior0 + 2.0);
	vector<double> t(l);
	double fApB;
	
	// Initial Point and Initial Fun Value
	double A = 0.0, B = log((prior0 + 1.0) / (prior1 + 1.0));
	double fval = 0.0;

	for (size_t i = 0; i < l; ++i)
	{
		if (labels[i] > 0) t[i] = hiTarget;
		else t[i] = loTarget;

		fApB = dec_value[i] * A + B;

		if (fApB >= 0)
			fval += t[i] * fApB + log(1 + exp(-fApB));
		else
			fval += (t[i] - 1) * fApB + log(1 + exp(fApB));
	}

	int iter;
	for (iter = 0; iter < max_iter; ++iter)
	{
		// Update Gradient and Hessian (use H' = H + sigma I)
		double h11 = sigma; // numerically ensures strict PD
		double h22 = sigma;
		double h21 = 0.0;
		double g1 = 0.0;
		double g2 = 0.0;

		for (size_t i = 0; i < l; ++i)
		{
			double p, q;
			
			fApB = dec_value[i] * A + B;

			if (fApB >= 0)
			{
				p = exp(-fApB) / (1.0 + exp(-fApB));
				q = 1.0 / (1.0 + exp(-fApB));
			}
			else
			{
				p = 1.0 / (1.0 + exp(fApB));
				q = exp(fApB) / (1.0 + exp(fApB));
			}

			double d2 = p * q;
			h11 += dec_value[i] * dec_value[i] * d2;
			h22 += d2;
			h21 += dec_value[i] * d2;
			double d1 = t[i] - p;
			g1 += dec_value[i] * d1;
			g2 +=d1;
		}

		// Stopping Criteria
		if (fabs(g1) < eps and fabs(g2) < eps)
			break;

		// Finding Newton direction: -inv(H') * g
		double det = h11 * h22 - h21 * h21;
		double dA = -( h22 * g1 - h21 * g2) / det;
		double dB = -(-h21 * g1 + h11 * g2) / det;
		double gd = g1 * dA + g2 * dB;

		double stepsize = 1;		// Line Search
		while (stepsize >= min_step)
		{
			double newA = A + stepsize * dA;
			double newB = B + stepsize * dB;

			// New function value
			double newf = 0.0;
			for (size_t i = 0; i < l; ++i)
			{
				fApB = dec_value[i] * newA + newB;
				if (fApB >= 0)
					newf += t[i] * fApB + log(1 + exp(-fApB));
				else
					newf += (t[i] - 1) * fApB + log(1 + exp(fApB));
			}

			// Check sufficient decrease
			if (newf < fval + 0.0001 * stepsize * gd)
			{
				A = newA; B = newB; fval = newf;
				break;
			}
			else
				stepsize = stepsize / 2.0;
		}

		if (stepsize < min_step)
		{
			cerr << "Line search fails in two-class probability estimates" << endl;
			break;
		}
	}

	if (iter >= max_iter)
		cerr << "Reaching maximal iterations in two-class probability estimates" << endl;

	return make_tuple(A, B);
}

double SigmoidPredict(double decision_value, double A, double B)
{
	double fApB = decision_value * A + B;
	// 1 - p used later; avoid catastrophic cancellation
	if (fApB >= 0)
		return exp(-fApB) / (1.0 + exp(-fApB));
	else
		return 1.0 / (1 + exp(fApB));
}

// Method 2 from the multiclass_prob paper by Wu, Lin, and Weng
vector<double> MultiClassProbability(size_t k, const vector<vector<double>>& r)
{
	size_t max_iter = 100;
	if (max_iter < k)
		max_iter = k;

	vector<vector<double>> Q(k, vector<double>(k));
	vector<double> Qp(k);
	vector<double> p(k);

	double eps = 0.005 / k;
	
	for (size_t t = 0; t < k; ++t)
	{
		p[t] = 1.0 / k;  // Valid if k = 1
		Q[t][t] = 0;

		for (size_t j = 0; j < t; ++j)
		{
			Q[t][t] += r[j][t] * r[j][t];
			Q[t][j] = Q[j][t];
		}

		for (size_t j = t + 1; j < k; ++j)
		{
			Q[t][t] += r[j][t] * r[j][t];
			Q[t][j] = -r[j][t] * r[t][j];
		}
	}

	size_t iter;
	for (iter = 0; iter < max_iter; ++iter)
	{
		// stopping condition, recalculate QP,pQP for numerical accuracy
		double pQp = 0;

		for (size_t t = 0; t < k; ++t)
		{
			Qp[t] = 0;
			for (size_t j = 0; j < k; ++j)
				Qp[t] += Q[t][j] * p[j];
			pQp += p[t] * Qp[t];
		}

		double max_error = 0;
		for (size_t t = 0; t < k; ++t)
		{
			double error = fabs(Qp[t] - pQp);
			if (error > max_error)
				max_error = error;
		}

		if (max_error < eps)
			break;
		
		for (size_t t = 0; t < k; ++t)
		{
			double diff = (-Qp[t] + pQp) / Q[t][t];
			p[t] += diff;
			pQp = (pQp + diff * (diff * Q[t][t] + 2 * Qp[t])) / (1 + diff) / (1 + diff);
			for (size_t j = 0; j < k; ++j)
			{
				Qp[j] = (Qp[j] + diff * Q[t][j]) / (1 + diff);
				p[j] /= (1 + diff);
			}
		}
	}

	if (iter >= max_iter)
		cerr << "Exceeds max_iter in multiclass_prob" << endl;

	return p;
}

template<KernelType K>
vector<double> Model<K, C_SVC>::PredictWithProbability(const Vector& x) const
{
	size_t nr_class = this->m_class.size();

	vector<double> sums;
	sums.reserve(nr_class * (nr_class - 1) / 2);
	
	Predict(x, sums);
	
	double min_prob = 1e-7;
	vector<vector<double>> pairwise_prob(nr_class, vector<double>(nr_class));
	
	for (size_t i = 0, k = 0; i < nr_class; ++i)
	{
		for (size_t j = i + 1; j < nr_class; ++j)
		{
			double p = SigmoidPredict(sums[k], m_probA[k], m_probB[k]);

			if (p < min_prob)
				p = min_prob;
			if (p > 1 - min_prob)
				p = 1 - min_prob;
			
			pairwise_prob[i][j] = p;
			pairwise_prob[j][i] = 1 - p;
			++k;
		}
	}
	
	return MultiClassProbability(nr_class, pairwise_prob);
//	vector<double> result = MultiClassProbability(nr_class, pairwise_prob);
//
//	auto voted = distance(result.begin(),
//		max_element(result.begin(), result.end()));
//	
//	return this->m_class[voted].label;
}

template<KernelType K>
ModelBase* SVM<K,C_SVC>::Train(const vector<label_type>& y, const Data& m)
{
	size_t l = y.size();
	
	vector<label_type> label;
	vector<size_t> perm, count, start;
	
	GroupClasses(y, label, count, start, perm);

	if (cif::VERBOSE > 2)
	{
		cout << "Nr of classes: " << label.size() << endl
			 << " => { ";
		for (auto l: label)
			cout << int32_t(l) << ' ';
		cout << '}' << endl;
	}
	
	size_t nr_class = label.size();
	size_t N = nr_class * (nr_class - 1) / 2;
	
	if (nr_class == 1)
		throw runtime_error("Training data in only one class");
	
	View x(m, perm);
	
	vector<double> weighted_C(nr_class, this->m_params.C);
	for (auto weight: m_weight)
	{
		auto i = find(label.begin(), label.end(), weight.label);
		if (i == label.end())
			throw runtime_error("Class label " + to_string(weight.label) + " specified in weight is not found");
		weighted_C[i - label.begin()] *= weight.weight;
	}
	
	vector<bool> nonzero(l, false);
	vector<decision_function> f(N);
	
	vector<double> probA, probB;
	
	if (this->m_params.probability)
	{
		probA = vector<double>(N);
		probB = vector<double>(N);
	}
	
	size_t p = 0;
	for (size_t i = 0; i < nr_class; ++i)
	{
		for (size_t j = i + 1; j < nr_class; ++j)
		{
			size_t si = start[i], sj = start[j];
			size_t ci = count[i], cj = count[j];
			
			size_t l = ci * cj;
			
			View sx(x, vector<size_t>());	sx.reserve(l);
			vector<label_type> sy;				sy.reserve(l);
			
			for (size_t k = 0; k < ci; ++k)
			{
				sx.push_back(x[si + k]);
				sy.push_back(+1);
			}
			
			for (size_t k = 0; k < cj; ++k)
			{
				sx.push_back(x[sj + k]);
				sy.push_back(-1);
			}
			
			if (this->m_params.probability)
				tie(probA[p], probB[p]) = BinarySSVCProbability(sy, sx, weighted_C[i], weighted_C[j]);
			
			f[p] = TrainOne(sy, sx, weighted_C[i], weighted_C[j]);
			
			for (size_t k = 0; k < ci; ++k)
			{
				if (not nonzero[si + k] and fabs(f[p].alpha[k]) > 0)
					nonzero[si + k] = true;
			}

			for (size_t k = 0; k < cj; ++k)
			{
				if (not nonzero[sj + k] and fabs(f[p].alpha[ci + k]) > 0)
					nonzero[sj + k] = true;
			}
			
			++p;
		}
	}
	
	vector<double> rho(N);
	transform(f.begin(), f.end(), rho.begin(),
		[](const decision_function& df) -> double { return df.rho; });

	// build model
	typedef ModelBase::ClassInfo class_info_type;
	vector<class_info_type> classes;
	
	size_t totalSV = 0;
	for (size_t i = 0; i < nr_class; ++i)
	{
		class_info_type ci = { label[i] };
		
		for (size_t j = 0; j < count[i]; ++j)
			if (nonzero[start[i] + j])
			{
				++ci.nr_sv;
				++totalSV;
			}
		
		classes.push_back(ci);
	}
	
//		if (cif::VERBOSE > 1)
//			cout << "Total nSV = " << totalSV << endl;

	Matrix sv;
	for (size_t i = 0; i < l; ++i)
	{
		if (nonzero[i])
			sv.push_back(x[i]);
	}
	
	vector<size_t> nz_start({ 0 });
	for (size_t i = 1; i < nr_class; ++i)
		nz_start.push_back(nz_start.back() + classes[i - 1].nr_sv);
	
	vector<vector<double>> sv_coef(totalSV, vector<double>(nr_class - 1));
	
	p = 0;
	for (size_t i = 0; i < nr_class; ++i)
		for (size_t j = i + 1; j < nr_class; ++j)
		{
			size_t si = start[i], sj = start[j];
			size_t ci = count[i], cj = count[j];
			
			size_t q = nz_start[i];
			for (size_t k = 0; k < ci; ++k)
				if (nonzero[si + k])
					sv_coef[q++][j - 1] = f[p].alpha[k];
			q = nz_start[j];
			for (size_t k = 0; k < cj; ++k)
				if (nonzero[sj + k])
					sv_coef[q++][i] = f[p].alpha[ci + k];
			++p;
		}
	
	model_type* result;

	if (this->m_params.probability)
		result = new model_type(this->m_params, sv, classes, rho, sv_coef, probA, probB);
	else
		result = new model_type(this->m_params, sv, classes, rho, sv_coef);
	
	return result;
}

// Cross-validation decision values for probability estimates

template<KernelType K>
tuple<double,double> SVM<K,C_SVC>::BinarySSVCProbability(
	const vector<label_type>& y, const Data& x, double Cp, double Cn)
{
	size_t l = y.size();
	const size_t nr_fold = 5;
	
	vector<size_t> perm(l);
	iota(perm.begin(), perm.end(), 0);
	
	mt19937 rng(m_rd());
	shuffle(perm.begin(), perm.end(), rng);
	
	vector<double> dec_value(l);
	
	for (size_t i = 0; i < nr_fold; ++i)
	{
		auto begin = i * l / nr_fold;
		auto end = (i + 1) * l / nr_fold;
		
		auto subl = l - (end - begin);
		View subx(x, vector<size_t>());
		vector<label_type> suby;
		suby.reserve(subl);
		
		for (size_t j = 0; j < begin; ++j)
		{
			subx.push_back(x[perm[j]]);
			suby.push_back(y[perm[j]]);
		}
		
		for (size_t j = end; j < l; ++j)
		{
			subx.push_back(x[perm[j]]);
			suby.push_back(y[perm[j]]);
		}
		
		size_t p_count = 0, n_count = 0;
		for (auto l: suby)
			if (l > 0)
				++p_count;
			else
				++n_count;
		
		if (p_count == 0 or n_count == 0)
		{
			int dec_v = 0;
			if (p_count > 0)
				dec_v = 1;
			if (n_count > 0)
				dec_v = -1;
			
			for (size_t j = begin; j < end; ++j)
				dec_value[perm[j]] = dec_v;
		}
		else
		{
			param_type params = this->m_params;
			params.C = 1;
			params.probability = false;
			
			SVM sub_svm(params);
			sub_svm.m_weight = vector<weight_type>(
				{
					{ +1, Cp },
					{ -1, Cn }
				});
			
			unique_ptr<model_type> submodel(
				static_cast<model_type*>(sub_svm.Train(suby, subx)));
			
			for (size_t j = begin; j < end; ++j)
			{
				vector<double> sdv;
				submodel->Predict(x[perm[j]], sdv);
				
				assert(sdv.size() == 1);
				dec_value[perm[j]] = sdv[0] * suby[0];
			}
		}
	}
	
	return SigmoidTrain(dec_value, y);
}

template<KernelType K>
vector<label_type> SVM<K,C_SVC>::CrossValidation(const vector<label_type>& y, 
	const Matrix& x, size_t nr_fold)
{
	size_t l = x.size();
	
	if (nr_fold > l)
	{
		nr_fold = l;
		cerr << "WARNING: # folds > # data. Will use # folds = # data instead (i.e., leave-one-out cross validation)" << endl;
	}
	
	vector<label_type> result(l);

	mt19937 rng(m_rd());
	
	vector<label_type> label;
	vector<size_t> perm, count, start;
	
	GroupClasses(y, label, count, start, perm);

	size_t nr_class = label.size();
	
	// we shuffle the indices for the various classes and take nr_fold samples out of them
	vector<size_t> index = perm;
	for (size_t c = 0; c < nr_class; ++c)
		shuffle(index.begin() + start[c], index.begin() + start[c] + count[c], rng);
	
	size_t last_fold_start = 0;
	
	for (size_t i = 0; i < nr_fold; ++i)
	{
		// calculate the size of this fold
		auto fold_count = accumulate(count.begin(), count.end(), 0UL,
			[i,nr_fold](size_t s, size_t c) -> size_t { return s + (i + 1) * c / nr_fold - i * c / nr_fold; });
		
		size_t fold_start = last_fold_start;
		last_fold_start += fold_count;

		// adjust the permutation index for this fold
		auto pi = perm.begin() + fold_start;
		
		for (size_t c = 0; c < nr_class; ++c)
		{
			auto begin = start[c] + i * count[c] / nr_fold;
			auto end = start[c] + (i + 1) * count[c] / nr_fold;

			for (auto j = begin; j < end; ++j)
				*pi++ = index[j];
		}

		// construct a model based on the 'rest' of the samples, the ones in the fold will be used for checking
		auto begin = fold_start;
		auto end = fold_start + fold_count;
		
		View sx(x, vector<size_t>());
		vector<label_type> sy;
		
		for (size_t j = 0; j < begin; ++j)
		{
			sx.push_back(x[perm[j]]);
			sy.push_back(y[perm[j]]);
		}
		for (size_t j = end; j < l; ++j)
		{
			sx.push_back(x[perm[j]]);
			sy.push_back(y[perm[j]]);
		}
		
		auto model = unique_ptr<ModelBase>(Train(sy, sx));

		// if probability
		// else
		for (size_t j = begin; j < end; ++j)
			result[perm[j]] = model->Predict(x[perm[j]]);
	}
	
	return result;
}

template<KernelType K>
void SVM<K,C_SVC>::GroupClasses(const vector<label_type>& y,
	vector<label_type>& label, vector<size_t>& count,
	vector<size_t>& start, vector<size_t>& perm)
{
	size_t l = y.size();
	
	vector<size_t> data_label;
	
	for (auto ly: y)
	{
		auto i = find(label.begin(), label.end(), ly);
		if (i != label.end())
		{
			data_label.push_back(i - label.begin());
			count[i - label.begin()] += 1;
		}
		else
		{
			data_label.push_back(label.size());
			label.push_back(ly);
			count.push_back(1);
		}
	}
	
	size_t nr_class = label.size();
	
	if (nr_class == 2 and label[0] == -1 and label[1] == 1)
	{
		swap(label[0], label[1]);
		swap(count[0], count[1]);
		for_each(data_label.begin(), data_label.end(),
			[](size_t& s) { s = (s == 0 ? 1 : 0); });
	}
	
	start = vector<size_t>(nr_class);
	start[0] = 0;
	for (size_t i = 1; i < nr_class; ++i)
		start[i] = start[i - 1] + count[i - 1];
	
	perm = vector<size_t>(l);
	for (size_t i = 0; i < l; ++i)
	{
		perm[start[data_label[i]]] = i;
		++start[data_label[i]];
	}
	
	start[0] = 0;
	for (size_t i = 1; i < nr_class; ++i)
		start[i] = start[i - 1] + count[i - 1];
}

template<KernelType K>
typename SVMBase<K>::decision_function SVM<K,C_SVC>::TrainOne(const vector<label_type>& y, const Data& x,
	double Cp, double Cn)
{
	size_t l = x.size();
	
	decision_function result = { vector<double>(l, 0) };
	vector<int8_t> ny(l);
	
	transform(y.begin(), y.end(), ny.begin(),
		[](label_type l) -> int8_t { return l > 0 ? +1 : -1; });

	typedef QMatrix<C_SVC,RBF> QM;

const double eps = 1e-3;
const double tau = 1e-12;

	QM Q(ny, x, this->m_params);
	Solver<QM> s(Q, eps, tau, this->m_params.C);
	
	double obj;
	tie(result.alpha, result.rho, obj) = s.Solve(y);
	
	if (cif::VERBOSE > 1 and Cp == Cn)
	{
		double sum_alpha = accumulate(
			result.alpha.begin(), result.alpha.end(), 0.0);
		cout << "nu = " << (sum_alpha / (Cp * l)) << endl;
	}
	
	size_t nSV = 0, nBSV = 0;
	for (size_t i = 0; i < l; ++i)
	{
		if (fabs(result.alpha[i]) > 0)
		{
			++nSV;
			if ((y[i] > 0 and fabs(result.alpha[i] >= Cp)) or
				(y[i] < 0 and fabs(result.alpha[i] >= Cn)))
			{
				++nBSV;
			}
		}
		
		result.alpha[i] *= y[i];
	}
	
	if (cif::VERBOSE > 1)
		cout << "obj = " << obj << ", rho = " << result.rho << endl
			 << "nSV = " << nSV << ", nBSV = " << nBSV << endl; 

	return result;
}

// explicitly instantiate the various classes here 
template class SVM<RBF,C_SVC>;

}
