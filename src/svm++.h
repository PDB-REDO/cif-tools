/*
	based on code from libsvm, which requires the following statement: 
	
	Copyright (c) 2000-2014 Chih-Chung Chang and Chih-Jen Lin
	All rights reserved.

*/

#pragma once

#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>
#include <random>

#include <cassert>

#include "zeep/xml/serialize.hpp"

namespace svm
{

// Code based on libsvm.

//	All label and value types are fixed. int8_t and float should be sufficient anyway.

typedef signed char label_type;		// The label, should be sufficient
typedef float value_type;			// The value type for the input data

// to avoid typing too much
using std::vector;
using std::array;
using std::tuple;

// forward declarations for the data classes
class Data;
class Matrix;
class View;

// --------------------------------------------------------------------
//	Data classes, to store the vectors and matrices
//
// 	Vector is the basis. It stores the observed and scaled values for a sample.
//	This is a candidate for improvement. Currently memory is wasted in sparse situations.

template<int N>
class Vector_
{
  public:

	struct reference
	{
		reference(Vector_& v, size_t index)
			: m_v(v), m_ix(index) {}
		
		reference& operator=(float v)
		{
			if (std::isnan(v))
				m_v.m_used[m_ix] = false;
			else
			{
//				if (m_ix >= m_v.m_used.size())
//				{
//					size_t d = m_ix - m_v.m_used.size() + 1;
//					m_v.m_used.insert(m_v.m_used.end(), d, false);
//					m_v.m_values.insert(m_v.m_values.end(), d, 0);
//				}

				m_v.m_values[m_ix] = v;
				m_v.m_used[m_ix] = true;
			}
			
			return *this;
		}
		
		operator value_type () const
		{
			float result;
			if (m_v.m_used[m_ix])
				result = m_v.m_values[m_ix];
			else
				result = nan("");
			return result;
		}
	
	  private:
		Vector_&	m_v;
		size_t		m_ix;
	};

	Vector_() { std::fill(m_used.begin(), m_used.end(), false); }
	Vector_(const Vector_& rhs)
		: m_used(rhs.m_used), m_values(rhs.m_values) {}
	Vector_& operator=(const Vector_& rhs)
	{
		m_used = rhs.m_used;
		m_values = rhs.m_values;
		return *this;
	}

	Vector_(Vector_&& rhs)
		: m_used(std::move(rhs.m_used)), m_values(std::move(rhs.m_values)) {}

	template<class InputIter>
	Vector_(InputIter b, InputIter e)
	{
		size_t l = std::distance(b, e);
		assert(l <= N);
		
		for (size_t i = 0; i < l; ++i)
		{
			m_values[i] = *b;
			m_used[i] = not std::isnan(*b);
			++b;
		}  
	}

	Vector_& operator=(Vector_&& rhs)
	{
		m_used = rhs.m_used;
		m_values = std::move(rhs.m_values);
		return *this;
	}
	
	size_t size() const
	{
		return m_values.size();
	}
	
	reference operator[](size_t index)
	{
		assert(index < N);
		return reference(*this, index);
	}
	
	value_type operator[](size_t index) const
	{
		assert(contains(index));
		return m_values[index];
	}
	
	void fill(value_type v)
	{
		m_used.fill(true);
		m_values.fill(v);
	}

	bool contains(size_t index) const
	{
		return index < m_used.size() and m_used[index];
	}
	
	double dot(const Vector_& rhs) const
	{
		double result = 0;
		for (size_t i = 0; i < N; ++i)
		{
			if (m_used[i] and rhs.m_used[i])
				result += m_values[i] * rhs.m_values[i];
		}
		
		return result;
	}

  private:
	array<bool,N>		m_used;
	array<value_type,N>	m_values;
};

typedef Vector_<9>	Vector;

inline double dot(const Vector& a, const Vector& b)
{
	return a.dot(b);
}

// --------------------------------------------------------------------
//	Data is an abstract base class for Matrix and View

class Data
{
  public:
	virtual ~Data() {}

	typedef Vector		vector_type;
	
	virtual size_t size() const = 0;
	virtual vector_type& operator[](size_t row) = 0;
	virtual const vector_type& operator[](size_t row) const = 0;
};

// --------------------------------------------------------------------
//	Matrix contains actual data.

class Matrix : public Data
{
  public:
	friend class View; 

	typedef typename vector<Vector>::iterator iterator;
	typedef typename vector<Vector>::const_iterator const_iterator;

	Matrix() {}
	
	virtual size_t size() const { return m_data.size(); }
	
	virtual Vector& operator[](size_t row)
	{
		if (row >= m_data.size())
			throw std::runtime_error("Row index out of range for matrix");
		return m_data[row];
	}
	
	virtual const Vector& operator[](size_t row) const
	{
		if (row >= m_data.size())
			throw std::runtime_error("Row index out of range for matrix");
		return m_data[row];
	}
	
	iterator begin()				{ return m_data.begin(); }
	iterator end()				 	{ return m_data.end(); }

	const_iterator begin() const	{ return m_data.begin(); }
	const_iterator end() const		{ return m_data.end(); }

	void push_back(const Vector& v)
	{
		m_data.push_back(v);
	}

	void push_back(Vector&& v)
	{
		m_data.push_back(std::move(v));
	}

	// Scale each column to the range [0..1] returning the
	// min and max values observed for each column.
	vector<std::pair<value_type,value_type>> Scale();

  private:
	vector<Vector>	m_data;
};

// --------------------------------------------------------------------
//	KernelTypes, there's only one at the moment

enum KernelType
{
	RBF
};

template<KernelType K>
struct KernelTypeTraits;

template<>
struct KernelTypeTraits<RBF>
{
	static const char* GetKernelTypeString() { return "rbf"; }
};

template<KernelType K> class Kernel {};

// --------------------------------------------------------------------
//	RBF Kernel definition

template<>
class Kernel<RBF>
{
  public:
	typedef KernelTypeTraits<RBF>		kernel_type_traits;

	struct param_type
	{
		double gamma;
		
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version)
		{
			ar & zeep::xml::make_element_nvp("gamma", gamma);
		}
	};

	Kernel(const Data& x, const param_type& params)
		: m_x(x), m_gamma(params.gamma) {}
	
	Kernel(const Kernel&) = delete;
	Kernel& operator=(const Kernel&) = delete;
	
	static double f(const Vector& x, const Vector& y, const param_type& params);
	double operator()(std::size_t i, std::size_t j) const;
	
	const Data&	m_x;
	double		m_gamma;
};

// --------------------------------------------------------------------
//	SVMTypes, only one at the moment: C_SVC

enum SVMType
{
	C_SVC
};

template<SVMType S>
struct SVMTypeTraits;

template<>
struct SVMTypeTraits<C_SVC>
{
	static const char* GetSVMTypeString() { return "c_svc"; }
};

template<KernelType K, SVMType S>
class SVM
{
};

// --------------------------------------------------------------------
//	ModelBase, the abstract base class for Models

class ModelBase
{
  public:
	virtual ~ModelBase() {}

	struct ClassInfo
	{
		label_type		label;
		size_t			nr_sv;
	};

	ModelBase(const ModelBase&) = delete;
	ModelBase& operator=(const ModelBase&) = delete;

	virtual label_type Predict(const Vector& x) const = 0;
	virtual vector<double> PredictWithProbability(const Vector& x) const;

	static ModelBase* Create(zeep::xml::element& n);
	
	virtual zeep::xml::element* GetConfig();

  protected:

	ModelBase(const Matrix& sv, const vector<ClassInfo>& classes,
		const vector<double>& rho, const vector<vector<double>>& sv_coef)
		: m_sv(sv), m_class(classes), m_rho(rho), m_sv_coef(sv_coef)
	{
	}
	
	ModelBase(zeep::xml::element& config);
	
	Matrix					m_sv;
	vector<ClassInfo>		m_class;
	vector<double>			m_rho;
	vector<vector<double>>	m_sv_coef;
};

template<KernelType K, SVMType S>
class Model
{
};

template<KernelType K>
class Model<K, C_SVC> : public ModelBase
{
  public:
	typedef ModelBase::ClassInfo					class_info_type;
	typedef Kernel<K>								kernel_type;
	typedef SVM<K,C_SVC>							svm_type;
	typedef typename svm_type::param_type			param_type;

	Model(const param_type& params,
		const Matrix& sv, const vector<class_info_type>& classes,
		const vector<double>& rho, const vector<vector<double>>& sv_coef,
		const vector<double>& probA, const vector<double>& probB)
		: ModelBase(sv, classes, rho, sv_coef), m_params(params)
		, m_probA(probA), m_probB(probB)
	{
	}

	Model(const param_type& params,
		const Matrix& sv, const vector<class_info_type>& classes,
		const vector<double>& rho, const vector<vector<double>>& sv_coef)
		: ModelBase(sv, classes, rho, sv_coef), m_params(params)
	{
	}

	Model(zeep::xml::element& config)
		: ModelBase(config)
	{
		zeep::xml::deserializer sr(config);
		sr.deserialize_element("params", m_params);
		sr.deserialize_element("probA", m_probA);
		sr.deserialize_element("probB", m_probB);
	}
	
	Model(const Model&) = delete;
	Model& operator=(const Model&) = delete;
	
	virtual label_type Predict(const Vector& x) const;
	virtual vector<double> PredictWithProbability(const Vector& x) const;
	void Predict(const Vector& x, vector<double>& sums) const;

	virtual zeep::xml::element* GetConfig();

  private:
	param_type			m_params;
	vector<double>		m_probA, m_probB;
};

// --------------------------------------------------------------------

template<KernelType K>
class SVMBase
{
  public:
	typedef Kernel<K>							kernel_type;
	
	struct decision_function
	{
		vector<double>	alpha;
		double			rho;
	}; 

	struct weight_type
	{
		label_type	label;
		double		weight;
	};

	struct param_type : public kernel_type::param_type
	{
		double C;
		bool probability;
	};

	SVMBase(const param_type& params)
		: m_params(params)
	{
	}

	virtual ModelBase* Train(const vector<label_type>& y, const Data& m) = 0;
	virtual vector<label_type> CrossValidation(const vector<label_type>& y, 
		const Matrix& m, size_t nr_fold) = 0;

  protected:
	param_type m_params;
};

// --------------------------------------------------------------------

template<KernelType K>
class SVM<K, C_SVC> : public SVMBase<K>
{
  public:
	typedef SVMBase<K>								base_type;
	typedef Kernel<K>								kernel_type;
	typedef typename base_type::param_type			param_type;
	typedef Model<K,C_SVC>							model_type;
	typedef typename base_type::decision_function	decision_function;
	typedef typename base_type::weight_type			weight_type;
	
	SVM(const param_type& params)
		: base_type(params)
	{
	}

	SVM(const SVM&) = delete;
	SVM& operator=(const SVM&) = delete;

	virtual ModelBase* Train(const vector<label_type>& y, const Data& m);
	virtual vector<label_type> CrossValidation(const vector<label_type>& y, 
		const Matrix& x, size_t nr_fold);

  private:

	void GroupClasses(const vector<label_type>& y,
		vector<label_type>& label, vector<std::size_t>& count,
		vector<std::size_t>& start, vector<size_t>& perm);
	decision_function TrainOne(const vector<label_type>& y, const Data& x,
		double Cp, double Cn);

	tuple<double,double> BinarySSVCProbability(const vector<label_type>& y, const Data& x,
		double Cp, double Cn); 

	vector<weight_type>	m_weight;
	std::random_device	m_rd;
};

typedef SVM<RBF, C_SVC> SVM_C_SVC_RBF;

}
