#include "libpr.h"

#include <fstream>

#include <zeep/xml/serialize.hpp>
#include <zeep/xml/document.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "svm++.h"

using namespace std;
namespace ba = boost::algorithm;

//struct svm_class
//{
//	int8_t label;
//	size_t nr_sv;
//	
//	template<class Archive>
//	void serialize(Archive& ar, const unsigned int version)
//	{
//		ar & zeep::xml::make_attribute_nvp("label", label)
//		   & zeep::xml::make_attribute_nvp("nr_sv", nr_sv);
//	}
//};
//
//struct svm_sv
//{
//	int index;
//	double value;
//
//	template<class Archive>
//	void serialize(Archive& ar, const unsigned int version)
//	{
//		ar & zeep::xml::make_attribute_nvp("index", index)
//		   & zeep::xml::make_element_nvp(".", value);
//	}
//};
//
//struct svm_node
//{
//	vector<double> sv_coef;
//	vector<svm_sv> sv;
//
//	template<class Archive>
//	void serialize(Archive& ar, const unsigned int version)
//	{
//		ar & zeep::xml::make_element_nvp("sv_coef", sv_coef)
//		   & zeep::xml::make_element_nvp("sv", sv);
//	}
//};
//
//struct svm_config
//{
//	string					svm_type;
//	string					kernel_type;
//	boost::optional<double>	gamma;
//	vector<double>			rho;
//	vector<svm_class>		classes;
//	vector<svm_node>		sv;
//
//	template<class Archive>
//	void serialize(Archive& ar, const unsigned int version)
//	{
//		ar & zeep::xml::make_element_nvp("svm-type", svm_type)
//		   & zeep::xml::make_element_nvp("kernel-type", kernel_type)
//		   & zeep::xml::make_element_nvp("gamma", gamma)
//		   & zeep::xml::make_element_nvp("rho", rho)
//		   & zeep::xml::make_element_nvp("class", classes)
//		   & zeep::xml::make_element_nvp("svm-node", sv);
//	}
//};

void predict_file(svm::ModelBase* model, const char* filename)
{
//	typedef svm::ModelBase	SVMModel;
	typedef svm::Vector		SVMVector;

	ifstream test(filename);
	string line;

	while (getline(test, line))
	{
		list<string> f;
		ba::split(f, line, ba::is_any_of(" "));
		f.pop_front();
		
		SVMVector v;
		for (string& fs: f)
		{
			auto c = fs.find(':');
			if (c == string::npos)
				continue;
			
			int index = boost::lexical_cast<int>(fs.substr(0, c));
			float value = boost::lexical_cast<float>(fs.substr(c + 1));
			
			v[index] = value;
		}

		char predictedClass = model->Predict(v);
	
		auto prob = model->PredictWithProbability(v);

		auto voted = distance(prob.begin(),
			max_element(prob.begin(), prob.end()));
		
		cout << predictedClass << '\t'
			 << (voted ? '0' : '1') << '\t'
			 << "probabilities: ";
		copy(prob.begin(), prob.end(), ostream_iterator<double>(cout, ", "));
		cout << endl;
	}
}

//int centrifuge_predict(int argc, char* argv[])
//{
//	ifstream file("1cbs-0.8.svm");
//	
////	svm_config config;
//	
//	zeep::xml::document doc(file);
////	doc.deserialize("svm-config", config);
//
//	auto config = doc.find_first("//svm-config");
//	if (not config)
//		throw runtime_error("invalid svm file");
//	
//	auto model = svm::ModelBase::Create(*config);
//
//	predict_file(model, "1ctn-0.8-scaled.txt");
//
//	return 0;
//}

int centrifuge_test(int argc, char* argv[])
{
	cif::VERBOSE = 1;
	
	ifstream data("1cbs-0.8.txt.scale");
	if (not data.is_open())
		throw runtime_error("no such file");

	typedef svm::SVM_C_SVC_RBF SVM;
	typedef SVM::param_type SVMParams;
	typedef svm::Matrix SVMMatrix;

	SVMParams params;
	params.gamma = 0.5;
	params.C = 8;
	params.probability = true;

	SVM svm(params);

	vector<int8_t> labels;
	SVMMatrix m;
	size_t row = 0;

	for (;;)
	{
		string line;
		getline(data, line);
		
		if (line.empty())
		{
			if (data.eof())
				break;
			continue;
		}
		
		vector<string> f;
		ba::split(f, line, ba::is_any_of(" "));
		
		if (f.size() < 1)
			throw runtime_error("invalid data");
		
		int8_t c = boost::lexical_cast<int8_t>(f.front());
		f.erase(f.begin(), f.begin() + 1);
		
		labels.push_back(c);
		
		for (auto& fs: f)
		{
			auto c = fs.find(':');
			if (c == string::npos)
				continue;
			
			int index = boost::lexical_cast<int>(fs.substr(0, c));
			float value = boost::lexical_cast<float>(fs.substr(c + 1));

			m[row][index - 1] = value;
		}
		
		++row;
	}

	auto model = svm.Train(labels, m);
	
//	model->Print();

	auto cv = svm.CrossValidation(labels, m, 5);
	
	size_t correct = 0;
	for (size_t i = 0; i < labels.size(); ++i)
		if (cv[i] == labels[i])
			++correct;
	
	cout << "Cross validation Accuracy = " << (100.0 * correct / labels.size()) << '%' << endl;

	predict_file(model, "1ctn-0.8-scaled.txt");
	
	zeep::xml::document doc;
	doc.root()->append(model->GetConfig());
	
	ofstream model_file("test.model");
	if (model_file.is_open())
		model_file << doc;

	return 0;
}
