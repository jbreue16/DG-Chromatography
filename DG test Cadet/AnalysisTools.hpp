#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <fstream>

#include "gnuplot-iostream.h"

#include "DGspecific.hpp"

/**
* @file includes functions for analysis and gnuplot API
*/


/**
* Plotter API for gnuplot
*/
class Plotter {
public:
	Gnuplot gp;
	Plotter() : gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"") {}

	Plotter& title(std::string_view title) {
		gp << "set title '" << title << "'\n";
		return *this;
	}
	Plotter& xlabel(std::string_view xlabel) {
		gp << "set xlabel '" << xlabel << "'\n";
		return *this;
	}
	Plotter& ylabel(std::string_view ylabel) {
		gp << "set ylabel '" << ylabel << "'\n";
		return *this;
	}
	void oneLine_plot(Eigen::VectorXd data, Eigen::VectorXd xAxis) {
		std::vector<boost::tuple<double, double> > pts_C;
		for (int i = 0; i < data.size(); i++) {
			pts_C.push_back(boost::make_tuple(xAxis[i], data[i]));
		}
		gp << "plot '-' with lines\n";
		gp.send(pts_C);
		std::cin.get();
	}
	void twoLine_plot(Eigen::VectorXd data1, Eigen::VectorXd data2, Eigen::VectorXd xAxis) {
		std::vector<boost::tuple<double, double, double> > pts_C;
		for (int i = 0; i < data1.size(); i++) {
			pts_C.push_back(boost::make_tuple(xAxis[i], data1[i], data2[i]));
		}
		gp << "plot '-' with lines\n";
		gp.send(pts_C);
		std::cin.get();
	}
	void oneLine_plot(Eigen::VectorXd data, Eigen::VectorXd xAxis, std::string_view title, std::string_view xlabel, std::string_view ylabel) {
		this->title(title);
		this->xlabel(xlabel);
		this->ylabel(ylabel);
		this->oneLine_plot(data, xAxis);
	}
	void twoLine_plot(Eigen::VectorXd data1, Eigen::VectorXd data2, Eigen::VectorXd xAxis, std::string_view title, std::string_view xlabel, std::string_view ylabel) {
		this->title(title);
		this->xlabel(xlabel);
		this->ylabel(ylabel);
		this->twoLine_plot(data1, data2, xAxis);
	}
};

void write_csv(Eigen::VectorXd vec, std::string name) {
	std::ofstream myFile(name);
	//myFile << "c" << "\n"; // column name
	for (int i = 0; i < vec.size(); ++i)
	{
		myFile << vec[i] << "\n";
	}
	myFile.close();
}
Eigen::VectorXd read_csv(std::string name, int size) {
	std::ifstream myFile(name);
	if (!myFile.is_open()) throw std::runtime_error("Could not open file");
	Eigen::VectorXd result(size);
	int count = 0;
	while (count < size) {
		myFile >> result[count];
		count++;
	}
	myFile.close();
	return result;
}

/**
* @brief returns state vector of one component
*/
VectorXd getComponent(VectorXd state, ParameterProvider para, int Comp) {
	VectorXd Component = VectorXd::Zero(state.size() / para.nComp);
	for (int i = 0; i < para.nCells * para.nNodes; i++) {
		Component[i] = state[i * para.strideNode() + Comp * para.strideComp()];
	}
	return Component;
}