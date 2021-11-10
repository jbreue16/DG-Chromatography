#include <iostream>
#include <Eigen/Dense>
#include <string>

#include "gnuplot-iostream.h"

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
	void line_plot(Eigen::VectorXd data, Eigen::VectorXd xAxis) {
		std::vector<boost::tuple<double, double> > pts_C;
		for (int i = 0; i < data.size(); i++) {
			pts_C.push_back(boost::make_tuple(xAxis[i], data[i]));
		}
		gp << "plot '-' with lines\n";// << gp.file1d(cache.c) << endl;
		gp.send(pts_C);
		std::cin.get();
	}
	void line_plot(Eigen::VectorXd data, Eigen::VectorXd xAxis, std::string_view title, std::string_view xlabel, std::string_view ylabel) {
		this->title(title);
		this->xlabel(xlabel);
		this->ylabel(ylabel);
		this->line_plot(data, xAxis);
	}
};