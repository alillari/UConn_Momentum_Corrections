#include "fitting/FitModel.h"
#include <TF1.h>

FitModel::FitModel(int bg_order) 
	: bg_order_(bg_order) {}

TF1* FitModel::makeTF1(const std::string& name, double xmin, double xmax) const
{
    std::string formula = "gaus(0)";
    formula += " + pol" + std::to_string(bg_order_) + "(3)";

    return new TF1(name.c_str(), formula.c_str(), xmin, xmax);
}
