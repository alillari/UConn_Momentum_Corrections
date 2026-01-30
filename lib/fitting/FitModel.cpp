#include "fitting/FitModel.h"
#include "fitting/FitResult.h"
#include "fitting/PeakFitter.h"

#include <TF1.h>
#include <string>

FitModel::FitModel(int bg_order) 
	: bg_order_(bg_order) {}

TF1* FitModel::makeTF1(const std::string& name, double xmin, double xmax) const {
    std::string formula = "gaus(0)";
    formula += " + pol" + std::to_string(bg_order_) + "(3)";

    return new TF1(name.c_str(), formula.c_str(), xmin, xmax);
}

void FitModel::setParameters(TF1* f, const PeakGuess& signal_guess) const {
	f->SetParameter(0, signal_guess.amplitude);
	f->SetParameter(1, signal_guess.mean);
	f->SetParameter(2, signal_guess.sigma);

	int bg_start = 3;
	for (int i = 0; i <= bg_order_; ++i) {
    		f->SetParameter(bg_start + i, 0.0);
	}
	
}

FitResult FitModel::extractParameters(TF1* f) const {

	FitResult result;

	result.bg_order = bg_order_;

	result.mean = f->GetParameter(1);	
	result.mean_err = f->GetParError(1);

	result.sigma = f->GetParameter(2);
	result.sigma_err = f->GetParError(2);

	int ndf = f->GetNDF();
	result.chi2_ndf = (ndf > 0) ? f->GetChisquare() / ndf : -1.0;

	return result;

}
