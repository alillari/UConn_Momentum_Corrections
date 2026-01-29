#pragma once

#include "fitting/FitConfig.h"
#include "fitting/FitResult.h"
#include "fitting/FitModel.h"

class TH1;

class PeakFitter {
public:
    explicit PeakFitter(const FitConfig& cfg);

    FitResult fit(TH1* h) const;

private:
    FitConfig cfg_;

    PeakGuess guessPeak(TH1* h) const;	//I am going to keep this here as it doesn't care about if this peak is gaussian or not
    FitResult tryModel(TH1* h, const PeakGuess& guess, int bg_order) const;
};
