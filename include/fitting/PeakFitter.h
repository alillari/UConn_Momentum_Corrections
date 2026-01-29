#pragma once

#include "fitting/FitConfig.h"
#include "fitting/FitResult.h"
#include "fitting/FitModel.h"

class TH1;

struct PeakGuess {
        double mean;
        double sigma;
        double amplitude;
};

class PeakFitter {
public:
    explicit PeakFitter(const FitConfig& cfg);

    FitResult fit(TH1* h) const;

private:
    FitConfig cfg_;

    PeakGuess guessPeak(TH1* h) const;
    FitResult tryModel(TH1* h, const PeakGuess& guess, int bg_order) const;
};
