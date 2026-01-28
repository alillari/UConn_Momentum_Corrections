#pragma once

#include "fitting/FitConfig.h"
#include "fitting/FitResult.h"

class TH1;

class PeakFitter {
public:
    explicit PeakFitter(const FitConfig& cfg);

    FitResult fit(TH1* h) const;

private:
    FitConfig cfg_;

    struct PeakGuess {
        double mean;
        double sigma;
        double amplitude;
    };

    PeakGuess guessPeak(TH1* h) const;
    FitResult tryModel(TH1* h, const PeakGuess& guess, int bg_order) const;
};
