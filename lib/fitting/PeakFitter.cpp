#include "fitting/PeakFitter.h"
#include "fitting/FitModel.h"

#include <TH1.h>
#include <TF1.h>

PeakFitter::PeakFitter(const FitConfig& cfg)
: cfg_(cfg) {}

FitResult PeakFitter::fit(TH1* h) const
{
    FitResult best;
    best.chi2_ndf = 1e9;

    if (h->GetEntries() < cfg_.min_entries) {
        best.skipped = true;
        return best;
    }

    PeakGuess guess = guessPeak(h);

    for (int order : cfg_.bg_orders) {
        FitResult r = tryModel(h, guess, order);
        if (!r.success) continue;

        if (r.chi2_ndf < best.chi2_ndf) {
            best = r;
        }
    }

    return best;
}
