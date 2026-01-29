#include "fitting/PeakFitter.h"
#include "fitting/FitModel.h"

#include <TH1.h>
#include <TF1.h>

//Setting global limits of fit sigma
constexpr double MIN_SIGMA_FRACTION = 0.3;
constexpr double MAX_SIGMA_FRACTION = 3.0;

PeakFitter::PeakFitter(const FitConfig& cfg)
: cfg_(cfg) {}

PeakFitter::PeakGuess PeakFitter::guessPeak(TH1* h) const {
    PeakGuess g;

    // fallback defaults
    g.mean  = h->GetMean();
    g.sigma = cfg_.expected_sigma;
    g.amplitude = h->GetMaximum();

    std::unique_ptr<TH1> h_smooth(static_cast<TH1*>(h->Clone("h_smooth")));
    h_smooth->Smooth(1);

    TSpectrum spec(3);
    int n_peaks = spec.Search(h_smooth.get(), 2, "", 0.05);

    if (n_peaks <= 0)
        return g;

    double* peaks = spec.GetPositionX();

    int best_bin = -1;
    double best_height = -1;

    for (int i = 0; i < n_peaks; ++i) {
        int bin = h->FindBin(peaks[i]);
        double height = h->GetBinContent(bin);
        if (height > best_height) {
            best_height = height;
            best_bin = bin;
        }
    }

    if (best_bin < 0)
        return g;

    g.mean = h->GetBinCenter(best_bin);
    g.amplitude = best_height;

    // local RMS estimate
    int lo = h->FindBin(g.mean - 2 * cfg_.expected_sigma);
    int hi = h->FindBin(g.mean + 2 * cfg_.expected_sigma);

    double sum = 0, sumw = 0;
    for (int i = lo; i <= hi; ++i) {
        double x = h->GetBinCenter(i);
        double w = h->GetBinContent(i);
        sum += w * (x - g.mean) * (x - g.mean);
        sumw += w;
    }

    if (sumw > 0)
        g.sigma = std::sqrt(sum / sumw);

    g.sigma = std::clamp(g.sigma, 0.5 * cfg_.expected_sigma, 2.0 * cfg_.expected_sigma);

    return g;
}

FitResult PeakFitter::tryModel(TH1* h, PeakGuess guess, int order) const {
    FitResult result;

    //TODO: determine how fits should be named
    TF1* f = makeTF1(...);

    //TODO: Figure out how to set initial parameters



    int status = h->Fit(f, "QNR");

    if (status != 0) {
        delete f;
        return result;
    }

    //TODO: Figure out how to extract the correct parameters


    //TODO: Apply guardrails?

    delete f;
    return result;	

}

FitResult PeakFitter::fit(TH1* h) const {
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
