#pragma once
#include <string>

class TF1;
struct FitResult;
struct PeakGuess;

class FitModel {
public:
    explicit FitModel(int bg_order);

    TF1* makeTF1(const std::string& name, double xmin, double xmax) const;

    int bgOrder() const { return bg_order_; }

    void setParameters(TF1* f, PeakGuess signal_guess) const;

    FitResult extractParameters(TF1* f) const;

private:
    int bg_order_;
};
