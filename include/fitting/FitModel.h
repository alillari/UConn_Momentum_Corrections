#pragma once
#include <string>

class TF1;

class FitModel {
public:
    explicit FitModel(int bg_order);

    TF1* makeTF1(const std::string& name,
                 double xmin,
                 double xmax) const;

    int bgOrder() const { return bg_order_; }

private:
    int bg_order_;
};
