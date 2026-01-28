#pragma once
#include <vector>
#include <string>

struct FitConfig {
    std::string name;

    double xmin;
    double xmax;
    double expected_sigma;

    int min_entries = 50;
    std::vector<int> bg_orders;

    bool allow_no_peak = true;
};
