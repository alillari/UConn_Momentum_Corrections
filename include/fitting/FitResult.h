#pragma once

struct FitResult {
    bool success = false;
    bool skipped = false;

    int bg_order = -1;

    double mean = 0;
    double mean_err = 0;
    double sigma = 0;

    double chi2_ndf = 0;
    bool suspicious = false;
};
