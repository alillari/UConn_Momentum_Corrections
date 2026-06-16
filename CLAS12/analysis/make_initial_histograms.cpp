#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TString.h>
#include <TList.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include "Math/Vector4D.h"

#include <iostream>
#include <vector>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstring>

#include <MomCorrParticle.h>
#include <MomCorrConfig.h>

#include <Kinematics.h>
#include <ExclusiveKinematics.h>
#include <MissingMassRDF.h>
#include <DeltaPRDF.h>

int main(int argc, char* argv[]) {
    ROOT::EnableImplicitMT();

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <directory/path/> <num_files | all> <config.json>\n";
        return 1;
    }

    TString pathPattern = argv[1];

    TSystemDirectory dir("inputDir", pathPattern);
    TList* files = dir.GetListOfFiles();

    if (!files) {
        std::cerr << "Could not list files in directory: " << pathPattern << "\n";
        return 1;
    }

    std::vector<TString> rootFiles;
    TIter next(files);
    TSystemFile* file;

    while ((file = static_cast<TSystemFile*>(next()))) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root")) {
            rootFiles.push_back(pathPattern + "/" + fname);
        }
    }

    if (rootFiles.empty()) {
        std::cerr << "No ROOT files found in directory: " << pathPattern << "\n";
        return 1;
    }

    int numFilesToProcess = 0;

    if (TString(argv[2]) == "all") {
        numFilesToProcess = rootFiles.size();
    } else {
        numFilesToProcess = std::stoi(argv[2]);

        if (numFilesToProcess <= 0 ||
            numFilesToProcess > static_cast<int>(rootFiles.size())) {
            std::cerr << "Invalid number of files requested.\n";
            return 1;
        }
    }

    std::string configFile = argv[3];
    MomCorrConfig config(configFile);

    TChain chain("h22");

    for (int i = 0; i < numFilesToProcess; ++i) {
        chain.Add(rootFiles[i]);
    }

    const double beam_energy = config.GetBeamEnergy();
    const TString channel_name = config.GetChannel();
    const TString bending = config.GetBending();
    const TString specifiers = config.GetSpecifier();

    const double missing_mass_low = config.GetMissingMassLow();
    const double missing_mass_high = config.GetMissingMassHigh();
    const double missing_mass_width = config.GetMissingMassBinWidth();

    const double dp_low = config.GetDpLow();
    const double dp_high = config.GetDpHigh();
    const double dp_bin_width = config.GetDpBinWidth();

    const int missing_mass_bins =
        static_cast<int>((missing_mass_high - missing_mass_low) / missing_mass_width);

    const int dp_bins =
        static_cast<int>((dp_high - dp_low) / dp_bin_width);

    TString dataOutLocation = "../analysis_out/" + specifiers + "/";
    gSystem->mkdir(dataOutLocation, true);

    TString outputFileName =
        dataOutLocation + specifiers + "_raw_histograms.root";

    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Could not create output file: " << outputFileName << "\n";
        return 1;
    }

    TTree config_tree("config", "Analysis configuration");

    std::ifstream in(configFile);
    std::ostringstream ss;
    ss << in.rdbuf();
    std::string configs_string = ss.str();

    config_tree.Branch("json", &configs_string);
    config_tree.Fill();
    config_tree.Write();

    TTree* meta = new TTree("hist_meta", "Histogram metadata");

    int hist_id;
    int meta_sector;
    int is_phi_binning;
    int meta_phi_bin;

    char meta_particle[16];
    char meta_detector[16];

    char x_var[32];

    int x_bins;
    double x_min;
    double x_max;

    char y_var[32];

    int y_bins;
    double y_min;
    double y_max;

    char meta_hist_name[128];

    meta->Branch("hist_id", &hist_id);
    meta->Branch("sector", &meta_sector);
    meta->Branch("is_phi_binning", &is_phi_binning);
    meta->Branch("phi_bin", &meta_phi_bin);

    meta->Branch("particle", meta_particle, "particle/C");
    meta->Branch("detector", meta_detector, "detector/C");

    meta->Branch("x_var", x_var, "x_var/C");
    meta->Branch("y_var", y_var, "y_var/C");

    meta->Branch("x_bins", &x_bins);
    meta->Branch("x_min", &x_min);
    meta->Branch("x_max", &x_max);

    meta->Branch("y_bins", &y_bins);
    meta->Branch("y_min", &y_min);
    meta->Branch("y_max", &y_max);

    meta->Branch("hist_name", meta_hist_name, "hist_name/C");

    int unique_hist_id = 0;

    auto FillMeta = [&](int sector,
                        int phi_binning,
                        int phi_bin,
                        const MomCorrParticle& particle,
                        const std::string& xv,
                        int xb,
                        double xmin,
                        double xmax,
                        const std::string& yv,
                        int yb,
                        double ymin,
                        double ymax,
                        const std::string& hist_name) {
        hist_id = unique_hist_id;
        meta_sector = sector;
        is_phi_binning = phi_binning;
        meta_phi_bin = phi_bin;

        std::strncpy(meta_particle, particle.GetName().c_str(), sizeof(meta_particle));
        meta_particle[sizeof(meta_particle) - 1] = '\0';

        std::strncpy(meta_detector, particle.GetDetector().c_str(), sizeof(meta_detector));
        meta_detector[sizeof(meta_detector) - 1] = '\0';

        std::strncpy(x_var, xv.c_str(), sizeof(x_var));
        x_var[sizeof(x_var) - 1] = '\0';

        x_bins = xb;
        x_min = xmin;
        x_max = xmax;

        std::strncpy(y_var, yv.c_str(), sizeof(y_var));
        y_var[sizeof(y_var) - 1] = '\0';

        y_bins = yb;
        y_min = ymin;
        y_max = ymax;

        std::strncpy(meta_hist_name, hist_name.c_str(), sizeof(meta_hist_name));
        meta_hist_name[sizeof(meta_hist_name) - 1] = '\0';

        meta->Fill();
    };

    ROOT::RDataFrame df_base(chain);
    ROOT::RDF::RNode df = df_base;

    std::vector<MomCorrParticle> particle_list = config.BuildParticles();

    for (const auto& particle : particle_list) {
        df = particle.AddBranches(df);
    }

    df = DefineMissingMass(df, particle_list, beam_energy);
    df = momcorr::DefineDeltaP_epipN(df, particle_list, beam_energy);

    for (const auto& particle : particle_list) {
        const std::string mom_branch = particle.GetMagBranch();
        const std::string sector_branch = particle.GetSectorBranch();

        const double mom_low = particle.GetMomentumMin();
        const double mom_high = particle.GetMomentumMax();
        const double mom_bin_width = particle.GetMomentumBinWidth();

        const int mom_bins =
            static_cast<int>((mom_high - mom_low) / mom_bin_width);

        const bool make_delta_p = (particle.GetPID() == 211);

        std::vector<std::pair<std::string, std::string>> dp_branches;

        if (make_delta_p) {
            dp_branches = {
                {"deltaP_fourvec", particle.GetName() + "_deltaP_fourvec"},
                {"deltaP_legacy",  particle.GetName() + "_deltaP_legacy"},
                {"deltaP_diff",    particle.GetName() + "_deltaP_diff"}
            };
        }

        for (int sector : particle.GetSectors()) {
            auto df_sector =
                df.Filter(sector_branch + " == " + std::to_string(sector));

            if (particle.IsPhiBinningEnabled()) {
                const std::string phi_bin_branch = particle.GetPhiBinBranch();

                for (const auto& [phi_bin, label] : particle.GetPhiBinningLabels()) {
                    auto df_phi =
                        df_sector.Filter(phi_bin_branch + " == " + std::to_string(phi_bin));

                    const std::string mm_hist_name =
                        "hMM_vs_" + particle.GetName() +
                        "_sec" + std::to_string(sector) +
                        "_phi" + label +
                        "_id" + std::to_string(unique_hist_id);

                    const std::string mm_hist_title =
                        "MM vs " + particle.GetName() +
                        " Momentum [Sector " + std::to_string(sector) +
                        ", Phi Bin: " + label +
                        "]; Momentum (GeV); MM (GeV/c^{2})";

                    auto h_mm = df_phi.Histo2D(
                        {mm_hist_name.c_str(),
                         mm_hist_title.c_str(),
                         mom_bins, mom_low, mom_high,
                         missing_mass_bins, missing_mass_low, missing_mass_high},
                        mom_branch,
                        "missing_mass"
                    );

                    h_mm->Write();

                    FillMeta(
                        sector,
                        1,
                        phi_bin,
                        particle,
                        "p",
                        mom_bins,
                        mom_low,
                        mom_high,
                        "missing_mass",
                        missing_mass_bins,
                        missing_mass_low,
                        missing_mass_high,
                        mm_hist_name
                    );

                    unique_hist_id++;

                    if (make_delta_p) {
                        for (const auto& [tag, branch] : dp_branches) {
                            const std::string dp_hist_name =
                                "h_" + tag + "_vs_" + particle.GetName() +
                                "_sec" + std::to_string(sector) +
                                "_phi" + label +
                                "_id" + std::to_string(unique_hist_id);

                            const std::string dp_hist_title =
                                "#Delta p " + tag +
                                " vs " + particle.GetName() +
                                " Momentum [Sector " + std::to_string(sector) +
                                ", Phi Bin: " + label +
                                "]; Momentum (GeV); #Delta p (GeV)";

                            auto h_dp = df_phi.Histo2D(
                                {dp_hist_name.c_str(),
                                 dp_hist_title.c_str(),
                                 mom_bins, mom_low, mom_high,
                                 dp_bins, dp_low, dp_high},
                                mom_branch,
                                branch
                            );

                            h_dp->Write();

                            FillMeta(
                                sector,
                                1,
                                phi_bin,
                                particle,
                                "p",
                                mom_bins,
                                mom_low,
                                mom_high,
                                tag,
                                dp_bins,
                                dp_low,
                                dp_high,
                                dp_hist_name
                            );

                            unique_hist_id++;
                        }
                    }
                }
            } else {
                const std::string mm_hist_name =
                    "hMM_vs_" + particle.GetName() +
                    "_sec" + std::to_string(sector) +
                    "_id" + std::to_string(unique_hist_id);

                const std::string mm_hist_title =
                    "MM vs " + particle.GetName() +
                    " Momentum [Sector " + std::to_string(sector) +
                    "]; Momentum (GeV); MM (GeV/c^{2})";

                auto h_mm = df_sector.Histo2D(
                    {mm_hist_name.c_str(),
                     mm_hist_title.c_str(),
                     mom_bins, mom_low, mom_high,
                     missing_mass_bins, missing_mass_low, missing_mass_high},
                    mom_branch,
                    "missing_mass"
                );

                h_mm->Write();

                FillMeta(
                    sector,
                    0,
                    -1,
                    particle,
                    "p",
                    mom_bins,
                    mom_low,
                    mom_high,
                    "missing_mass",
                    missing_mass_bins,
                    missing_mass_low,
                    missing_mass_high,
                    mm_hist_name
                );

                unique_hist_id++;

                if (make_delta_p) {
                    for (const auto& [tag, branch] : dp_branches) {
                        const std::string dp_hist_name =
                            "h_" + tag + "_vs_" + particle.GetName() +
                            "_sec" + std::to_string(sector) +
                            "_id" + std::to_string(unique_hist_id);

                        const std::string dp_hist_title =
                            "#Delta p " + tag +
                            " vs " + particle.GetName() +
                            " Momentum [Sector " + std::to_string(sector) +
                            "]; Momentum (GeV); #Delta p (GeV)";

                        auto h_dp = df_sector.Histo2D(
                            {dp_hist_name.c_str(),
                             dp_hist_title.c_str(),
                             mom_bins, mom_low, mom_high,
                             dp_bins, dp_low, dp_high},
                            mom_branch,
                            branch
                        );

                        h_dp->Write();

                        FillMeta(
                            sector,
                            0,
                            -1,
                            particle,
                            "p",
                            mom_bins,
                            mom_low,
                            mom_high,
                            tag,
                            dp_bins,
                            dp_low,
                            dp_high,
                            dp_hist_name
                        );

                        unique_hist_id++;
                    }
                }
            }
        }
    }

    meta->Write();
    outputFile->Close();

    std::cout << "Wrote raw histograms to: " << outputFileName << "\n";

    return 0;
}
