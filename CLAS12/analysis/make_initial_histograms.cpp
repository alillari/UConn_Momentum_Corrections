#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystemDirectory.h>
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

#include <MomCorrParticle.h>
#include <MomCorrConfig.h>

#include <Kinematics.h>
#include <ExclusiveKinematics.h>
#include <MissingMassRDF.h>

int main(int argc, char* argv[]){
    ROOT::EnableImplicitMT();

    //Logic for correct usage of command

    if (argc < 4) {
    	std::cerr << "Usage: " << argv[0] << " <directory/path/> <num_files | all> <config.json>\n";
    	return 1;
    }

    TString pathPattern = argv[1];
    
    TSystemDirectory dir("inputDir", pathPattern);
    TList* files = dir.GetListOfFiles();

    std::vector<TString> rootFiles;
    TIter next(files);
    TSystemFile* file;

    while ((file = (TSystemFile*)next())) {
    	TString fname = file->GetName();
    	if (!file->IsDirectory() && fname.EndsWith(".root")) {
        	rootFiles.push_back(
            	pathPattern + "/" + fname
        	);
    	}
    }
    
    if (rootFiles.empty()) {
        std::cerr << "No ROOT files found matching pattern: " << pathPattern << "\n";
        return 1;
    }

    int numFilesToProcess = 0;
    if (TString(argv[2]) == "all") {
        numFilesToProcess = rootFiles.size();
    } else {
        numFilesToProcess = std::stoi(argv[2]);
        if (numFilesToProcess <= 0 || numFilesToProcess > static_cast<int>(rootFiles.size())) {
            std::cerr << "Invalid number of files requested.\n";
            return 1;
        }
    }

    //TODO: Add logic for checking if the last input is actually a JSON file
	
    std::string configFile = argv[3];
    MomCorrConfig config(configFile);

    TChain chain("h22");  // TTree name is "h22", change as needed
			  // h22 is the default TTree name for the UConn group, for some reason
    for (int i = 0; i < numFilesToProcess; ++i) {
        chain.Add(rootFiles[i]);
    } 

    const double beam_energy = config.GetBeamEnergy();
    const TString channel_name = config.GetChannel();
    const TString bending = config.GetBending();
    const TString specifiers = config.GetSpecifier();

    const double mom_bin = config.GetDefaultMomentumBin();

    const double missing_mass_low  = config.GetMissingMassLow();
    const double missing_mass_high = config.GetMissingMassHigh();
    const double missing_mass_width = config.GetMissingMassBinWidth();
    
    const double dp_low = config.GetDpLow();
    const double dp_high = config.GetDpHigh();
    const double dp_bin_width = config.GetDpBinWidth();

    TString dataOutLocation = "../analysis_out/" + specifiers + "/";
    TString outputFileName = dataOutLocation + specifiers + "_raw_histograms.root";
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    //Storing JSON metadata inside the ROOT file so one doesn't need to keep track of the JSON files seperately from their output
    TTree config_tree("config", "Analysis configuration");
    
    std::ifstream in(configFile);
    std::ostringstream ss;
    ss << in.rdbuf();
    std::string configs_string = ss.str();

    config_tree.Branch("json", &configs_string);
    config_tree.Fill();
    config_tree.Write();

    //Storing histogram metadata
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

    ROOT::RDataFrame df_base(chain);
    ROOT::RDF::RNode df = df_base;

    //Fundamental particle constant
    const double El_mass = 0.000511;
    const double Mu_mass = 0.105658;
    const double Pro_mass = 0.938272088;
    const double Pip_mass = 0.140;

    std::vector<MomCorrParticle> particle_list = config.BuildParticles();

    std::vector<std::string> momentum_columns;
    //Main logic loop for creating histograms
    for(auto& particle: particle_list){
    	df = particle.AddBranches(df);
    }

    df = DefineMissingMass(df, particle_list, beam_energy);
    //df = DefineDP(df, particle_list, beam_energy);

    int unique_hist_id = 0;

    for (const auto& particle : particle_list) {
    	std::string mom_branch = particle.GetName() + "_mag";
    	std::string sector_branch = particle.GetSectorBranch();

    	double mom_low = particle.GetMomentumMin();
    	double mom_high = particle.GetMomentumMax();
    	double mom_bin_width = particle.GetMomentumBinWidth();
    	int mom_bins = static_cast<int>((mom_high - mom_low) / mom_bin_width);

    	int missing_mass_bins = static_cast<int>((missing_mass_high - missing_mass_low) / missing_mass_width);

    	for (int sector : particle.GetSectors()) {
    	    // Filter the dataframe for the current sector
    	    auto df_sector = df.Filter(sector_branch + " == " + std::to_string(sector));

    	    if (particle.IsPhiBinningEnabled()) {
    	        std::string phi_bin_branch = particle.GetName() + "_phiBin";
    	        for (const auto& [phi_bin, label] : particle.GetPhiBinningLabels()) {
    	            auto df_phi = df_sector.Filter(phi_bin_branch + " == " + std::to_string(phi_bin));

    	            auto h = df_phi.Histo2D(
    	                {("hMM_vs_" + particle.GetName() + "_sec" + std::to_string(sector) + "_phi" + label + "_id" + std::to_string(unique_hist_id)).c_str(),
    	                 ("MM vs " + particle.GetName() + " Momentum [Sector " + std::to_string(sector) + ", Phi Bin: " + label + "]; Momentum (GeV); MM (GeV/c^2)").c_str(),
    	                 mom_bins, mom_low, mom_high, missing_mass_bins, missing_mass_low, missing_mass_high},
    	                mom_branch, "missing_mass"
    	            );
		    h->Write();

                    hist_id  = unique_hist_id;
		    meta_sector   = sector;
		    is_phi_binning  = 1;
		    meta_phi_bin = phi_bin;

		    strncpy(meta_particle, particle.GetName().c_str(), sizeof(meta_particle));
		    meta_particle[sizeof(meta_particle)-1] = '\0';

		    strncpy(meta_detector, particle.GetDetector().c_str(), sizeof(meta_detector));
		    meta_detector[sizeof(meta_detector)-1] = '\0';

		    strncpy(x_var, "p", sizeof(x_var));
		    x_var[sizeof(x_var)-1] = '\0';

		    x_bins = mom_bins;
		    x_min = mom_low;
		    x_max = mom_high;

		    strncpy(y_var, "missing_mass", sizeof(y_var));
		    y_var[sizeof(y_var)-1] = '\0';

		    y_bins = missing_mass_bins;
		    y_min = missing_mass_low;
		    y_max = missing_mass_high;

		    strncpy(meta_hist_name, ("hMM_vs_" + particle.GetName() + "_sec" + std::to_string(sector) + "_phi" + label + "_id" + std::to_string(unique_hist_id)).c_str(), sizeof(meta_hist_name));
		    meta_hist_name[sizeof(meta_hist_name)-1] = '\0';

		    meta->Fill();

		    unique_hist_id++;

    	        }
    	    } else {
    	        auto h = df_sector.Histo2D(
    	            {("hMM_vs_" + particle.GetName() + "_sec" + std::to_string(sector) + "_id" + std::to_string(unique_hist_id)).c_str(),
    	             ("MM vs " + particle.GetName() + " Momentum [Sector " + std::to_string(sector) + "]; Momentum (GeV); MM (GeV/c^2)").c_str(),
    	             mom_bins, mom_low, mom_high, missing_mass_bins, missing_mass_low, missing_mass_high},
    	            mom_branch, "missing_mass"
    	        );
		h->Write();

		hist_id  = unique_hist_id;
                meta_sector   = sector;
                is_phi_binning  = 0;
                meta_phi_bin = -1;

                strncpy(meta_particle, particle.GetName().c_str(), sizeof(meta_particle));
                meta_particle[sizeof(meta_particle)-1] = '\0';

                strncpy(meta_detector, particle.GetDetector().c_str(), sizeof(meta_detector));
                meta_detector[sizeof(meta_detector)-1] = '\0';

                strncpy(x_var, "p", sizeof(x_var));
                x_var[sizeof(x_var)-1] = '\0';

                x_bins = mom_bins;
                x_min = mom_low;
                x_max = mom_high;

                strncpy(y_var, "missing_mass", sizeof(y_var));
                y_var[sizeof(y_var)-1] = '\0';

                y_bins = missing_mass_bins;
                y_min = missing_mass_low;
                y_max = missing_mass_high;

                strncpy(meta_hist_name, ("hMM_vs_" + particle.GetName() + "_sec" + std::to_string(sector) + "_id" + std::to_string(unique_hist_id)).c_str(), sizeof(meta_hist_name));
		meta_hist_name[sizeof(meta_hist_name)-1] = '\0';

                meta->Fill();

		unique_hist_id++;
    	    }
    	}
    }

    meta->Write();

    outputFile->Close();

    return 0;
}
