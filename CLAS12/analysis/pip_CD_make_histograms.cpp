#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
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

#include <MomCorrParticle.h>
#include <ComputeMissingMassHelper.h>
#include <MomCorrConfig.h>

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
    const TString specifiers = config.GetSpecifiers();

    const double mom_bin = config.GetDefaultMomentumBin();

    const double missing_mass_low  = config.GetMissingMassMin();
    const double missing_mass_high = config.GetMissingMassMax();
    const double missing_mass_width = config.GetMissingMassBinWidth();
    
    const double dp_low = config.GetDpMin();
    const double dp_high = config.GetDpHigh();
    const double dp_bin_width = config.GetDpBinWidth();

    TString dataOutLocation = "../analysis_out/" + specifiers + "/";
    TString outputFileName = dataOutLocation + specifiers + "_raw_histograms.root";
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

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
    	                {("hMM_vs_" + particle.GetName() + "_sec" + std::to_string(sector) + "_phi" + label).c_str(),
    	                 ("MM vs " + particle.GetName() + " Momentum [Sector " + std::to_string(sector) + ", Phi Bin: " + label + "]; Momentum (GeV); MM (GeV/c^2)").c_str(),
    	                 mom_bins, mom_low, mom_high, missing_mass_bins, missing_mass_low, missing_mass_high},
    	                mom_branch, "missing_mass"
    	            );
		    h->Write();

		    //auto h_DP = df_phi.Histo2D(
		    //	{().c_str()},

    	        }
    	    } else {
    	        auto h = df_sector.Histo2D(
    	            {("hMM_vs_" + particle.GetName() + "_sec" + std::to_string(sector)).c_str(),
    	             ("MM vs " + particle.GetName() + " Momentum [Sector " + std::to_string(sector) + "]; Momentum (GeV); MM (GeV/c^2)").c_str(),
    	             mom_bins, mom_low, mom_high, missing_mass_bins, missing_mass_low, missing_mass_high},
    	            mom_branch, "missing_mass"
    	        );
		h->Write();
    	    }
    	}
    }

    outputFile->Close();

    return 0;
}
