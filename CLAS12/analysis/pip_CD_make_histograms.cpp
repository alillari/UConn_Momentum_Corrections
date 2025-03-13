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


namespace fs = std::filesystem;

int main(int argc, char* argv[]){
    ROOT::EnableImplicitMT();

    //Logic for correct usage of command

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <directory/path/*.root> <num_files | all>\n";
        return 1;
    }

    TString pathPattern = argv[1];
    TString fileList = gSystem->GetFromPipe(Form("ls %s 2>/dev/null", pathPattern.Data()));

    std::vector<TString> rootFiles;
    std::istringstream stream(fileList.Data());
    TString file;
    while (stream >> file) {
        rootFiles.push_back(file);
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

    TChain chain("h22");  // TTree name is "h22", change as needed
    for (int i = 0; i < numFilesToProcess; ++i) {
        chain.Add(rootFiles[i]);
    } 

    //Specify all the information you feel is necessary to recognize your root file, I recommend DETECTOR_SUBSYSTEM_RUNGROUP_RUNPERIOD_MAGNET_PARTICLE
    //If you're not doing CLAS12, make your own specifier
    //This is just for your own book-keeping, change it as required or to taste :)
    //Change dataOutLocation to your heart's content, just make sure you know where it is and run the enviromental sh script with the correct location

    TString specifiers = "CLAS12_CD_RGA_Sp19_In_pip";
    TString dataOutLocation = "../analysis_out/" + specifiers + "/";
    TString outputFileName = dataOutLocation + specifiers + "_raw_histograms.root";
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    //RDataFrame and particle lines below
    
    ROOT::RDataFrame df(chain);

    //Beam and particle information, adjust for your data and particles
    //Current beam energy is for RGA Sp19
    const double beam_energy = 10.1998;
    const double proton_mass = 0.938;
    const double El_mass = 0.000511;
    const double Pro_mass = 0.938272088;
    const double Pip_mass = 0.140;

    //Making TLorentzVectors from beam and target information. Don't touch.
    ROOT::Math::PxPyPzMVector beam_vector(0,0,sqrt(pow(beam_energy,2) - pow(El_mass,2)),El_mass);
    ROOT::Math::PxPyPzMVector target_vector(0,0,0,Pro_mass);

    //Make your the particles in your dataset here!
    //Particles are made with the constructor as seen below
    //It goes:
    //	name, mass (GeV), branch names for momentums, branch name for sector, vector for sectors, detector (2 for FD, 3 for CD), min mommentum (GeV), max mommentum (GeV), mommentum bin width (GeV)
    //NOTE: Be generous with range and number of bins. You can always decrease the range and increase bin size in the fitting code
    //	    But you can't increase the number of bins after
    

    double mom_bin = .05;

    double El_mom_low = 2;
    double El_mom_high = 9;
    bool El_phi_flag = true;


    double Pip_mom_low = 0;
    double Pip_mom_high = 10;
    bool Pip_phi_flag = false;

    const std::vector<int> six_sector = {1,2,3,4,5,6};
    int El_detector = 2;
    int Pip_detector = 3;

    //Define your phi shift. I pulled this from Richard's code. Good luck.
    //Variables passed in need to be phi, momentum, sector. Otherwise code can't call them correctly.
    auto El_compute_local_phi  = [](double ElPhi, double El, int esec) {
    	double localPhi = ElPhi - (esec - 1) * 60;
    	return localPhi - (30 / El);
    };

    //Determine what your phi bins are
    double El_phi_divider = 5;
    auto El_phi_binning = [El_phi_divider](double localPhi) {
	if(localPhi <= -El_phi_divider){ return 1; }
	else if(-El_phi_divider < localPhi && localPhi <= El_phi_divider){ return 2; }
	else { return 3; }
    };

    std::unordered_map<int, std::string> El_phi_bin_map = {{1, "negative"}, {2, "neutral"}, {3,"positive"}};


    MomCorrParticle Electron("El", El_mass, "e_px", "e_py", "e_pz", "esec", El_detector, six_sector, El_mom_low, El_mom_high, mom_bin, El_phi_flag, El_compute_local_phi, El_phi_binning, El_phi_bin_map);
    //MomCorrParticle Pip("Pip", Pip_mass, "pip_px", "pip_py", "pip_pz", "pipsec", six_sector, Pip_detector, Pip_mom_low, Pip_mom_high, mom_bin, Pip_phi_flag);

    //std::vector<MomCorrParticle> particle_list = {electron, pip};

    return 0;
}
