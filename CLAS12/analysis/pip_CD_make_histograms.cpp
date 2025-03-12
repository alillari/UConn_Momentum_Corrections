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

#include <MomCorrParticle.h>

int main(int argc, char* argv[]){
   
    //Logic for correct usage of command

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <directory/path/*.root>\n";
        return 1;
    }

    std::string pathPattern = argv[1];
    std::string directory = fs::path(pathPattern).parent_path().string();
    std::string filenamePattern = fs::path(pathPattern).filename().string();

    if (!fs::exists(directory) || !fs::is_directory(directory)) {
        std::cerr << "Error: Directory does not exist.\n";
        return 1;
    }

    std::vector<std::string> rootFiles;
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".root" && 
            entry.path().filename().string().find(filenamePattern.substr(0, filenamePattern.find('*'))) == 0) {
            rootFiles.push_back(entry.path().string());
        }
    }

    if (rootFiles.empty()) {
        std::cerr << "No ROOT files found matching the pattern.\n";
        return 1;
    }

    std::cout << "Found " << rootFiles.size() << " ROOT file(s):\n";
    for (const auto& file : rootFiles) {
        std::cout << "  " << file << "\n";
    }

    //Check if the ROOT file has the right tree name, change if relevant
    //h22 is standard in this group for some reason? I suggest you conform to it
    TTree* tree = nullptr;
    TString treeName("h22");
    file.GetObject(treeName, tree);
    if (!tree) {
        std::cerr << "Did not find specified TTree in file.\n";
        return 1;
    }

    //Define TChain and output file
    TChain* inputChain = new TChain(treeName);
    inputChain->Add(pathPattern.c_str());

    //Specify all the information you feel is necessary to recognize your root file, I recommend DETECTOR_SUBSYSTEM_RUNGROUP_RUNPERIOD_MAGNET_PARTICLE
    //If you're not doing CLAS12, make your own specifier
    //This is just for your own book-keeping, change it as required or to taste :)
    //Change dataOutLocation to your heart's content, just make sure you know where it is and run the enviromental sh script with the correct location

    TString specifiers = "CLAS12_CD_RGA_Sp19_In_pip";
    TString dataOutLocation = "../analysis_out/" + specifiers + "/";
    TString outputFileName = dataOutLocation + specifiers + "_raw_histograms.root";
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    //RDataFrame and particle lines below
    
    ROOT::RDataFrame df(inputChain);

    //Beam and particle information, adjust for your data
    //Current beam energy is for RGA Sp19
    const double beam_energy = 10.1998;
    const double proton_mass = 0.938;
    const double el_mass = 0.000511;
    const double pip_mass = 0.140;

    //Make your the particles in your dataset here!
    //Particles are made with the constructor as seen below
    //It goes 
    //	name, mass (GeV), branch names for momentums, branch name for sector, vector for sectors, detector (2 for FD, 3 for CD), min mommentum (GeV), max mommentum (GeV), mommentum bin width (GeV)
    //NOTE: Be generous with range and number of bins. You can always decrease the range and increase bin size in the fitting code
    //	    But you can't increase the number of bins after
    

    double mom_bin = .05;

    double el_mom_low = 2;
    double el_mom_high = 9;
    bool el_phi_flag = true;

    double pip_mom_low = 0;
    double pip_mom_high = 10;
    bool pip_phi_flag = false;

    std::vector<int> six_sector = {1,2,3,4,5,6};
    int el_detector = 2;
    int pip_detector = 3;

    MomCorrParticle electron("el", el_mass, "e_px", "e_py", "e_pz", "esec", six_sector, el_detector, el_mom_low, el_mom_high, mom_bin, el_phi_flag, );
    MomCorrParticle pip("pip", pip_mass, "pip_px", "pip_py", "pip_pz", "pipsec", six_sector, pip_detector, pip_mom_low, pip_mom_high, mom_bin);

    std::vector<MomCorrParticle> particleList = {electron, pip};

    
}
