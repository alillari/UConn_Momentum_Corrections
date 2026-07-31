#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <TFile.h>
#include <TInterpreter.h>
#include <TList.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <TString.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <MomCorrConfig.h>
#include <MomCorrParticle.h>
#include <ObservableRegistry.h>

namespace {

std::vector<TString> FindRootFiles(const TString& directory) {
    TSystemDirectory dir("inputDir", directory);
    TList* files = dir.GetListOfFiles();

    if (!files) {
        throw std::runtime_error(
            "Could not list files in directory: " + std::string(directory.Data()));
    }

    std::vector<TString> rootFiles;
    TIter next(files);
    TSystemFile* file = nullptr;

    while ((file = static_cast<TSystemFile*>(next()))) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root")) {
            rootFiles.push_back(directory + "/" + fname);
        }
    }

    std::sort(rootFiles.begin(), rootFiles.end());
    if (rootFiles.empty()) {
        throw std::runtime_error(
            "No ROOT files found in directory: " + std::string(directory.Data()));
    }
    return rootFiles;
}

int ParseFileCount(const char* arg, std::size_t available) {
    if (TString(arg) == "all") {
        return static_cast<int>(available);
    }

    std::size_t consumed = 0;
    const int requested = std::stoi(arg, &consumed);
    if (consumed != std::strlen(arg) || requested <= 0 ||
        requested > static_cast<int>(available)) {
        throw std::runtime_error("Invalid number of files requested");
    }
    return requested;
}

std::set<std::string> BranchNames(TChain& chain) {
    std::set<std::string> names;
    TObjArray* branches = chain.GetListOfBranches();
    if (!branches) {
        throw std::runtime_error("Could not inspect input tree branches");
    }

    for (int i = 0; i < branches->GetEntries(); ++i) {
        names.insert(branches->At(i)->GetName());
    }
    return names;
}

void RequireBranch(const std::set<std::string>& branches,
                   const std::string& name,
                   const std::string& context) {
    if (branches.find(name) == branches.end()) {
        throw std::runtime_error(
            "Missing input branch '" + name + "' required by " + context);
    }
}

void ValidateInputBranches(TChain& chain,
                           const std::vector<MomCorrParticle>& particles) {
    const auto branches = BranchNames(chain);
    for (const auto& particle : particles) {
        const std::string context =
            "particle role '" + particle.GetRole() + "'";
        RequireBranch(branches, particle.GetInputPxBranch(), context);
        RequireBranch(branches, particle.GetInputPyBranch(), context);
        RequireBranch(branches, particle.GetInputPzBranch(), context);
        RequireBranch(branches, particle.GetInputSectorBranch(), context);
    }
}

std::string SafeLabel(std::string label) {
    for (char& c : label) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_')) {
            c = '_';
        }
    }
    return label;
}

std::string AxisTitle(const std::string& variable) {
    if (variable == "p") return "Momentum (GeV)";
    if (variable == "missing_mass") return "MM (GeV/c^{2})";
    if (variable.find("deltaP") != std::string::npos) return "#Delta p (GeV)";
    return variable;
}

} // namespace

int main(int argc, char* argv[]) {
    ROOT::EnableImplicitMT();

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <directory/path/> <num_files | all> <config.json>\n";
        return 1;
    }

    try {
        const TString inputDirectory = argv[1];
        const auto rootFiles = FindRootFiles(inputDirectory);
        const int numFilesToProcess = ParseFileCount(argv[2], rootFiles.size());

        const std::string configFile = argv[3];
        MomCorrConfig config(configFile);
        std::vector<MomCorrParticle> particles = config.BuildParticles();

        TChain chain(config.GetTreeName().c_str());
        for (int i = 0; i < numFilesToProcess; ++i) {
            chain.Add(rootFiles[i]);
        }

        ValidateInputBranches(chain, particles);

        // Make compiled observable helpers visible to RDF string expressions.
        gInterpreter->AddIncludePath(MOMCORR_INCLUDE_DIR);
        gInterpreter->Declare("#include \"Kinematics.h\"\n#include \"DeltaP.h\"");

        TString outputDirectory = config.GetOutputDirectory();
        gSystem->mkdir(outputDirectory, true);

        TString outputFileName =
            outputDirectory + "/" + config.GetSpecifier() + "_raw_histograms.root";
        TFile outputFile(outputFileName, "RECREATE");
        if (outputFile.IsZombie()) {
            throw std::runtime_error(
                "Could not create output file: " + std::string(outputFileName.Data()));
        }

        TTree configTree("config", "Analysis configuration");
        std::ifstream in(configFile);
        std::ostringstream ss;
        ss << in.rdbuf();
        std::string configString = ss.str();
        configTree.Branch("json", &configString);
        configTree.Fill();
        configTree.Write();

        TTree meta("hist_meta", "Histogram metadata");
        int histId = 0;
        int sector = 0;
        int isPhiBinning = 0;
        int phiBin = 0;
        int xBins = 0;
        int yBins = 0;
        double xMin = 0.0;
        double xMax = 0.0;
        double yMin = 0.0;
        double yMax = 0.0;
        std::string histogramRequest;
        std::string particleRole;
        std::string particleName;
        std::string detector;
        std::string xVar;
        std::string yVar;
        std::string histName;

        meta.Branch("hist_id", &histId);
        meta.Branch("sector", &sector);
        meta.Branch("is_phi_binning", &isPhiBinning);
        meta.Branch("phi_bin", &phiBin);
        meta.Branch("histogram_request", &histogramRequest);
        meta.Branch("particle_role", &particleRole);
        meta.Branch("particle", &particleName);
        meta.Branch("detector", &detector);
        meta.Branch("x_var", &xVar);
        meta.Branch("y_var", &yVar);
        meta.Branch("x_bins", &xBins);
        meta.Branch("x_min", &xMin);
        meta.Branch("x_max", &xMax);
        meta.Branch("y_bins", &yBins);
        meta.Branch("y_min", &yMin);
        meta.Branch("y_max", &yMax);
        meta.Branch("hist_name", &histName);

        ROOT::RDataFrame dfBase(chain);
        ROOT::RDF::RNode df = dfBase;

        for (const auto& particle : particles) {
            df = particle.AddBranches(df);
        }
        df = momcorr::DefineObservables(
            df, config.GetObservables(), particles, config);

        for (const auto& hist : config.GetHistograms()) {
            const auto& particle =
                config.ParticleByRole(particles, hist.particleRole);

            std::vector<int> sectors = {0};
            if (hist.splitSector) {
                sectors = particle.GetSectors();
            }

            std::vector<std::pair<int, std::string>> phiBins = {{0, "all"}};
            if (hist.splitPhi && particle.IsPhiBinningEnabled()) {
                phiBins = particle.GetPhiBinningLabels();
            }

            for (const int currentSector : sectors) {
                ROOT::RDF::RNode dfSector = df;
                if (hist.splitSector) {
                    dfSector = df.Filter(
                        particle.GetSectorBranch() + " == " +
                        std::to_string(currentSector));
                }

                for (const auto& [currentPhiBin, phiLabel] : phiBins) {
                    ROOT::RDF::RNode dfBinned = dfSector;
                    if (hist.splitPhi && particle.IsPhiBinningEnabled()) {
                        dfBinned = dfSector.Filter(
                            particle.GetPhiBinBranch() + " == " +
                            std::to_string(currentPhiBin));
                    }

                    histId++;
                    const std::string sectorLabel =
                        hist.splitSector
                            ? "_sec" + std::to_string(currentSector)
                            : "_allsec";
                    const std::string phiPart =
                        (hist.splitPhi && particle.IsPhiBinningEnabled())
                            ? "_phi" + SafeLabel(phiLabel)
                            : "";

                    histName =
                        "h_" + hist.name + "_" + particle.GetName() +
                        sectorLabel + phiPart + "_id" +
                        std::to_string(histId);

                    const std::string title =
                        hist.name + " for " + particle.GetName() +
                        sectorLabel + phiPart + ";" +
                        AxisTitle(hist.x) + ";" + AxisTitle(hist.y);

                    auto h2 = dfBinned.Histo2D(
                        {histName.c_str(),
                         title.c_str(),
                         hist.xAxis.bins, hist.xAxis.low, hist.xAxis.high,
                         hist.yAxis.bins, hist.yAxis.low, hist.yAxis.high},
                        particle.GetMagBranch(),
                        hist.y
                    );

                    h2->Write();

                    sector = hist.splitSector ? currentSector : 0;
                    isPhiBinning =
                        (hist.splitPhi && particle.IsPhiBinningEnabled()) ? 1 : 0;
                    phiBin = isPhiBinning ? currentPhiBin : 0;
                    xBins = hist.xAxis.bins;
                    yBins = hist.yAxis.bins;
                    xMin = hist.xAxis.low;
                    xMax = hist.xAxis.high;
                    yMin = hist.yAxis.low;
                    yMax = hist.yAxis.high;
                    histogramRequest = hist.name;
                    particleRole = particle.GetRole();
                    particleName = particle.GetName();
                    detector = particle.GetDetector();
                    xVar = hist.x;
                    yVar = hist.y;
                    meta.Fill();
                }
            }
        }

        meta.Write();
        outputFile.Close();

        std::cout << "Wrote raw histograms to: " << outputFileName << "\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "make_initial_histograms: " << ex.what() << "\n";
        return 1;
    }
}
