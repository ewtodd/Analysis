#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>
#include <future>
#include <mutex>
#include <thread>

static std::mutex print_mutex;

Int_t GetCrystalIndex(Float_t x, Float_t y) {
  if (x < 0 && y < 0)
    return 0;
  if (x > 0 && y < 0)
    return 1;
  if (x < 0 && y > 0)
    return 2;
  if (x > 0 && y > 0)
    return 3;
  return -1; // on boundary
}

Bool_t FilterFile(TString filename) {
  TString filepath = "root_files/" + filename + ".root";
  TFile *file = new TFile(filepath, "UPDATE");
  TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));

  Float_t energy = 0;
  Float_t x = 0, y = 0, z = 0;
  UInt_t eventTime = 0;
  Int_t liveTime = 0;
  Int_t nInteractions = 0;
  Int_t interaction = 0;

  tree->SetBranchAddress("energykeV", &energy);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);
  tree->SetBranchAddress("zmm", &z);
  tree->SetBranchAddress("eventTime", &eventTime);
  tree->SetBranchAddress("liveTime", &liveTime);
  tree->SetBranchAddress("nInteractions", &nInteractions);
  tree->SetBranchAddress("interaction", &interaction);

  Int_t n_entries = tree->GetEntries();

  Float_t outEnergy = 0;
  Float_t outX = 0, outY = 0, outZ = 0;
  UInt_t outEventTime = 0;
  Int_t outLiveTime = 0;
  Int_t outNInteractions = 0;
  Int_t outInteraction = 0;

  Float_t crystalLiveTime_s[Constants::N_CRYSTALS];
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++)
    crystalLiveTime_s[c] = 0.0;

  TTree *filteredTrees[Constants::N_CRYSTALS];
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    filteredTrees[c] = new TTree(Form("crystal%d_filtered_tree", c),
                                 Form("All interactions for crystal %d", c));
    filteredTrees[c]->Branch("energykeV", &outEnergy, "energykeV/F");
    filteredTrees[c]->Branch("xmm", &outX, "xmm/F");
    filteredTrees[c]->Branch("ymm", &outY, "ymm/F");
    filteredTrees[c]->Branch("zmm", &outZ, "zmm/F");
    filteredTrees[c]->Branch("eventTime", &outEventTime, "eventTime/i");
    filteredTrees[c]->Branch("liveTime", &outLiveTime, "liveTime/I");
    filteredTrees[c]->Branch("nInteractions", &outNInteractions,
                             "nInteractions/I");
    filteredTrees[c]->Branch("interaction", &outInteraction, "interaction/I");
  }

  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);

    Int_t crystal = GetCrystalIndex(x, y);
    if (crystal < 0)
      continue;

    outEnergy = energy;
    outX = x;
    outY = y;
    outZ = z;
    outEventTime = eventTime;
    outLiveTime = liveTime;
    outNInteractions = nInteractions;
    outInteraction = interaction;

    filteredTrees[crystal]->Fill();

    // Accumulate live time only for the first interaction of each event
    if (interaction == 0)
      crystalLiveTime_s[crystal] += liveTime * Constants::TENS_OF_NS_TO_S;
  }

  {
    std::lock_guard<std::mutex> lock(print_mutex);
    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      std::cout << filename << " crystal " << c << ": "
                << filteredTrees[c]->GetEntries() << " interactions"
                << std::endl;
    }
  }

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    filteredTrees[c]->Write("", TObject::kOverwrite);
  }

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TParameter<Float_t> *ltParam = new TParameter<Float_t>(
        Form("LiveTime_Filtered_Crystal%d_s", c), crystalLiveTime_s[c]);
    ltParam->Write("", TObject::kOverwrite);
  }

  file->Close();
  return kTRUE;
}

void Filter() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);
  ROOT::EnableThreadSafety();

  std::vector<TString> filenames = Constants::ALL_DATASETS;
  Int_t n_files = Int_t(filenames.size());
  Int_t n_workers = Int_t(std::thread::hardware_concurrency());
  n_workers = TMath::Min(n_workers, n_files);

  std::cout << "Filtering " << n_files << " files with " << n_workers
            << " workers." << std::endl;

  for (Int_t i = 0; i < n_files; i += n_workers) {
    std::vector<std::future<Bool_t>> futures;
    Int_t batch_end = TMath::Min(i + n_workers, n_files);

    for (Int_t j = i; j < batch_end; ++j) {
      futures.push_back(
          std::async(std::launch::async, FilterFile, filenames[j]));
    }

    for (size_t j = 0; j < futures.size(); ++j) {
      Bool_t result = futures[j].get();
      if (!result) {
        std::cerr << "FAILED: " << filenames[i + j] << std::endl;
      }
    }
  }

  std::cout << "All filtering complete." << std::endl;
}
