#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>

void AddHistogram(TString filename) {
  TString filepath = "root_files/" + filename + ".root";
  TFile *file = new TFile(filepath, "UPDATE");

  Bool_t isFiltered = filename.Contains("_filtered");

  if (isFiltered && !Constants::FILTERED) {
    std::cerr << "File name contains filtered but constant is not set."
              << std::endl;
  }

  if ((isFiltered || Constants::FILTERED) && Constants::NORMALIZE_BY_TIME) {
    std::cout << "WARNING: normalization by time on filtered data may not be "
                 "reliable."
              << std::endl;
  }

  TTree *tree = nullptr;
  TString treeName;
  TString energyBranchName;

  if (isFiltered) {
    treeName = "bef_tree";
    energyBranchName = "energykeV";
    tree = static_cast<TTree *>(file->Get(treeName));
  } else {
    treeName = "bef_tree_event_summary";
    energyBranchName = "totalEnergykeV";
    tree = static_cast<TTree *>(file->Get(treeName));
  }

  if (!tree) {
    std::cerr << "ERROR: Could not find tree '" << treeName << "' in file "
              << filepath << std::endl;
    file->Close();
    return;
  }

  TParameter<Double_t> *param =
      (TParameter<Double_t> *)file->Get("N42_RealTime_Total");

  if (!param) {
    std::cerr << "WARNING: Could not find N42_RealTime_Total parameter in "
              << filepath << std::endl;
  }

  Double_t n42_time = param ? param->GetVal() : 1.0;

  Float_t energy;
  tree->SetBranchAddress(energyBranchName, &energy);

  Int_t n_entries = tree->GetEntries();

  TString perSecond = Constants::NORMALIZE_BY_TIME && param ? " / s" : "";

  TH1F *hist = new TH1F(
      PlottingUtils::GetRandomName(),
      Form("%s; Energy [keV]; Counts / %d eV%s", filename.Data(),
           Constants::BIN_WIDTH_EV, perSecond.Data()),
      Constants::HIST_NBINS, Constants::HIST_XMIN, Constants::HIST_XMAX);

  TH1F *zoomedHist = new TH1F(
      PlottingUtils::GetRandomName(),
      Form("%s; Energy [keV]; Counts / %d eV%s", filename.Data(),
           Constants::BIN_WIDTH_EV, perSecond.Data()),
      Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN, Constants::ZOOMED_XMAX);

  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    hist->Fill(energy);
    if (Constants::ZOOMED_XMIN < energy && energy < Constants::ZOOMED_XMAX)
      zoomedHist->Fill(energy);
  }

  if (Constants::NORMALIZE_BY_TIME && param) {
    Float_t normalization = 1.0 / n42_time;
    hist->Scale(normalization);
    zoomedHist->Scale(normalization);
  }

  std::cout << "Created histograms for " << filename
            << " (using tree: " << treeName << ", branch: " << energyBranchName
            << ")" << std::endl;

  PlottingUtils::ConfigureHistogram(hist, kP10Violet);

  PlottingUtils::ConfigureHistogram(zoomedHist, kP10Violet);

  hist->GetYaxis()->SetTitleOffset(1.2);
  zoomedHist->GetYaxis()->SetTitleOffset(1.2);

  hist->Write("hist", TObject::kOverwrite);
  zoomedHist->Write("zoomedHist", TObject::kOverwrite);

  std::cout << "Wrote histograms for " << filename << std::endl;

  file->Close();
}

void Histogram() {
  InitUtils::SetROOTPreferences();

  TString suffix = Constants::FILTERED ? "_filtered" : "";

  AddHistogram("01122026-PassiveBackground" + suffix);
  AddHistogram("01122026-Calibration" + suffix);

  AddHistogram("01132026-ActiveBackground-5Percent" + suffix);
  AddHistogram("01132026-ActiveBackground-25Percent" + suffix);
  AddHistogram("01132026-ActiveBackground-90Percent" + suffix);

  AddHistogram("01132026-CdShield-GeSamplesIn-10Percent" + suffix);
  AddHistogram("01132026-CdShield-GeSamplesIn-25Percent" + suffix);

  AddHistogram("01132026-CdShield-ActiveBackground-10Percent" + suffix);

  AddHistogram("01132026-CuShield-GeSamplesIn-10Percent" + suffix);

  AddHistogram("01132026-CuShield-ActiveBackground-Am241-10Percent" + suffix);

  AddHistogram("01132026-PostReactor-Calibration" + suffix);

  AddHistogram("01142026-CuShield-GeSamplesIn-10Percent" + suffix);

  AddHistogram("01142026-CuShield-ActiveBackground-10Percent" + suffix);

  AddHistogram("01142026-CuShield-GeSamplesIn-MovedBack-90Percent" + suffix);

  AddHistogram("01152026-NewSetup-GeSamplesIn-5Percent" + suffix);
  AddHistogram("01152026-NewSetup-ActiveBackground-5Percent" + suffix);
  AddHistogram("01152026-NewSetup-PostReactor-Am241" + suffix);
  AddHistogram("01152026-NewSetup-PostReactor-Ba133" + suffix);
  AddHistogram("01152026-NewSetup-ShutterClosed-5Percent" + suffix);

  AddHistogram("01162026-NoShield-GeOnCZT-0_5Percent" + suffix);
  AddHistogram("01162026-NoShield-ActiveBackground-0_5Percent" + suffix);
  AddHistogram("01162026-NoShield-GeSamplesIn-GraphiteCastle-10Percent" +
               suffix);
  AddHistogram("01162026-NoShield-ActiveBackground-GraphiteCastle-10Percent" +
               suffix);
  AddHistogram("01162026-NoShield-PostReactor-Am241-Ba133" + suffix);
}
