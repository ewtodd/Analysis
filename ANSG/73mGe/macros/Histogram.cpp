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

  TH1F *peakHist = new TH1F(
      PlottingUtils::GetRandomName(),
      Form("%s; Energy [keV]; Counts / %d eV%s", filename.Data(),
           Constants::BIN_WIDTH_EV, perSecond.Data()),
      Constants::PEAK_NBINS, Constants::PEAK_XMIN, Constants::PEAK_XMAX);

  tree->LoadBaskets();
  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    hist->Fill(energy);
    if (Constants::ZOOMED_XMIN < energy && energy < Constants::ZOOMED_XMAX)
      zoomedHist->Fill(energy);
    if (Constants::PEAK_XMIN < energy && energy < Constants::PEAK_XMAX)
      peakHist->Fill(energy);
  }

  if (Constants::NORMALIZE_BY_TIME && param) {
    Float_t normalization = 1.0 / n42_time;
    hist->Scale(normalization);
    zoomedHist->Scale(normalization);
    peakHist->Scale(normalization);
  }

  std::cout << "Created histograms for " << filename
            << " (using tree: " << treeName << ", branch: " << energyBranchName
            << ")" << std::endl;

  PlottingUtils::ConfigureHistogram(hist, kP10Violet);
  PlottingUtils::ConfigureHistogram(zoomedHist, kP10Violet);
  PlottingUtils::ConfigureHistogram(peakHist, kP10Violet);

  hist->GetYaxis()->SetTitleOffset(1.2);
  zoomedHist->GetYaxis()->SetTitleOffset(1.2);
  peakHist->GetYaxis()->SetTitleOffset(1.2);

  hist->Write("hist", TObject::kOverwrite);
  zoomedHist->Write("zoomedHist", TObject::kOverwrite);
  peakHist->Write("peakHist", TObject::kOverwrite);

  std::cout << "Wrote histograms for " << filename << std::endl;

  file->Close();
}

void AddAllHistograms(std::vector<TString> filenames) {
  Int_t n_files = filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString filename = filenames.at(i);
    AddHistogram(filename);
  }
}

void Histogram() {
  InitUtils::SetROOTPreferences();
  std::vector<TString> filenames;

  // January 12, 2026
  filenames.push_back(Constants::PASSIVEBACKGROUND_01122026);
  filenames.push_back(Constants::CALIBRATION_01122026);

  //  // January 13, 2026
  filenames.push_back(Constants::ACTIVEBACKGROUND_TEST_5PERCENT_01132026);
  filenames.push_back(Constants::ACTIVEBACKGROUND_TEST_90PERCENT_01132026);
  filenames.push_back(Constants::CDSHIELDSIGNAL_10PERCENT_01132026);
  filenames.push_back(Constants::CDSHIELDBACKGROUND_10PERCENT_01132026);
  filenames.push_back(Constants::CDSHIELDSIGNAL_25PERCENT_01132026);
  filenames.push_back(Constants::CDSHIELDBACKGROUND_25PERCENT_01132026);
  filenames.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_01132026);
  filenames.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_01132026);
  filenames.push_back(Constants::POSTREACTOR_AM241_01132026);

  // January 14, 2026
  filenames.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_01142026);
  filenames.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_01142026);
  filenames.push_back(Constants::CUSHIELDSIGNAL_90PERCENT_01142026);

  // January 15, 2026
  filenames.push_back(Constants::NOSHIELDSIGNAL_5PERCENT_01152026);
  filenames.push_back(Constants::NOSHIELDBACKGROUND_5PERCENT_01152026);
  filenames.push_back(Constants::POSTREACTOR_AM241_01152026);
  filenames.push_back(Constants::POSTREACTOR_BA133_01152026);
  filenames.push_back(Constants::SHUTTERCLOSED_01152026);

  // January 16, 2026
  filenames.push_back(Constants::NOSHIELD_GEONCZT_0_5PERCENT_01162026);
  filenames.push_back(Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_01162026);
  filenames.push_back(
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_01162026);
  filenames.push_back(
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_01162026);
  filenames.push_back(Constants::POSTREACTOR_AM241_BA133_01162026);

  AddAllHistograms(filenames);
}
