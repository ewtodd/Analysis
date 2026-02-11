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

  TTree *tree = nullptr;
  TString treeName;
  TString energyBranchName;

  treeName = "CZT";
  energyBranchName = "fEDep";
  tree = static_cast<TTree *>(file->Get(treeName));

  if (!tree) {
    std::cerr << "ERROR: Could not find tree '" << treeName << "' in file "
              << filepath << std::endl;
    file->Close();
    return;
  }

  Double_t energy;
  tree->SetBranchAddress(energyBranchName, &energy);

  Int_t n_entries = tree->GetEntries();
  TString perSecond = "";

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

void SimHistogram() {
  InitUtils::SetROOTPreferences();
  std::vector<TString> filenames;

  filenames.push_back(Constants::SIM_1E6);
  filenames.push_back(Constants::SIM_5E7);

  AddAllHistograms(filenames);
}
