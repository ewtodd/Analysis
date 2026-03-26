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

  treeName = "bef_tree";
  energyBranchName = "energykeV";
  tree = static_cast<TTree *>(file->Get(treeName));

  if (!tree) {
    std::cerr << "ERROR: Could not find tree '" << treeName << "' in file "
              << filepath << std::endl;
    file->Close();
    return;
  }

  Float_t energy;
  tree->SetBranchAddress(energyBranchName, &energy);

  Float_t x = 0, y = 0, z = 0;

  Int_t nInteractions = 0;

  tree->SetBranchAddress("xum", &x);
  tree->SetBranchAddress("yum", &y);
  tree->SetBranchAddress("zum", &z);
  tree->SetBranchAddress("nInteractions", &nInteractions);

  Int_t n_entries = tree->GetEntries();

  TH1F *hist = new TH1F(PlottingUtils::GetRandomName(),
                        Form("%s; Energy [keV]; Counts / %d eV",
                             filename.Data(), Constants::BIN_WIDTH_EV),
                        Constants::HIST_NBINS, Constants::HIST_XMIN,
                        Constants::HIST_XMAX);

  TH1F *zoomedHist = new TH1F(PlottingUtils::GetRandomName(),
                              Form("%s; Energy [keV]; Counts / %d eV",
                                   filename.Data(), Constants::BIN_WIDTH_EV),
                              Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                              Constants::ZOOMED_XMAX);

  TH1F *peakHist = new TH1F(PlottingUtils::GetRandomName(),
                            Form("%s; Energy [keV]; Counts / %d eV",
                                 filename.Data(), Constants::BIN_WIDTH_EV),
                            Constants::PEAK_NBINS, Constants::PEAK_XMIN,
                            Constants::PEAK_XMAX);

  tree->LoadBaskets();
  Bool_t in_excluded_region;

  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);

    in_excluded_region = kFALSE;
    if (nInteractions != 1)
      in_excluded_region = kTRUE;
    if (z < Constants::FILTER_DEPTH_UM)
      in_excluded_region = kTRUE;

    if (!in_excluded_region) {
      for (size_t r = 0; r < Constants::FILTER_REGIONS_EXCLUDE_XY_UM.size();
           r++) {
        if (x >= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmin &&
            x <= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmax &&
            y >= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymin &&
            y <= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymax) {
          in_excluded_region = kTRUE;
          break;
        }
      }
    }

    if (!in_excluded_region) {
      hist->Fill(energy);
      if (Constants::ZOOMED_XMIN < energy && energy < Constants::ZOOMED_XMAX)
        zoomedHist->Fill(energy);
      if (Constants::PEAK_XMIN < energy && energy < Constants::PEAK_XMAX)
        peakHist->Fill(energy);
    }
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
  AddAllHistograms(Constants::ALL_DATASETS);
}
