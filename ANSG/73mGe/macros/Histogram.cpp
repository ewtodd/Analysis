#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>

void AddHistogram(TString filename, Bool_t reprocess = kTRUE) {
  if (!reprocess)
    return;
  TString filepath = "root_files/" + filename + ".root";
  TFile *file = new TFile(filepath, "UPDATE");

  Bool_t isFiltered = filename.Contains("_filtered");

  if (isFiltered && !Constants::FILTERED) {
    std::cerr << "File name contains filtered but constant is not set."
              << std::endl;
  }

  if ((isFiltered || Constants::FILTERED) && Constants::NORMALIZE_BY_TIME) {
    std::cerr << "Cannot normalize by time on filtered data." << std::endl;
  }

  TTree *tree = nullptr;
  TString treeName;
  TString energyBranchName;

  if (isFiltered) {
    treeName = "bef_tree";
    energyBranchName = "energy";
    tree = static_cast<TTree *>(file->Get(treeName));
  } else {
    treeName = "bef_tree_event_summary";
    energyBranchName = "totalEnergy";
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

  UInt_t energy_uint = 0;
  tree->SetBranchAddress(energyBranchName, &energy_uint);

  Int_t n_entries = tree->GetEntries();

  double binWidth_keV = Constants::BIN_WIDTH_EV / 1000.0;

  int hist_xmin = 0, hist_xmax = 1500;
  int zoom_xmin = 50, zoom_xmax = 90;
  int hist_nbins = (hist_xmax - hist_xmin) / binWidth_keV;
  int zoomed_nbins = (zoom_xmax - zoom_xmin) / binWidth_keV;

  TString perSecond = Constants::NORMALIZE_BY_TIME && param ? " / s" : "";

  TH1F *hist =
      new TH1F(PlottingUtils::GetRandomName(),
               Form("%s; Energy [keV]; Counts / %d eV%s", filename.Data(),
                    Constants::BIN_WIDTH_EV, perSecond.Data()),
               hist_nbins, hist_xmin, hist_xmax);

  TH1F *zoomedHist =
      new TH1F(PlottingUtils::GetRandomName(),
               Form("%s; Energy [keV]; Counts / %d eV%s", filename.Data(),
                    Constants::BIN_WIDTH_EV, perSecond.Data()),
               zoomed_nbins, zoom_xmin, zoom_xmax);

  Float_t energy = 0;

  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    energy = energy_uint / 1000.0;
    hist->Fill(energy);
    if (zoom_xmin < energy && energy < zoom_xmax)
      zoomedHist->Fill(energy);
  }

  if (Constants::NORMALIZE_BY_TIME && param) {
    hist->Scale(1.0 / n42_time);
    zoomedHist->Scale(1.0 / n42_time);
  }

  std::cout << "Created histograms for " << filename
            << " (using tree: " << treeName << ", branch: " << energyBranchName
            << ")" << std::endl;

  PlottingUtils::ConfigureHistogram(hist, kP10Violet);

  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas, kFALSE);
  PlottingUtils::ConfigureAndDrawHistogram(zoomedHist, kP10Violet);
  PlottingUtils::SaveFigure(canvas, filename + "_zoomedHist.png", kFALSE);

  hist->Write("hist", TObject::kOverwrite);
  zoomedHist->Write("zoomedHist", TObject::kOverwrite);

  std::cout << "Wrote histograms for " << filename << std::endl;

  file->Close();
}

void Histogram() {
  Bool_t reprocess = kTRUE;

  InitUtils::SetROOTPreferences();

  TString suffix = Constants::FILTERED ? "_filtered" : "";

  AddHistogram("01122026-Calibration" + suffix, reprocess);

  AddHistogram("01132026-CdShield-GeSamplesIn-25Percent" + suffix, reprocess);

  AddHistogram("01132026-CdShield-GeSamplesIn-10Percent" + suffix, reprocess);

  AddHistogram("01132026-CdShield-ActiveBackground-10Percent" + suffix,
               reprocess);

  AddHistogram("01132026-CuShield-GeSamplesIn-10Percent" + suffix, reprocess);

  AddHistogram("01132026-CuShield-ActiveBackground-Am241-10Percent" + suffix,
               reprocess);

  AddHistogram("01132026-PostReactor-Calibration" + suffix, reprocess);

  AddHistogram("01142026-CuShield-GeSamplesIn-10Percent" + suffix, reprocess);

  AddHistogram("01142026-CuShield-ActiveBackground-10Percent" + suffix,
               reprocess);

  AddHistogram("01142026-CuShield-GeSamplesIn-MovedBack-90Percent" + suffix,
               reprocess);
}
