#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>
#include <iomanip>

struct DatasetPair {
  TString name;
  TString signalFile;
  TString backgroundFile;
};

Int_t GetCrystalIndex(Double_t x, Double_t y) {
  if (x < 0 && y < 0)
    return 0;
  if (x > 0 && y < 0)
    return 1;
  if (x < 0 && y > 0)
    return 2;
  if (x > 0 && y > 0)
    return 3;
  return -1;
}

Double_t ComputeLiveTimePerCrystal(TFile *file, Int_t crystal) {
  TString paramName = Form("LiveTime_Unfiltered_Crystal%d_s", crystal);
  TParameter<Double_t> *param =
      (TParameter<Double_t> *)file->Get(paramName);
  if (!param) {
    std::cerr << "WARNING: " << paramName << " not found" << std::endl;
    return 0.0;
  }
  return param->GetVal();
}

struct CrystalHistograms {
  TH1D *wide;
};

struct MultiIntHistograms {
  TH1D *wide;
  std::vector<CrystalHistograms> crystalHists;
};

MultiIntHistograms BuildMultiIntHistograms(TString filename) {
  MultiIntHistograms result;
  TString filepath = "root_files/" + filename + ".root";
  TFile *file = new TFile(filepath, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << filepath << std::endl;
    return {};
  }

  TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));
  if (!tree) {
    std::cerr << "ERROR: No bef_tree in " << filepath << std::endl;
    file->Close();
    delete file;
    return {};
  }

  Double_t totalEnergy = 0;
  Double_t x = 0, y = 0;
  Int_t nInteractions = 0;
  Int_t interaction = 0;

  tree->SetBranchAddress("totalEnergykeV", &totalEnergy);
  tree->SetBranchAddress("xmm", &x);
  tree->SetBranchAddress("ymm", &y);
  tree->SetBranchAddress("nInteractions", &nInteractions);
  tree->SetBranchAddress("interaction", &interaction);

  // Per-crystal histograms
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    CrystalHistograms ch;
    ch.wide = new TH1D(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        Constants::HIST_NBINS, Constants::HIST_XMIN, Constants::HIST_XMAX);
    ch.wide->SetDirectory(0);
    result.crystalHists.push_back(ch);
  }

  Int_t n_entries = tree->GetEntries();
  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);

    if (nInteractions == 1)
      continue;
    if (interaction != 0)
      continue;

    Int_t crystal = GetCrystalIndex(x, y);
    if (crystal >= 0 && crystal < Constants::N_CRYSTALS)
      result.crystalHists[crystal].wide->Fill(totalEnergy);
  }

  // Get per-crystal live times and normalize
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    Double_t lt = ComputeLiveTimePerCrystal(file, c);
    if (lt > 0) {
      result.crystalHists[c].wide->Scale(1.0 / lt);
    }
  }

  file->Close();
  delete file;

  // Sum crystals
  result.wide = static_cast<TH1D *>(
      result.crystalHists[0].wide->Clone(PlottingUtils::GetRandomName()));
  result.wide->SetDirectory(0);

  for (Int_t c = 1; c < Constants::N_CRYSTALS; c++) {
    result.wide->Add(result.crystalHists[c].wide);
  }

  std::cout << filename << ": built multi-interaction histograms" << std::endl;
  return result;
}

void SubtractBackgroundMultiInt() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  std::vector<DatasetPair> datasets;

  datasets.push_back({"CdShield_10_01132026",
                      Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
                      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113});

  datasets.push_back({"CdShield_25_01132026",
                      Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
                      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113});

  datasets.push_back({"CuShield_10_01132026",
                      Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113});

  datasets.push_back({"CuShield_10_01142026",
                      Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114});

  datasets.push_back({"NoShield_5_01152026",
                      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115});

  datasets.push_back(
      {"NoShield_GraphiteCastle_10_01162026",
       Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
       Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116});

  TH1D *wideCombined = new TH1D(
      "multiint_wide_combined",
      Form("Multi-Int BkgSub; Energy [keV]; Counts / %d eV / s",
           Constants::BIN_WIDTH_EV),
      Constants::HIST_NBINS, Constants::HIST_XMIN, Constants::HIST_XMAX);
  wideCombined->SetDirectory(0);

  for (Int_t i = 0; i < (Int_t)datasets.size(); i++) {
    std::cout << "Processing: " << datasets[i].name << std::endl;

    MultiIntHistograms signal = BuildMultiIntHistograms(datasets[i].signalFile);
    MultiIntHistograms background =
        BuildMultiIntHistograms(datasets[i].backgroundFile);

    // Subtract: signal - background (already normalized by live time)
    signal.wide->Add(background.wide, -1);

    wideCombined->Add(signal.wide);
  }

  TCanvas *canvasWide = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(wideCombined, kP10Violet);
  wideCombined->SetFillStyle(0);
  wideCombined->SetLineWidth(2);
  wideCombined->Draw("HIST");
  PlottingUtils::SaveFigure(canvasWide, "multiint_background_subtracted",
                            "backgroundSubtraction", PlotSaveOptions::kLOG);

  TFile *outputFile =
      new TFile("root_files/Combined_MultiInt_BkgSubtracted.root", "RECREATE");
  wideCombined->Write("wideCombined", TObject::kOverwrite);
  outputFile->Close();
  delete outputFile;

  std::cout << "Multi-interaction background subtraction complete!"
            << std::endl;
}
