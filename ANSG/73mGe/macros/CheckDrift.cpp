#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TParameter.h>
#include <TTree.h>
#include <iostream>
#include <vector>

const Float_t E_PB_KA1 = 72.8042;
const Float_t PB_FIT_LOW = 66;
const Float_t PB_FIT_HIGH = 81;
const Int_t N_POINTS = 100000;

void ProcessDrift(TString filename) {
  TString filepath = "root_files/" + filename + ".root";
  TFile *file = new TFile(filepath, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << filepath << std::endl;
    return;
  }

  TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));
  if (!tree) {
    std::cerr << "ERROR: Cannot find bef_tree in " << filepath << std::endl;
    file->Close();
    delete file;
    return;
  }

  UInt_t eventTime;
  tree->SetBranchAddress("eventTime", &eventTime);
  Int_t liveTime;
  tree->SetBranchAddress("liveTime", &liveTime);

  Int_t n_entries = tree->GetEntries();

  TGraph *eventPlot = new TGraph(n_entries);
  Float_t event;
  Float_t t0;

  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    event = eventTime / 1e3;

    if (i == 0)
      t0 = event;

    eventPlot->SetPoint(i, i, event - t0);
  }

  TCanvas *canvasEvent = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawGraph(eventPlot, kBlue, "; Sample; Time [s]");
  PlottingUtils::SaveFigure(canvasEvent, filename + "_Event",
                            PlotSaveOptions::kLINEAR);

  file->Close();

  delete file;
  delete canvasEvent;

  std::cout << "Processed drift check for " << filename << std::endl;
}

void ProcessPeakDrift(TString filename) {
  TString filepath = "root_files/" + filename + ".root";
  TFile *file = new TFile(filepath, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << filepath << std::endl;
    return;
  }

  TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));
  if (!tree) {
    std::cerr << "ERROR: Cannot find bef_tree in " << filepath << std::endl;
    file->Close();
    delete file;
    return;
  }

  Float_t energy;
  tree->SetBranchAddress("energykeV", &energy);

  Int_t n_entries = tree->GetEntries();

  // Fill one histogram with all events in the Pb Ka1 window
  TH1F *peakHist = new TH1F(PlottingUtils::GetRandomName(), "", 100, PB_FIT_LOW,
                            PB_FIT_HIGH);
  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    if (energy >= PB_FIT_LOW && energy <= PB_FIT_HIGH)
      peakHist->Fill(energy);
  }
  peakHist->SetDirectory(0);

  // Fit the whole dataset once
  FittingUtils *fitter = new FittingUtils(peakHist, PB_FIT_LOW, PB_FIT_HIGH,
                                          kTRUE, kTRUE, kTRUE, kTRUE, kTRUE);
  FitResult result = fitter->FitPeak(filename, "PbKa1");

  if (!result.valid) {
    std::cerr << "ERROR: Fit failed for " << filename << std::endl;
    delete fitter;
    delete peakHist;
    file->Close();
    delete file;
    return;
  }

  TF1 *fitFunc = fitter->GetFitFunction();

  // Sample N_POINTS from the fitted function and plot value vs sample index
  TGraph *samplePlot = new TGraph(N_POINTS);
  for (Int_t i = 0; i < N_POINTS; i++) {
    Float_t val = fitFunc->GetRandom(PB_FIT_LOW, PB_FIT_HIGH);
    samplePlot->SetPoint(i, i, val);
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureAndDrawGraph(samplePlot, kRed,
                                       "; Sample; Pb K#alpha_{1} Energy [keV]");
  PlottingUtils::SaveFigure(canvas, filename + "_PbKa1Drift",
                            PlotSaveOptions::kLINEAR);

  file->Close();

  delete fitter;
  delete peakHist;
  delete file;
  delete canvas;

  std::cout << "Processed Pb Ka1 peak drift for " << filename << std::endl;
}

void CheckDrift() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  std::vector<TString> filenames;

  filenames.push_back(Constants::PASSIVEBACKGROUND_20260112);
  filenames.push_back(Constants::CALIBRATION_20260112);
  filenames.push_back(Constants::ACTIVEBACKGROUND_TEST_5PERCENT_20260113);
  filenames.push_back(Constants::ACTIVEBACKGROUND_TEST_90PERCENT_20260113);
  filenames.push_back(Constants::CDSHIELDSIGNAL_10PERCENT_20260113);
  filenames.push_back(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113);
  filenames.push_back(Constants::CDSHIELDSIGNAL_25PERCENT_20260113);
  filenames.push_back(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113);
  filenames.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_20260113);
  filenames.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113);
  filenames.push_back(Constants::POSTREACTOR_AM241_20260113);
  filenames.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_20260114);
  filenames.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114);
  filenames.push_back(Constants::CUSHIELDSIGNAL_90PERCENT_20260114);
  filenames.push_back(Constants::NOSHIELDSIGNAL_5PERCENT_20260115);
  filenames.push_back(Constants::NOSHIELDBACKGROUND_5PERCENT_20260115);
  filenames.push_back(Constants::POSTREACTOR_AM241_20260115);
  filenames.push_back(Constants::POSTREACTOR_BA133_20260115);
  filenames.push_back(Constants::SHUTTERCLOSED_20260115);
  filenames.push_back(Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116);
  filenames.push_back(Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116);
  filenames.push_back(
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116);
  filenames.push_back(
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116);
  filenames.push_back(Constants::POSTREACTOR_AM241_BA133_20260116);

  Int_t n_files = filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    ProcessDrift(filenames.at(i));
  }

  std::vector<TString> pb_backgrounds;
  pb_backgrounds.push_back(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113);
  pb_backgrounds.push_back(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113);
  pb_backgrounds.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113);

  Int_t n_pb = pb_backgrounds.size();
  for (Int_t i = 0; i < n_pb; i++) {
    ProcessPeakDrift(pb_backgrounds.at(i));
  }
}
