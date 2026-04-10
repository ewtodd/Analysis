#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>

void MakePlot(TString signal_filename, TString background_filename) {

  TString signal_filepath = "root_files/" + signal_filename + ".root";
  TFile *signal_file = new TFile(signal_filepath, "UPDATE");
  TTree *signal_tree_with_pos =
      static_cast<TTree *>(signal_file->Get("bef_tree"));

  TString background_filepath = "root_files/" + background_filename + ".root";
  TFile *background_file = new TFile(background_filepath, "UPDATE");
  TTree *background_tree_with_pos =
      static_cast<TTree *>(background_file->Get("bef_tree"));

  Float_t signal_energy = 0;
  Float_t signal_x = 0, signal_y = 0, signal_z = 0;
  Int_t signal_nInteractions = 0;

  signal_tree_with_pos->SetBranchAddress("energykeV", &signal_energy);
  signal_tree_with_pos->SetBranchAddress("xmm", &signal_x);
  signal_tree_with_pos->SetBranchAddress("ymm", &signal_y);
  signal_tree_with_pos->SetBranchAddress("zmm", &signal_z);
  signal_tree_with_pos->SetBranchAddress("nInteractions",
                                         &signal_nInteractions);
  Float_t background_energy = 0;
  Float_t background_x = 0, background_y = 0, background_z = 0;
  Int_t background_nInteractions = 0;

  background_tree_with_pos->SetBranchAddress("energykeV", &background_energy);
  background_tree_with_pos->SetBranchAddress("xmm", &background_x);
  background_tree_with_pos->SetBranchAddress("ymm", &background_y);
  background_tree_with_pos->SetBranchAddress("zmm", &background_z);
  background_tree_with_pos->SetBranchAddress("nInteractions",
                                             &background_nInteractions);

  Int_t signal_n_entries = signal_tree_with_pos->GetEntries();
  Int_t background_n_entries = background_tree_with_pos->GetEntries();

  TH1D *signalSpectrum = new TH1D(
      PlottingUtils::GetRandomName(),
      Form("%s; Energy [keV]; Counts / %d eV", signal_filename.Data(),
           Constants::BIN_WIDTH_EV),
      Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN, Constants::ZOOMED_XMAX);

  TH1D *backgroundSpectrum = new TH1D(
      PlottingUtils::GetRandomName(),
      Form("%s; Energy [keV]; Counts / %d eV", signal_filename.Data(),
           Constants::BIN_WIDTH_EV),
      Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN, Constants::ZOOMED_XMAX);

  Bool_t signal_in_excluded_region;
  Bool_t background_in_excluded_region;

  for (Int_t i = 0; i < signal_n_entries; i++) {
    signal_in_excluded_region = kFALSE;

    signal_tree_with_pos->GetEntry(i);

    if (signal_nInteractions != 1)
      signal_in_excluded_region = kTRUE;

    if (signal_z < Constants::FILTER_DEPTH_MM)
      signal_in_excluded_region = kTRUE;

    if (signal_energy > Constants::ZOOMED_XMIN &&
        signal_energy < Constants::ZOOMED_XMAX) {
      if (!signal_in_excluded_region) {
        signalSpectrum->Fill(signal_energy);
      }
    }
  }

  for (Int_t i = 0; i < background_n_entries; i++) {
    background_in_excluded_region = kFALSE;

    background_tree_with_pos->GetEntry(i);

    if (background_nInteractions != 1)
      background_in_excluded_region = kTRUE;

    if (background_z < Constants::FILTER_DEPTH_MM)
      background_in_excluded_region = kTRUE;

    if (background_energy > Constants::ZOOMED_XMIN &&
        background_energy < Constants::ZOOMED_XMAX) {
      if (!background_in_excluded_region) {
        backgroundSpectrum->Fill(background_energy);
      }
    }
  }

  Float_t maxY = TMath::Max(signalSpectrum->GetMaximum(),
                            backgroundSpectrum->GetMaximum());

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();

  PlottingUtils::ConfigureHistogram(signalSpectrum, kRed);
  signalSpectrum->GetYaxis()->SetRangeUser(0, maxY * 1.1);
  signalSpectrum->GetYaxis()->SetTitleOffset(1.5);
  signalSpectrum->Draw("HIST");

  PlottingUtils::ConfigureHistogram(backgroundSpectrum, kBlack);
  backgroundSpectrum->Draw("HIST SAME");

  TLegend *leg = PlottingUtils::AddLegend(0.72, 0.9, 0.75, 0.88);
  leg->AddEntry(signalSpectrum, "Sample In", "l");
  leg->AddEntry(backgroundSpectrum, "Sample Out", "l");
  leg->Draw();

  canvas->SetLeftMargin(0.2);
  PlottingUtils::SaveFigure(canvas, "CZTSpecExample", "AlexPaper",
                            PlotSaveOptions::kLINEAR);
  signal_file->Close();
  background_file->Close();
}

void SignalBackground() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  MakePlot(Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
           Constants::CDSHIELDBACKGROUND_10PERCENT_20260113);
}
