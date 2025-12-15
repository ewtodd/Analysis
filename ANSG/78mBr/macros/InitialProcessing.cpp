#include "PlottingUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

void InitialWaveformProcessing(const TString filepath,
                               const TString output_name, const Int_t color,
                               const Bool_t reprocess = kFALSE) {
  if (!(gSystem->AccessPathName(output_name + ".root")) &&
      reprocess == kFALSE) {
    std::cout << "Output file already exists and reprocess is not set to true, "
                 "skipping..."
              << std::endl;
    return;
  }
  PlottingUtils::SetROOTPreferences();
  WaveformProcessingUtils *processor = new WaveformProcessingUtils();
  processor->SetPolarity(-1);
  processor->SetTriggerThreshold(0.15);
  processor->SetSampleWindows(15, 125);
  processor->SetGates(5, 10, 200);
  processor->SetMaxEvents(-1);
  processor->SetVerbose(kTRUE);
  processor->ProcessFile(filepath, output_name);

  std::cout << "MACRO: Processed raw waveforms." << std::endl;
  std::cout << "MACRO: Now reading results file." << std::endl;

  TFile *output = new TFile(output_name + ".root", "UPDATE");
  TTree *features_tree = static_cast<TTree *>(output->Get("features"));

  Float_t long_integral;
  features_tree->SetBranchAddress("long_integral", &long_integral);
  Float_t pulse_height;
  features_tree->SetBranchAddress("pulse_height", &pulse_height);

  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas, kFALSE);

  TH1F *long_integral_hist =
      new TH1F("", "; Pulse Integral [a.u.]; Counts", 500, 0, 4e5);

  TH1F *pulse_height_hist =
      new TH1F("", "; Pulse Height [ADC]; Counts", 500, 0, 16384);

  Int_t num_entries = features_tree->GetEntries();

  for (Int_t i = 0; i < num_entries; i++) {
    features_tree->GetEntry(i);
    long_integral_hist->Fill(long_integral);
    pulse_height_hist->Fill(pulse_height);
  }

  PlottingUtils::ConfigureAndDrawHistogram(long_integral_hist, color);
  PlottingUtils::SaveFigure(canvas, output_name + "_long_integral.png");

  PlottingUtils::ConfigureAndDrawHistogram(pulse_height_hist, color);
  canvas->Update();

  PlottingUtils::SaveFigure(canvas, output_name + "_pulse_height.png");

  output->cd();
  long_integral_hist->Write("Pulse Integral", TObject::kOverwrite);
  pulse_height_hist->Write("Pulse Height", TObject::kOverwrite);
  output->Close();

  delete canvas;
  delete processor;
}
void InitialProcessing() {
  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();

  TString filepath_Am241 =
      "/home/e-work/LABDATA/ANSG/78mBr/half_life_2/DAQ/"
      "59_5keV_calibration_300s/RAW/"
      "DataR_CH0@DT5730S_31017_59_5keV_calibration_300s.root";
  TString output_name_Am241 = "calibration_Am241";
  InitialWaveformProcessing(filepath_Am241, output_name_Am241, colors[0],
                            kTRUE);
  TString filepath_Eu152 =
      "/home/e-work/LABDATA/ANSG/78mBr/half_life_2/DAQ/"
      "Europium_calibration_300s/RAW/"
      "DataR_CH0@DT5730S_31017_Europium_calibration_300s.root";
  TString output_name_Eu152 = "calibration_Eu152";
  InitialWaveformProcessing(filepath_Eu152, output_name_Eu152, colors[1],
                            kTRUE);
  TString filepath_bkg = "/home/e-work/LABDATA/ANSG/78mBr/day2/"
                         "bkg_day2/RAW/DataR_CH0@DT5730B_969_bkg_day2.root";
  TString output_name_bkg = "background";
  InitialWaveformProcessing(filepath_bkg, output_name_bkg, colors[2], kTRUE);
}
