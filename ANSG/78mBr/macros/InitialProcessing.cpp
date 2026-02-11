#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

void InitialWaveformProcessing(const TString filepath,
                               const TString output_name, const Int_t color) {

  WaveformProcessingUtils *processor = new WaveformProcessingUtils();
  processor->SetPolarity(-1);
  processor->SetTriggerThreshold(0.15);
  processor->SetNumberOfSamplesForBaseline(10);
  processor->SetSampleWindows(18, 135);
  processor->SetGates(5, 40, 220);
  processor->SetSaveSampleWaveforms(5);
  processor->SetMaxEvents(-1);
  processor->SetVerbose(kTRUE);
  processor->ProcessFile(filepath, output_name);
  std::cout << "MACRO: Processed raw waveforms." << std::endl;
  std::cout << "MACRO: Now reading results file." << std::endl;

  TString output_filepath = "root_files/" + output_name + ".root";

  TFile *output = new TFile(output_filepath, "UPDATE");
  TTree *features_tree = static_cast<TTree *>(output->Get("features"));

  Float_t long_integral;
  features_tree->SetBranchAddress("long_integral", &long_integral);
  Float_t pulse_height;
  features_tree->SetBranchAddress("pulse_height", &pulse_height);

  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas, kFALSE);

  TH1F *long_integral_hist =
      new TH1F("",
               Form("; Pulse Integral [a.u.]; Counts / %d a.u.",
                    Constants::PI_BIN_WIDTH),
               Constants::PI_HIST_NBINS, Constants::PI_HIST_XMIN,
               Constants::PI_HIST_XMAX);

  TH1F *pulse_height_hist = new TH1F(
      "",
      Form("; Pulse Height [ADC]; Counts / %d ADC", Constants::PH_BIN_WIDTH),
      Constants::PH_HIST_NBINS, Constants::PH_HIST_XMIN,
      Constants::PH_HIST_XMAX);

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

void ProcessFiles(
    std::vector<TString> filepaths, std::vector<TString> output_names,
    std::vector<Int_t> colors = PlottingUtils::GetDefaultColors()) {
  Int_t entries = filepaths.size();
  for (Int_t i = 0; i < entries; i++) {
    TString filepath = filepaths[i];
    TString output_name = output_names[i];
    Int_t color = colors[i];
    InitialWaveformProcessing(filepath, output_name, color);
  }
}

void InitialProcessing() {
  InitUtils::SetROOTPreferences();

  std::vector<TString> filepaths = {
      "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/"
      "59_5keV_calibration_300s/RAW/"
      "DataR_CH0@DT5730S_31017_59_5keV_calibration_300s.root",
      "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/"
      "Europium_calibration_300s/RAW/"
      "DataR_CH0@DT5730S_31017_Europium_calibration_300s.root",
      "/home/e-work/LabData/ANSG/78mBr/day2/bkg_day2/RAW/"
      "DataR_CH0@DT5730B_969_bkg_day2.root",
      "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/irradiation_1/RAW/"
      "DataR_CH0@DT5730S_31017_irradiation_1.root",
      "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/irradiation_2/RAW/"
      "DataR_CH0@DT5730S_31017_irradiation_2.root",
      "/home/e-work/LabData/ANSG/78mBr/half_life_2/DAQ/irradiation_3/RAW/"
      "DataR_CH0@DT5730S_31017_irradiation_3.root",
      "/home/e-work/LabData/ANSG/78mBr/day2/irradiation_day2/RAW/"
      "DataR_CH0@DT5730B_969_irradiation_day2.root"};

  std::vector<TString> output_names = {
      Constants::CALIBRATION_AM241, Constants::CALIBRATION_EU152,
      Constants::BACKGROUND,        Constants::IRRADIATION_ONE,
      Constants::IRRADIATION_TWO,   Constants::IRRADIATION_THREE,
      Constants::IRRADIATION_FOUR};

  ProcessFiles(filepaths, output_names);
}
