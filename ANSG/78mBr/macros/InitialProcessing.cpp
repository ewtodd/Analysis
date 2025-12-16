#include "PlottingUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <TROOT.h>
#include <TSystem.h>
#include <algorithm>
#include <iostream>
#include <vector>

void InitialWaveformProcessing(const TString filepath,
                               const TString output_name, const Int_t color) {

  PlottingUtils::SetROOTPreferences();
  WaveformProcessingUtils *processor = new WaveformProcessingUtils();
  processor->SetPolarity(-1);
  processor->SetTriggerThreshold(0.15);
  processor->SetNumberOfSamplesForBaseline(10);
  processor->SetSampleWindows(18, 125);
  processor->SetGates(5, 40, 200);
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

  TString filepath_Am241 =
      "/home/e-work/LABDATA/ANSG/78mBr/half_life_2/DAQ/"
      "59_5keV_calibration_300s/RAW/"
      "DataR_CH0@DT5730S_31017_59_5keV_calibration_300s.root";
  TString output_name_Am241 = "calibration_Am241";
  TString filepath_Eu152 =
      "/home/e-work/LABDATA/ANSG/78mBr/half_life_2/DAQ/"
      "Europium_calibration_300s/RAW/"
      "DataR_CH0@DT5730S_31017_Europium_calibration_300s.root";
  TString output_name_Eu152 = "calibration_Eu152";
  TString filepath_bkg = "/home/e-work/LABDATA/ANSG/78mBr/day2/"
                         "bkg_day2/RAW/DataR_CH0@DT5730B_969_bkg_day2.root";
  TString output_name_bkg = "background";

  TString filepath_irradiation_one =
      "/home/e-work/LABDATA/ANSG/78mBr/half_life_2/DAQ/irradiation_1/RAW/"
      "DataR_CH0@DT5730S_31017_irradiation_1.root";
  TString output_name_irradiation_one = "irradiation_one";

  TString filepath_irradiation_two =
      "/home/e-work/LABDATA/ANSG/78mBr/half_life_2/DAQ/irradiation_2/RAW/"
      "DataR_CH0@DT5730S_31017_irradiation_2.root";
  TString output_name_irradiation_two = "irradiation_two";

  TString filepath_irradiation_three =
      "/home/e-work/LABDATA/ANSG/78mBr/half_life_2/DAQ/irradiation_3/RAW/"
      "DataR_CH0@DT5730S_31017_irradiation_3.root";
  TString output_name_irradiation_three = "irradiation_three";

  TString filepath_irradiation_four =
      "/home/e-work/LABDATA/ANSG/78mBr/day2/irradiation_day2/RAW/"
      "DataR_CH0@DT5730B_969_irradiation_day2.root";
  TString output_name_irradiation_four = "irradiation_four";

  std::vector<TString> filepaths = {filepath_Am241,
                                    filepath_Eu152,
                                    filepath_bkg,
                                    filepath_irradiation_one,
                                    filepath_irradiation_two,
                                    filepath_irradiation_three,
                                    filepath_irradiation_four};
  std::vector<TString> output_names = {output_name_Am241,
                                       output_name_Eu152,
                                       output_name_bkg,
                                       output_name_irradiation_one,
                                       output_name_irradiation_two,
                                       output_name_irradiation_three,
                                       output_name_irradiation_four};

  ProcessFiles(filepaths, output_names);
}
