#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <vector>

FitResultDetailed FitSinglePeak(const TString input_name,
                                const TString peak_name,
                                const Float_t expected_mu) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TH1F *zoomedHist = static_cast<TH1F *>(file->Get("zoomedHist"));

  zoomedHist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResultDetailed result;

  if (peak_name == "Am_59.5keV") {
    fitter =
        new FittingUtils(zoomedHist, 50, 64, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE);
  }
  if (peak_name == "Ba_80.98keV") {
    fitter = new FittingUtils(zoomedHist, 74, 85, kFALSE, kTRUE, kTRUE, kTRUE,
                              kTRUE);
    std::vector<Double_t> params = {
        81.346,   // Mu
        0.755432, // Sigma
        8091.01,  // GausAmplitude
        376.894,  // BkgConst
        -1.46859, // BkgSlope
        1361.7,   // StepAmplitude
        14490,    // LowTailAmplitude
        1.31752,  // LowTailSlope
        9507.4,   // HighTailAmplitude (was fixed)
        0.710555  // HighTailSlope (was fixed)
    };
    fitter->SetManualParameters(params);
  }
  result = fitter->FitPeakDetailed(input_name, peak_name);
  delete zoomedHist;
  delete fitter;
  return result;
}

FitResultDoublePeakDetailed FitDoublePeak(const TString input_name,
                                          const TString peak_name,
                                          const Float_t mu1_init,
                                          const Float_t mu2_init) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TH1F *zoomedHist = static_cast<TH1F *>(file->Get("zoomedHist"));

  zoomedHist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResultDoublePeakDetailed result;

  if (peak_name == "Am_59.5_Ba_53.16keV") {
    fitter = new FittingUtils(zoomedHist, 50, 67, kFALSE, kTRUE, kTRUE, kTRUE,
                              kTRUE);
  }
  if (peak_name == "Pb_KAlpha") {
    if (input_name == Constants::CDSHIELDBACKGROUND_25PERCENT_01132026) {
      fitter = new FittingUtils(zoomedHist, 63, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_01142026) {
      fitter = new FittingUtils(zoomedHist, 60, 80, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);

    } else
      fitter = new FittingUtils(zoomedHist, 65, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
  }

  result =
      fitter->FitDoublePeakDetailed(input_name, peak_name, mu1_init, mu2_init);
  delete zoomedHist;
  delete fitter;
  return result;
}

std::vector<FitResultDetailed>
FitMultiplePeaks(std::vector<TString> input_names,
                 std::vector<TString> peak_names,
                 std::vector<Float_t> mu_guesses) {

  std::vector<FitResultDetailed> results;
  Int_t entries = input_names.size();
  for (Int_t i = 0; i < entries; i++) {
    FitResultDetailed result;
    if (input_names[i] != "zero") {

      Bool_t isFiltered = input_names[i].Contains("_filtered");

      if (isFiltered && !Constants::FILTERED) {
        std::cerr << "File name contains filtered but constant is not set."
                  << std::endl;
      }

      TString tree_name;
      TString energyBranchName;

      if (isFiltered) {
        tree_name = "bef_tree";
        energyBranchName = "energykeV";
      } else {
        tree_name = "bef_tree_event_summary";
        energyBranchName = "totalEnergy";
      }

      result = FitSinglePeak(input_names[i], peak_names[i], mu_guesses[i]);
    } else {
      result.mu = 0;
      result.mu_error = 0;
      result.sigma = 0;
      result.sigma_error = 0;
    }
    results.push_back(result);
  }
  return results;
}

TF1 *CreateAndSaveCalibration(std::vector<Float_t> mu,
                              std::vector<Float_t> calibration_values_keV,
                              std::vector<Float_t> mu_errors) {

  Int_t size = calibration_values_keV.size();

  TGraph *calibration_curve =
      new TGraph(size, mu.data(), calibration_values_keV.data());
  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas);
  PlottingUtils::ConfigureGraph(
      calibration_curve, kBlue,
      "; Precalibrated Energy [keV]; Deposited Energy [keV]");

  calibration_curve->GetXaxis()->SetRangeUser(-5, 100);
  calibration_curve->GetYaxis()->SetRangeUser(-5, 100);
  calibration_curve->GetXaxis()->SetNdivisions(506);
  calibration_curve->SetMarkerStyle(5);
  calibration_curve->SetMarkerSize(2);
  calibration_curve->Draw("AP");

  TF1 *calibration_fit = new TF1("linear", "pol1", -10, 100);
  calibration_fit->SetParameter(0, 0);
  calibration_fit->SetParameter(1, 1);
  calibration_fit->SetNpx(1000);

  TFitResultPtr fit_result = calibration_curve->Fit(calibration_fit, "LRE");

  calibration_fit->Draw("SAME");

  PlottingUtils::SaveFigure(canvas, "calibration.png", kFALSE);
  return calibration_fit;
}

void PulseHeightToDepositedEnergy(
    std::vector<TString> input_names, TF1 *calibration_function,
    std::vector<Int_t> colors = PlottingUtils::GetDefaultColors()) {
  Int_t entries = input_names.size();

  TString calibration_function_filepath =
      "root_files/calibration_function.root";
  TFile *calibration_file =
      new TFile(calibration_function_filepath, "RECREATE");
  calibration_function->Write("calibration", TObject::kOverwrite);

  for (Int_t i = 0; i < entries; i++) {
    TString input_name = input_names[i];
    TH1F *light_output_hist =
        new TH1F("", "; Light Output [keVee]; Counts", 500, 0,
                 calibration_function->Eval(16384));
    TCanvas *canvas = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvas);

    TString output_filepath = "root_files/" + input_name + ".root";

    TFile *output = new TFile(output_filepath, "UPDATE");

    output->cd();

    TTree *features_tree = static_cast<TTree *>(output->Get("features"));

    Float_t pulse_height, light_output_keVee;
    features_tree->SetBranchAddress("pulse_height", &pulse_height);
    features_tree->Branch("light_output", &light_output_keVee,
                          "light_output/F");

    Int_t num_entries = features_tree->GetEntries();

    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      light_output_keVee = calibration_function->Eval(pulse_height);
      features_tree->GetBranch("light_output")->Fill();
      light_output_hist->Fill(light_output_keVee);
    }

    PlottingUtils::ConfigureAndDrawHistogram(light_output_hist, colors[i]);
    PlottingUtils::SaveFigure(canvas, input_name + "_light_output.png");

    features_tree->Write("", TObject::kOverwrite);
    light_output_hist->Write("Light Output", TObject::kOverwrite);
    output->Close();
    delete canvas;
    delete light_output_hist;
  }
  calibration_file->Close();
}

void Calibration() {
  InitUtils::SetROOTPreferences();

  std::vector<TString> input_names_calibrations = {
      "zero", Constants::CALIBRATION_01122026,
      Constants::POSTREACTOR_AM241_01132026,
      Constants::POSTREACTOR_AM241_01152026,
      Constants::POSTREACTOR_BA133_01152026};

  std::vector<TString> peak_names = {"zero", "Am_59.5keV", "Am_59.5keV",
                                     "Am_59.5keV", "Ba_80.98keV"};

  std::vector<Float_t> mu_guesses_single = {0, 59.5409, 59.5409, 59.5409,
                                            80.98};

  std::vector<FitResultDetailed> fit_results =
      FitMultiplePeaks(input_names_calibrations, peak_names, mu_guesses_single);

  FitResultDoublePeakDetailed am_ba =
      FitDoublePeak(Constants::POSTREACTOR_AM241_BA133_01162026,
                    "Am_59.5_Ba_53.16keV", 53.1, 59.5);

  FitResultDoublePeakDetailed pb_kalpha_results_1 =
      FitDoublePeak(Constants::CDSHIELDBACKGROUND_10PERCENT_01132026,
                    "Pb_KAlpha", 72.8042, 74.9694);

  FitResultDoublePeakDetailed pb_kalpha_results_2 =
      FitDoublePeak(Constants::CDSHIELDBACKGROUND_25PERCENT_01132026,
                    "Pb_KAlpha", 72.8042, 74.9694);

  FitResultDoublePeakDetailed pb_kalpha_results_3 =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_01132026,
                    "Pb_KAlpha", 72.8042, 74.9694);

  FitResultDoublePeakDetailed pb_kalpha_results_4 =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_01142026,
                    "Pb_KAlpha", 72.8042, 74.9694);

  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> calibration_values_keV;

  for (Int_t i = 0; i < fit_results.size(); i++) {
    mu.push_back(fit_results[i].mu);
    mu_errors.push_back(fit_results[i].mu_error);
    if (i == 0) {
      calibration_values_keV.push_back(0);
    } else if (i == 1 || i == 2 | i == 3) {
      calibration_values_keV.push_back(59.5409);
    } else if (i == 4) {
      calibration_values_keV.push_back(80.98);
    } else {
      std::cerr << "Warning: Unexpected index " << i << " in fit_results!"
                << std::endl;
    }
  }

  mu.push_back(am_ba.peak1.mu);
  mu_errors.push_back(am_ba.peak1.mu_error);
  calibration_values_keV.push_back(53.16);

  mu.push_back(am_ba.peak2.mu);
  mu_errors.push_back(am_ba.peak2.mu_error);
  calibration_values_keV.push_back(59.5409);

  mu.push_back(pb_kalpha_results_1.peak1.mu);
  mu_errors.push_back(pb_kalpha_results_1.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);

  mu.push_back(pb_kalpha_results_1.peak2.mu);
  mu_errors.push_back(pb_kalpha_results_1.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);

  mu.push_back(pb_kalpha_results_2.peak1.mu);
  mu_errors.push_back(pb_kalpha_results_2.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);

  mu.push_back(pb_kalpha_results_2.peak2.mu);
  mu_errors.push_back(pb_kalpha_results_2.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);

  mu.push_back(pb_kalpha_results_3.peak1.mu);
  mu_errors.push_back(pb_kalpha_results_3.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);

  mu.push_back(pb_kalpha_results_3.peak2.mu);
  mu_errors.push_back(pb_kalpha_results_3.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);

  mu.push_back(pb_kalpha_results_4.peak1.mu);
  mu_errors.push_back(pb_kalpha_results_4.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);

  mu.push_back(pb_kalpha_results_4.peak2.mu);
  mu_errors.push_back(pb_kalpha_results_4.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);

  TF1 *calibration_function =
      CreateAndSaveCalibration(mu, calibration_values_keV, mu_errors);
}
