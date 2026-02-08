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
                                const Float_t expected_mu,
                                const Bool_t use_calibrated = kFALSE) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TString histName = use_calibrated ? "calibrated_zoomedHist" : "zoomedHist";
  TH1F *zoomedHist = static_cast<TH1F *>(file->Get(histName));

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

FitResultDoublePeakDetailed
FitDoublePeak(const TString input_name, const TString peak_name,
              const Float_t mu1_init, const Float_t mu2_init,
              const Bool_t use_calibrated = kFALSE) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TString histName = use_calibrated ? "calibrated_zoomedHist" : "zoomedHist";
  TH1F *zoomedHist = static_cast<TH1F *>(file->Get(histName));

  zoomedHist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResultDoublePeakDetailed result;

  if (peak_name == "Am_59.5_Ba_53.16keV") {
    fitter = new FittingUtils(zoomedHist, 47, 67, kFALSE, kTRUE, kTRUE, kTRUE,
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

FitResultTriplePeakDetailed
FitTriplePeak(const TString input_name, const TString peak_name,
              const FitResultDoublePeakDetailed &constrained_peaks,
              const Float_t mu3_init, const Bool_t use_calibrated = kFALSE) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TString histName = use_calibrated ? "calibrated_zoomedHist" : "zoomedHist";
  TH1F *zoomedHist = static_cast<TH1F *>(file->Get(histName));

  zoomedHist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResultTriplePeakDetailed result;

  if (peak_name == "Ge") {
    if (input_name == Constants::CDSHIELDSIGNAL_25PERCENT_01132026) {
      fitter = new FittingUtils(zoomedHist, 63, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name == Constants::CUSHIELDSIGNAL_10PERCENT_01142026) {
      fitter = new FittingUtils(zoomedHist, 60, 80, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);

    } else
      fitter = new FittingUtils(zoomedHist, 65, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
  }

  result = fitter->FitTriplePeakDetailed(input_name, peak_name,
                                         constrained_peaks, mu3_init);
  delete zoomedHist;
  delete fitter;
  return result;
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
  calibration_fit->FixParameter(0, 0);
  calibration_fit->SetParameter(1, 1);
  calibration_fit->SetNpx(1000);

  TFitResultPtr fit_result = calibration_curve->Fit(calibration_fit, "LRE");

  calibration_fit->Draw("SAME");

  PlottingUtils::SaveFigure(canvas, "calibration.png", kFALSE);
  return calibration_fit;
}

void PulseHeightToDepositedEnergy(std::vector<TString> input_names,
                                  TF1 *calibration_function) {
  Int_t entries = input_names.size();

  TString calibration_function_filepath =
      "root_files/calibration_function.root";
  TFile *calibration_file =
      new TFile(calibration_function_filepath, "RECREATE");
  calibration_function->Write("calibration", TObject::kOverwrite);

  for (Int_t i = 0; i < entries; i++) {
    TString input_name = input_names[i];

    TH1F *hist = new TH1F(PlottingUtils::GetRandomName(),
                          Form("; Deposited Energy [keV]; Counts / %d eV",
                               Constants::BIN_WIDTH_EV),
                          Constants::HIST_NBINS, Constants::HIST_XMIN,
                          Constants::HIST_XMAX);

    TH1F *zoomedHist = new TH1F(PlottingUtils::GetRandomName(),
                                Form("; Deposited Energy [keV]; Counts / %d eV",
                                     Constants::BIN_WIDTH_EV),
                                Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                                Constants::ZOOMED_XMAX);

    TH1F *peakHist = new TH1F(PlottingUtils::GetRandomName(),
                              Form("; Deposited Energy [keV]; Counts / %d eV",
                                   Constants::BIN_WIDTH_EV),
                              Constants::PEAK_NBINS, Constants::PEAK_XMIN,
                              Constants::PEAK_XMAX);

    TString output_filepath = "root_files/" + input_name + ".root";

    TFile *output = new TFile(output_filepath, "UPDATE");

    output->cd();

    Bool_t isFiltered = input_name.Contains("_filtered");

    TTree *tree = nullptr;
    TString treeName;
    TString energyBranchName;

    if (isFiltered) {
      treeName = "bef_tree";
      energyBranchName = "energykeV";
      tree = static_cast<TTree *>(output->Get(treeName));
    } else {
      treeName = "bef_tree_event_summary";
      energyBranchName = "totalEnergykeV";
      tree = static_cast<TTree *>(output->Get(treeName));
    }

    if (!tree) {
      std::cerr << "ERROR: Could not find tree '" << treeName << "' in file "
                << output_filepath << std::endl;
      output->Close();
      delete hist;
      delete zoomedHist;
      delete peakHist;
      continue;
    }

    Float_t energy, deposited_energy_keV;
    tree->SetBranchAddress(energyBranchName, &energy);
    tree->Branch("deposited_energy", &deposited_energy_keV,
                 "deposited_energy/F");

    Int_t num_entries = tree->GetEntries();

    tree->LoadBaskets();
    for (Int_t j = 0; j < num_entries; j++) {
      tree->GetEntry(j);
      deposited_energy_keV = calibration_function->Eval(energy);
      tree->GetBranch("deposited_energy")->Fill();
      hist->Fill(deposited_energy_keV);
      zoomedHist->Fill(deposited_energy_keV);
      peakHist->Fill(deposited_energy_keV);
    }

    PlottingUtils::ConfigureHistogram(hist, kViolet);

    tree->Write("", TObject::kOverwrite);
    hist->Write("calibrated_hist", TObject::kOverwrite);
    zoomedHist->Write("calibrated_zoomedHist", TObject::kOverwrite);
    peakHist->Write("calibrated_peakHist", TObject::kOverwrite);
    output->Close();
    delete hist;
    delete zoomedHist;
    delete peakHist;
  }
  calibration_file->Close();
}

void Calibration() {
  InitUtils::SetROOTPreferences();

  FitResultDetailed zero_result;
  zero_result.mu = 0;
  zero_result.mu_error = 0;
  zero_result.sigma = 0;
  zero_result.sigma_error = 0;

  FitResultDetailed calibration_01122026_result =
      FitSinglePeak(Constants::CALIBRATION_01122026, "Am_59.5keV", 59.5409);

  FitResultDetailed postreactor_am241_01132026_result = FitSinglePeak(
      Constants::POSTREACTOR_AM241_01132026, "Am_59.5keV", 59.5409);

  FitResultDetailed postreactor_am241_01152026_result = FitSinglePeak(
      Constants::POSTREACTOR_AM241_01152026, "Am_59.5keV", 59.5409);

  FitResultDetailed postreactor_ba133_01152026_result = FitSinglePeak(
      Constants::POSTREACTOR_BA133_01152026, "Ba_80.98keV", 80.98);

  FitResultDoublePeakDetailed am_ba =
      FitDoublePeak(Constants::POSTREACTOR_AM241_BA133_01162026,
                    "Am_59.5_Ba_53.16keV", 53.1, 59.5);

  FitResultDoublePeakDetailed cd_shield_background_10_percent_01132026 =
      FitDoublePeak(Constants::CDSHIELDBACKGROUND_10PERCENT_01132026,
                    "Pb_KAlpha", 72.8042, 74.9694);

  FitResultDoublePeakDetailed cd_shield_background_25_percent_01132026 =
      FitDoublePeak(Constants::CDSHIELDBACKGROUND_25PERCENT_01132026,
                    "Pb_KAlpha", 72.8042, 74.9694);

  FitResultDoublePeakDetailed cu_shield_background_10_percent_01132026 =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_01132026,
                    "Pb_KAlpha", 72.8042, 74.9694);

  FitResultDoublePeakDetailed cu_shield_background_10_percent_01142026 =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_01142026,
                    "Pb_KAlpha", 72.8042, 74.9694);

  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> calibration_values_keV;

  mu.push_back(zero_result.mu);
  mu_errors.push_back(zero_result.mu_error);
  calibration_values_keV.push_back(0);

  mu.push_back(calibration_01122026_result.mu);
  mu_errors.push_back(calibration_01122026_result.mu_error);
  calibration_values_keV.push_back(59.5409);

  mu.push_back(postreactor_am241_01132026_result.mu);
  mu_errors.push_back(postreactor_am241_01132026_result.mu_error);
  calibration_values_keV.push_back(59.5409);

  //  mu.push_back(postreactor_am241_01152026_result.mu);
  //  mu_errors.push_back(postreactor_am241_01152026_result.mu_error);
  //  calibration_values_keV.push_back(59.5409);
  //
  //  mu.push_back(postreactor_ba133_01152026_result.mu);
  //  mu_errors.push_back(postreactor_ba133_01152026_result.mu_error);
  //  calibration_values_keV.push_back(80.98);

  //  mu.push_back(am_ba.peak1.mu);
  //  mu_errors.push_back(am_ba.peak1.mu_error);
  //  calibration_values_keV.push_back(53.16);
  //
  //  mu.push_back(am_ba.peak2.mu);
  //  mu_errors.push_back(am_ba.peak2.mu_error);
  //  calibration_values_keV.push_back(59.5409);

  mu.push_back(cd_shield_background_10_percent_01132026.peak1.mu);
  mu_errors.push_back(cd_shield_background_10_percent_01132026.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);

  mu.push_back(cd_shield_background_10_percent_01132026.peak2.mu);
  mu_errors.push_back(cd_shield_background_10_percent_01132026.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);

  mu.push_back(cd_shield_background_25_percent_01132026.peak1.mu);
  mu_errors.push_back(cd_shield_background_25_percent_01132026.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);

  mu.push_back(cd_shield_background_25_percent_01132026.peak2.mu);
  mu_errors.push_back(cd_shield_background_25_percent_01132026.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);

  mu.push_back(cu_shield_background_10_percent_01132026.peak1.mu);
  mu_errors.push_back(cu_shield_background_10_percent_01132026.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);

  mu.push_back(cu_shield_background_10_percent_01132026.peak2.mu);
  mu_errors.push_back(cu_shield_background_10_percent_01132026.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);

  mu.push_back(cu_shield_background_10_percent_01142026.peak1.mu);
  mu_errors.push_back(cu_shield_background_10_percent_01142026.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);

  mu.push_back(cu_shield_background_10_percent_01142026.peak2.mu);
  mu_errors.push_back(cu_shield_background_10_percent_01142026.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);

  TF1 *calibration_function =
      CreateAndSaveCalibration(mu, calibration_values_keV, mu_errors);
  PulseHeightToDepositedEnergy(Constants::ALL_DATASETS, calibration_function);
}
