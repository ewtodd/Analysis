#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iomanip>
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

  TString histName = use_calibrated ? "calibrated_hist" : "hist";
  TH1F *zoomedHist = static_cast<TH1F *>(file->Get(histName));

  zoomedHist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResultDetailed result;

  if (peak_name == "Am_59.5keV") {
    fitter =
        new FittingUtils(zoomedHist, 50, 65, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE);
  }
  if (peak_name == "Ba_80.98keV") {
    fitter = new FittingUtils(zoomedHist, 70, 90, kFALSE, kTRUE, kTRUE, kTRUE,
                              kTRUE);
  }
  if (peak_name == "Cd114m_95.9keV") {
    fitter = new FittingUtils(zoomedHist, 92, 100, kTRUE, kTRUE, kTRUE, kTRUE,
                              kTRUE);
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
    fitter = new FittingUtils(zoomedHist, 45, 70, kFALSE, kTRUE, kTRUE, kTRUE,
                              kTRUE);
  }
  if (peak_name == "Pb_KAlpha") {
    if (input_name == Constants::CDSHIELDBACKGROUND_25PERCENT_20260113) {
      fitter = new FittingUtils(zoomedHist, 63, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260114) {
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
    if (input_name == Constants::CDSHIELDSIGNAL_25PERCENT_20260113) {
      fitter = new FittingUtils(zoomedHist, 63, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name == Constants::CUSHIELDSIGNAL_10PERCENT_20260114) {
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

void CheckDrift() {
  InitUtils::SetROOTPreferences();

  FitResultDetailed zero_result;
  zero_result.mu = 0;
  zero_result.mu_error = 0;
  zero_result.sigma = 0;
  zero_result.sigma_error = 0;

  FitResultDetailed calibration_20260112_result =
      FitSinglePeak(Constants::CALIBRATION_20260112, "Am_59.5keV", 59.5409);

  FitResultDetailed postreactor_am241_20260113_result = FitSinglePeak(
      Constants::POSTREACTOR_AM241_20260113, "Am_59.5keV", 59.5409);

  FitResultDetailed postreactor_am241_20260115_result = FitSinglePeak(
      Constants::POSTREACTOR_AM241_20260115, "Am_59.5keV", 59.5409);

  FitResultDetailed postreactor_ba133_20260115_result = FitSinglePeak(
      Constants::POSTREACTOR_BA133_20260115, "Ba_80.98keV", 80.98);

  FitResultDoublePeakDetailed am_ba =
      FitDoublePeak(Constants::POSTREACTOR_AM241_BA133_20260116,
                    "Am_59.5_Ba_53.16keV", 53.1, 59.5);

  FitResultDoublePeakDetailed cd_shield_background_10_percent_20260113 =
      FitDoublePeak(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
                    "Pb_KAlpha", 72.8042, 74.9694);

  FitResultDoublePeakDetailed cd_shield_background_25_percent_20260113 =
      FitDoublePeak(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
                    "Pb_KAlpha", 72.8042, 74.9694);

  FitResultDoublePeakDetailed cu_shield_background_10_percent_20260113 =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                    "Pb_KAlpha", 72.8042, 74.9694);

  FitResultDoublePeakDetailed cu_shield_background_10_percent_20260114 =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                    "Pb_KAlpha", 72.8042, 74.9694);

  // Cd-114m activation peak fits for gain monitoring
  FitResultDetailed cd_shield_signal_25_percent_20260113_cd114m = FitSinglePeak(
      Constants::CDSHIELDSIGNAL_25PERCENT_20260113, "Cd114m_95.9keV", 95.9023);

  FitResultDetailed cd_shield_signal_10_percent_20260113_cd114m = FitSinglePeak(
      Constants::CDSHIELDSIGNAL_10PERCENT_20260113, "Cd114m_95.9keV", 95.9023);

  FitResultDetailed cu_shield_signal_10_percent_20260113_cd114m = FitSinglePeak(
      Constants::CUSHIELDSIGNAL_10PERCENT_20260113, "Cd114m_95.9keV", 95.9023);

  FitResultDetailed cu_shield_signal_10_percent_20260114_cd114m = FitSinglePeak(
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114, "Cd114m_95.9keV", 95.9023);

  FitResultDetailed noshieldsignal_5percent_20260115_cd114m = FitSinglePeak(
      Constants::NOSHIELDSIGNAL_5PERCENT_20260115, "Cd114m_95.9keV", 95.9023);

  FitResultDetailed noshield_graphitecastlesignal_10percent_20260116_cd114m =
      FitSinglePeak(Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
                    "Cd114m_95.9keV", 95.9023);

  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> calibration_values_keV;
  std::vector<Float_t> reduced_chi2;
  std::vector<TString> run_names;

  run_names.push_back("Zero Point");
  mu.push_back(zero_result.mu);
  mu_errors.push_back(zero_result.mu_error);
  calibration_values_keV.push_back(0);
  reduced_chi2.push_back(0);

  run_names.push_back("Calibration Am-241 (01/12/2026)");
  mu.push_back(calibration_20260112_result.mu);
  mu_errors.push_back(calibration_20260112_result.mu_error);
  calibration_values_keV.push_back(59.5409);
  reduced_chi2.push_back(calibration_20260112_result.reduced_chi2);

  run_names.push_back("Post-reactor Am-241 (01/13/2026)");
  mu.push_back(postreactor_am241_20260113_result.mu);
  mu_errors.push_back(postreactor_am241_20260113_result.mu_error);
  calibration_values_keV.push_back(59.5409);
  reduced_chi2.push_back(postreactor_am241_20260113_result.reduced_chi2);

  run_names.push_back("Post-reactor Am-241 (01/15/2026)");
  mu.push_back(postreactor_am241_20260115_result.mu);
  mu_errors.push_back(postreactor_am241_20260115_result.mu_error);
  calibration_values_keV.push_back(59.5409);
  reduced_chi2.push_back(postreactor_am241_20260115_result.reduced_chi2);

  run_names.push_back("Post-reactor Ba-133 (01/15/2026)");
  mu.push_back(postreactor_ba133_20260115_result.mu);
  mu_errors.push_back(postreactor_ba133_20260115_result.mu_error);
  calibration_values_keV.push_back(80.98);
  reduced_chi2.push_back(postreactor_ba133_20260115_result.reduced_chi2);

  run_names.push_back("Am-241+Ba-133 Ba-53.16keV (01/16/2026)");
  mu.push_back(am_ba.peak1.mu);
  mu_errors.push_back(am_ba.peak1.mu_error);
  calibration_values_keV.push_back(53.16);
  reduced_chi2.push_back(am_ba.reduced_chi2);

  run_names.push_back("Am-241+Ba-133 Am-59.5keV (01/16/2026)");
  mu.push_back(am_ba.peak2.mu);
  mu_errors.push_back(am_ba.peak2.mu_error);
  calibration_values_keV.push_back(59.5409);
  reduced_chi2.push_back(-1);

  run_names.push_back("Am-241+Ba-133 Ba-80.98keV (01/16/2026)");
  mu.push_back(postreactor_ba133_20260115_result.mu);
  mu_errors.push_back(postreactor_ba133_20260115_result.mu_error);
  calibration_values_keV.push_back(80.98);
  reduced_chi2.push_back(postreactor_ba133_20260115_result.reduced_chi2);

  run_names.push_back("Cd Shield Bkg 10% Pb-Ka1 (01/13/2026)");
  mu.push_back(cd_shield_background_10_percent_20260113.peak1.mu);
  mu_errors.push_back(cd_shield_background_10_percent_20260113.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);
  reduced_chi2.push_back(cd_shield_background_10_percent_20260113.reduced_chi2);

  run_names.push_back("Cd Shield Bkg 10% Pb-Ka2 (01/13/2026)");
  mu.push_back(cd_shield_background_10_percent_20260113.peak2.mu);
  mu_errors.push_back(cd_shield_background_10_percent_20260113.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);
  reduced_chi2.push_back(-1);

  run_names.push_back("Cd Shield Bkg 25% Pb-Ka1 (01/13/2026)");
  mu.push_back(cd_shield_background_25_percent_20260113.peak1.mu);
  mu_errors.push_back(cd_shield_background_25_percent_20260113.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);
  reduced_chi2.push_back(cd_shield_background_25_percent_20260113.reduced_chi2);

  run_names.push_back("Cd Shield Bkg 25% Pb-Ka2 (01/13/2026)");
  mu.push_back(cd_shield_background_25_percent_20260113.peak2.mu);
  mu_errors.push_back(cd_shield_background_25_percent_20260113.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);
  reduced_chi2.push_back(-1);

  run_names.push_back("Cu Shield Bkg 10% Pb-Ka1 (01/13/2026)");
  mu.push_back(cu_shield_background_10_percent_20260113.peak1.mu);
  mu_errors.push_back(cu_shield_background_10_percent_20260113.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);
  reduced_chi2.push_back(cu_shield_background_10_percent_20260113.reduced_chi2);

  run_names.push_back("Cu Shield Bkg 10% Pb-Ka2 (01/13/2026)");
  mu.push_back(cu_shield_background_10_percent_20260113.peak2.mu);
  mu_errors.push_back(cu_shield_background_10_percent_20260113.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);
  reduced_chi2.push_back(-1);

  run_names.push_back("Cu Shield Bkg 10% Pb-Ka1 (01/14/2026)");
  mu.push_back(cu_shield_background_10_percent_20260114.peak1.mu);
  mu_errors.push_back(cu_shield_background_10_percent_20260114.peak1.mu_error);
  calibration_values_keV.push_back(72.8042);
  reduced_chi2.push_back(cu_shield_background_10_percent_20260114.reduced_chi2);

  run_names.push_back("Cu Shield Bkg 10% Pb-Ka2 (01/14/2026)");
  mu.push_back(cu_shield_background_10_percent_20260114.peak2.mu);
  mu_errors.push_back(cu_shield_background_10_percent_20260114.peak2.mu_error);
  calibration_values_keV.push_back(74.9694);
  reduced_chi2.push_back(-1);

  run_names.push_back("Cd Shield Signal 25% Cd-114m (01/13/2026)");
  mu.push_back(cd_shield_signal_25_percent_20260113_cd114m.mu);
  mu_errors.push_back(cd_shield_signal_25_percent_20260113_cd114m.mu_error);
  calibration_values_keV.push_back(95.9023);
  reduced_chi2.push_back(
      cd_shield_signal_25_percent_20260113_cd114m.reduced_chi2);

  run_names.push_back("Cd Shield Signal 10% Cd-114m (01/13/2026)");
  mu.push_back(cd_shield_signal_10_percent_20260113_cd114m.mu);
  mu_errors.push_back(cd_shield_signal_10_percent_20260113_cd114m.mu_error);
  calibration_values_keV.push_back(95.9023);
  reduced_chi2.push_back(
      cd_shield_signal_10_percent_20260113_cd114m.reduced_chi2);

  run_names.push_back("Cu Shield Signal 10% Cd-114m (01/13/2026)");
  mu.push_back(cu_shield_signal_10_percent_20260113_cd114m.mu);
  mu_errors.push_back(cu_shield_signal_10_percent_20260113_cd114m.mu_error);
  calibration_values_keV.push_back(95.9023);
  reduced_chi2.push_back(
      cu_shield_signal_10_percent_20260113_cd114m.reduced_chi2);

  run_names.push_back("Cu Shield Signal 10% Cd-114m (01/14/2026)");
  mu.push_back(cu_shield_signal_10_percent_20260114_cd114m.mu);
  mu_errors.push_back(cu_shield_signal_10_percent_20260114_cd114m.mu_error);
  calibration_values_keV.push_back(95.9023);
  reduced_chi2.push_back(
      cu_shield_signal_10_percent_20260114_cd114m.reduced_chi2);

  run_names.push_back("No Shield Signal 5% Cd-114m (01/15/2026)");
  mu.push_back(noshieldsignal_5percent_20260115_cd114m.mu);
  mu_errors.push_back(noshieldsignal_5percent_20260115_cd114m.mu_error);
  calibration_values_keV.push_back(95.9023);
  reduced_chi2.push_back(noshieldsignal_5percent_20260115_cd114m.reduced_chi2);

  run_names.push_back(
      "No Shield Graphite Castle Signal 10% Cd-114m (01/16/2026)");
  mu.push_back(noshield_graphitecastlesignal_10percent_20260116_cd114m.mu);
  mu_errors.push_back(
      noshield_graphitecastlesignal_10percent_20260116_cd114m.mu_error);
  calibration_values_keV.push_back(95.9023);
  reduced_chi2.push_back(
      noshield_graphitecastlesignal_10percent_20260116_cd114m.reduced_chi2);

  std::cout << "\nCalibration Points Summary:" << std::endl;
  for (size_t i = 0; i < mu.size(); ++i) {
    std::cout << std::left << std::setw(55) << run_names[i] << ": "
              << std::fixed << std::setprecision(4) << mu[i] << " +/- "
              << mu_errors[i] << " keV";
    if (reduced_chi2[i] >= 0) {
      std::cout << " (χ²/ndf = " << std::setprecision(3) << reduced_chi2[i]
                << ")";
    }
    std::cout << std::endl;
  }

  TF1 *calibration_function =
      CreateAndSaveCalibration(mu, calibration_values_keV, mu_errors);
}
