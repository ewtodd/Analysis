#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <cmath>
#include <iomanip>
#include <vector>

const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_53 = 53.16;
const Float_t E_BA133_81 = 80.98;
const Float_t E_PB_KA1 = 72.8042;
const Float_t E_PB_KA2 = 74.9694;
const Float_t E_CD114M = 95.9023;

struct CalibrationData {
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> calibration_values_keV;
  std::vector<Float_t> reduced_chi2;
  std::vector<TString> run_names;
};

void PrintCalibrationSummary(const CalibrationData &cal_data,
                             TString date_label);

static const Int_t N_CRYSTALS = 4;

TH1F *LoadHistogram(const TString input_name, Int_t crystal,
                    const Bool_t use_calibrated = kFALSE) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TString histName = use_calibrated
      ? Form("calibrated_zoomedHist_crystal%d", crystal)
      : Form("zoomedHist_crystal%d", crystal);
  TH1F *hist = static_cast<TH1F *>(file->Get(histName));
  if (!hist) {
    std::cerr << "ERROR: Cannot find " << histName << " in " << input_name
              << std::endl;
    file->Close();
    delete file;
    return nullptr;
  }
  hist->SetDirectory(0);
  file->Close();
  delete file;
  return hist;
}

FitResult FitCalibrationPeak(const TString input_name, const TString peak_name,
                             Int_t crystal, const Bool_t use_calibrated,
                             const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name, crystal, use_calibrated);
  if (!hist)
    return {};

  Bool_t use_flat_background = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter = nullptr;

  if (peak_name == "Am_59.5keV") {
    if (input_name == Constants::POSTREACTOR_AM241_20260113)
      fitter = new FittingUtils(hist, 50, 70, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
    else if (input_name == Constants::POSTREACTOR_AM241_BA133_20260116)
      fitter = new FittingUtils(hist, 55, 70, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
    else
      fitter = new FittingUtils(hist, 51, 71, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
  }
  if (peak_name == "Ba_80.98keV") {
    fitter =
        new FittingUtils(hist, 75, 90, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  }
  if (peak_name == "Cd114m_95.9keV") {
    if (input_name == Constants::NOSHIELDBACKGROUND_5PERCENT_20260115)
      fitter = new FittingUtils(hist, 91, 100, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
    else
      fitter = new FittingUtils(hist, 91, 103, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
  }

  if (interactive)
    fitter->SetInteractive();
  TString label = Form("%s_crystal%d", input_name.Data(), crystal);
  FitResult result = fitter->FitSinglePeak(label, peak_name);
  delete hist;
  delete fitter;
  return result;
}

FitResult FitPbKAlpha(const TString input_name, Int_t crystal,
                      const Bool_t use_calibrated,
                      const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name, crystal, use_calibrated);
  if (!hist)
    return {};

  Bool_t use_flat_background = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter = nullptr;

  if (input_name == Constants::CDSHIELDBACKGROUND_25PERCENT_20260113)
    fitter =
        new FittingUtils(hist, 66, 81, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260114)
    fitter =
        new FittingUtils(hist, 66, 82, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260113) {
    use_flat_background = kFALSE;
    fitter =
        new FittingUtils(hist, 65, 82, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  } else {
    use_flat_background = kFALSE;
    fitter =
        new FittingUtils(hist, 65, 81, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  }

  if (interactive)
    fitter->SetInteractive();
  TString label = Form("%s_crystal%d", input_name.Data(), crystal);
  FitResult result =
      fitter->FitDoublePeak(label, "Pb_KAlpha", E_PB_KA1, E_PB_KA2);
  delete hist;
  delete fitter;
  return result;
}

TF1 *CreateAndSaveCalibration(
    const std::vector<Float_t> &mu,
    const std::vector<Float_t> &calibration_values_keV,
    const std::vector<Float_t> &mu_errors, const TString date_label) {

  Int_t size = calibration_values_keV.size();

  TGraph *calibration_curve =
      new TGraph(size, mu.data(), calibration_values_keV.data());
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(
      calibration_curve, kBlue,
      "; Precalibrated Energy [keV]; Deposited Energy [keV]");

  calibration_curve->GetXaxis()->SetRangeUser(-5, 100);
  calibration_curve->GetYaxis()->SetRangeUser(-5, 100);
  calibration_curve->GetXaxis()->SetNdivisions(506);
  calibration_curve->SetMarkerStyle(5);
  calibration_curve->SetMarkerSize(2);
  calibration_curve->Draw("AP");

  TF1 *calibration_fit = new TF1("linear_" + date_label, "pol1", -10, 100);
  calibration_fit->FixParameter(0, 0);
  calibration_fit->SetParameter(1, 1);
  calibration_fit->SetNpx(1000);

  TFitResultPtr fit_result = calibration_curve->Fit(calibration_fit, "LRE");

  calibration_fit->Draw("SAME");

  PlottingUtils::SaveFigure(canvas, "calibration_" + date_label, "",
                            PlotSaveOptions::kLINEAR);
  return calibration_fit;
}

void AddZeroPoint(CalibrationData &cal_data) {
  cal_data.run_names.push_back("Zero Point");
  cal_data.mu.push_back(0);
  cal_data.mu_errors.push_back(0);
  cal_data.calibration_values_keV.push_back(0);
  cal_data.reduced_chi2.push_back(0);
}

void AddCalibrationPoint(CalibrationData &cal_data, const TString run_name,
                         const Float_t mu, const Float_t mu_error,
                         const Float_t true_energy,
                         const Float_t reduced_chi2) {
  cal_data.run_names.push_back(run_name);
  cal_data.mu.push_back(mu);
  cal_data.mu_errors.push_back(mu_error);
  cal_data.calibration_values_keV.push_back(true_energy);
  cal_data.reduced_chi2.push_back(reduced_chi2);
}

void PrintCalibrationSummary(const CalibrationData &cal_data,
                             TString date_label) {
  std::cout << "Calibration Points for " << date_label << std::endl;
  for (size_t i = 0; i < cal_data.mu.size(); ++i) {
    std::cout << std::left << std::setw(45) << cal_data.run_names[i] << ": "
              << std::fixed << std::setprecision(4) << cal_data.mu[i] << " +/- "
              << cal_data.mu_errors[i] << " keV";
    if (cal_data.reduced_chi2[i] > 0) {
      std::cout << " (chi2/ndf = " << std::setprecision(3)
                << cal_data.reduced_chi2[i] << ")";
    }
    std::cout << std::endl;
  }
}

TF1 *BuildCalibration_20260113(Int_t crystal, const Bool_t use_calibrated,
                               const Bool_t interactive) {
  TString crystalLabel = Form("crystal%d", crystal);
  std::cout << "Building Calibration for Jan 13 " << crystalLabel << std::endl;

  CalibrationData cal_data;
  AddZeroPoint(cal_data);

  FitResult am241_result =
      FitCalibrationPeak(Constants::POSTREACTOR_AM241_20260113, "Am_59.5keV",
                         crystal, use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Am-241", am241_result.peaks.at(0).mu,
      am241_result.peaks.at(0).mu_error, E_AM241, am241_result.reduced_chi2);

  FitResult cu_bkg_pb =
      FitPbKAlpha(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                  crystal, use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka1",
                      cu_bkg_pb.peaks.at(0).mu, cu_bkg_pb.peaks.at(0).mu_error,
                      E_PB_KA1, cu_bkg_pb.reduced_chi2);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka2",
                      cu_bkg_pb.peaks.at(1).mu, cu_bkg_pb.peaks.at(1).mu_error,
                      E_PB_KA2, -1);

  FitResult cu_bkg_cd =
      FitCalibrationPeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                         "Cd114m_95.9keV", crystal, use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Cd-114m",
                      cu_bkg_cd.peaks.at(0).mu, cu_bkg_cd.peaks.at(0).mu_error,
                      E_CD114M, cu_bkg_cd.reduced_chi2);

  TString dateLabel = Form("20260113_%s", crystalLabel.Data());
  PrintCalibrationSummary(cal_data, dateLabel);

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, dateLabel);
  return cal_func;
}

TF1 *BuildCalibration_20260114(Int_t crystal, const Bool_t use_calibrated,
                               const Bool_t interactive) {
  TString crystalLabel = Form("crystal%d", crystal);
  std::cout << "Building Calibration for Jan 14 " << crystalLabel << std::endl;

  CalibrationData cal_data;
  AddZeroPoint(cal_data);

  FitResult cu_bkg_pb =
      FitPbKAlpha(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                  crystal, use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka1",
                      cu_bkg_pb.peaks.at(0).mu, cu_bkg_pb.peaks.at(0).mu_error,
                      E_PB_KA1, cu_bkg_pb.reduced_chi2);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka2",
                      cu_bkg_pb.peaks.at(1).mu, cu_bkg_pb.peaks.at(1).mu_error,
                      E_PB_KA2, -1);

  FitResult cu_bkg_cd =
      FitCalibrationPeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                         "Cd114m_95.9keV", crystal, use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Cd-114m",
                      cu_bkg_cd.peaks.at(0).mu, cu_bkg_cd.peaks.at(0).mu_error,
                      E_CD114M, cu_bkg_cd.reduced_chi2);

  TString dateLabel = Form("20260114_%s", crystalLabel.Data());
  PrintCalibrationSummary(cal_data, dateLabel);

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, dateLabel);
  return cal_func;
}

TF1 *BuildCalibration_20260115(Int_t crystal, const Bool_t use_calibrated,
                               const Bool_t interactive) {
  TString crystalLabel = Form("crystal%d", crystal);
  std::cout << "Building Calibration for Jan 15 " << crystalLabel << std::endl;

  CalibrationData cal_data;
  AddZeroPoint(cal_data);

  FitResult am241_result =
      FitCalibrationPeak(Constants::POSTREACTOR_AM241_20260115, "Am_59.5keV",
                         crystal, use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Am-241", am241_result.peaks.at(0).mu,
      am241_result.peaks.at(0).mu_error, E_AM241, am241_result.reduced_chi2);

  FitResult ba133_result =
      FitCalibrationPeak(Constants::POSTREACTOR_BA133_20260115, "Ba_80.98keV",
                         crystal, use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Ba-133 80.98keV", ba133_result.peaks.at(0).mu,
      ba133_result.peaks.at(0).mu_error, E_BA133_81, ba133_result.reduced_chi2);

  FitResult cd_fit =
      FitCalibrationPeak(Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
                         "Cd114m_95.9keV", crystal, use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "No Shield Bkg Cd-114m", cd_fit.peaks.at(0).mu,
                      cd_fit.peaks.at(0).mu_error, E_CD114M,
                      cd_fit.reduced_chi2);

  TString dateLabel = Form("20260115_%s", crystalLabel.Data());
  PrintCalibrationSummary(cal_data, dateLabel);

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, dateLabel);
  return cal_func;
}

TF1 *BuildCalibration_20260116(Int_t crystal, const Bool_t use_calibrated,
                               const Bool_t interactive) {
  TString crystalLabel = Form("crystal%d", crystal);
  std::cout << "Building Calibration for Jan 16 " << crystalLabel << std::endl;

  CalibrationData cal_data;
  AddZeroPoint(cal_data);

  FitResult am241_result =
      FitCalibrationPeak(Constants::POSTREACTOR_AM241_BA133_20260116,
                         "Am_59.5keV", crystal, use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Am-241", am241_result.peaks.at(0).mu,
      am241_result.peaks.at(0).mu_error, E_AM241, am241_result.reduced_chi2);

  FitResult ba133_result =
      FitCalibrationPeak(Constants::POSTREACTOR_AM241_BA133_20260116,
                         "Ba_80.98keV", crystal, use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Ba-133 80.98keV", ba133_result.peaks.at(0).mu,
      ba133_result.peaks.at(0).mu_error, E_BA133_81, ba133_result.reduced_chi2);

  FitResult cd_fit = FitCalibrationPeak(
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      "Cd114m_95.9keV", crystal, use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "GraphiteCastle Bkg Cd-114m",
                      cd_fit.peaks.at(0).mu, cd_fit.peaks.at(0).mu_error,
                      E_CD114M, cd_fit.reduced_chi2);

  TString dateLabel = Form("20260116_%s", crystalLabel.Data());
  PrintCalibrationSummary(cal_data, dateLabel);

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, dateLabel);
  return cal_func;
}

void PulseHeightToDepositedEnergy(const std::vector<TString> &input_names,
                                  TF1 *calibration_functions[N_CRYSTALS],
                                  const TString date_label) {
  for (Int_t c = 0; c < N_CRYSTALS; c++) {
    TString cal_filepath = Form("root_files/calibration_function_%s_crystal%d.root",
                                date_label.Data(), c);
    TFile *cal_file = new TFile(cal_filepath, "RECREATE");
    calibration_functions[c]->Write("calibration", TObject::kOverwrite);
    cal_file->Close();
    delete cal_file;
  }

  TString treeType = Constants::USE_FILTERED ? "filtered" : "unfiltered";

  Int_t n_files = input_names.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_name = input_names[i];
    TString filepath = "root_files/" + input_name + ".root";
    TFile *file = new TFile(filepath, "UPDATE");

    for (Int_t c = 0; c < N_CRYSTALS; c++) {
      TString treeName = Form("crystal%d_%s_tree", c, treeType.Data());
      TTree *tree = static_cast<TTree *>(file->Get(treeName));
      if (!tree) {
        std::cerr << "ERROR: Could not find " << treeName << " in " << filepath
                  << std::endl;
        continue;
      }

      Float_t energy = 0;
      tree->SetBranchAddress("energykeV", &energy);

      Float_t deposited_energy = 0;
      TString calTreeName = Form("calibrated_crystal%d_tree", c);
      TTree *cal_tree = new TTree(calTreeName,
          Form("Calibrated events for crystal %d", c));
      cal_tree->SetDirectory(file);
      cal_tree->Branch("deposited_energy", &deposited_energy,
                       "deposited_energy/F");

      TH1F *hist = new TH1F(PlottingUtils::GetRandomName(),
                             Form("; Deposited Energy [keV]; Counts / %d eV",
                                  Constants::BIN_WIDTH_EV),
                             Constants::HIST_NBINS, Constants::HIST_XMIN,
                             Constants::HIST_XMAX);
      hist->SetDirectory(0);

      TH1F *zoomedHist = new TH1F(PlottingUtils::GetRandomName(),
                                   Form("; Deposited Energy [keV]; Counts / %d eV",
                                        Constants::BIN_WIDTH_EV),
                                   Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                                   Constants::ZOOMED_XMAX);
      zoomedHist->SetDirectory(0);

      TH1F *peakHist = new TH1F(PlottingUtils::GetRandomName(),
                                 Form("; Deposited Energy [keV]; Counts / %d eV",
                                      Constants::BIN_WIDTH_EV),
                                 Constants::PEAK_NBINS, Constants::PEAK_XMIN,
                                 Constants::PEAK_XMAX);
      peakHist->SetDirectory(0);

      Int_t n_entries = tree->GetEntries();
      for (Int_t j = 0; j < n_entries; j++) {
        tree->GetEntry(j);
        deposited_energy = calibration_functions[c]->Eval(energy);
        cal_tree->Fill();
        hist->Fill(deposited_energy);
        zoomedHist->Fill(deposited_energy);
        peakHist->Fill(deposited_energy);
      }

      file->cd();
      cal_tree->Write(calTreeName, TObject::kOverwrite);
      hist->Write(Form("calibrated_hist_crystal%d", c), TObject::kOverwrite);
      zoomedHist->Write(Form("calibrated_zoomedHist_crystal%d", c), TObject::kOverwrite);
      peakHist->Write(Form("calibrated_peakHist_crystal%d", c), TObject::kOverwrite);

      delete hist;
      delete zoomedHist;
      delete peakHist;
    }

    std::cout << "Calibrated " << input_name << std::endl;
    file->Close();
    delete file;
  }
}

void Calibration() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t use_calibrated = kFALSE;
  Bool_t interactive = kTRUE;

  TF1 *cal_jan13[N_CRYSTALS];
  TF1 *cal_jan14[N_CRYSTALS];
  TF1 *cal_jan15[N_CRYSTALS];
  TF1 *cal_jan16[N_CRYSTALS];

  for (Int_t c = 0; c < N_CRYSTALS; c++) {
    cal_jan13[c] = BuildCalibration_20260113(c, use_calibrated, interactive);
    cal_jan14[c] = BuildCalibration_20260114(c, use_calibrated, interactive);
    cal_jan15[c] = BuildCalibration_20260115(c, use_calibrated, interactive);
    cal_jan16[c] = BuildCalibration_20260116(c, use_calibrated, interactive);
  }

  std::vector<TString> datasets_jan13 = {
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
      Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
      Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
      Constants::POSTREACTOR_AM241_20260113,
      Constants::ACTIVEBACKGROUND_TEST_5PERCENT_20260113,
      Constants::ACTIVEBACKGROUND_TEST_90PERCENT_20260113};
  PulseHeightToDepositedEnergy(datasets_jan13, cal_jan13, "20260113");

  std::vector<TString> datasets_jan14 = {
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_90PERCENT_20260114};
  PulseHeightToDepositedEnergy(datasets_jan14, cal_jan14, "20260114");

  std::vector<TString> datasets_jan15 = {
      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
      Constants::POSTREACTOR_AM241_20260115,
      Constants::POSTREACTOR_BA133_20260115, Constants::SHUTTERCLOSED_20260115};
  PulseHeightToDepositedEnergy(datasets_jan15, cal_jan15, "20260115");

  std::vector<TString> datasets_jan16 = {
      Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116,
      Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116};
  PulseHeightToDepositedEnergy(datasets_jan16, cal_jan16, "20260116");
}
