#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "HyperEMGFitHelpers.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <iomanip>
#include <vector>

const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_53 = 53.16;
const Float_t E_BA133_81 = 80.98;
const Float_t E_BA133_276 = 276.3989;
const Float_t E_BA133_303 = 302.8508;
const Float_t E_BA133_356 = 356.0129;
const Float_t E_BA133_384 = 383.8485;
const Float_t E_PB_KA1 = 72.8042;
const Float_t E_PB_KA2 = 74.9694;
const Float_t E_CD114M = 95.9023;
const Float_t E_ANNIHILATION = 510.999;

struct CalibrationData {
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> calibration_values_keV;
  std::vector<Float_t> reduced_chi2;
  std::vector<TString> run_names;
};

void PrintCalibrationSummary(const CalibrationData &cal_data,
                             TString date_label);

TH1F *LoadHistogram(const TString input_name) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TH1F *hist = static_cast<TH1F *>(file->Get("snip_zoomedHist"));
  if (!hist) {
    std::cerr << "ERROR: No snip_zoomedHist in " << input_name << std::endl;
    file->Close();
    delete file;
    return nullptr;
  }
  hist->SetDirectory(0);
  hist->Rebin(Constants::REBIN_FACTOR);
  hist->SetTitle(Form("; Energy [keV]; Counts / %d eV",
                      Constants::REBIN_FACTOR * Constants::BIN_WIDTH_EV));

  file->Close();
  delete file;
  return hist;
}

TH1F *LoadFullHistogram(const TString input_name) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TH1F *hist = static_cast<TH1F *>(file->Get("snip_hist"));
  if (!hist) {
    std::cerr << "ERROR: No snip_hist in " << input_name << std::endl;
    file->Close();
    delete file;
    return nullptr;
  }
  hist->SetDirectory(0);
  hist->Rebin(Constants::REBIN_FACTOR);
  hist->SetTitle(Form("; Energy [keV]; Counts / %d eV",
                      Constants::REBIN_FACTOR * Constants::BIN_WIDTH_EV));

  file->Close();
  delete file;
  return hist;
}

HyperEMGFitResult FitCalibrationPeakFull(const TString input_name,
                                         const TString peak_name,
                                         const Bool_t interactive) {
  TH1F *hist = LoadFullHistogram(input_name);
  if (!hist)
    return {};

  Double_t fit_low = 0, fit_high = 0;
  if (peak_name == "Ba_276.40keV") {
    fit_low = 270;
    fit_high = 283;
  } else if (peak_name == "Ba_302.85keV") {
    fit_low = 296;
    fit_high = 310;
  } else if (peak_name == "Ba_356.01keV") {
    fit_low = 349;
    fit_high = 363;
  } else if (peak_name == "Ba_383.85keV") {
    fit_low = 377;
    fit_high = 391;
  } else if (peak_name == "Annihilation_511keV") {
    fit_low = 500;
    fit_high = 525;
  } else {
    delete hist;
    return {};
  }

  HyperEMGFitResult result = FitSingleHyperEMG(
      hist, fit_low, fit_high, kTRUE, input_name, peak_name, interactive);
  delete hist;
  return result;
}

HyperEMGFitResult FitCalibrationPeak(const TString input_name,
                                     const TString peak_name,
                                     const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Double_t fit_low = 0, fit_high = 0;

  if (peak_name == "Am_59.5keV") {
    if (input_name == Constants::POSTREACTOR_AM241_20260113) {
      fit_low = 50;
      fit_high = 70;
    } else if (input_name == Constants::POSTREACTOR_AM241_BA133_20260116) {
      fit_low = 55;
      fit_high = 70;
    } else {
      fit_low = 51;
      fit_high = 71;
    }
  } else if (peak_name == "Ba_53.16keV") {
    fit_low = 48;
    fit_high = 58;
  } else if (peak_name == "Ba_80.98keV") {
    fit_low = 75;
    fit_high = 90;
  } else if (peak_name == "Cd114m_95.9keV") {
    if (input_name == Constants::NOSHIELDBACKGROUND_5PERCENT_20260115) {
      fit_low = 91;
      fit_high = 100;
    } else {
      fit_low = 91;
      fit_high = 103;
    }
  }

  HyperEMGFitResult result = FitSingleHyperEMG(
      hist, fit_low, fit_high, kTRUE, input_name, peak_name, interactive);
  delete hist;
  return result;
}

HyperEMGFitResult FitPbKAlpha(const TString input_name,
                              const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Double_t fit_low = 0, fit_high = 0;

  if (input_name == Constants::CDSHIELDBACKGROUND_25PERCENT_20260113) {
    fit_low = 66;
    fit_high = 81;
  } else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260114) {
    fit_low = 66;
    fit_high = 82;
  } else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260113) {
    fit_low = 65;
    fit_high = 82;
  } else {
    fit_low = 65;
    fit_high = 81;
  }

  HyperEMGFitResult result =
      FitDoubleHyperEMG(hist, fit_low, fit_high, E_PB_KA1, E_PB_KA2, kTRUE,
                        input_name, "Pb_KAlpha", interactive);
  delete hist;
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

  calibration_curve->GetXaxis()->SetRangeUser(-5, 600);
  calibration_curve->GetYaxis()->SetRangeUser(-5, 600);
  calibration_curve->GetXaxis()->SetNdivisions(506);
  calibration_curve->SetMarkerStyle(5);
  calibration_curve->SetMarkerSize(2);
  calibration_curve->Draw("AP");

  TF1 *calibration_fit = new TF1("cal_" + date_label, "pol2", -10, 600);
  calibration_fit->SetParameter(0, 0);
  calibration_fit->SetParameter(1, 1);
  calibration_fit->SetParameter(2, 0);
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

// Build master calibration from Jan 15 (most calibration sources available)
// Uses Am-241, Ba-133 (53, 81, 303, 356, 384), 511 keV → pol2
TF1 *BuildMasterCalibration(const Bool_t interactive) {
  std::cout << "Building Master Calibration from Jan 15" << std::endl;

  CalibrationData cal_data;

  TString am_dataset = Constants::POSTREACTOR_AM241_20260115;
  TString ba_dataset = Constants::POSTREACTOR_BA133_20260115;
  TString sig_dataset = Constants::NOSHIELDSIGNAL_5PERCENT_20260115;

  // Am-241 59.5 keV
  HyperEMGFitResult am =
      FitCalibrationPeak(am_dataset, "Am_59.5keV", interactive);
  if (am.valid)
    AddCalibrationPoint(cal_data, "Am-241 59.5keV", am.peaks.at(0).mu,
                        am.peaks.at(0).mu_error, E_AM241, am.reduced_chi2);

  // Ba-133 53.16 keV
  HyperEMGFitResult ba53 =
      FitCalibrationPeak(ba_dataset, "Ba_53.16keV", interactive);
  if (ba53.valid)
    AddCalibrationPoint(cal_data, "Ba-133 53.16keV", ba53.peaks.at(0).mu,
                        ba53.peaks.at(0).mu_error, E_BA133_53,
                        ba53.reduced_chi2);

  // Ba-133 80.98 keV
  HyperEMGFitResult ba81 =
      FitCalibrationPeak(ba_dataset, "Ba_80.98keV", interactive);
  if (ba81.valid)
    AddCalibrationPoint(cal_data, "Ba-133 80.98keV", ba81.peaks.at(0).mu,
                        ba81.peaks.at(0).mu_error, E_BA133_81,
                        ba81.reduced_chi2);

  // Ba-133 276.40 keV
  HyperEMGFitResult ba276 =
      FitCalibrationPeakFull(ba_dataset, "Ba_276.40keV", interactive);
  if (ba276.valid)
    AddCalibrationPoint(cal_data, "Ba-133 276.40keV", ba276.peaks.at(0).mu,
                        ba276.peaks.at(0).mu_error, E_BA133_276,
                        ba276.reduced_chi2);

  // Ba-133 302.85 keV
  HyperEMGFitResult ba303 =
      FitCalibrationPeakFull(ba_dataset, "Ba_302.85keV", interactive);
  if (ba303.valid)
    AddCalibrationPoint(cal_data, "Ba-133 302.85keV", ba303.peaks.at(0).mu,
                        ba303.peaks.at(0).mu_error, E_BA133_303,
                        ba303.reduced_chi2);

  // Ba-133 356.01 keV
  HyperEMGFitResult ba356 =
      FitCalibrationPeakFull(ba_dataset, "Ba_356.01keV", interactive);
  if (ba356.valid)
    AddCalibrationPoint(cal_data, "Ba-133 356.01keV", ba356.peaks.at(0).mu,
                        ba356.peaks.at(0).mu_error, E_BA133_356,
                        ba356.reduced_chi2);

  // Ba-133 383.85 keV
  HyperEMGFitResult ba384 =
      FitCalibrationPeakFull(ba_dataset, "Ba_383.85keV", interactive);
  if (ba384.valid)
    AddCalibrationPoint(cal_data, "Ba-133 383.85keV", ba384.peaks.at(0).mu,
                        ba384.peaks.at(0).mu_error, E_BA133_384,
                        ba384.reduced_chi2);

  // 511 keV annihilation
  HyperEMGFitResult ann =
      FitCalibrationPeakFull(sig_dataset, "Annihilation_511keV", interactive);
  if (ann.valid)
    AddCalibrationPoint(cal_data, "511 keV Annihilation", ann.peaks.at(0).mu,
                        ann.peaks.at(0).mu_error, E_ANNIHILATION,
                        ann.reduced_chi2);

  PrintCalibrationSummary(cal_data, "Master_20260115");

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, "Master_20260115");
  return cal_func;
}

// For other days: fit available peaks, find linear drift mapping
// E_master = alpha + beta * E_dayX, then compose with master calibration.
// Returns a composed TF1: master_cal(alpha + beta * x)
TF1 *BuildDriftCalibration(TF1 *master, const TString date_label,
                           const std::vector<TString> &peak_datasets,
                           const std::vector<TString> &peak_names,
                           const std::vector<Float_t> &true_energies,
                           const std::vector<Bool_t> &use_full,
                           const Bool_t interactive) {
  std::cout << "Building Drift Calibration for " << date_label << std::endl;

  // Fit peaks for this day
  std::vector<Float_t> measured;
  std::vector<Float_t> master_equiv; // where this peak falls in master coords
  std::vector<Float_t> measured_errors;

  for (Int_t p = 0; p < (Int_t)peak_names.size(); p++) {
    HyperEMGFitResult fit;
    if (use_full[p])
      fit =
          FitCalibrationPeakFull(peak_datasets[p], peak_names[p], interactive);
    else if (peak_names[p].Contains("Pb_KAlpha")) {
      fit = FitPbKAlpha(peak_datasets[p], interactive);
      if (fit.valid) {
        // Pb Ka double peak: add both
        for (Int_t k = 0; k < (Int_t)fit.peaks.size(); k++) {
          measured.push_back(fit.peaks[k].mu);
          measured_errors.push_back(fit.peaks[k].mu_error);
          master_equiv.push_back(true_energies[p + k]);
        }
        p++; // skip the second Pb entry
        continue;
      }
    } else
      fit = FitCalibrationPeak(peak_datasets[p], peak_names[p], interactive);

    if (fit.valid) {
      measured.push_back(fit.peaks.at(0).mu);
      measured_errors.push_back(fit.peaks.at(0).mu_error);
      master_equiv.push_back(true_energies[p]);
    }
  }

  if (measured.size() < 2) {
    std::cerr << "WARNING: Only " << measured.size()
              << " valid peaks for drift, need >= 2" << std::endl;
    return nullptr;
  }

  // For each measured peak, find what gain-matched energy in the master
  // would give the same true energy. Invert master: E_true = master(E_gm)
  // => E_gm_master = master^-1(E_true). Then the drift maps
  // E_measured_dayX -> E_gm_master via linear fit.
  //
  // For pol2: E_true = p0 + p1*x + p2*x^2, solve for x given E_true
  // using quadratic formula: x = (-p1 + sqrt(p1^2 - 4*p2*(p0 - E_true))) /
  // (2*p2)
  Double_t p0 = master->GetParameter(0);
  Double_t p1 = master->GetParameter(1);
  Double_t p2 = master->GetParameter(2);

  std::vector<Float_t> master_gm; // gain-matched energy in master coordinates
  for (Int_t i = 0; i < (Int_t)master_equiv.size(); i++) {
    Double_t E_true = master_equiv[i];
    Double_t e_gm = 0;
    if (TMath::Abs(p2) > 1e-12) {
      Double_t disc = p1 * p1 - 4.0 * p2 * (p0 - E_true);
      if (disc >= 0)
        e_gm = (-p1 + TMath::Sqrt(disc)) / (2.0 * p2);
      else
        e_gm = (E_true - p0) / p1; // fallback to linear
    } else {
      e_gm = (p1 != 0) ? (E_true - p0) / p1 : E_true;
    }
    master_gm.push_back(e_gm);
  }

  // Fit linear drift: E_gm_master = alpha + beta * E_measured_dayX
  Int_t n = measured.size();
  TGraph *drift_graph = new TGraph(n, measured.data(), master_gm.data());
  TF1 *drift_fit = new TF1("drift_" + date_label, "pol1", -10, 600);
  drift_fit->SetParameter(0, 0);
  drift_fit->SetParameter(1, 1);
  drift_graph->Fit(drift_fit, "RQ");

  Double_t alpha = drift_fit->GetParameter(0);
  Double_t beta = drift_fit->GetParameter(1);

  std::cout << "  Drift: alpha = " << std::fixed << std::setprecision(6)
            << alpha << ", beta = " << beta << " (" << n << " peaks)"
            << std::endl;

  // Print per-peak residuals
  for (Int_t i = 0; i < n; i++) {
    Double_t mapped = alpha + beta * measured[i];
    Double_t cal_true = master->Eval(mapped);
    std::cout << "    E_meas = " << std::setprecision(4) << measured[i]
              << ", E_true = " << master_equiv[i] << ", cal = " << cal_true
              << ", residual = " << (cal_true - master_equiv[i]) << " keV"
              << std::endl;
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(drift_graph, kBlue,
                                "; Measured [keV]; Master GM [keV]");
  drift_graph->SetMarkerStyle(5);
  drift_graph->SetMarkerSize(2);
  drift_graph->Draw("AP");
  drift_fit->Draw("SAME");
  PlottingUtils::SaveFigure(canvas, "drift_" + date_label, "",
                            PlotSaveOptions::kLINEAR);
  delete drift_graph;

  // Compose: E_true = master(alpha + beta * x) = p0 + p1*(a+b*x) +
  // p2*(a+b*x)^2
  TF1 *cal_func =
      new TF1("cal_" + date_label,
              Form("[0] + [1]*(%f + %f*x) + [2]*(%f + %f*x)*(%f + %f*x)", alpha,
                   beta, alpha, beta, alpha, beta),
              -10, 600);
  cal_func->SetParameter(0, p0);
  cal_func->SetParameter(1, p1);
  cal_func->SetParameter(2, p2);
  cal_func->SetNpx(1000);

  // Save
  TString cal_filepath =
      "root_files/calibration_function_" + date_label + ".root";
  TFile *cal_file = new TFile(cal_filepath, "RECREATE");
  cal_func->Write("calibration", TObject::kOverwrite);
  drift_fit->Write("drift", TObject::kOverwrite);
  cal_file->Close();
  delete cal_file;
  delete drift_fit;

  return cal_func;
}

// Invert calibration function at a given E_true value.
// For the composed drift+master function, use Newton's method.
Double_t InvertCalibration(TF1 *cal_func, Double_t e_true) {
  // Start with linear estimate
  Double_t x = e_true;
  for (Int_t iter = 0; iter < 20; iter++) {
    Double_t f = cal_func->Eval(x) - e_true;
    Double_t df = cal_func->Derivative(x);
    if (TMath::Abs(df) < 1e-15)
      break;
    Double_t dx = f / df;
    x -= dx;
    if (TMath::Abs(dx) < 1e-6)
      break;
  }
  return x;
}

// Remap a histogram's x-axis through a calibration function.
// For each output bin, numerically invert the calibration to find the source.
TH1F *RemapHistogram(TH1F *input, TF1 *cal_func, Int_t nbins, Float_t xmin,
                     Float_t xmax) {
  TH1F *output = new TH1F(
      PlottingUtils::GetRandomName(),
      Form("; Deposited Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
      nbins, xmin, xmax);
  output->SetDirectory(0);

  for (Int_t bin = 1; bin <= nbins; bin++) {
    Double_t e_cal = output->GetBinCenter(bin);
    Double_t e_gm = InvertCalibration(cal_func, e_cal);
    Int_t src_bin = input->FindBin(e_gm);
    if (src_bin >= 1 && src_bin <= input->GetNbinsX())
      output->SetBinContent(bin, input->GetBinContent(src_bin));
  }

  return output;
}

void PulseHeightToDepositedEnergy(const std::vector<TString> &input_names,
                                  TF1 *calibration_function,
                                  const TString date_label) {
  TString cal_filepath =
      "root_files/calibration_function_" + date_label + ".root";
  TFile *cal_file = new TFile(cal_filepath, "RECREATE");
  calibration_function->Write("calibration", TObject::kOverwrite);
  cal_file->Close();
  delete cal_file;

  Int_t n_files = input_names.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_name = input_names[i];
    TString filepath = "root_files/" + input_name + ".root";
    TFile *file = new TFile(filepath, "UPDATE");

    // Load SNIP-subtracted gain-matched histograms
    TH1F *snipHist = static_cast<TH1F *>(file->Get("snip_hist"));
    TH1F *snipZoomed = static_cast<TH1F *>(file->Get("snip_zoomedHist"));
    TH1F *snipPeak = static_cast<TH1F *>(file->Get("snip_peakHist"));

    if (!snipHist || !snipZoomed || !snipPeak) {
      std::cerr << "ERROR: Missing snip histograms in " << filepath
                << std::endl;
      file->Close();
      delete file;
      continue;
    }

    // Remap through calibration function
    TH1F *hist =
        RemapHistogram(snipHist, calibration_function, Constants::HIST_NBINS,
                       Constants::HIST_XMIN, Constants::HIST_XMAX);
    TH1F *zoomedHist = RemapHistogram(
        snipZoomed, calibration_function, Constants::ZOOMED_NBINS,
        Constants::ZOOMED_XMIN, Constants::ZOOMED_XMAX);
    TH1F *peakHist =
        RemapHistogram(snipPeak, calibration_function, Constants::PEAK_NBINS,
                       Constants::PEAK_XMIN, Constants::PEAK_XMAX);

    file->cd();
    hist->Write("calibrated_hist", TObject::kOverwrite);
    zoomedHist->Write("calibrated_zoomedHist", TObject::kOverwrite);
    peakHist->Write("calibrated_peakHist", TObject::kOverwrite);

    delete hist;
    delete zoomedHist;
    delete peakHist;

    std::cout << "Calibrated " << input_name << std::endl;
    file->Close();
    delete file;
  }
}

void Calibration() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t interactive = kTRUE;

  // Build master calibration from Jan 15 (pol2, all Ba + Am + 511)
  TF1 *master = BuildMasterCalibration(interactive);

  std::cout << std::endl;
  std::cout << "Master Calibration (pol2): p0 = " << std::fixed
            << std::setprecision(6) << master->GetParameter(0)
            << ", p1 = " << master->GetParameter(1)
            << ", p2 = " << std::scientific << std::setprecision(4)
            << master->GetParameter(2) << std::endl;
  std::cout << std::endl;

  // Jan 13: drift from master using Am-241, Pb Ka, 511
  std::vector<TString> ds13 = {Constants::POSTREACTOR_AM241_20260113,
                               Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                               Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                               Constants::CUSHIELDSIGNAL_10PERCENT_20260113};
  std::vector<TString> pn13 = {"Am_59.5keV", "Pb_KAlpha", "Pb_KAlpha_2",
                               "Annihilation_511keV"};
  std::vector<Float_t> te13 = {E_AM241, E_PB_KA1, E_PB_KA2, E_ANNIHILATION};
  std::vector<Bool_t> uf13 = {kFALSE, kFALSE, kFALSE, kTRUE};
  TF1 *cal_jan13 = BuildDriftCalibration(master, "20260113", ds13, pn13, te13,
                                         uf13, interactive);

  // Jan 14: drift from master using Pb Ka, 511
  std::vector<TString> ds14 = {Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                               Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                               Constants::CUSHIELDSIGNAL_10PERCENT_20260114};
  std::vector<TString> pn14 = {"Pb_KAlpha", "Pb_KAlpha_2",
                               "Annihilation_511keV"};
  std::vector<Float_t> te14 = {E_PB_KA1, E_PB_KA2, E_ANNIHILATION};
  std::vector<Bool_t> uf14 = {kFALSE, kFALSE, kTRUE};
  TF1 *cal_jan14 = BuildDriftCalibration(master, "20260114", ds14, pn14, te14,
                                         uf14, interactive);

  // Jan 15: master is already the calibration
  TF1 *cal_jan15 = master;

  // Jan 16: drift from master using Am-241, Ba-133 53/81, 511
  std::vector<TString> ds16 = {
      Constants::POSTREACTOR_AM241_BA133_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116};
  std::vector<TString> pn16 = {"Am_59.5keV", "Ba_53.16keV", "Ba_80.98keV",
                               "Annihilation_511keV"};
  std::vector<Float_t> te16 = {E_AM241, E_BA133_53, E_BA133_81, E_ANNIHILATION};
  std::vector<Bool_t> uf16 = {kFALSE, kFALSE, kFALSE, kTRUE};
  TF1 *cal_jan16 = BuildDriftCalibration(master, "20260116", ds16, pn16, te16,
                                         uf16, interactive);

  std::cout << std::endl;

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
  PulseHeightToDepositedEnergy(datasets_jan15, cal_jan15, "Master_20260115");

  std::vector<TString> datasets_jan16 = {
      Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116,
      Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116};
  PulseHeightToDepositedEnergy(datasets_jan16, cal_jan16, "20260116");
}
