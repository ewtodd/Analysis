#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TROOT.h>
#include <iomanip>

const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_81 = 80.9979;
const Float_t E_BA133_303 = 302.8508;
const Float_t E_BA133_356 = 356.0129;
const Float_t E_BA133_384 = 383.8485;
const Float_t E_ANNIHILATION = 510.999;

TH1F *LoadFullHistogram(const TString input_name, Int_t crystal) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TString histName = Form("hist_crystal%d", crystal);
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

TH1F *LoadZoomedHistogram(const TString input_name, Int_t crystal) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TString histName = Form("zoomedHist_crystal%d", crystal);
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

struct GainMatchResult {
  Float_t centroid;
  Float_t centroid_error;
  Float_t sigma;
  Float_t sigma_error;
  Float_t fwhm;
  Float_t resolution_pct;
  Float_t reduced_chi2;
  Bool_t valid;
};

GainMatchResult FitPeak(TH1F *hist, Float_t fit_low, Float_t fit_high,
                        Bool_t use_flat_background, const TString label,
                        const TString peak_name, Bool_t interactive) {
  GainMatchResult result;
  result.valid = kFALSE;
  if (!hist)
    return result;

  FittingUtils *fitter =
      new FittingUtils(hist, fit_low, fit_high, use_flat_background, kFALSE,
                       kTRUE, kTRUE, kTRUE);
  if (interactive)
    fitter->SetInteractive();
  FitResult fit = fitter->FitSinglePeak(label, peak_name);
  if (fit.valid) {
    result.valid = kTRUE;
    result.centroid = fit.peaks[0].mu;
    result.centroid_error = fit.peaks[0].mu_error;
    result.sigma = fit.peaks[0].sigma;
    result.sigma_error = fit.peaks[0].sigma_error;
    result.fwhm = 2.3548 * fit.peaks[0].sigma;
    result.resolution_pct = 100.0 * result.fwhm / fit.peaks[0].mu;
    result.reduced_chi2 = fit.reduced_chi2;
  }
  delete fitter;
  return result;
}

void PrintResults(const TString peak_name,
                  GainMatchResult results[Constants::N_CRYSTALS]) {
  std::cout << std::endl;
  std::cout << "  " << peak_name << std::endl;
  std::cout
      << "  Crystal | Centroid [keV]              | FWHM [keV] | Resolution "
         "[%] | chi2/ndf"
      << std::endl;
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    if (results[c].valid) {
      std::cout << "        " << c << " | " << std::fixed
                << std::setprecision(4) << results[c].centroid << " +/- "
                << std::setw(8) << results[c].centroid_error << " | "
                << std::setprecision(4) << std::setw(10) << results[c].fwhm
                << " | " << std::setprecision(2) << std::setw(13)
                << results[c].resolution_pct << "% | " << std::setprecision(3)
                << results[c].reduced_chi2 << std::endl;
    } else {
      std::cout << "        " << c << " | FAILED" << std::endl;
    }
  }
  std::cout << std::endl;
}

void GainMatch() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t interactive = kTRUE;
  TString am_dataset = Constants::POSTREACTOR_AM241_20260115;
  TString ba_dataset = Constants::POSTREACTOR_BA133_20260115;
  TString anni_dataset = Constants::NOSHIELDSIGNAL_5PERCENT_20260115;

  GainMatchResult am59_results[Constants::N_CRYSTALS];
  GainMatchResult ba81_results[Constants::N_CRYSTALS];
  GainMatchResult ba303_results[Constants::N_CRYSTALS];
  GainMatchResult ba356_results[Constants::N_CRYSTALS];
  GainMatchResult ba384_results[Constants::N_CRYSTALS];
  GainMatchResult ann511_results[Constants::N_CRYSTALS];

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TString crystal_suffix = Form("_crystal%d", c);

    // Ba-133 81 keV
    TH1F *ba_hist = LoadZoomedHistogram(ba_dataset, c);
    ba81_results[c] =
        FitPeak(ba_hist, 78, 84, kTRUE, ba_dataset + crystal_suffix,
                "Ba133_81keV", interactive);
    delete ba_hist;

    // Am-241 59.5 keV
    TH1F *am_hist = LoadZoomedHistogram(am_dataset, c);
    am59_results[c] =
        FitPeak(am_hist, 45, 65, kTRUE, am_dataset + crystal_suffix,
                "Am241_59keV", interactive);
    delete am_hist;

    // Ba-133 higher-energy lines (from full histogram)
    TH1F *ba_full = LoadFullHistogram(ba_dataset, c);

    ba303_results[c] =
        FitPeak(ba_full, 296, 310, kFALSE, ba_dataset + crystal_suffix,
                "Ba133_303keV", interactive);

    ba356_results[c] =
        FitPeak(ba_full, 349, 363, kFALSE, ba_dataset + crystal_suffix,
                "Ba133_356keV", interactive);

    ba384_results[c] =
        FitPeak(ba_full, 377, 391, kFALSE, ba_dataset + crystal_suffix,
                "Ba133_384keV", interactive);

    delete ba_full;

    // 511 keV annihilation
    TH1F *bkg_full = LoadFullHistogram(anni_dataset, c);
    ann511_results[c] =
        FitPeak(bkg_full, 500, 525, kFALSE, anni_dataset + crystal_suffix,
                "Annihilation_511keV", interactive);
    delete bkg_full;
  }

  PrintResults("Am-241 59.5 keV", am59_results);
  PrintResults("Ba-133 81 keV", ba81_results);
  PrintResults("Ba-133 303 keV", ba303_results);
  PrintResults("Ba-133 356 keV", ba356_results);
  PrintResults("Ba-133 384 keV", ba384_results);
  PrintResults("511 keV Annihilation", ann511_results);

  // Build per-crystal gain corrections
  TString output_path = "root_files/gain_match.root";
  TFile *outfile = new TFile(output_path, "RECREATE");

  std::cout << std::endl;
  std::cout << "  Per-Crystal Gain Factors:" << std::endl;

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    std::vector<Float_t> measured;
    std::vector<Float_t> true_energies;

    if (am59_results[c].valid) {
      measured.push_back(am59_results[c].centroid);
      true_energies.push_back(E_AM241);
    }
    if (ba81_results[c].valid) {
      measured.push_back(ba81_results[c].centroid);
      true_energies.push_back(E_BA133_81);
    }
    if (ba303_results[c].valid) {
      measured.push_back(ba303_results[c].centroid);
      true_energies.push_back(E_BA133_303);
    }
    if (ba356_results[c].valid) {
      measured.push_back(ba356_results[c].centroid);
      true_energies.push_back(E_BA133_356);
    }
    if (ba384_results[c].valid) {
      measured.push_back(ba384_results[c].centroid);
      true_energies.push_back(E_BA133_384);
    }
    if (ann511_results[c].valid) {
      measured.push_back(ann511_results[c].centroid);
      true_energies.push_back(E_ANNIHILATION);
    }

    Int_t n_points = measured.size();
    TGraph *graph = new TGraph(n_points, measured.data(), true_energies.data());
    TString func_name = Form("gain_crystal%d", c);
    TF1 *fit = new TF1(func_name, "pol2", -10, 600);
    fit->SetParameter(0, 0);
    fit->SetParameter(1, 1);
    fit->SetParameter(2, 0);
    fit->SetNpx(1000);
    graph->Fit(fit, "RQ");

    std::cout << "    Crystal " << c << ": p0 = " << std::fixed
              << std::setprecision(6) << fit->GetParameter(0) << " +/- "
              << fit->GetParError(0) << ", p1 = " << fit->GetParameter(1)
              << " +/- " << fit->GetParError(1) << ", p2 = " << std::scientific
              << std::setprecision(4) << fit->GetParameter(2) << " +/- "
              << fit->GetParError(2) << " (" << n_points << " peaks)"
              << std::endl;

    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
    PlottingUtils::ConfigureGraph(graph, kBlue,
                                  "; Measured Energy [keV]; True Energy [keV]");
    graph->SetMarkerStyle(5);
    graph->SetMarkerSize(2);
    graph->Draw("AP");
    fit->Draw("SAME");
    PlottingUtils::SaveFigure(canvas, "gain_match_" + func_name, "",
                              PlotSaveOptions::kLINEAR);

    fit->Write(func_name);
    delete graph;
  }

  outfile->Close();
  delete outfile;

  std::cout << std::endl;
  std::cout << "  Gain match factors saved to " << output_path << std::endl;

  // Per-crystal FWHM calibration: FWHM(E) = sqrt(c + b*E + a*E^2)
  TString fwhm_output_path = "root_files/fwhm_calibration.root";
  TFile *fwhmFile = new TFile(fwhm_output_path, "RECREATE");

  std::cout << std::endl;
  std::cout << "  Per-Crystal FWHM Calibration:" << std::endl;

  const Int_t N_FWHM_PEAKS = 6;
  GainMatchResult *all_results[N_FWHM_PEAKS] = {
      am59_results,  ba81_results,  ba303_results,
      ba356_results, ba384_results, ann511_results};
  Float_t all_energies[N_FWHM_PEAKS] = {E_AM241,     E_BA133_81,  E_BA133_303,
                                         E_BA133_356, E_BA133_384, E_ANNIHILATION};

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    std::vector<Double_t> energies;
    std::vector<Double_t> fwhms;
    std::vector<Double_t> energy_errors;
    std::vector<Double_t> fwhm_errors;

    for (Int_t p = 0; p < N_FWHM_PEAKS; p++) {
      if (all_results[p][c].valid) {
        energies.push_back(all_energies[p]);
        fwhms.push_back(all_results[p][c].fwhm);
        energy_errors.push_back(0.0);
        fwhm_errors.push_back(2.3548 * all_results[p][c].sigma_error);
      }
    }

    Int_t n_fwhm = energies.size();
    if (n_fwhm < 3) {
      std::cerr << "WARNING: Crystal " << c << " has only " << n_fwhm
                << " valid FWHM points, skipping" << std::endl;
      continue;
    }

    TGraphErrors *fwhm_graph =
        new TGraphErrors(n_fwhm, energies.data(), fwhms.data(),
                         energy_errors.data(), fwhm_errors.data());

    TString fwhm_func_name = Form("fwhm_calibration_crystal%d", c);
    TF1 *fwhm_fit =
        new TF1(fwhm_func_name, "sqrt([0] + [1]*x + [2]*x*x)", 0, 600);
    fwhm_fit->SetParNames("c", "b", "a");
    fwhm_fit->SetParameter(0, 1.0);
    fwhm_fit->SetParameter(1, 0.01);
    fwhm_fit->SetParameter(2, 1e-6);
    fwhm_fit->SetNpx(1000);

    fwhm_graph->Fit(fwhm_fit, "R");

    std::cout << "    Crystal " << c << ": c = " << std::scientific
              << std::setprecision(4) << fwhm_fit->GetParameter(0)
              << ", b = " << fwhm_fit->GetParameter(1)
              << ", a = " << fwhm_fit->GetParameter(2) << " (" << n_fwhm
              << " points)" << std::endl;

    TCanvas *fwhm_canvas = PlottingUtils::GetConfiguredCanvas();
    PlottingUtils::ConfigureGraph(fwhm_graph, kBlue,
                                  "; Energy [keV]; FWHM [keV]");
    fwhm_graph->SetMarkerStyle(5);
    fwhm_graph->SetMarkerSize(2);
    fwhm_graph->Draw("AP");
    fwhm_fit->Draw("SAME");
    PlottingUtils::SaveFigure(fwhm_canvas, "fwhm_" + fwhm_func_name, "",
                              PlotSaveOptions::kLINEAR);

    fwhm_fit->Write(fwhm_func_name);
    fwhm_graph->Write(Form("fwhm_graph_crystal%d", c));
    delete fwhm_graph;
  }

  fwhmFile->Close();
  delete fwhmFile;

  std::cout << std::endl;
  std::cout << "  FWHM calibration saved to " << fwhm_output_path << std::endl;
}
