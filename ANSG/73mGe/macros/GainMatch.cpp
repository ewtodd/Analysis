#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "HyperEMGFitHelpers.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>
#include <iomanip>

Int_t GetCrystalIndex(Float_t x, Float_t y) {
  if (x < 0 && y < 0)
    return 0;
  if (x > 0 && y < 0)
    return 1;
  if (x < 0 && y > 0)
    return 2;
  if (x > 0 && y > 0)
    return 3;
  return -1;
}

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

  TString histName = Form("hist_crystal%d_smoothed", crystal);
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

  TString histName = Form("zoomedHist_crystal%d_smoothed", crystal);
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

  HyperEMGFitResult fit = FitSingleHyperEMG(hist, fit_low, fit_high,
                                             use_flat_background, label,
                                             peak_name, interactive);
  if (fit.valid) {
    result.valid = kTRUE;
    result.centroid = fit.peaks[0].mu;
    result.centroid_error = fit.peaks[0].mu_error;
    result.sigma = fit.peaks[0].sigma;
    result.sigma_error = fit.peaks[0].sigma_error;
    result.reduced_chi2 = fit.reduced_chi2;
  }
  return result;
}

void PrintResults(const TString peak_name,
                  GainMatchResult results[Constants::N_CRYSTALS]) {
  std::cout << std::endl;
  std::cout << "  " << peak_name << std::endl;
  std::cout << "  Crystal | Centroid [keV]              | chi2/ndf"
            << std::endl;
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    if (results[c].valid) {
      std::cout << "        " << c << " | " << std::fixed
                << std::setprecision(4) << results[c].centroid << " +/- "
                << std::setw(8) << results[c].centroid_error << " | "
                << std::setprecision(3) << results[c].reduced_chi2 << std::endl;
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
        FitPeak(ba_hist, 78, 84, kFALSE, ba_dataset + crystal_suffix,
                "Ba133_81keV", interactive);
    delete ba_hist;

    // Am-241 59.5 keV
    TH1F *am_hist = LoadZoomedHistogram(am_dataset, c);
    am59_results[c] =
        FitPeak(am_hist, 45, 65, kFALSE, am_dataset + crystal_suffix,
                "Am241_59keV", interactive);
    delete am_hist;

    // Ba-133 higher-energy lines (from full histogram, fresh copy each)
    TH1F *ba_full_303 = LoadFullHistogram(ba_dataset, c);
    ba303_results[c] =
        FitPeak(ba_full_303, 296, 310, kFALSE, ba_dataset + crystal_suffix,
                "Ba133_303keV", interactive);
    delete ba_full_303;

    TH1F *ba_full_356 = LoadFullHistogram(ba_dataset, c);
    ba356_results[c] =
        FitPeak(ba_full_356, 349, 363, kFALSE, ba_dataset + crystal_suffix,
                "Ba133_356keV", interactive);
    delete ba_full_356;

    TH1F *ba_full_384 = LoadFullHistogram(ba_dataset, c);
    ba384_results[c] =
        FitPeak(ba_full_384, 377, 391, kFALSE, ba_dataset + crystal_suffix,
                "Ba133_384keV", interactive);
    delete ba_full_384;

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

  // ---- Phase 2: Build gain-matched totalEnergy trees for all datasets ----
  std::cout << std::endl;
  std::cout << "  Building gain-matched totalEnergy trees..." << std::endl;

  TF1 *gain_functions[Constants::N_CRYSTALS];
  {
    TFile *gf = new TFile(output_path, "READ");
    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      TString name = Form("gain_crystal%d", c);
      TF1 *func = static_cast<TF1 *>(gf->Get(name));
      gain_functions[c] = static_cast<TF1 *>(func->Clone());
    }
    gf->Close();
    delete gf;
  }

  for (Int_t d = 0; d < (Int_t)Constants::ALL_DATASETS.size(); d++) {
    TString dataset = Constants::ALL_DATASETS[d];
    TString filepath = "root_files/" + dataset + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    if (!file || file->IsZombie()) {
      std::cerr << "ERROR: Cannot open " << filepath << std::endl;
      continue;
    }

    TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));
    if (!tree) {
      std::cerr << "ERROR: No bef_tree in " << filepath << std::endl;
      file->Close();
      delete file;
      continue;
    }

    Float_t energy = 0, x = 0, y = 0;
    Int_t nInteractions = 0, interaction = 0;
    Int_t liveTime = 0;
    UInt_t eventTime = 0;

    tree->SetBranchAddress("energykeV", &energy);
    tree->SetBranchAddress("xmm", &x);
    tree->SetBranchAddress("ymm", &y);
    tree->SetBranchAddress("nInteractions", &nInteractions);
    tree->SetBranchAddress("interaction", &interaction);
    tree->SetBranchAddress("liveTime", &liveTime);
    tree->SetBranchAddress("eventTime", &eventTime);

    Float_t totalEnergyGM = 0;
    UInt_t outEventTime = 0;
    Int_t outLiveTime = 0;

    TTree *gmTree =
        new TTree("gain_matched_tree", "Gain-matched total energy per event");
    gmTree->Branch("totalEnergykeV", &totalEnergyGM, "totalEnergykeV/F");
    gmTree->Branch("eventTime", &outEventTime, "eventTime/i");
    gmTree->Branch("liveTime", &outLiveTime, "liveTime/I");

    Int_t n_entries = tree->GetEntries();
    Float_t accumE = 0;

    for (Int_t i = 0; i < n_entries; i++) {
      tree->GetEntry(i);

      Int_t crystal = GetCrystalIndex(x, y);
      Float_t gainMatched =
          (crystal >= 0) ? gain_functions[crystal]->Eval(energy) : energy;

      if (interaction == 0) {
        accumE = gainMatched;
        outEventTime = eventTime;
        outLiveTime = liveTime;
      } else {
        accumE += gainMatched;
      }

      if (interaction == nInteractions - 1) {
        totalEnergyGM = accumE;
        gmTree->Fill();
      }
    }

    gmTree->Write("gain_matched_tree", TObject::kOverwrite);
    std::cout << "    " << dataset << ": " << gmTree->GetEntries()
              << " gain-matched events" << std::endl;

    file->Close();
    delete file;
  }

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++)
    delete gain_functions[c];

  std::cout << std::endl;
  std::cout << "  Gain matching complete." << std::endl;
}
