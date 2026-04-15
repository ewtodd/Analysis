#include "Constants.hpp"
#include "InitUtils.hpp"
#include "InteractiveSNIPEditor.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TTree.h>
#include <iomanip>

const Int_t MARKOV_WINDOW = 3;
const Double_t SNIP_M_MAX = 30.0;
const TString SNIP_PARAM_FILE = "plots/snipBackground/snip_parameters.snip";

void SmoothMarkov(TH1F *hist) {
  TSpectrum spec;
  Int_t nbins = hist->GetNbinsX();
  Double_t *buf = new Double_t[nbins];
  for (Int_t i = 0; i < nbins; i++)
    buf[i] = hist->GetBinContent(i + 1);
  spec.SmoothMarkov(buf, nbins, MARKOV_WINDOW);
  for (Int_t i = 0; i < nbins; i++)
    hist->SetBinContent(i + 1, buf[i]);
  delete[] buf;
}

TH1F *BuildHistFromTree(TFile *file, Int_t nbins, Float_t xmin, Float_t xmax) {
  TTree *tree = static_cast<TTree *>(file->Get("gain_matched_tree"));
  if (!tree)
    return nullptr;

  TH1F *hist =
      new TH1F(PlottingUtils::GetRandomName(),
               Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
               nbins, xmin, xmax);
  hist->SetDirectory(0);

  Float_t totalEnergy = 0;
  tree->SetBranchAddress("totalEnergykeV", &totalEnergy);
  Int_t n = tree->GetEntries();
  for (Int_t i = 0; i < n; i++) {
    tree->GetEntry(i);
    hist->Fill(totalEnergy);
  }

  return hist;
}

void PlotSNIPWindow(const Double_t *m) {
  const Int_t n_points = 300;
  Double_t energies[n_points];
  Double_t windows[n_points];

  for (Int_t i = 0; i < n_points; i++) {
    Double_t energy =
        Constants::HIST_XMIN +
        (Constants::HIST_XMAX - Constants::HIST_XMIN) * i / (n_points - 1.0);
    energies[i] = energy;
    windows[i] = EvalSNIPWindow(energy, m);
  }

  TGraph *graph = new TGraph(n_points, energies, windows);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureGraph(graph, kBlue,
                                "; Energy [keV]; SNIP Window [keV]");
  graph->SetLineWidth(2);
  graph->Draw("AL");

  Double_t ctrl_e[N_SNIP_REGIONS];
  Double_t ctrl_m[N_SNIP_REGIONS];
  for (Int_t r = 0; r < N_SNIP_REGIONS; r++) {
    ctrl_e[r] = (r + 0.5) * SNIP_REGION_WIDTH;
    ctrl_m[r] = m[r];
  }
  TGraph *points = new TGraph(N_SNIP_REGIONS, ctrl_e, ctrl_m);
  points->SetMarkerStyle(20);
  points->SetMarkerSize(1.5);
  points->SetMarkerColor(kRed);
  points->Draw("P SAME");

  PlottingUtils::SaveFigure(canvas, "snip_window_vs_energy", "snipBackground",
                            PlotSaveOptions::kLOG);
  delete graph;
  delete points;
}

void SNIPBackground() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t interactive = kTRUE;

  Double_t m[N_SNIP_REGIONS] = {2.0, 3.0, 5.0, 7.0, 10.0};
  Bool_t skip_interactive = kFALSE;

  if (LoadSNIPParameters(SNIP_PARAM_FILE, m)) {
    std::cout << "Loaded SNIP parameters from " << SNIP_PARAM_FILE << ":";
    for (Int_t r = 0; r < N_SNIP_REGIONS; r++)
      std::cout << " " << m[r];
    std::cout << " keV" << std::endl;
  }

  std::vector<TString> datasets = Constants::ALL_DATASETS;

  for (Int_t d = 0; d < (Int_t)datasets.size(); d++) {
    TString dataset = datasets[d];
    std::cout << "SNIP Background: " << dataset << std::endl;

    TString filepath = "root_files/" + dataset + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    if (!file || file->IsZombie()) {
      std::cerr << "ERROR: Cannot open " << filepath << std::endl;
      continue;
    }

    TH1F *hist = BuildHistFromTree(file, Constants::HIST_NBINS,
                                   Constants::HIST_XMIN, Constants::HIST_XMAX);
    TH1F *zoomedHist =
        BuildHistFromTree(file, Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                          Constants::ZOOMED_XMAX);
    TH1F *peakHist =
        BuildHistFromTree(file, Constants::PEAK_NBINS, Constants::PEAK_XMIN,
                          Constants::PEAK_XMAX);

    if (!hist || !zoomedHist || !peakHist) {
      std::cerr << "ERROR: Cannot build histograms for " << dataset
                << std::endl;
      file->Close();
      delete file;
      continue;
    }

    if (interactive && !skip_interactive) {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);

      SNIPParameters params = LaunchSNIPEditor(hist, m, SNIP_M_MAX, dataset);

      gROOT->SetBatch(was_batch);

      if (params.accepted) {
        for (Int_t r = 0; r < N_SNIP_REGIONS; r++)
          m[r] = params.m[r];

        SaveSNIPParameters(SNIP_PARAM_FILE, m);

        if (params.use_for_all)
          skip_interactive = kTRUE;
      }
    }

    std::cout << "  SNIP windows:";
    for (Int_t r = 0; r < N_SNIP_REGIONS; r++)
      std::cout << " " << std::fixed << std::setprecision(1) << m[r];
    std::cout << " keV" << std::endl;

    TH1F *background = ComputeAdaptiveSNIPBackground(hist, m);
    hist->Add(background, -1.0);
    delete background;

    TH1F *fullForSub =
        BuildHistFromTree(file, Constants::HIST_NBINS, Constants::HIST_XMIN,
                          Constants::HIST_XMAX);
    TH1F *bgFull = ComputeAdaptiveSNIPBackground(fullForSub, m);

    for (Int_t bin = 1; bin <= zoomedHist->GetNbinsX(); bin++) {
      Double_t x = zoomedHist->GetBinCenter(bin);
      Int_t bg_bin = bgFull->FindBin(x);
      zoomedHist->SetBinContent(bin, zoomedHist->GetBinContent(bin) -
                                         bgFull->GetBinContent(bg_bin));
    }

    for (Int_t bin = 1; bin <= peakHist->GetNbinsX(); bin++) {
      Double_t x = peakHist->GetBinCenter(bin);
      Int_t bg_bin = bgFull->FindBin(x);
      peakHist->SetBinContent(bin, peakHist->GetBinContent(bin) -
                                       bgFull->GetBinContent(bg_bin));
    }

    delete fullForSub;
    delete bgFull;

    SmoothMarkov(hist);
    SmoothMarkov(zoomedHist);
    SmoothMarkov(peakHist);

    hist->Write("snip_hist", TObject::kOverwrite);
    zoomedHist->Write("snip_zoomedHist", TObject::kOverwrite);
    peakHist->Write("snip_peakHist", TObject::kOverwrite);

    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
    PlottingUtils::ConfigureHistogram(zoomedHist, kP10Violet);
    zoomedHist->Draw("HIST");
    PlottingUtils::SaveFigure(canvas, dataset + "_snip_zoomed",
                              "snipBackground", PlotSaveOptions::kLOG);

    delete hist;
    delete zoomedHist;
    delete peakHist;

    file->Close();
    delete file;

    std::cout << "  Saved SNIP-subtracted histograms for " << dataset
              << std::endl;
  }

  PlotSNIPWindow(m);
  std::cout << "SNIP background subtraction complete." << std::endl;
}
