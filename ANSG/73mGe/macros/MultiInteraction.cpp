#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TTree.h>

const Int_t WIDE_XMAX = 6000;
const Int_t WIDE_NBINS =
    (WIDE_XMAX - Constants::HIST_XMIN) / Constants::BIN_WIDTH_KEV;

void MultiInteraction() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  std::vector<TString> filenames = Constants::ALL_DATASETS;

  for (Int_t f = 0; f < (Int_t)filenames.size(); f++) {
    TString filename = filenames[f];
    TString filepath = "root_files/" + filename + ".root";
    TFile *file = new TFile(filepath, "READ");
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

    Double_t totalEnergy = 0;
    Int_t nInteractions = 0;
    Int_t interaction = 0;

    tree->SetBranchAddress("totalEnergykeV", &totalEnergy);
    tree->SetBranchAddress("nInteractions", &nInteractions);
    tree->SetBranchAddress("interaction", &interaction);

    TH1D *hSummedWide = new TH1D(
        Form("multiint_summed_wide_%s", filename.Data()),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        WIDE_NBINS, Constants::HIST_XMIN, WIDE_XMAX);
    TH1D *hSummedZoomed = new TH1D(
        Form("multiint_summed_zoomed_%s", filename.Data()),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
        Constants::ZOOMED_XMAX);
    TH1D *hSummedPeak = new TH1D(
        Form("multiint_summed_peak_%s", filename.Data()),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        Constants::PEAK_NBINS, Constants::PEAK_XMIN, Constants::PEAK_XMAX);

    hSummedWide->SetDirectory(0);
    hSummedZoomed->SetDirectory(0);
    hSummedPeak->SetDirectory(0);

    Int_t n_entries = tree->GetEntries();
    Int_t nMultiEvents = 0;

    for (Int_t i = 0; i < n_entries; i++) {
      tree->GetEntry(i);

      if (nInteractions == 1)
        continue;
      if (interaction != 0)
        continue;

      hSummedWide->Fill(totalEnergy);
      if (totalEnergy > Constants::ZOOMED_XMIN &&
          totalEnergy < Constants::ZOOMED_XMAX)
        hSummedZoomed->Fill(totalEnergy);
      if (totalEnergy > Constants::PEAK_XMIN &&
          totalEnergy < Constants::PEAK_XMAX)
        hSummedPeak->Fill(totalEnergy);
      nMultiEvents++;
    }

    std::cout << filename << ": " << nMultiEvents << " multi-interaction events"
              << std::endl;

    file->Close();
    delete file;

    // Plot wide histograms
    TCanvas *cSummed = PlottingUtils::GetConfiguredCanvas();
    PlottingUtils::ConfigureHistogram(hSummedWide, kP10Violet);
    hSummedWide->SetFillStyle(0);
    hSummedWide->SetLineWidth(2);
    hSummedWide->Draw("HIST");
    PlottingUtils::SaveFigure(cSummed, filename + "_multiint_summed",
                              "multiInteraction", PlotSaveOptions::kLOG);

    // Save histograms to existing file
    TFile *outFile = new TFile(filepath, "UPDATE");
    hSummedWide->Write("multiint_summed_wide", TObject::kOverwrite);
    hSummedZoomed->Write("multiint_summed_zoomed", TObject::kOverwrite);
    hSummedPeak->Write("multiint_summed_peak", TObject::kOverwrite);
    outFile->Close();
    delete outFile;

    delete hSummedWide;
    delete hSummedZoomed;
    delete hSummedPeak;
  }
}
