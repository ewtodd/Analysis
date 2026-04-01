#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>
#include <iomanip>

static const Int_t N_CRYSTALS = 4;


void LiveTimeSums(std::vector<TString> filenames) {
  Double_t const TENS_OF_NS_TO_S = 1e-8;
  Double_t const MS_TO_S = 1e-3;

  Int_t n_files = Int_t(filenames.size());
  for (Int_t j = 0; j < n_files; j++) {
    TString filename = filenames[j];
    TString filepath = "root_files/" + filename + ".root";
    TFile *file = new TFile(filepath, "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Could not open " << filepath << std::endl;
      continue;
    }

    TParameter<Double_t> *param =
        (TParameter<Double_t> *)file->Get("N42_RealTime_Total");
    Double_t daqLiveTime_s = param ? param->GetVal() : -1;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "--- " << filename << " ---" << std::endl;
    std::cout << "  DAQ live time: " << daqLiveTime_s << " s" << std::endl;

    for (Int_t c = 0; c < N_CRYSTALS; c++) {
      TString treeName = Form("crystal%d_unfiltered_tree", c);
      TTree *tree = static_cast<TTree *>(file->Get(treeName));
      if (!tree) {
        std::cerr << filename << ": missing " << treeName << std::endl;
        continue;
      }

      Int_t liveTime = 0;
      UInt_t eventTime = 0;
      Int_t interaction = 0;
      tree->SetBranchAddress("liveTime", &liveTime);
      tree->SetBranchAddress("eventTime", &eventTime);
      tree->SetBranchAddress("interaction", &interaction);

      Int_t n_entries = tree->GetEntries();
      if (n_entries == 0)
        continue;

      Double_t sumLiveTime_s = 0;
      UInt_t firstEventTime = 0;
      UInt_t lastEventTime = 0;
      Int_t nPhysicalEvents = 0;

      for (Int_t i = 0; i < n_entries; i++) {
        tree->GetEntry(i);
        if (interaction != 0)
          continue;
        sumLiveTime_s += liveTime * TENS_OF_NS_TO_S;
        if (nPhysicalEvents == 0)
          firstEventTime = eventTime;
        lastEventTime = eventTime;
        nPhysicalEvents++;
      }

      Double_t eventSpan_s = (Double_t)(lastEventTime - firstEventTime) * MS_TO_S;

      std::cout << "  Crystal " << c << ":"
                << " physEvents=" << nPhysicalEvents
                << " (entries=" << n_entries << ")"
                << " sumLT=" << sumLiveTime_s << " s"
                << " span=" << eventSpan_s << " s"
                << " sum/span=" << (eventSpan_s > 0 ? sumLiveTime_s / eventSpan_s : 0)
                << " sum/DAQ=" << (daqLiveTime_s > 0 ? sumLiveTime_s / daqLiveTime_s : 0)
                << std::endl;
    }

    std::cout << std::endl;

    file->Close();
  }
}

void LiveTimeDiagnostic() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  std::vector<TString> filenames;
  filenames.push_back(Constants::CDSHIELDSIGNAL_10PERCENT_20260113);
  filenames.push_back(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113);
  filenames.push_back(Constants::CDSHIELDSIGNAL_25PERCENT_20260113);
  filenames.push_back(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113);

  LiveTimeSums(filenames);
}
