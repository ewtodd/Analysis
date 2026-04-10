#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TGraph.h>
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

void LiveTimeSums(std::vector<TString> filenames) {
  Float_t const MS_TO_S = 1e-3;

  Int_t n_files = Int_t(filenames.size());
  for (Int_t j = 0; j < n_files; j++) {
    TString filename = filenames[j];
    TString filepath = "root_files/" + filename + ".root";
    TFile *file = new TFile(filepath, "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Could not open " << filepath << std::endl;
      continue;
    }

    TParameter<Float_t> *param =
        (TParameter<Float_t> *)file->Get("N42_RealTime_Total");
    Float_t daqLiveTime_s = param ? param->GetVal() : -1;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "--- " << filename << " ---" << std::endl;
    std::cout << "  DAQ live time: " << daqLiveTime_s << " s" << std::endl;

    TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));
    if (!tree) {
      std::cerr << filename << ": missing bef_tree" << std::endl;
      file->Close();
      continue;
    }

    Float_t x = 0, y = 0;
    UInt_t eventTime = 0;
    Int_t interaction = 0;
    tree->SetBranchAddress("xmm", &x);
    tree->SetBranchAddress("ymm", &y);
    tree->SetBranchAddress("eventTime", &eventTime);
    tree->SetBranchAddress("interaction", &interaction);

    Int_t n_entries = tree->GetEntries();

    UInt_t firstEventTime[Constants::N_CRYSTALS];
    UInt_t lastEventTime[Constants::N_CRYSTALS];
    Int_t nPhysicalEvents[Constants::N_CRYSTALS];
    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      firstEventTime[c] = 0;
      lastEventTime[c] = 0;
      nPhysicalEvents[c] = 0;
    }

    for (Int_t i = 0; i < n_entries; i++) {
      tree->GetEntry(i);
      if (interaction != 0)
        continue;
      Int_t c = GetCrystalIndex(x, y);
      if (c < 0)
        continue;
      if (nPhysicalEvents[c] == 0)
        firstEventTime[c] = eventTime;
      lastEventTime[c] = eventTime;
      nPhysicalEvents[c]++;
    }

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      TString paramName = Form("LiveTime_Unfiltered_Crystal%d_s", c);
      TParameter<Float_t> *ltParam =
          (TParameter<Float_t> *)file->Get(paramName);
      Float_t sumLiveTime_s = ltParam ? ltParam->GetVal() : -1;

      Float_t eventSpan_s =
          (Float_t)(lastEventTime[c] - firstEventTime[c]) * MS_TO_S;

      std::cout << "  Crystal " << c << ":"
                << " physEvents=" << nPhysicalEvents[c]
                << " sumLT=" << sumLiveTime_s << " s"
                << " span=" << eventSpan_s << " s"
                << " sum/span="
                << (eventSpan_s > 0 ? sumLiveTime_s / eventSpan_s : 0)
                << " sum/FAQ="
                << (daqLiveTime_s > 0 ? sumLiveTime_s / daqLiveTime_s : 0)
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
