#include "Constants.hpp"
#include "InitUtils.hpp"
#include <TROOT.h>
#include <TSystem.h>
#include <algorithm>
#include <vector>

void CalculateChargeComparison(const std::vector<TString> output_names) {
  Int_t n_files = output_names.size();

  for (Int_t entry = 0; entry < n_files; entry++) {
    TString output_name = output_names.at(entry);
    TString filepath = "root_files/" + output_name + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree = static_cast<TTree *>(file->Get("features"));

    TArrayF *samples = nullptr;
    Float_t light_output_keVee;
    Int_t trigger_position;

    tree->SetBranchAddress("Samples", &samples);
    tree->SetBranchAddress("light_output", &light_output_keVee);
    tree->SetBranchAddress("trigger_position", &trigger_position);

    Float_t charge_comparison;
    tree->Branch("charge_comparison", &charge_comparison,
                 "charge_comparison/F");

    Int_t n_entries = tree->GetEntries();

    Float_t short_calc, long_calc;

    for (Int_t i = 0; i < n_entries; i++) {
      short_calc = 0;
      long_calc = 0;
      tree->GetEntry(i);

      Int_t n_samples = samples->GetSize();
      Int_t start = std::max<Int_t>(
          trigger_position - Constants::DEFAULT_PROCESSING_CONFIG.pre_gate, 0);
      Int_t short_end =
          std::min(start + Constants::OPTIMAL_SHORT_GATE, n_samples);
      Int_t long_end =
          std::min(start + Constants::OPTIMAL_LONG_GATE, n_samples);

      for (Int_t j = 0; j < n_samples; j++) {
        Float_t value = samples->GetAt(j);
        if (j < short_end)
          short_calc += value;
        long_calc += value;
      }

      if (long_calc > short_calc) {
        charge_comparison = 1 - (short_calc / long_calc);
      } else {
        charge_comparison = -1;
      }

      tree->GetBranch("charge_comparison")->Fill();
    }

    tree->Write("features", TObject::kOverwrite);
    file->Close();
    delete file;
    std::cout << "Calculated charge comparison for " << output_name
              << std::endl;
  }
}

void ChargeComparison() {
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT);
  CalculateChargeComparison(Constants::ALL_OUTPUT_NAMES);
}
