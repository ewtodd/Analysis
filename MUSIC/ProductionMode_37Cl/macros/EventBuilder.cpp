#include "Constants.hpp"
#include "InitUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <ostream>
#include <vector>

void ResetEvent(Int_t leftdE[16], Int_t rightdE[16], Int_t totaldE[16],
                ULong64_t timestamp_distribution[36], UInt_t all_flags[36],
                Int_t hits[36], Int_t &cathode, Int_t &strip0, Int_t &strip17,
                Int_t &grid) {

  for (Int_t k = 0; k < 16; k++) {
    leftdE[k] = rightdE[k] = totaldE[k] = 0;
  }
  for (Int_t k = 0; k < 36; k++) {
    timestamp_distribution[k] = 0;
    all_flags[k] = 0;
    hits[k] = 0;
  }
  cathode = strip0 = strip17 = grid = 0;
}

Int_t GetStripNumber(TString map_name) {
  TString num_str = map_name;
  num_str.Remove(0, 1);
  return num_str.Atoi();
}

void BuildEvents(std::vector<TString> input_filenames,
                 std::vector<TString> output_names, Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  const ULong64_t TIME_WINDOW = Constants::COINCIDENCE_WINDOW;
  std::cout << "Using coincidence window of " << TIME_WINDOW * 1e-6
            << " microseconds..." << std::endl;

  Int_t n_files = input_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filename = input_filenames[i];
    TString input_filepath = "root_files/" + input_filename + ".root";
    TFile *input_file = new TFile(input_filepath, "READ");

    if (!input_file || input_file->IsZombie()) {
      std::cerr << "Error opening file: " << input_filepath << std::endl;
      continue;
    }

    TTree *input_tree = static_cast<TTree *>(input_file->Get("Data_R"));
    if (!input_tree) {
      std::cerr << "Error getting tree from: " << input_filepath << std::endl;
      input_file->Close();
      continue;
    }

    UShort_t board, channel, energy;
    UInt_t flags;
    ULong64_t timestamp;
    input_tree->SetBranchAddress("Board", &board);
    input_tree->SetBranchAddress("Channel", &channel);
    input_tree->SetBranchAddress("Energy", &energy);
    input_tree->SetBranchAddress("ShiftedTimestamp", &timestamp);
    input_tree->SetBranchAddress("Flags", &flags);

    TString output_name = output_names[i];
    TString output_filepath = "root_files/" + output_name + ".root";
    TFile *output_file = new TFile(output_filepath, "RECREATE");

    Int_t leftdE[16], rightdE[16], totaldE[16];
    ULong64_t timestamp_distribution[36];
    UInt_t all_flags[36];
    Int_t hits[36];
    Int_t cathode, strip0, strip17, grid;

    TTree *output_tree = new TTree("event", "MUSIC events");
    output_tree->Branch("LeftdE", leftdE, "LeftdE[16]/I");
    output_tree->Branch("RightdE", rightdE, "RightdE[16]/I");
    output_tree->Branch("TotaldE", totaldE, "TotaldE[16]/I");
    output_tree->Branch("TimestampDistribution", timestamp_distribution,
                        "TimestampDistribution[36]/l");
    output_tree->Branch("AllFlags", all_flags, "AllFlags[36]/i");
    output_tree->Branch("Hits", hits, "Hits[36]/I");
    output_tree->Branch("Cathode", &cathode, "Cathode/I");
    output_tree->Branch("Strip0", &strip0, "Strip0/I");
    output_tree->Branch("Strip17", &strip17, "Strip17/I");
    output_tree->Branch("Grid", &grid, "Grid/I");

    ULong64_t event_time_zero = 0;
    Bool_t in_event = kFALSE;

    ResetEvent(leftdE, rightdE, totaldE, timestamp_distribution, all_flags,
               hits, cathode, strip0, strip17, grid);

    Long64_t n_entries = input_tree->GetEntries();
    std::cout << "Building events from " << n_entries << " entries..."
              << std::endl;

    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      TString map_name = Constants::channelMap.at({board, channel});

      if (map_name == "")
        continue;

      if (map_name == "Grid") {
        if (in_event) {
          for (Int_t k = 0; k < 36; k++) {
            if (timestamp_distribution[k] > 0) {
              timestamp_distribution[k] -= event_time_zero;
            }
          }
          output_tree->Fill();
        }

        ResetEvent(leftdE, rightdE, totaldE, timestamp_distribution, all_flags,
                   hits, cathode, strip0, strip17, grid);

        grid = energy;
        event_time_zero = timestamp;
        timestamp_distribution[35] = timestamp;
        all_flags[35] = flags;
        hits[35]++;
        in_event = kTRUE;
      }

      else {
        if (!in_event)
          continue;

        Long64_t time_diff = (timestamp > event_time_zero)
                                 ? (timestamp - event_time_zero)
                                 : (event_time_zero - timestamp);

        if (TMath::Abs(time_diff) > TIME_WINDOW) {
          for (Int_t k = 0; k < 36; k++) {
            if (timestamp_distribution[k] > 0) {
              timestamp_distribution[k] -= event_time_zero;
            }
          }
          output_tree->Fill();

          ResetEvent(leftdE, rightdE, totaldE, timestamp_distribution,
                     all_flags, hits, cathode, strip0, strip17, grid);
          in_event = kFALSE;
          continue;
        }
        std::cout << time_diff * 1e-6 << std::endl;
        if (map_name == "Strip0") {
          strip0 = energy;
          timestamp_distribution[0] = timestamp;
          all_flags[0] = flags;
          hits[0]++;

        } else if (map_name.BeginsWith("L")) {
          Int_t strip_num = GetStripNumber(map_name);
          leftdE[strip_num] = energy;
          totaldE[strip_num] += energy;
          timestamp_distribution[strip_num] = timestamp;
          all_flags[strip_num] = flags;
          hits[strip_num]++;

        } else if (map_name.BeginsWith("R")) {
          Int_t strip_num = GetStripNumber(map_name);
          rightdE[strip_num] = energy;
          totaldE[strip_num] += energy;
          timestamp_distribution[strip_num + 16] = timestamp;
          all_flags[strip_num + 16] = flags;
          hits[strip_num + 16]++;

        } else if (map_name == "Strip17") {
          strip17 = energy;
          timestamp_distribution[33] = timestamp;
          all_flags[33] = flags;
          hits[33]++;

        } else if (map_name == "Cathode") {
          cathode = energy;
          timestamp_distribution[34] = timestamp;
          all_flags[34] = flags;
          hits[34]++;
        }
      }

      if (j % 1000000 == 0) {
        std::cout << "  Progress: " << j << "/" << n_entries << std::endl;
      }
    }

    if (in_event && grid > 0) {
      for (Int_t k = 0; k < 36; k++) {
        if (timestamp_distribution[k] > 0) {
          timestamp_distribution[k] -= event_time_zero;
        }
      }
      output_tree->Fill();
    }
    output_file->cd();
    output_tree->Write();
    output_file->Close();
    input_file->Close();

    std::cout << "Events saved to: " << output_filepath << std::endl;
  }
}

void EventBuilder() {
  Bool_t reprocess_initial = kTRUE;

  std::vector<TString> filenames, output_names;

  for (Int_t i = 0; i < Constants::N_FILES; i++) {
    filenames.push_back(Form("Sorted_Run37_%d", i));
    std::cout << "Processing file: " << std::endl;
    std::cout << filenames[i] << std::endl;
    output_names.push_back(Form("Events_Run37_%d", i));
  }

  InitUtils::SetROOTPreferences();
  BuildEvents(filenames, output_names, reprocess_initial);
}
