#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TMarker.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TSystem.h>
#include <TTree.h>
#include <vector>

void CalculatePSPvsLO(std::vector<TString> input_names,
                      Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;
  Int_t entries = input_names.size();
  for (Int_t i = 0; i < entries; i++) {
    TString input_name = input_names[i];

    TString output_filepath = "root_files/" + input_name + ".root";

    TFile *output = new TFile(output_filepath, "UPDATE");
    TTree *features_tree = static_cast<TTree *>(output->Get("features"));

    Float_t short_integral, long_integral, light_output_keVee, psp;
    features_tree->SetBranchAddress("short_integral", &short_integral);
    features_tree->SetBranchAddress("long_integral", &long_integral);
    features_tree->SetBranchAddress("light_output", &light_output_keVee);
    features_tree->Branch("psp", &psp, "psp/F");

    Int_t num_entries = features_tree->GetEntries();
    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      psp = long_integral > 1e-3
                ? (long_integral - short_integral) / long_integral
                : 0;
      features_tree->GetBranch("psp")->Fill();
    }

    features_tree->Write("", TObject::kOverwrite);
    output->Close();
  }
}

void PlotPSPvsLO(std::vector<TString> input_names) {
  Int_t entries = input_names.size();
  for (Int_t i = 0; i < entries; i++) {
    TString input_name = input_names[i];

    TH2F *PSPvsLO = new TH2F("", "; Light Output [keVee]; PSP; Counts", 250, 0,
                             1200, 100, 0, 1);

    TCanvas *canvas = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvas);

    TString output_filepath = "root_files/" + input_name + ".root";
    TFile *output = new TFile(output_filepath, "UPDATE");
    TTree *features_tree = static_cast<TTree *>(output->Get("features"));
    Float_t light_output_keVee, psp;
    features_tree->SetBranchAddress("light_output", &light_output_keVee);
    features_tree->SetBranchAddress("psp", &psp);

    Int_t num_entries = features_tree->GetEntries();
    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      PSPvsLO->Fill(light_output_keVee, psp);
    }

    PlottingUtils::ConfigureAndDraw2DHistogram(PSPvsLO, canvas);
    PlottingUtils::SaveFigure(canvas, input_name + "_psp_vs_light_output.png",
                              kFALSE);

    PSPvsLO->Write("PSP vs. Light Output", TObject::kOverwrite);
    output->Close();
    delete canvas;
    delete PSPvsLO;
  }
}

void FindCandidateWaveforms(std::vector<TString> input_names,
                            Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t entries = input_names.size();

  TString candidates_filepath = "root_files/candidates.root";

  TFile *output = new TFile(candidates_filepath, "RECREATE");
  TTree *output_tree = new TTree("features", "Candidate waveform features.");
  TArrayS *output_samples = nullptr;
  Float_t output_psp, output_light_output_keVee;

  output_tree->Branch("light_output", &output_light_output_keVee,
                      "light_output/F");
  output_tree->Branch("psp", &output_psp, "output_psp/F");
  output_tree->Branch("Samples", &output_samples);

  for (Int_t i = 0; i < entries; i++) {

    TString input_name = input_names[i];

    TString input_filepath = "root_files/" + input_name + ".root";

    TFile *input = new TFile(input_filepath, "READ");
    TTree *features_tree = static_cast<TTree *>(input->Get("features"));

    Float_t light_output_keVee, psp;
    features_tree->SetBranchAddress("light_output", &light_output_keVee);
    features_tree->SetBranchAddress("psp", &psp);

    TArrayS *samples = nullptr;
    features_tree->SetBranchAddress("Samples", &samples);

    Int_t num_entries = features_tree->GetEntries();
    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      if ((0.08 < psp && psp < 0.25) &&
          (125 < light_output_keVee && light_output_keVee < 170)) {
        output_psp = psp;
        output_light_output_keVee = light_output_keVee;
        if (output_samples)
          delete output_samples;
        output_samples = new TArrayS(*samples);
        output_tree->Fill();
      }
    }
    input->Close();
  }

  output->cd();
  output_tree->Write("", TObject::kOverwrite);
  output->Close();
}

void CreateTemplateWaveforms(Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  TString template_148keV_filepath = "root_files/calibration_Eu152.root";
  TFile *input_148keV = new TFile(template_148keV_filepath, "READ");
  TTree *features_tree_148keV =
      static_cast<TTree *>(input_148keV->Get("features"));

  TArrayS *samples_148keV = nullptr;
  Float_t pulse_height_148keV;
  features_tree_148keV->SetBranchAddress("Samples", &samples_148keV);
  features_tree_148keV->SetBranchAddress("pulse_height", &pulse_height_148keV);

  Int_t num_entries_148keV = features_tree_148keV->GetEntries();
  std::vector<std::vector<Double_t>> selected_waveforms_148keV;
  Int_t nsamples_148keV = 0;

  for (Int_t i = 0; i < num_entries_148keV; i++) {
    features_tree_148keV->GetEntry(i);

    if (pulse_height_148keV > 1500 && pulse_height_148keV < 1800) {
      nsamples_148keV = samples_148keV->GetSize();
      std::vector<Double_t> waveform(nsamples_148keV);
      for (Int_t j = 0; j < nsamples_148keV; j++) {
        waveform[j] = (Double_t)samples_148keV->At(j);
      }
      selected_waveforms_148keV.push_back(waveform);
    }
  }

  input_148keV->Close();

  std::cout << "Found " << selected_waveforms_148keV.size()
            << " waveforms for template 148 keV" << std::endl;

  std::vector<Double_t> template_waveform_148keV(nsamples_148keV, 0.0);
  for (const auto &wf : selected_waveforms_148keV) {
    for (Int_t j = 0; j < nsamples_148keV; j++) {
      template_waveform_148keV[j] += wf[j];
    }
  }

  for (Int_t j = 0; j < nsamples_148keV; j++) {
    template_waveform_148keV[j] /= selected_waveforms_148keV.size();
  }

  std::vector<Double_t> x_values(nsamples_148keV);
  for (Int_t j = 0; j < nsamples_148keV; j++) {
    x_values[j] = 2 * j;
  }

  TString output_filepath = "root_files/template_waveforms.root";
  TFile *output = new TFile(output_filepath, "RECREATE");
  TGraph *template_graph_148keV = new TGraph(nsamples_148keV, x_values.data(),
                                             template_waveform_148keV.data());
  template_graph_148keV->SetName("template_148keV");
  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas);
  PlottingUtils::ConfigureGraph(
      template_graph_148keV, kRed,
      "122 keV 152Eu Template Waveform;Time [ns];ADC");
  template_graph_148keV->SetLineWidth(2);
  template_graph_148keV->Draw();
  PlottingUtils::SaveFigure(canvas, "template_148keV.png", kFALSE);
  template_graph_148keV->Write("", TObject::kOverwrite);

  TString template_32keV_filepath = "root_files/calibration_Eu152.root";
  TFile *input_32keV = new TFile(template_32keV_filepath, "READ");
  TTree *features_tree_32keV =
      static_cast<TTree *>(input_32keV->Get("features"));

  TArrayS *samples_32keV = nullptr;
  Float_t pulse_height_32keV;
  features_tree_32keV->SetBranchAddress("Samples", &samples_32keV);
  features_tree_32keV->SetBranchAddress("pulse_height", &pulse_height_32keV);

  Int_t num_entries_32keV = features_tree_32keV->GetEntries();
  std::vector<std::vector<Double_t>> selected_waveforms_32keV;
  Int_t nsamples_32keV = 0;

  for (Int_t i = 0; i < num_entries_32keV; i++) {
    features_tree_32keV->GetEntry(i);

    if (pulse_height_32keV > 450 && pulse_height_32keV < 600) {
      nsamples_32keV = samples_32keV->GetSize();
      std::vector<Double_t> waveform(nsamples_32keV);
      for (Int_t j = 0; j < nsamples_32keV; j++) {
        waveform[j] = (Double_t)samples_32keV->At(j);
      }
      selected_waveforms_32keV.push_back(waveform);
    }
  }

  input_32keV->Close();

  std::cout << "Found " << selected_waveforms_32keV.size()
            << " waveforms for template 32 keV" << std::endl;

  std::vector<Double_t> template_waveform_32keV(nsamples_32keV, 0.0);
  for (const auto &wf : selected_waveforms_32keV) {
    for (Int_t j = 0; j < nsamples_32keV; j++) {
      template_waveform_32keV[j] += wf[j];
    }
  }

  for (Int_t j = 0; j < nsamples_32keV; j++) {
    template_waveform_32keV[j] /= selected_waveforms_32keV.size();
  }

  TGraph *template_graph_32keV = new TGraph(nsamples_32keV, x_values.data(),
                                            template_waveform_32keV.data());
  template_graph_32keV->SetName("template_32keV");
  canvas->Clear();
  PlottingUtils::ConfigureCanvas(canvas);
  PlottingUtils::ConfigureGraph(
      template_graph_32keV, kBlue,
      "33.4 keV La K-alpha Template Waveform;Time [ns];ADC");
  template_graph_32keV->SetLineWidth(2);
  template_graph_32keV->Draw();
  PlottingUtils::SaveFigure(canvas, "template_32keV.png", kFALSE);
  output->cd();
  template_graph_32keV->Write("", TObject::kOverwrite);

  output->Close();
}

void FitDoubleWaveforms(Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  TString calibration_filepath = "root_files/calibration_function.root";
  TFile *calib_file = new TFile(calibration_filepath, "READ");
  TF1 *calibration_function =
      static_cast<TF1 *>(calib_file->Get("calibration"));

  Double_t calibration_slope = calibration_function->GetParameter(1);

  TString template_filepath = "root_files/template_waveforms.root";
  TFile *template_file = new TFile(template_filepath, "READ");
  TGraph *template_graph_148 =
      static_cast<TGraph *>(template_file->Get("template_148keV"));
  TGraph *template_graph_32 =
      static_cast<TGraph *>(template_file->Get("template_32keV"));

  Double_t *template_148_y = template_graph_148->GetY();
  Double_t *template_32_y = template_graph_32->GetY();
  Int_t template_npoints = template_graph_148->GetN();

  Int_t template_148_peak_pos = TMath::LocMax(template_npoints, template_148_y);
  Int_t template_32_peak_pos = TMath::LocMax(template_npoints, template_32_y);
  Double_t template_148_peak_height = template_148_y[template_148_peak_pos];
  Double_t template_32_peak_height = template_32_y[template_32_peak_pos];

  TString input_filepath = "root_files/candidates.root";
  TFile *input = new TFile(input_filepath, "READ");
  TTree *features_tree = static_cast<TTree *>(input->Get("features"));

  TArrayS *samples = nullptr;
  Float_t light_output, psp;
  features_tree->SetBranchAddress("Samples", &samples);
  features_tree->SetBranchAddress("light_output", &light_output);
  features_tree->SetBranchAddress("psp", &psp);

  TString output_filepath = "root_files/fitted_doubles.root";
  TFile *output = new TFile(output_filepath, "RECREATE");
  TTree *output_tree =
      new TTree("features", "Fitted double waveform features.");
  Float_t peak1_height, peak2_height, time_difference, chi_squared;

  output_tree->Branch("peak1_height", &peak1_height, "peak1_height/F");
  output_tree->Branch("peak2_height", &peak2_height, "peak2_height/F");
  output_tree->Branch("time_difference", &time_difference, "time_difference/F");
  output_tree->Branch("chi_squared", &chi_squared, "chi_squared/F");

  Int_t feature_tree_entries = features_tree->GetEntries();
  Int_t num_entries = TMath::Min(10000, feature_tree_entries);
  Int_t plot_count = TMath::Min(10, num_entries);

  std::cout << "Fitting " << num_entries << " candidate waveforms..."
            << std::endl;

  for (Int_t i = 0; i < num_entries; i++) {
    features_tree->GetEntry(i);

    if (i % 10000 == 0) {
      std::cout << "Processing waveform " << i << "/" << num_entries
                << std::endl;
    }

    Int_t nsamples = samples->GetSize();
    std::vector<Double_t> waveform(nsamples);
    std::vector<Double_t> x_values(nsamples);
    for (Int_t j = 0; j < nsamples; j++) {
      waveform[j] = (Double_t)samples->At(j);
      x_values[j] = j * 2.0;
    }

    Double_t *waveform_ptr = waveform.data();
    Int_t candidate_peak_pos = TMath::LocMax(nsamples, waveform_ptr);
    Double_t candidate_peak_height = waveform[candidate_peak_pos];

    Double_t initial_energy_148 =
        calibration_function->Eval(candidate_peak_height);

    Double_t energy_148_min = initial_energy_148 * 0.9;
    Double_t energy_148_max = initial_energy_148 * 1.1;
    Int_t n_energy_148_steps = 10;

    Double_t energy_32_min = 10.0;
    Double_t energy_32_max = 40.0;
    Int_t n_energy_32_steps = 10;

    Double_t time_diff_min = 10.0;
    Double_t time_diff_max = 200.0;
    Int_t n_time_steps = (time_diff_max - time_diff_min) / 2;

    Double_t best_chi2 = 1e20;
    Double_t best_energy_148 = initial_energy_148;
    Double_t best_energy_32 = 32.0;
    Double_t best_time_diff = 30;

    for (Int_t ie1 = 0; ie1 < n_energy_148_steps; ie1++) {
      Double_t energy_148 =
          energy_148_min +
          ie1 * (energy_148_max - energy_148_min) / (n_energy_148_steps - 1);

      Double_t ph_148 = energy_148 / calibration_slope;
      Double_t scale_148 = ph_148 / template_148_peak_height;

      for (Int_t ie2 = 0; ie2 < n_energy_32_steps; ie2++) {
        Double_t energy_32 =
            energy_32_min +
            ie2 * (energy_32_max - energy_32_min) / (n_energy_32_steps - 1);

        Double_t ph_32 = energy_32 / calibration_slope;
        Double_t scale_32 = ph_32 / template_32_peak_height;

        for (Int_t it = 0; it < n_time_steps; it++) {
          Double_t time_diff =
              time_diff_min +
              it * (time_diff_max - time_diff_min) / (n_time_steps - 1);

          Int_t time_diff_samples = TMath::Nint(time_diff / 2.0);

          std::vector<Double_t> double_template(nsamples, 0.0);

          for (Int_t j = 0; j < nsamples; j++) {
            Int_t template_index =
                j - candidate_peak_pos + template_148_peak_pos;
            if (template_index >= 0 && template_index < template_npoints) {
              double_template[j] += template_148_y[template_index] * scale_148;
            }
          }

          Int_t peak_32_pos = candidate_peak_pos + time_diff_samples;
          for (Int_t j = 0; j < nsamples; j++) {
            Int_t template_index = j - peak_32_pos + template_32_peak_pos;
            if (template_index >= 0 && template_index < template_npoints) {
              double_template[j] += template_32_y[template_index] * scale_32;
            }
          }

          Double_t chi2 = 0.0;
          for (Int_t j = 0; j < nsamples; j++) {
            Double_t diff = waveform[j] - double_template[j];
            chi2 += diff * diff;
          }

          if (chi2 < best_chi2) {
            best_chi2 = chi2;
            best_energy_148 = energy_148;
            best_energy_32 = energy_32;
            best_time_diff = time_diff;
          }
        }
      }
    }

    peak1_height = best_energy_148;
    peak2_height = best_energy_32;
    time_difference = best_time_diff;
    chi_squared = best_chi2;

    output_tree->Fill();

    if (i < plot_count) {
      Double_t ph_148_best = best_energy_148 / calibration_slope;
      Double_t scale_148_best = ph_148_best / template_148_peak_height;

      Double_t ph_32_best = best_energy_32 / calibration_slope;
      Double_t scale_32_best = ph_32_best / template_32_peak_height;

      Int_t time_diff_samples_best = TMath::Nint(best_time_diff / 2.0);

      std::vector<Double_t> best_template(nsamples, 0.0);
      std::vector<Double_t> residuals(nsamples, 0.0);

      for (Int_t j = 0; j < nsamples; j++) {
        Int_t template_index = j - candidate_peak_pos + template_148_peak_pos;
        if (template_index >= 0 && template_index < template_npoints) {
          best_template[j] += template_148_y[template_index] * scale_148_best;
        }
      }

      Int_t peak_32_pos_best = candidate_peak_pos + time_diff_samples_best;
      for (Int_t j = 0; j < nsamples; j++) {
        Int_t template_index = j - peak_32_pos_best + template_32_peak_pos;
        if (template_index >= 0 && template_index < template_npoints) {
          best_template[j] += template_32_y[template_index] * scale_32_best;
        }
      }

      for (Int_t j = 0; j < nsamples; j++) {
        residuals[j] = waveform[j] - best_template[j];
      }

      TGraph *original = new TGraph(nsamples, x_values.data(), waveform.data());
      TGraph *fitted =
          new TGraph(nsamples, x_values.data(), best_template.data());
      TGraph *residual_graph =
          new TGraph(nsamples, x_values.data(), residuals.data());

      TCanvas *canvas = new TCanvas("", "", 1200, 800);
      canvas->Divide(1, 2);

      canvas->cd(1);
      PlottingUtils::ConfigureAndDrawGraph(original, kBlue + 1,
                                           "Original Waveform;Time [ns];ADC");
      fitted->SetLineColor(kRed);
      fitted->SetLineWidth(2);
      fitted->Draw("L SAME");

      TLegend *leg1 = new TLegend(0.65, 0.62, 0.88, 0.88);
      leg1->SetBorderSize(1);
      leg1->SetFillColor(kWhite);
      leg1->SetTextSize(0.045);
      leg1->SetMargin(0.1);
      leg1->SetTextFont(132);
      leg1->AddEntry(original, "Original", "lp");
      leg1->AddEntry(fitted, "Best Fit", "l");
      leg1->AddEntry((TObject *)0, Form("E_{1} = %.1f keV", best_energy_148),
                     "");
      leg1->AddEntry((TObject *)0, Form("E_{2} = %.1f keV", best_energy_32),
                     "");
      leg1->AddEntry((TObject *)0, Form("#Deltat = %.1f ns", best_time_diff),
                     "");
      leg1->AddEntry((TObject *)0, Form("#chi^{2} = %.0f", best_chi2), "");
      leg1->Draw();

      canvas->cd(2);
      PlottingUtils::ConfigureAndDrawGraph(residual_graph, kGreen + 2,
                                           "Residuals;Time [ns];ADC");

      TLine *zero_line = new TLine(x_values[0], 0, x_values[nsamples - 1], 0);
      zero_line->SetLineStyle(2);
      zero_line->SetLineColor(kBlack);
      zero_line->Draw();

      PlottingUtils::SaveFigure(canvas, Form("double_waveform_fit_%d.png", i),
                                kFALSE);

      delete canvas;
      delete original;
      delete fitted;
      delete residual_graph;
    }
  }

  input->Close();
  template_file->Close();
  calib_file->Close();

  output->cd();
  output_tree->Write();
  output->Close();
}

void ConvertAndPlotDoublePeaks(Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;
  TString calibration_filepath = "root_files/calibration_function.root";
  TFile *calib_file = new TFile(calibration_filepath, "READ");
  TF1 *calibration_function =
      static_cast<TF1 *>(calib_file->Get("calibration"));

  TString output_filepath = "root_files/fitted_doubles.root";
  TFile *output = new TFile(output_filepath, "UPDATE");
  TTree *features_tree = static_cast<TTree *>(output->Get("features"));

  Float_t peak1_height, peak2_height, time_difference;
  Float_t peak1_light_output, peak2_light_output;

  features_tree->SetBranchAddress("peak1_height", &peak1_height);
  features_tree->SetBranchAddress("peak2_height", &peak2_height);
  features_tree->SetBranchAddress("time_difference", &time_difference);

  features_tree->Branch("peak1_light_output", &peak1_light_output,
                        "peak1_light_output/F");
  features_tree->Branch("peak2_light_output", &peak2_light_output,
                        "peak2_light_output/F");
  TH2F *LO1vsLO2 =
      new TH2F("", "; Peak 1 Light Output [keVee]; Peak 2 Light Output [keVee]",
               150, 100, 200, 150, 0, 50);
  TH1F *peak1_hist =
      new TH1F("", "; Peak 1 Light Output [keVee]; Counts", 150, 100, 200);
  TH1F *peak2_hist =
      new TH1F("", "; Peak 2 Light Output [keVee]; Counts", 150, 0, 50);

  Int_t num_entries = features_tree->GetEntries();
  for (Int_t j = 0; j < num_entries; j++) {
    features_tree->GetEntry(j);
    peak1_light_output = calibration_function->Eval(peak1_height);
    peak2_light_output = calibration_function->Eval(peak2_height);
    features_tree->GetBranch("peak1_light_output")->Fill();
    features_tree->GetBranch("peak2_light_output")->Fill();
    LO1vsLO2->Fill(peak1_light_output, peak2_light_output);
    peak1_hist->Fill(peak1_light_output);
    peak2_hist->Fill(peak2_light_output);
  }

  output->cd();
  features_tree->Write("", TObject::kOverwrite);

  features_tree->SetBranchAddress("peak1_light_output", &peak1_light_output);
  features_tree->SetBranchAddress("peak2_light_output", &peak2_light_output);

  TCanvas *canvas2d = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas2d);
  PlottingUtils::ConfigureAndDraw2DHistogram(LO1vsLO2, canvas2d);
  PlottingUtils::SaveFigure(canvas2d, "peak1_vs_peak2_light_output.png",
                            kFALSE);

  TCanvas *canvas1 = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas1);
  PlottingUtils::ConfigureAndDrawHistogram(peak1_hist, kBlue + 1);
  PlottingUtils::SaveFigure(canvas1, "peak1_light_output.png");

  TCanvas *canvas2 = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas2);
  PlottingUtils::ConfigureAndDrawHistogram(peak2_hist, kRed + 1);
  PlottingUtils::SaveFigure(canvas2, "peak2_light_output.png");

  output->cd();
  LO1vsLO2->Write("Peak1 vs Peak2 Light Output", TObject::kOverwrite);
  peak1_hist->Write("Peak1 Light Output", TObject::kOverwrite);
  peak2_hist->Write("Peak2 Light Output", TObject::kOverwrite);

  delete canvas2d;
  delete canvas1;
  delete canvas2;
  delete LO1vsLO2;
  delete peak1_hist;
  delete peak2_hist;

  output->Close();
  delete output;
  calib_file->Close();
  delete calib_file;
}

void FitAndExtractHalfLife(Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;
  TString input_filepath = "root_files/fitted_doubles.root";
  TFile *input = new TFile(input_filepath, "READ");
  TTree *features_tree = static_cast<TTree *>(input->Get("features"));
  TH2F *LO1vsLO2 =
      static_cast<TH2F *>(input->Get("Peak1 vs Peak2 Light Output"));

  Float_t peak1_light_output, peak2_light_output, time_difference;
  features_tree->SetBranchAddress("peak1_light_output", &peak1_light_output);
  features_tree->SetBranchAddress("peak2_light_output", &peak2_light_output);
  features_tree->SetBranchAddress("time_difference", &time_difference);

  TF2 *gaussian2d =
      new TF2("gaussian2d", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",
              0, 200, 0, 200);
  gaussian2d->SetParameters(LO1vsLO2->GetMaximum(), LO1vsLO2->GetMean(1),
                            LO1vsLO2->GetRMS(1), LO1vsLO2->GetMean(2),
                            LO1vsLO2->GetRMS(2));

  LO1vsLO2->Fit(gaussian2d, "Q");

  Double_t mean_x = gaussian2d->GetParameter(1);
  Double_t sigma_x = gaussian2d->GetParameter(2);
  Double_t mean_y = gaussian2d->GetParameter(3);
  Double_t sigma_y = gaussian2d->GetParameter(4);

  std::cout << "2D Gaussian fit: " << std::endl;
  std::cout << "Peak 1 mean: " << mean_x << " +/- " << sigma_x << std::endl;
  std::cout << "Peak 2 mean: " << mean_y << " +/- " << sigma_y << std::endl;

  Int_t lower = 0;
  Int_t upper = 300;

  Int_t bins = (upper - lower) / 2;

  TH1F *time_diff_hist =
      new TH1F("", "; Time Difference [ns]; Counts / 2 ns", bins, lower, upper);
  time_diff_hist->Sumw2();

  Int_t num_entries = features_tree->GetEntries();
  for (Int_t j = 0; j < num_entries; j++) {
    features_tree->GetEntry(j);

    Double_t dist_x = (peak1_light_output - mean_x) / sigma_x;
    Double_t dist_y = (peak2_light_output - mean_y) / sigma_y;
    Double_t distance = TMath::Sqrt(dist_x * dist_x + dist_y * dist_y);

    if (distance < 2.0) {
      time_diff_hist->Fill(time_difference);
    }
  }

  Int_t fit_lower = 55;
  Int_t fit_upper = 220;
  TF1 *exponential =
      new TF1("exponential", "[0]*TMath::Exp(-x/[1])", fit_lower, fit_upper);
  exponential->SetParameters(time_diff_hist->GetMaximum(), 100);
  exponential->SetParNames("Amplitude", "Lifetime");

  time_diff_hist->Fit(exponential, "R");

  Double_t lifetime = exponential->GetParameter(1);
  Double_t lifetime_error = exponential->GetParError(1);
  Double_t half_life = lifetime * TMath::Log(2);
  Double_t half_life_error = lifetime_error * TMath::Log(2);

  std::cout << "Lifetime: " << lifetime << " +/- " << lifetime_error << " ns"
            << std::endl;
  std::cout << "Half-life: " << half_life << " +/- " << half_life_error << " ns"
            << std::endl;

  TCanvas *canvas_2d = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas_2d);
  PlottingUtils::ConfigureAndDraw2DHistogram(LO1vsLO2, canvas_2d);
  gaussian2d->Draw("CONT3 SAME");
  PlottingUtils::SaveFigure(canvas_2d, "fitted_2d_gaussian.png", kFALSE);

  TCanvas *canvas_time = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas_time);
  PlottingUtils::ConfigureHistogram(time_diff_hist, kBlue + 1);
  time_diff_hist->SetMarkerSize(2);
  time_diff_hist->Draw("E SAME");
  TMarker *marker_lower =
      new TMarker(fit_lower, exponential->Eval(fit_lower), 21);
  marker_lower->SetMarkerColor(kRed);
  marker_lower->SetMarkerSize(1);
  marker_lower->Draw();

  Float_t marker_upper_y =
      exponential->Eval(fit_upper) > time_diff_hist->GetMinimum()
          ? exponential->Eval(fit_upper)
          : time_diff_hist->GetMinimum();
  TMarker *marker_upper = new TMarker(fit_upper, marker_upper_y, 21);
  marker_upper->SetMarkerColor(kRed);
  marker_upper->SetMarkerSize(1);
  marker_upper->Draw();
  exponential->SetRange(lower, upper);
  exponential->Draw("SAME");

  TLegend *leg = new TLegend(0.6, 0.75, 0.88, 0.88);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.05);
  leg->SetTextFont(132);
  leg->AddEntry((TObject *)0,
                Form("t_{1/2} = %.1f #pm %.1f ns", half_life, half_life_error),
                "");
  leg->AddEntry(marker_upper, "Fit region");
  leg->SetMargin(0.1);
  leg->Draw();

  PlottingUtils::SaveFigure(canvas_time, "time_difference_fit.png");

  TString output_filepath = "root_files/half_life_results.root";
  TFile *output = new TFile(output_filepath, "RECREATE");
  LO1vsLO2->Write("2D_Gaussian_Fit");
  time_diff_hist->Write("Time_Difference");
  gaussian2d->Write("gaussian2d_fit");
  exponential->Write("exponential_fit");
  output->Close();

  delete canvas_2d;
  delete canvas_time;
  delete time_diff_hist;
  delete gaussian2d;
  delete exponential;
  delete output;

  input->Close();
  delete input;
}

void HalfLife() {
  Bool_t reprocess_initial = kTRUE;
  Bool_t reprocess_candidate = kTRUE;
  Bool_t reprocess_waveform = kFALSE;
  Bool_t reprocess_halflife = kFALSE;
  PlottingUtils::SetROOTPreferences();

  TString input_name_Am241 = "calibration_Am241";
  TString input_name_Eu152 = "calibration_Eu152";
  TString input_name_bkg = "background";
  TString input_name_irradiation_one = "irradiation_one";
  TString input_name_irradiation_two = "irradiation_two";
  TString input_name_irradiation_three = "irradiation_three";
  TString input_name_irradiation_four = "irradiation_four";

  std::vector<Int_t> defaultColors = PlottingUtils::GetDefaultColors();
  std::vector<TString> input_names = {input_name_Am241,
                                      input_name_Eu152,
                                      input_name_bkg,
                                      input_name_irradiation_one,
                                      input_name_irradiation_two,
                                      input_name_irradiation_three,
                                      input_name_irradiation_four};

  CalculatePSPvsLO(input_names, reprocess_initial);
  PlotPSPvsLO(input_names);

  std::vector<TString> irradiation_names = {
      input_name_irradiation_one, input_name_irradiation_two,
      input_name_irradiation_three, input_name_irradiation_four};

  FindCandidateWaveforms(irradiation_names, reprocess_candidate);
  if (reprocess_candidate)
    PlotPSPvsLO({"candidates"});
  CreateTemplateWaveforms(reprocess_waveform);
  FitDoubleWaveforms(reprocess_waveform);
  if (reprocess_waveform)
    PlotPSPvsLO({"fitted_doubles"});
  ConvertAndPlotDoublePeaks(reprocess_halflife);
  FitAndExtractHalfLife(reprocess_halflife);
}
