#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
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

    TH2F *PSPvsLO =
        new TH2F("", "; Light Output [keVee]; Counts", 250, 0, 1200, 100, 0, 1);

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
      if ((0.07 < psp && psp < 0.25) &&
          (130 < light_output_keVee && light_output_keVee < 170)) {
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

void AnalyzeDoubleWaveforms(Bool_t reprocess = kFALSE) {
  if (!reprocess) {
    return;
  }

  TString candidates_filepath = "root_files/candidates.root";

  TFile *input = new TFile(candidates_filepath, "READ");
  TTree *features_tree = static_cast<TTree *>(input->Get("features"));

  TArrayS *samples = nullptr;
  Float_t light_output, psp;
  features_tree->SetBranchAddress("Samples", &samples);
  features_tree->SetBranchAddress("light_output", &light_output);
  features_tree->SetBranchAddress("psp", &psp);

  TString output_filepath = "root_files/double_waveforms.root";

  TFile *output = new TFile(output_filepath, "RECREATE");
  TTree *output_tree = new TTree("features", "Double waveform features.");
  TArrayS *output_samples = nullptr;
  Float_t output_psp, output_light_output;
  Int_t peak1_position, peak2_position;

  output_tree->Branch("light_output", &output_light_output, "light_output/F");
  output_tree->Branch("psp", &output_psp, "psp/F");
  output_tree->Branch("peak1_position", &peak1_position, "peak1_position/I");
  output_tree->Branch("peak2_position", &peak2_position, "peak2_position/I");
  output_tree->Branch("Samples", &output_samples);

  Int_t num_entries = features_tree->GetEntries();

  for (Int_t i = 0; i < num_entries; i++) {
    features_tree->GetEntry(i);

    Int_t nsamples = samples->GetSize();
    std::vector<Double_t> waveform_array(nsamples);
    for (Int_t j = 0; j < nsamples; j++) {
      waveform_array[j] = (Double_t)samples->At(j);
    }

    TSpectrum *spectrum = new TSpectrum(10);
    Int_t npeaks =
        spectrum->SearchHighRes(waveform_array.data(), waveform_array.data(),
                                nsamples, 2, 5, kFALSE, 3, kTRUE, 3);

    if (npeaks == 2) {
      Double_t *peak_positions = spectrum->GetPositionX();

      output_psp = psp;
      output_light_output = light_output;
      peak1_position = (Int_t)peak_positions[0];
      peak2_position = (Int_t)peak_positions[1];

      if (output_samples)
        delete output_samples;
      output_samples = new TArrayS(*samples);

      output_tree->Fill();
    }

    delete spectrum;
  }

  input->Close();
  output->cd();
  output_tree->Write("", TObject::kOverwrite);
  output->Close();
}

void CreateTemplateWaveform(Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  TString template_148keV_filepath = "root_files/calibration_Eu152.root";
  TFile *input = new TFile(template_148keV_filepath, "READ");
  TTree *features_tree = static_cast<TTree *>(input->Get("features"));

  TArrayS *samples = nullptr;
  Float_t pulse_height;
  features_tree->SetBranchAddress("Samples", &samples);
  features_tree->SetBranchAddress("pulse_height", &pulse_height);

  Int_t num_entries = features_tree->GetEntries();
  std::vector<std::vector<Double_t>> selected_waveforms;
  Int_t nsamples = 0;

  for (Int_t i = 0; i < num_entries; i++) {
    features_tree->GetEntry(i);

    if (pulse_height > 1500 && pulse_height < 1800) {
      nsamples = samples->GetSize();
      std::vector<Double_t> waveform(nsamples);
      for (Int_t j = 0; j < nsamples; j++) {
        waveform[j] = (Double_t)samples->At(j);
      }
      selected_waveforms.push_back(waveform);
    }
  }

  input->Close();

  std::cout << "Found " << selected_waveforms.size()
            << " waveforms for template" << std::endl;

  std::vector<Double_t> template_waveform(nsamples, 0.0);
  for (const auto &wf : selected_waveforms) {
    for (Int_t j = 0; j < nsamples; j++) {
      template_waveform[j] += wf[j];
    }
  }

  for (Int_t j = 0; j < nsamples; j++) {
    template_waveform[j] /= selected_waveforms.size();
  }

  std::vector<Double_t> x_values(nsamples);
  for (Int_t j = 0; j < nsamples; j++) {
    x_values[j] = j;
  }

  TString output_filepath = "root_files/template_waveform.root";
  TFile *output = new TFile(output_filepath, "RECREATE");
  TGraph *template_graph =
      new TGraph(nsamples, x_values.data(), template_waveform.data());
  template_graph->SetName("template");
  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas);
  PlottingUtils::ConfigureAndDrawGraph(
      template_graph, kRed, "122 keV 152Eu Template Waveform;Sample;ADC");
  PlottingUtils::SaveFigure(canvas, "template.png", kFALSE);
  template_graph->Write("", TObject::kOverwrite);
  output->Close();
}

void FitDoubleWaveforms(Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  TString template_148keV_filepath = "root_files/calibration_Eu152.root";

  TFile *template_file = new TFile(template_148keV_filepath, "READ");
  TGraph *template_graph =
      static_cast<TGraph *>(template_file->Get("template"));

  Double_t *template_y = template_graph->GetY();
  Int_t template_npoints = template_graph->GetN();
  Int_t template_peak_pos = TMath::LocMax(template_npoints, template_y);
  Double_t template_peak_height = template_y[template_peak_pos];

  TFile *input = new TFile("double_waveforms.root", "READ");
  TTree *features_tree = static_cast<TTree *>(input->Get("features"));

  TArrayS *samples = nullptr;
  Int_t peak1_position, peak2_position;
  Float_t light_output, psp;
  features_tree->SetBranchAddress("Samples", &samples);
  features_tree->SetBranchAddress("peak1_position", &peak1_position);
  features_tree->SetBranchAddress("peak2_position", &peak2_position);
  features_tree->SetBranchAddress("light_output", &light_output);
  features_tree->SetBranchAddress("psp", &psp);

  TString output_filepath = "root_files/fitted_doubles.root";
  TFile *output = new TFile(output_filepath, "RECREATE");
  TTree *output_tree =
      new TTree("features", "Fitted double waveform features.");
  Float_t peak1_height, peak2_height, time_difference;

  output_tree->Branch("peak1_height", &peak1_height, "peak1_height/F");
  output_tree->Branch("peak2_height", &peak2_height, "peak2_height/F");
  output_tree->Branch("time_difference", &time_difference, "time_difference/F");

  Int_t num_entries = features_tree->GetEntries();
  Int_t plot_count = TMath::Min(5, num_entries);

  for (Int_t i = 0; i < num_entries; i++) {
    features_tree->GetEntry(i);

    Int_t nsamples = samples->GetSize();
    std::vector<Double_t> waveform(nsamples);
    std::vector<Double_t> x_values(nsamples);
    for (Int_t j = 0; j < nsamples; j++) {
      waveform[j] = (Double_t)samples->At(j);
      x_values[j] = j;
    }

    TGraph *original = nullptr;
    if (i < plot_count) {
      original = new TGraph(nsamples, x_values.data(), waveform.data());
    }

    peak1_height = waveform[peak1_position];
    Double_t scale_factor = peak1_height / template_peak_height;

    for (Int_t j = 0; j < nsamples; j++) {
      Int_t template_index = j - peak1_position + template_peak_pos;
      Double_t template_val = 0;
      if (template_index >= 0 && template_index < template_npoints) {
        template_val = template_y[template_index] * scale_factor;
      }
      waveform[j] -= template_val;
    }

    peak2_height = waveform[peak2_position];
    time_difference = (peak2_position - peak1_position) * 2.0;

    output_tree->Fill();

    if (i < plot_count) {
      TGraph *subtracted =
          new TGraph(nsamples, x_values.data(), waveform.data());

      TCanvas *canvas = new TCanvas("", "", 1200, 800);
      canvas->Divide(1, 2);

      canvas->cd(1);
      PlottingUtils::ConfigureAndDrawGraph(
          original, kBlue + 1, "Original Double Waveform;Sample;ADC");

      canvas->cd(2);
      PlottingUtils::ConfigureAndDrawGraph(subtracted, kRed + 1,
                                           "Subtracted Waveform;Sample;ADC");

      PlottingUtils::SaveFigure(canvas, Form("double_waveform_%d.png", i),
                                kFALSE);

      delete canvas;
      delete original;
      delete subtracted;
    }
  }

  input->Close();
  template_file->Close();
  output->cd();
  output_tree->Write();
  output->Close();
}

void ConvertAndPlotDoublePeaks(Bool_t reprocess = kFALSE) {
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

  if (reprocess) {
    features_tree->Branch("peak1_light_output", &peak1_light_output,
                          "peak1_light_output/F");
    features_tree->Branch("peak2_light_output", &peak2_light_output,
                          "peak2_light_output/F");

    Int_t num_entries = features_tree->GetEntries();
    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      peak1_light_output = calibration_function->Eval(peak1_height);
      peak2_light_output = calibration_function->Eval(peak2_height);
      features_tree->GetBranch("peak1_light_output")->Fill();
      features_tree->GetBranch("peak2_light_output")->Fill();
    }

    output->cd();
    features_tree->Write("", TObject::kOverwrite);
  }

  features_tree->SetBranchAddress("peak1_light_output", &peak1_light_output);
  features_tree->SetBranchAddress("peak2_light_output", &peak2_light_output);

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
    LO1vsLO2->Fill(peak1_light_output, peak2_light_output);
    peak1_hist->Fill(peak1_light_output);
    peak2_hist->Fill(peak2_light_output);
  }

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

void FitAndExtractHalfLife() {
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

  TH1F *time_diff_hist =
      new TH1F("", "; Time Difference [ns]; Counts", 62, 0, 250);

  Int_t num_entries = features_tree->GetEntries();
  for (Int_t j = 0; j < num_entries; j++) {
    features_tree->GetEntry(j);

    Double_t dist_x = (peak1_light_output - mean_x) / sigma_x;
    Double_t dist_y = (peak2_light_output - mean_y) / sigma_y;
    Double_t distance = TMath::Sqrt(dist_x * dist_x + dist_y * dist_y);

    if (distance < 1.0) {
      time_diff_hist->Fill(time_difference);
    }
  }

  TF1 *exponential = new TF1("exponential", "[0]*TMath::Exp(-x/[1])", 49, 200);
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
  PlottingUtils::ConfigureAndDrawHistogram(time_diff_hist, kBlue + 1);
  exponential->Draw("SAME");

  TLegend *leg = new TLegend(0.6, 0.75, 0.88, 0.88);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.05);
  leg->SetTextFont(132);
  leg->AddEntry((TObject *)0,
                Form("t_{1/2} = %.1f #pm %.1f ns", half_life, half_life_error),
                "");
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

  CalculatePSPvsLO(input_names, kTRUE);
  PlotPSPvsLO(input_names);

  std::vector<TString> irradiation_names = {
      input_name_irradiation_one, input_name_irradiation_two,
      input_name_irradiation_three, input_name_irradiation_four};
  FindCandidateWaveforms(irradiation_names);
  PlotPSPvsLO({"candidates"});
  AnalyzeDoubleWaveforms(kTRUE);
  PlotPSPvsLO({"double_waveforms"});
  CreateTemplateWaveform(kTRUE);
  FitDoubleWaveforms(kTRUE);
  ConvertAndPlotDoublePeaks(kTRUE);
  FitAndExtractHalfLife();
}
