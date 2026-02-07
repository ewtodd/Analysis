#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TLine.h>

void DrawTailFunctions() {
  InitUtils::SetROOTPreferences();

  Double_t mu = 59.5;
  Double_t sigma = 1;
  Double_t tail_amplitude = 1;
  Double_t tail_range = 0.9;

  TF1 *low_tail = new TF1(
      "low_tail",
      [](double *x, double *p) { return FittingFunctions::LowTail(x, p); }, 40,
      80, 4);

  low_tail->SetParameter(0, mu);
  low_tail->SetParameter(1, sigma);
  low_tail->SetParameter(2, tail_amplitude);
  low_tail->SetParameter(3, tail_range);
  low_tail->SetNpx(1000);

  TF1 *high_tail = new TF1(
      "high_tail",
      [](double *x, double *p) { return FittingFunctions::HighTail(x, p); }, 40,
      80, 4);

  high_tail->SetParameter(0, mu);
  high_tail->SetParameter(1, sigma);
  high_tail->SetParameter(2, tail_amplitude);
  high_tail->SetParameter(3, tail_range);
  high_tail->SetNpx(1000);

  TCanvas *c1 = new TCanvas("c1", "Low Tail Test", 800, 600);
  PlottingUtils::ConfigureCanvas(c1, kFALSE);
  low_tail->SetTitle("Low Tail Function;x;Amplitude");
  low_tail->SetLineColor(kBlue);
  low_tail->SetLineWidth(2);
  low_tail->Draw();

  TLine *line1 = new TLine(mu, 0, mu, low_tail->GetMaximum());
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw();

  PlottingUtils::SaveFigure(c1, "low_tail_test.png", kTRUE);

  TCanvas *c2 = new TCanvas("c2", "High Tail Test", 800, 600);
  PlottingUtils::ConfigureCanvas(c2, kFALSE);
  high_tail->SetTitle("High Tail Function;x;Amplitude");
  high_tail->SetLineColor(kRed);
  high_tail->SetLineWidth(2);
  high_tail->Draw();

  TLine *line2 = new TLine(mu, 0, mu, high_tail->GetMaximum());
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(2);
  line2->SetLineWidth(2);
  line2->Draw();

  PlottingUtils::SaveFigure(c2, "high_tail_test.png", kTRUE);

  delete low_tail;
  delete high_tail;
}
