#ifndef INTERACTIVESNIPEDITOR_HPP
#define INTERACTIVESNIPEDITOR_HPP

#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TGraph.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGSlider.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TRootEmbeddedCanvas.h>
#include <TSpline.h>
#include <TSystem.h>
#include <TTimer.h>
#include <fstream>
#include <iostream>

const Int_t N_SNIP_REGIONS = 5;
const Double_t SNIP_REGION_WIDTH = 300.0; // keV per region

struct SNIPParameters {
  Double_t m[N_SNIP_REGIONS]; // window [keV] for each region
  Bool_t accepted;
  Bool_t use_for_all;
};

// Evaluate the SNIP window at a given energy using cubic spline interpolation
// through (0, 0) and the N_SNIP_REGIONS control points.
inline Double_t EvalSNIPWindow(Double_t energy, const Double_t *m) {
  const Int_t n_pts = N_SNIP_REGIONS + 1;
  Double_t x[n_pts], y[n_pts];

  // Anchor at origin
  x[0] = 0.0;
  y[0] = 0.0;

  for (Int_t r = 0; r < N_SNIP_REGIONS; r++) {
    x[r + 1] = (r + 0.5) * SNIP_REGION_WIDTH;
    y[r + 1] = m[r];
  }

  TGraph graph(n_pts, x, y);
  TSpline3 spline("snip_spline", &graph);

  Double_t w = spline.Eval(energy);
  if (w < 0)
    w = 0;
  return w;
}

// Compute adaptive SNIP background with cubic-spline interpolated window
// sizes from N_SNIP_REGIONS control points, anchored at (0, 0).
// Windows are monotonically increasing by construction of the control points.
inline TH1F *ComputeAdaptiveSNIPBackground(TH1F *hist, const Double_t *m) {
  Int_t nbins = hist->GetNbinsX();
  Double_t bin_width = hist->GetBinWidth(1);

  // Compute per-bin window sizes via cubic spline
  Int_t *window = new Int_t[nbins];
  Int_t max_window = 0;
  for (Int_t i = 0; i < nbins; i++) {
    Double_t energy = hist->GetBinCenter(i + 1);
    Double_t w = EvalSNIPWindow(energy, m);
    window[i] = TMath::Max(1, (Int_t)TMath::Nint(w / bin_width));
    if (window[i] > max_window)
      max_window = window[i];
  }

  // LLS transform
  Double_t *y = new Double_t[nbins];
  for (Int_t i = 0; i < nbins; i++) {
    Double_t val = hist->GetBinContent(i + 1);
    if (val < 0)
      val = 0;
    y[i] = TMath::Log(TMath::Log(TMath::Sqrt(val + 1.0) + 1.0) + 1.0);
  }

  // SNIP clipping: iterate p from max_window down to 1
  for (Int_t p = max_window; p >= 1; p--) {
    Double_t *y_new = new Double_t[nbins];
    for (Int_t i = 0; i < nbins; i++) {
      if (p <= window[i] && i - p >= 0 && i + p < nbins) {
        Double_t avg = (y[i - p] + y[i + p]) / 2.0;
        y_new[i] = TMath::Min(y[i], avg);
      } else {
        y_new[i] = y[i];
      }
    }
    for (Int_t i = 0; i < nbins; i++)
      y[i] = y_new[i];
    delete[] y_new;
  }

  // Inverse LLS transform
  TH1F *background =
      static_cast<TH1F *>(hist->Clone(PlottingUtils::GetRandomName()));
  background->SetDirectory(0);
  for (Int_t i = 0; i < nbins; i++) {
    Double_t v = TMath::Exp(TMath::Exp(y[i]) - 1.0) - 1.0;
    Double_t val = (v - 1.0) * (v - 1.0) - 1.0;
    if (val < 0)
      val = 0;
    background->SetBinContent(i + 1, val);
  }

  delete[] y;
  delete[] window;
  return background;
}

// Save SNIP parameters to a .snip text file
inline void SaveSNIPParameters(const TString &path, const Double_t *m) {
  std::ofstream out(path.Data());
  if (!out.is_open()) {
    std::cerr << "WARNING: Cannot save SNIP parameters to " << path
              << std::endl;
    return;
  }
  for (Int_t r = 0; r < N_SNIP_REGIONS; r++) {
    Double_t e_center = (r + 0.5) * SNIP_REGION_WIDTH;
    out << e_center << " " << m[r] << "\n";
  }
  out.close();
}

// Load SNIP parameters from a .snip text file. Returns true if loaded.
inline Bool_t LoadSNIPParameters(const TString &path, Double_t *m) {
  std::ifstream in(path.Data());
  if (!in.is_open())
    return kFALSE;

  for (Int_t r = 0; r < N_SNIP_REGIONS; r++) {
    Double_t e_center = 0;
    if (!(in >> e_center >> m[r])) {
      in.close();
      return kFALSE;
    }
  }
  in.close();
  return kTRUE;
}

class InteractiveSNIPEditor : public TGMainFrame {
private:
  static const Int_t kSliderRes = 1000;
  static const Int_t kBtnAccept = 1000;
  static const Int_t kBtnCancel = 1001;
  static const Int_t kBtnUseForAll = 1002;
  static const Int_t kSliderBase = 2000;
  static const Int_t kEntryBase = 3000;

  TH1F *hist_;
  TString label_;

  Double_t m_[N_SNIP_REGIONS];
  Double_t m_max_;

  TRootEmbeddedCanvas *embedded_canvas_;
  TPad *spectrum_pad_;
  TPad *subtracted_pad_;

  TH1F *hist_draw_;
  TH1F *background_draw_;
  TH1F *subtracted_draw_;
  TLatex *param_label_;

  TGHSlider *sliders_[N_SNIP_REGIONS];
  TGNumberEntry *entries_[N_SNIP_REGIONS];

  Bool_t needs_redraw_;
  Bool_t accepted_;
  Bool_t use_for_all_;
  Bool_t done_;
  Bool_t syncing_;
  TTimer *redraw_timer_;

  Double_t last_xmin_;
  Double_t last_xmax_;

  void BuildGUI() {
    TGHorizontalFrame *main_frame = new TGHorizontalFrame(this, 1200, 800);
    AddFrame(main_frame,
             new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    embedded_canvas_ =
        new TRootEmbeddedCanvas("SNIPCanvas", main_frame, 800, 750);
    main_frame->AddFrame(
        embedded_canvas_,
        new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    TGVerticalFrame *controls = new TGVerticalFrame(main_frame, 350, 750);
    main_frame->AddFrame(
        controls,
        new TGLayoutHints(kLHintsExpandY | kLHintsRight, 2, 2, 2, 2));

    if (label_.Length() > 0) {
      TGGroupFrame *info_grp =
          new TGGroupFrame(controls, "Dataset", kVerticalFrame);
      controls->AddFrame(info_grp,
                         new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
      TGLabel *info_label = new TGLabel(info_grp, label_);
      info_grp->AddFrame(info_label,
                         new TGLayoutHints(kLHintsCenterX, 4, 4, 4, 4));
    }

    for (Int_t r = 0; r < N_SNIP_REGIONS; r++) {
      Int_t lo = (Int_t)(r * SNIP_REGION_WIDTH);
      Int_t hi = (Int_t)((r + 1) * SNIP_REGION_WIDTH);
      TString title = Form("%d-%d keV", lo, hi);

      TGGroupFrame *grp =
          new TGGroupFrame(controls, title, kVerticalFrame);
      controls->AddFrame(grp,
                         new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 1));

      sliders_[r] =
          new TGHSlider(grp, 200, kSlider1 | kScaleBoth, kSliderBase + r);
      sliders_[r]->SetRange(0, kSliderRes);
      sliders_[r]->SetPosition(ValToSlider(m_[r]));
      sliders_[r]->Associate(this);
      grp->AddFrame(sliders_[r],
                    new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

      TGHorizontalFrame *w_row = new TGHorizontalFrame(grp, 200, 24);
      grp->AddFrame(w_row,
                    new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
      TGLabel *w_label = new TGLabel(w_row, "Window [keV]:");
      w_row->AddFrame(w_label,
                      new TGLayoutHints(kLHintsCenterY, 2, 4, 1, 1));
      entries_[r] = new TGNumberEntry(
          w_row, m_[r], 8, kEntryBase + r, TGNumberFormat::kNESReal,
          TGNumberFormat::kNEAPositive, TGNumberFormat::kNELNoLimits);
      entries_[r]->GetNumberEntry()->Associate(this);
      w_row->AddFrame(
          entries_[r],
          new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 2, 1, 1));
    }

    // Buttons
    TGHorizontalFrame *btn_frame = new TGHorizontalFrame(controls, 270, 40);
    controls->AddFrame(
        btn_frame,
        new TGLayoutHints(kLHintsExpandX | kLHintsBottom, 2, 2, 5, 5));

    TGTextButton *accept_btn =
        new TGTextButton(btn_frame, "Accept", kBtnAccept);
    accept_btn->Associate(this);
    TGTextButton *use_all_btn =
        new TGTextButton(btn_frame, "Use for all", kBtnUseForAll);
    use_all_btn->Associate(this);
    TGTextButton *cancel_btn =
        new TGTextButton(btn_frame, "Cancel", kBtnCancel);
    cancel_btn->Associate(this);

    TGLayoutHints *btn_hints =
        new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2);
    btn_frame->AddFrame(accept_btn, btn_hints);
    btn_frame->AddFrame(use_all_btn, btn_hints);
    btn_frame->AddFrame(cancel_btn, btn_hints);
  }

  void InitDrawing() {
    TCanvas *canvas = embedded_canvas_->GetCanvas();
    canvas->cd();

    spectrum_pad_ = new TPad("snip_spectrum", "snip_spectrum", 0, 0.45, 1, 1.0);
    spectrum_pad_->SetBottomMargin(0.04);
    spectrum_pad_->SetTopMargin(0.08);
    spectrum_pad_->SetGridx(1);
    spectrum_pad_->SetGridy(1);
    spectrum_pad_->Draw();

    subtracted_pad_ =
        new TPad("snip_subtracted", "snip_subtracted", 0, 0, 1, 0.45);
    subtracted_pad_->SetTopMargin(0.04);
    subtracted_pad_->SetBottomMargin(0.18);
    subtracted_pad_->SetGridx(1);
    subtracted_pad_->SetGridy(1);
    subtracted_pad_->Draw();

    spectrum_pad_->cd();
    hist_draw_ = static_cast<TH1F *>(hist_->Clone("snip_editor_hist"));
    hist_draw_->SetDirectory(0);
    hist_draw_->SetMarkerColor(kAzure);
    hist_draw_->SetMarkerStyle(20);
    hist_draw_->SetMarkerSize(0.5);
    hist_draw_->SetLineColor(kAzure);
    hist_draw_->GetXaxis()->SetLabelSize(0);
    hist_draw_->GetXaxis()->SetTitleSize(0);
    hist_draw_->Draw("HIST");

    background_draw_ = ComputeAdaptiveSNIPBackground(hist_, m_);
    background_draw_->SetLineColor(kRed);
    background_draw_->SetLineWidth(2);
    background_draw_->Draw("HIST SAME");

    param_label_ = new TLatex();
    param_label_->SetNDC();
    param_label_->SetTextSize(0.04);
    param_label_->SetTextAlign(31);
    UpdateParamLabel();
    param_label_->Draw();

    subtracted_pad_->cd();
    subtracted_draw_ =
        static_cast<TH1F *>(hist_->Clone("snip_editor_subtracted"));
    subtracted_draw_->SetDirectory(0);
    subtracted_draw_->Add(background_draw_, -1.0);
    subtracted_draw_->SetMarkerColor(kViolet);
    subtracted_draw_->SetMarkerStyle(20);
    subtracted_draw_->SetMarkerSize(0.5);
    subtracted_draw_->SetLineColor(kViolet);
    subtracted_draw_->GetXaxis()->SetTitleSize(0.08);
    subtracted_draw_->GetXaxis()->SetLabelSize(0.07);
    subtracted_draw_->GetYaxis()->SetTitleSize(0.08);
    subtracted_draw_->GetYaxis()->SetLabelSize(0.07);
    subtracted_draw_->GetYaxis()->SetTitleOffset(0.5);
    subtracted_draw_->Draw("HIST");
  }

  void UpdateParamLabel() {
    TString text = "m =";
    for (Int_t r = 0; r < N_SNIP_REGIONS; r++)
      text += Form(" %.1f", m_[r]);
    text += " keV";
    param_label_->SetText(0.88, 0.85, text);
  }

  void UpdateCanvas() {
    delete background_draw_;
    background_draw_ = ComputeAdaptiveSNIPBackground(hist_, m_);
    background_draw_->SetLineColor(kRed);
    background_draw_->SetLineWidth(2);

    spectrum_pad_->cd();
    hist_draw_->Draw("HIST");
    background_draw_->Draw("HIST SAME");
    UpdateParamLabel();
    param_label_->Draw();

    subtracted_pad_->cd();
    subtracted_draw_->Reset();
    subtracted_draw_->Add(hist_);
    subtracted_draw_->Add(background_draw_, -1.0);
    subtracted_draw_->Draw("HIST");

    TCanvas *canvas = embedded_canvas_->GetCanvas();
    spectrum_pad_->Modified();
    subtracted_pad_->Modified();
    canvas->Modified();
    canvas->Update();

    needs_redraw_ = kFALSE;
  }

  // Enforce monotonically increasing windows by clamping the changed region.
  // Can't go below the region below or above the region above.
  void EnforceMonotonic(Int_t r) {
    // Clamp to floor from region below
    if (r > 0 && m_[r] < m_[r - 1])
      m_[r] = m_[r - 1];
    // Clamp to ceiling from region above
    if (r < N_SNIP_REGIONS - 1 && m_[r] > m_[r + 1])
      m_[r] = m_[r + 1];

    syncing_ = kTRUE;
    sliders_[r]->SetPosition(ValToSlider(m_[r]));
    entries_[r]->SetNumber(m_[r]);
    syncing_ = kFALSE;
  }

  void SyncZoom() {
    Double_t spec_xmin = hist_draw_->GetXaxis()->GetFirst();
    Double_t spec_xmax = hist_draw_->GetXaxis()->GetLast();
    Double_t spec_lo = hist_draw_->GetXaxis()->GetBinLowEdge((Int_t)spec_xmin);
    Double_t spec_hi =
        hist_draw_->GetXaxis()->GetBinUpEdge((Int_t)spec_xmax);

    Double_t sub_xmin = subtracted_draw_->GetXaxis()->GetFirst();
    Double_t sub_xmax = subtracted_draw_->GetXaxis()->GetLast();
    Double_t sub_lo =
        subtracted_draw_->GetXaxis()->GetBinLowEdge((Int_t)sub_xmin);
    Double_t sub_hi =
        subtracted_draw_->GetXaxis()->GetBinUpEdge((Int_t)sub_xmax);

    Bool_t spec_changed =
        (TMath::Abs(spec_lo - last_xmin_) > 0.01 ||
         TMath::Abs(spec_hi - last_xmax_) > 0.01);
    Bool_t sub_changed =
        (TMath::Abs(sub_lo - last_xmin_) > 0.01 ||
         TMath::Abs(sub_hi - last_xmax_) > 0.01);

    if (spec_changed) {
      last_xmin_ = spec_lo;
      last_xmax_ = spec_hi;
      subtracted_pad_->cd();
      subtracted_draw_->GetXaxis()->SetRangeUser(spec_lo, spec_hi);
      subtracted_pad_->Modified();
    } else if (sub_changed) {
      last_xmin_ = sub_lo;
      last_xmax_ = sub_hi;
      spectrum_pad_->cd();
      hist_draw_->GetXaxis()->SetRangeUser(sub_lo, sub_hi);
      background_draw_->GetXaxis()->SetRangeUser(sub_lo, sub_hi);
      spectrum_pad_->Modified();
    }

    if (spec_changed || sub_changed) {
      TCanvas *canvas = embedded_canvas_->GetCanvas();
      canvas->Modified();
      canvas->Update();
    }
  }

  Int_t ValToSlider(Double_t val) {
    if (m_max_ <= 0)
      return 0;
    Double_t frac = val / m_max_;
    if (frac < 0.0)
      frac = 0.0;
    if (frac > 1.0)
      frac = 1.0;
    return static_cast<Int_t>(frac * kSliderRes);
  }

  Double_t SliderToVal(Int_t pos) {
    return (static_cast<Double_t>(pos) / kSliderRes) * m_max_;
  }

public:
  InteractiveSNIPEditor(const TGWindow *parent, TH1F *hist, const Double_t *m,
                        Double_t m_max, const TString &label = "")
      : TGMainFrame(parent, 1200, 800) {

    hist_ = hist;
    label_ = label;
    for (Int_t r = 0; r < N_SNIP_REGIONS; r++)
      m_[r] = m[r];
    m_max_ = m_max;
    accepted_ = kFALSE;
    use_for_all_ = kFALSE;
    done_ = kFALSE;
    needs_redraw_ = kFALSE;
    syncing_ = kFALSE;
    last_xmin_ = hist->GetXaxis()->GetXmin();
    last_xmax_ = hist->GetXaxis()->GetXmax();

    BuildGUI();
    InitDrawing();

    redraw_timer_ = new TTimer(this, 100);
    redraw_timer_->TurnOn();

    SetWindowName("SNIP Background Editor");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

    UpdateCanvas();
  }

  virtual ~InteractiveSNIPEditor() {
    if (redraw_timer_) {
      redraw_timer_->TurnOff();
      delete redraw_timer_;
    }
  }

  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t /*parm2*/) {
    if (syncing_)
      return kTRUE;

    switch (GET_MSG(msg)) {
    case kC_COMMAND:
      switch (GET_SUBMSG(msg)) {
      case kCM_BUTTON:
        if (parm1 == kBtnAccept) {
          accepted_ = kTRUE;
          done_ = kTRUE;
        } else if (parm1 == kBtnUseForAll) {
          accepted_ = kTRUE;
          use_for_all_ = kTRUE;
          done_ = kTRUE;
        } else if (parm1 == kBtnCancel) {
          accepted_ = kFALSE;
          done_ = kTRUE;
        }
        break;
      }
      break;

    case kC_HSLIDER: {
      Int_t idx = parm1 - kSliderBase;
      if (idx >= 0 && idx < N_SNIP_REGIONS) {
        m_[idx] = SliderToVal(sliders_[idx]->GetPosition());
        EnforceMonotonic(idx);
        needs_redraw_ = kTRUE;
      }
    } break;

    case kC_TEXTENTRY:
      if (GET_SUBMSG(msg) == kTE_ENTER || GET_SUBMSG(msg) == kTE_TAB) {
        Int_t idx = parm1 - kEntryBase;
        if (idx >= 0 && idx < N_SNIP_REGIONS) {
          m_[idx] = entries_[idx]->GetNumber();
          EnforceMonotonic(idx);
          needs_redraw_ = kTRUE;
        }
      }
      break;
    }

    return kTRUE;
  }

  virtual Bool_t HandleTimer(TTimer *timer) {
    if (timer == redraw_timer_) {
      SyncZoom();
      if (needs_redraw_)
        UpdateCanvas();
    }
    return kTRUE;
  }

  virtual void CloseWindow() {
    accepted_ = kFALSE;
    done_ = kTRUE;
  }

  Bool_t IsDone() const { return done_; }

  SNIPParameters GetResult() const {
    SNIPParameters p;
    for (Int_t r = 0; r < N_SNIP_REGIONS; r++)
      p.m[r] = m_[r];
    p.accepted = accepted_;
    p.use_for_all = use_for_all_;
    return p;
  }

  TTimer *GetRedrawTimer() { return redraw_timer_; }
  TRootEmbeddedCanvas *GetEmbeddedCanvas() { return embedded_canvas_; }
};

// Blocking launcher
inline SNIPParameters LaunchSNIPEditor(TH1F *hist, const Double_t *m,
                                       Double_t m_max,
                                       const TString &label = "") {
  SNIPParameters result;
  result.accepted = kFALSE;
  result.use_for_all = kFALSE;

  if (!gClient) {
    std::cerr << "SNIPEditor: GUI not available (gClient is null)."
              << std::endl;
    return result;
  }

  InteractiveSNIPEditor *editor =
      new InteractiveSNIPEditor(gClient->GetRoot(), hist, m, m_max, label);

  while (!editor->IsDone()) {
    gSystem->ProcessEvents();
    gSystem->Sleep(10);
  }

  result = editor->GetResult();

  editor->GetRedrawTimer()->TurnOff();

  gSystem->ProcessEvents();
  gSystem->Sleep(50);
  gSystem->ProcessEvents();

  TCanvas *ecanvas = editor->GetEmbeddedCanvas()->GetCanvas();
  if (ecanvas) {
    ecanvas->Clear();
    ecanvas->SetBatch(kTRUE);
    gROOT->GetListOfCanvases()->Remove(ecanvas);
  }

  editor->DontCallClose();
  editor->UnmapWindow();

  gSystem->ProcessEvents();
  gSystem->Sleep(50);
  gSystem->ProcessEvents();

  delete editor;

  return result;
}

#endif
