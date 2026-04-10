#ifndef INTERACTIVESUBTRACTIONEDITOR_HPP
#define INTERACTIVESUBTRACTIONEDITOR_HPP

#include <TBox.h>
#include <TCanvas.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGDoubleSlider.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGSlider.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMinuit.h>
#include <TPad.h>
#include <TROOT.h>
#include <TRootEmbeddedCanvas.h>
#include <TSystem.h>
#include <TTimer.h>
#include <cmath>
#include <iostream>

// Globals for TMinuit FCN
static TH1F *gSubEditorSignal = nullptr;
static TH1F *gSubEditorBackground = nullptr;
static Double_t gSubEditorFitLo = 0;
static Double_t gSubEditorFitHi = 0;

static void SubtractionFCN(Int_t & /*npar*/, Double_t * /*grad*/,
                           Double_t &fval, Double_t *par, Int_t /*iflag*/) {
  Double_t A = par[0];
  Double_t sum = 0;
  Int_t nbins = gSubEditorSignal->GetNbinsX();
  for (Int_t bin = 1; bin <= nbins; bin++) {
    Double_t x = gSubEditorSignal->GetBinCenter(bin);
    if (x < gSubEditorFitLo || x > gSubEditorFitHi)
      continue;
    Double_t s = gSubEditorSignal->GetBinContent(bin);
    Double_t b = gSubEditorBackground->GetBinContent(bin);
    Double_t residual = s - A * b;
    sum += residual * residual;
  }
  fval = sum;
}

class InteractiveSubtractionEditor : public TGMainFrame {
private:
  static const Int_t kSliderRes = 10000;

  static const Int_t kBtnFit = 1000;
  static const Int_t kBtnAccept = 1001;
  static const Int_t kBtnCancel = 1002;
  static const Int_t kLockCheck = 1003;
  static const Int_t kScaleSlider = 2000;
  static const Int_t kScaleEntry = 3000;
  static const Int_t kRangeSlider = 4000;
  static const Int_t kRangeLoEntry = 4001;
  static const Int_t kRangeHiEntry = 4002;

  TH1F *signal_hist_;
  TH1F *background_hist_;
  TString label_;

  Double_t scale_;
  Double_t scale_lo_;
  Double_t scale_hi_;
  Double_t initial_scale_;

  Double_t fit_range_lo_;
  Double_t fit_range_hi_;
  Double_t hist_x_min_;
  Double_t hist_x_max_;

  Double_t overlay_ymax_;
  Double_t subtraction_ymin_;
  Double_t subtraction_ymax_;

  TRootEmbeddedCanvas *embedded_canvas_;
  TPad *overlay_pad_;
  TPad *subtraction_pad_;

  TH1F *signal_draw_;
  TH1F *background_draw_;
  TH1F *subtraction_draw_;
  TLatex *scale_label_;
  TLine *fit_lo_line_overlay_;
  TLine *fit_hi_line_overlay_;
  TLine *fit_lo_line_sub_;
  TLine *fit_hi_line_sub_;

  TGHSlider *slider_;
  TGNumberEntry *entry_;
  TGCheckButton *lock_check_;
  TGDoubleHSlider *range_slider_;
  TGNumberEntry *range_lo_entry_;
  TGNumberEntry *range_hi_entry_;

  Bool_t needs_redraw_;
  Bool_t accepted_;
  Bool_t done_;
  Bool_t syncing_;
  TTimer *redraw_timer_;

  void BuildGUI() {
    TGHorizontalFrame *main_frame = new TGHorizontalFrame(this, 1000, 700);
    AddFrame(main_frame,
             new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    embedded_canvas_ =
        new TRootEmbeddedCanvas("SubEditorCanvas", main_frame, 750, 650);
    main_frame->AddFrame(
        embedded_canvas_,
        new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    TGVerticalFrame *controls = new TGVerticalFrame(main_frame, 280, 650);
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

    // Scale factor controls
    TGGroupFrame *scale_grp =
        new TGGroupFrame(controls, "Scale Factor", kVerticalFrame);
    controls->AddFrame(scale_grp,
                       new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));

    slider_ = new TGHSlider(scale_grp, 200, kSlider1, kScaleSlider);
    slider_->SetRange(0, kSliderRes);
    slider_->SetPosition(ValToSlider(scale_));
    slider_->Associate(this);
    scale_grp->AddFrame(slider_,
                        new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    TGHorizontalFrame *entry_row = new TGHorizontalFrame(scale_grp, 200, 28);
    scale_grp->AddFrame(entry_row,
                        new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    TGLabel *val_label = new TGLabel(entry_row, "Value:");
    entry_row->AddFrame(val_label,
                        new TGLayoutHints(kLHintsCenterY, 2, 4, 2, 2));
    entry_ = new TGNumberEntry(entry_row, scale_, 10, kScaleEntry,
                               TGNumberFormat::kNESReal,
                               TGNumberFormat::kNEAAnyNumber,
                               TGNumberFormat::kNELNoLimits);
    entry_->GetNumberEntry()->Associate(this);
    entry_row->AddFrame(
        entry_,
        new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 2, 2, 2));

    lock_check_ = new TGCheckButton(scale_grp, "Lock (skip fit)", kLockCheck);
    lock_check_->Associate(this);
    scale_grp->AddFrame(lock_check_,
                        new TGLayoutHints(kLHintsLeft, 4, 2, 2, 2));

    // Fit range controls
    TGGroupFrame *range_grp =
        new TGGroupFrame(controls, "Fit Range", kVerticalFrame);
    controls->AddFrame(range_grp,
                       new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));

    range_slider_ =
        new TGDoubleHSlider(range_grp, 200, kDoubleScaleNo, kRangeSlider);
    range_slider_->SetRange(hist_x_min_, hist_x_max_);
    range_slider_->SetPosition(fit_range_lo_, fit_range_hi_);
    range_slider_->Associate(this);
    range_grp->AddFrame(range_slider_,
                        new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    TGHorizontalFrame *range_entries =
        new TGHorizontalFrame(range_grp, 200, 28);
    range_grp->AddFrame(range_entries,
                        new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    TGLabel *range_lo_label = new TGLabel(range_entries, "Low:");
    range_entries->AddFrame(range_lo_label,
                            new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
    range_lo_entry_ = new TGNumberEntry(
        range_entries, fit_range_lo_, 8, kRangeLoEntry,
        TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
        TGNumberFormat::kNELNoLimits);
    range_lo_entry_->GetNumberEntry()->Associate(this);
    range_entries->AddFrame(
        range_lo_entry_,
        new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 2, 2, 2));

    TGLabel *range_hi_label = new TGLabel(range_entries, "High:");
    range_entries->AddFrame(range_hi_label,
                            new TGLayoutHints(kLHintsCenterY, 8, 2, 2, 2));
    range_hi_entry_ = new TGNumberEntry(
        range_entries, fit_range_hi_, 8, kRangeHiEntry,
        TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
        TGNumberFormat::kNELNoLimits);
    range_hi_entry_->GetNumberEntry()->Associate(this);
    range_entries->AddFrame(
        range_hi_entry_,
        new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 2, 2, 2));

    // Buttons
    TGHorizontalFrame *btn_frame = new TGHorizontalFrame(controls, 270, 40);
    controls->AddFrame(
        btn_frame,
        new TGLayoutHints(kLHintsExpandX | kLHintsBottom, 2, 2, 5, 5));

    TGTextButton *fit_btn = new TGTextButton(btn_frame, "Fit", kBtnFit);
    fit_btn->Associate(this);
    TGTextButton *accept_btn =
        new TGTextButton(btn_frame, "Accept", kBtnAccept);
    accept_btn->Associate(this);
    TGTextButton *cancel_btn =
        new TGTextButton(btn_frame, "Cancel", kBtnCancel);
    cancel_btn->Associate(this);

    TGLayoutHints *btn_hints =
        new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2);
    btn_frame->AddFrame(fit_btn, btn_hints);
    btn_frame->AddFrame(accept_btn, btn_hints);
    btn_frame->AddFrame(cancel_btn, btn_hints);
  }

  void ComputeFixedYRanges() {
    Double_t max_sig = signal_hist_->GetMaximum();
    Double_t max_bkg_scaled = background_hist_->GetMaximum() * scale_hi_;
    overlay_ymax_ =
        1.15 * (max_sig > max_bkg_scaled ? max_sig : max_bkg_scaled);

    Double_t sub_min = 1e30;
    Double_t sub_max = -1e30;
    Int_t nbins = signal_hist_->GetNbinsX();
    for (Int_t bin = 1; bin <= nbins; bin++) {
      Double_t s = signal_hist_->GetBinContent(bin);
      Double_t b = background_hist_->GetBinContent(bin);
      Double_t val_lo = s - scale_hi_ * b;
      Double_t val_hi = s - scale_lo_ * b;
      if (val_lo < sub_min)
        sub_min = val_lo;
      if (val_hi > sub_max)
        sub_max = val_hi;
      if (val_lo > sub_max)
        sub_max = val_lo;
      if (val_hi < sub_min)
        sub_min = val_hi;
    }
    Double_t sub_range = sub_max - sub_min;
    subtraction_ymin_ = sub_min - 0.1 * sub_range;
    subtraction_ymax_ = sub_max + 0.1 * sub_range;
  }

  void InitDrawing() {
    ComputeFixedYRanges();

    TCanvas *canvas = embedded_canvas_->GetCanvas();
    canvas->cd();

    overlay_pad_ = new TPad("sub_overlay", "sub_overlay", 0, 0.45, 1, 1.0);
    overlay_pad_->SetBottomMargin(0.04);
    overlay_pad_->SetTopMargin(0.08);
    overlay_pad_->SetGridx(1);
    overlay_pad_->SetGridy(1);
    overlay_pad_->Draw();

    subtraction_pad_ = new TPad("sub_result", "sub_result", 0, 0, 1, 0.45);
    subtraction_pad_->SetTopMargin(0.04);
    subtraction_pad_->SetBottomMargin(0.18);
    subtraction_pad_->SetGridx(1);
    subtraction_pad_->SetGridy(1);
    subtraction_pad_->Draw();

    // Overlay pad
    overlay_pad_->cd();

    signal_draw_ =
        static_cast<TH1F *>(signal_hist_->Clone("sub_editor_signal"));
    signal_draw_->SetDirectory(0);
    signal_draw_->SetMarkerColor(kAzure);
    signal_draw_->SetMarkerStyle(20);
    signal_draw_->SetMarkerSize(0.65);
    signal_draw_->SetLineColor(kAzure);
    signal_draw_->GetXaxis()->SetLabelSize(0);
    signal_draw_->GetXaxis()->SetTitleSize(0);
    signal_draw_->SetMaximum(overlay_ymax_);
    signal_draw_->SetMinimum(0);
    signal_draw_->Draw("P");

    background_draw_ =
        static_cast<TH1F *>(background_hist_->Clone("sub_editor_bkg"));
    background_draw_->SetDirectory(0);
    background_draw_->Scale(scale_);
    background_draw_->SetMarkerColor(kRed);
    background_draw_->SetMarkerStyle(20);
    background_draw_->SetMarkerSize(0.65);
    background_draw_->SetLineColor(kRed);
    background_draw_->Draw("P SAME");

    fit_lo_line_overlay_ =
        new TLine(fit_range_lo_, 0, fit_range_lo_, overlay_ymax_);
    fit_lo_line_overlay_->SetLineColor(kGreen + 2);
    fit_lo_line_overlay_->SetLineStyle(2);
    fit_lo_line_overlay_->SetLineWidth(2);
    fit_lo_line_overlay_->Draw();

    fit_hi_line_overlay_ =
        new TLine(fit_range_hi_, 0, fit_range_hi_, overlay_ymax_);
    fit_hi_line_overlay_->SetLineColor(kGreen + 2);
    fit_hi_line_overlay_->SetLineStyle(2);
    fit_hi_line_overlay_->SetLineWidth(2);
    fit_hi_line_overlay_->Draw();

    scale_label_ = new TLatex();
    scale_label_->SetNDC();
    scale_label_->SetTextSize(0.055);
    scale_label_->SetTextAlign(31);
    scale_label_->SetText(0.88, 0.85, Form("Scale = %.4f", scale_));
    scale_label_->Draw();

    // Subtraction pad
    subtraction_pad_->cd();

    subtraction_draw_ =
        static_cast<TH1F *>(signal_hist_->Clone("sub_editor_result"));
    subtraction_draw_->SetDirectory(0);
    subtraction_draw_->Add(background_hist_, -scale_);
    subtraction_draw_->SetMarkerColor(kViolet);
    subtraction_draw_->SetMarkerStyle(20);
    subtraction_draw_->SetMarkerSize(0.65);
    subtraction_draw_->SetLineColor(kViolet);
    subtraction_draw_->GetXaxis()->SetTitleSize(0.08);
    subtraction_draw_->GetXaxis()->SetLabelSize(0.07);
    subtraction_draw_->GetYaxis()->SetTitleSize(0.08);
    subtraction_draw_->GetYaxis()->SetLabelSize(0.07);
    subtraction_draw_->GetYaxis()->SetTitleOffset(0.5);
    subtraction_draw_->SetMaximum(subtraction_ymax_);
    subtraction_draw_->SetMinimum(subtraction_ymin_);
    subtraction_draw_->Draw("P");

    fit_lo_line_sub_ =
        new TLine(fit_range_lo_, subtraction_ymin_, fit_range_lo_,
                  subtraction_ymax_);
    fit_lo_line_sub_->SetLineColor(kGreen + 2);
    fit_lo_line_sub_->SetLineStyle(2);
    fit_lo_line_sub_->SetLineWidth(2);
    fit_lo_line_sub_->Draw();

    fit_hi_line_sub_ =
        new TLine(fit_range_hi_, subtraction_ymin_, fit_range_hi_,
                  subtraction_ymax_);
    fit_hi_line_sub_->SetLineColor(kGreen + 2);
    fit_hi_line_sub_->SetLineStyle(2);
    fit_hi_line_sub_->SetLineWidth(2);
    fit_hi_line_sub_->Draw();
  }

  void UpdateCanvas() {
    overlay_pad_->cd();
    background_draw_->Reset();
    background_draw_->Add(background_hist_, scale_);
    signal_draw_->SetMaximum(overlay_ymax_);
    signal_draw_->SetMinimum(0);

    fit_lo_line_overlay_->SetX1(fit_range_lo_);
    fit_lo_line_overlay_->SetX2(fit_range_lo_);
    fit_lo_line_overlay_->SetY2(overlay_ymax_);
    fit_hi_line_overlay_->SetX1(fit_range_hi_);
    fit_hi_line_overlay_->SetX2(fit_range_hi_);
    fit_hi_line_overlay_->SetY2(overlay_ymax_);

    scale_label_->SetText(0.88, 0.85, Form("Scale = %.4f", scale_));

    subtraction_pad_->cd();
    subtraction_draw_->Reset();
    subtraction_draw_->Add(signal_hist_);
    subtraction_draw_->Add(background_hist_, -scale_);
    subtraction_draw_->SetMaximum(subtraction_ymax_);
    subtraction_draw_->SetMinimum(subtraction_ymin_);

    fit_lo_line_sub_->SetX1(fit_range_lo_);
    fit_lo_line_sub_->SetX2(fit_range_lo_);
    fit_lo_line_sub_->SetY1(subtraction_ymin_);
    fit_lo_line_sub_->SetY2(subtraction_ymax_);
    fit_hi_line_sub_->SetX1(fit_range_hi_);
    fit_hi_line_sub_->SetX2(fit_range_hi_);
    fit_hi_line_sub_->SetY1(subtraction_ymin_);
    fit_hi_line_sub_->SetY2(subtraction_ymax_);

    TCanvas *canvas = embedded_canvas_->GetCanvas();
    overlay_pad_->Modified();
    subtraction_pad_->Modified();
    canvas->Modified();
    canvas->Update();

    needs_redraw_ = kFALSE;
  }

  void DoFit() {
    if (lock_check_->GetState() == kButtonDown) {
      std::cout << "  Scale locked at " << scale_ << ", skipping fit"
                << std::endl;
      return;
    }

    gSubEditorSignal = signal_hist_;
    gSubEditorBackground = background_hist_;
    gSubEditorFitLo = fit_range_lo_;
    gSubEditorFitHi = fit_range_hi_;

    TMinuit minuit(1);
    minuit.SetPrintLevel(-1);
    minuit.SetFCN(SubtractionFCN);
    minuit.DefineParameter(0, "A", scale_, scale_ * 0.01, 0, 0);
    minuit.Migrad();

    Double_t fitted_scale = 0;
    Double_t fitted_error = 0;
    minuit.GetParameter(0, fitted_scale, fitted_error);

    gSubEditorSignal = nullptr;
    gSubEditorBackground = nullptr;

    std::cout << "  Fit result: scale = " << fitted_scale << " +/- "
              << fitted_error << " (range: [" << fit_range_lo_ << ", "
              << fit_range_hi_ << "] keV)" << std::endl;

    scale_ = fitted_scale;

    if (scale_ < scale_lo_)
      scale_lo_ = scale_;
    if (scale_ > scale_hi_)
      scale_hi_ = scale_;

    syncing_ = kTRUE;
    slider_->SetPosition(ValToSlider(scale_));
    entry_->SetNumber(scale_);
    syncing_ = kFALSE;

    needs_redraw_ = kTRUE;
    UpdateCanvas();
  }

  void OnRangeChanged() {
    if (syncing_)
      return;

    Double_t new_lo = range_slider_->GetMinPosition();
    Double_t new_hi = range_slider_->GetMaxPosition();
    if (new_lo >= new_hi)
      return;

    fit_range_lo_ = new_lo;
    fit_range_hi_ = new_hi;

    syncing_ = kTRUE;
    range_lo_entry_->SetNumber(fit_range_lo_);
    range_hi_entry_->SetNumber(fit_range_hi_);
    syncing_ = kFALSE;

    needs_redraw_ = kTRUE;
  }

  Int_t ValToSlider(Double_t val) {
    if (scale_hi_ <= scale_lo_)
      return 0;
    Double_t frac = (val - scale_lo_) / (scale_hi_ - scale_lo_);
    if (frac < 0.0)
      frac = 0.0;
    if (frac > 1.0)
      frac = 1.0;
    return static_cast<Int_t>(frac * kSliderRes);
  }

  Double_t SliderToVal(Int_t pos) {
    Double_t frac = static_cast<Double_t>(pos) / kSliderRes;
    return scale_lo_ + frac * (scale_hi_ - scale_lo_);
  }

public:
  InteractiveSubtractionEditor(const TGWindow *parent, TH1F *signal,
                               TH1F *background, Double_t initial_scale,
                               Double_t fit_lo, Double_t fit_hi,
                               const TString &label = "")
      : TGMainFrame(parent, 1000, 700) {

    signal_hist_ = signal;
    background_hist_ = background;
    label_ = label;

    scale_ = initial_scale;
    initial_scale_ = initial_scale;
    scale_lo_ = 0.0;
    scale_hi_ = initial_scale * 2.0;

    fit_range_lo_ = fit_lo;
    fit_range_hi_ = fit_hi;
    hist_x_min_ = signal->GetXaxis()->GetXmin();
    hist_x_max_ = signal->GetXaxis()->GetXmax();

    accepted_ = kFALSE;
    done_ = kFALSE;
    needs_redraw_ = kFALSE;
    syncing_ = kFALSE;

    BuildGUI();
    InitDrawing();

    redraw_timer_ = new TTimer(this, 50);
    redraw_timer_->TurnOn();

    SetWindowName("Background Subtraction Editor");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

    UpdateCanvas();
  }

  virtual ~InteractiveSubtractionEditor() {
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
        if (parm1 == kBtnFit)
          DoFit();
        else if (parm1 == kBtnAccept) {
          accepted_ = kTRUE;
          done_ = kTRUE;
        } else if (parm1 == kBtnCancel) {
          scale_ = initial_scale_;
          accepted_ = kFALSE;
          done_ = kTRUE;
        }
        break;
      }
      break;

    case kC_HSLIDER:
      if (parm1 == kScaleSlider) {
        Int_t pos = slider_->GetPosition();
        scale_ = SliderToVal(pos);
        syncing_ = kTRUE;
        entry_->SetNumber(scale_);
        syncing_ = kFALSE;
        needs_redraw_ = kTRUE;
      } else if (parm1 == kRangeSlider) {
        OnRangeChanged();
      }
      break;

    case kC_TEXTENTRY:
      if (GET_SUBMSG(msg) == kTE_ENTER || GET_SUBMSG(msg) == kTE_TAB) {
        if (parm1 == kScaleEntry) {
          scale_ = entry_->GetNumber();
          if (scale_ < scale_lo_)
            scale_lo_ = scale_;
          if (scale_ > scale_hi_)
            scale_hi_ = scale_;
          syncing_ = kTRUE;
          slider_->SetPosition(ValToSlider(scale_));
          syncing_ = kFALSE;
          needs_redraw_ = kTRUE;
        } else if (parm1 == kRangeLoEntry || parm1 == kRangeHiEntry) {
          Double_t new_lo = range_lo_entry_->GetNumber();
          Double_t new_hi = range_hi_entry_->GetNumber();
          if (new_lo < new_hi) {
            fit_range_lo_ = new_lo;
            fit_range_hi_ = new_hi;
            syncing_ = kTRUE;
            range_slider_->SetPosition(fit_range_lo_, fit_range_hi_);
            syncing_ = kFALSE;
            needs_redraw_ = kTRUE;
          }
        }
      }
      break;
    }

    return kTRUE;
  }

  virtual Bool_t HandleTimer(TTimer *timer) {
    if (timer == redraw_timer_ && needs_redraw_) {
      UpdateCanvas();
    }
    return kTRUE;
  }

  virtual void CloseWindow() {
    scale_ = initial_scale_;
    accepted_ = kFALSE;
    done_ = kTRUE;
  }

  Bool_t WasAccepted() const { return accepted_; }
  Bool_t IsDone() const { return done_; }
  Double_t GetScale() const { return scale_; }
  TTimer *GetRedrawTimer() { return redraw_timer_; }
  TRootEmbeddedCanvas *GetEmbeddedCanvas() { return embedded_canvas_; }
};

// Blocking launcher — returns accepted scale, or -1 on cancel
inline Double_t LaunchSubtractionEditor(TH1F *signal, TH1F *background,
                                        Double_t initial_scale,
                                        Double_t fit_lo, Double_t fit_hi,
                                        const TString &label = "") {
  if (!gClient) {
    std::cerr << "SubtractionEditor: GUI not available (gClient is null)."
              << std::endl;
    return -1;
  }

  InteractiveSubtractionEditor *editor = new InteractiveSubtractionEditor(
      gClient->GetRoot(), signal, background, initial_scale, fit_lo, fit_hi,
      label);

  while (!editor->IsDone()) {
    gSystem->ProcessEvents();
    gSystem->Sleep(10);
  }

  Double_t result_scale = editor->WasAccepted() ? editor->GetScale() : -1.0;

  editor->GetRedrawTimer()->TurnOff();

  TCanvas *ecanvas = editor->GetEmbeddedCanvas()->GetCanvas();
  if (ecanvas) {
    gROOT->GetListOfCanvases()->Remove(ecanvas);
    ecanvas->SetBatch(kTRUE);
  }

  editor->DontCallClose();
  editor->UnmapWindow();
  gSystem->ProcessEvents();
  delete editor;

  return result_scale;
}

#endif
