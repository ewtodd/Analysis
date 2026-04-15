#ifndef HYPEREMGFITHELPERS_HPP
#define HYPEREMGFITHELPERS_HPP

// All hyper-EMG fit helpers are now in the FittingFunctions namespace upstream.
// This header provides unqualified names for backward compatibility with macros.

#include "FittingUtils.hpp"

using FittingFunctions::ExtractHyperEMGPeak;
using FittingFunctions::FitDoubleHyperEMG;
using FittingFunctions::FitSingleHyperEMG;
using FittingFunctions::FitTripleHyperEMG;
using FittingFunctions::LoadHyperEMGParams;
using FittingFunctions::PlotHyperEMGFit;
using FittingFunctions::SaveHyperEMGParams;
using FittingFunctions::SetupHyperEMGBackground;

#endif
