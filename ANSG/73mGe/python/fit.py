from analysis_utils import load_cpp_library
from analysis_utils.init import set_root_preferences
from analysis_utils.io import load_tree_data
from analysis_utils.fitting import single_peak_pdf, estimate_peak_params
from iminuit import Minuit, cost
import numpy as np
import constants as C
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = load_cpp_library()


def test_fit(filename):
    filepath = f"{C.ROOT_FILES_DIR}/{filename}.root"

    tree_name = "bef_tree_event_summary"
    energy_branch = "totalEnergykeV"

    df = load_tree_data(filepath, tree_name=tree_name)
    energy = df[energy_branch]

    pdf = single_peak_pdf(45, 75)
    c = cost.UnbinnedNLL(energy, pdf)
    values, limits, fixed = estimate_peak_params(energy, 45, 75, mu=59.5)

    m = Minuit(c, **values)
    for k, v in limits.items():
        m.limits[k] = v
    for k, v in fixed.items():
        m.fixed[k] = v

    # Initial fit: Gaussian + background only
    m.migrad()
    best_nll = m.fval
    gaus_amp = m.values["gaus_amp"]
    print(f"Initial fit: NLL = {best_nll:.2f}")

    # Low-side group test: enable step + low tails together
    print("Testing low-side group (step + low exp tail + low lin tail)...")
    m.fixed["step_amp"] = False
    m.limits["step_amp"] = (0, gaus_amp * 2)
    m.values["step_amp"] = gaus_amp * 0.1

    m.fixed["low_exp_amp"] = False
    m.fixed["low_exp_decay"] = False
    m.limits["low_exp_amp"] = (0, gaus_amp * 2)
    m.limits["low_exp_decay"] = (0.1, 50)
    m.values["low_exp_amp"] = gaus_amp * 0.15
    m.values["low_exp_decay"] = 1.0

    m.fixed["low_lin_amp"] = False
    m.fixed["low_lin_slope"] = False
    m.limits["low_lin_amp"] = (0, gaus_amp * 2)
    m.limits["low_lin_slope"] = (-1, 1)
    m.values["low_lin_amp"] = gaus_amp * 0.15
    m.values["low_lin_slope"] = 0.0

    m.migrad()
    print(
        f"  Low-side group NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
    )

    if m.fval < best_nll - 1:
        print("  Low-side group ACCEPTED — pruning individual components...")
        best_nll = m.fval

        # Try without step
        saved = m.values["step_amp"]
        m.fixed["step_amp"] = True
        m.values["step_amp"] = 0.0
        m.migrad()
        print(
            f"  Without step: NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
        )
        if m.fval < best_nll - 1:
            best_nll = m.fval
            print("  Step PRUNED")
        else:
            m.fixed["step_amp"] = False
            m.values["step_amp"] = saved
            m.migrad()
            print("  Step KEPT")

        # Try without low exp tail
        saved_amp = m.values["low_exp_amp"]
        saved_decay = m.values["low_exp_decay"]
        m.fixed["low_exp_amp"] = True
        m.fixed["low_exp_decay"] = True
        m.values["low_exp_amp"] = 0.0
        m.values["low_exp_decay"] = 1.0
        m.migrad()
        print(
            f"  Without low exp tail: NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
        )
        if m.fval < best_nll - 1:
            best_nll = m.fval
            print("  Low exp tail PRUNED")
        else:
            m.fixed["low_exp_amp"] = False
            m.fixed["low_exp_decay"] = False
            m.values["low_exp_amp"] = saved_amp
            m.values["low_exp_decay"] = saved_decay
            m.migrad()
            print("  Low exp tail KEPT")

        # Try without low lin tail
        saved_amp = m.values["low_lin_amp"]
        saved_slope = m.values["low_lin_slope"]
        m.fixed["low_lin_amp"] = True
        m.fixed["low_lin_slope"] = True
        m.values["low_lin_amp"] = 0.0
        m.values["low_lin_slope"] = 0.0
        m.migrad()
        print(
            f"  Without low lin tail: NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})"
        )
        if m.fval < best_nll - 1:
            best_nll = m.fval
            print("  Low lin tail PRUNED")
        else:
            m.fixed["low_lin_amp"] = False
            m.fixed["low_lin_slope"] = False
            m.values["low_lin_amp"] = saved_amp
            m.values["low_lin_slope"] = saved_slope
            m.migrad()
            print("  Low lin tail KEPT")
    else:
        print("  Low-side group REJECTED — re-fixing all")
        m.fixed["step_amp"] = True
        m.values["step_amp"] = 0.0
        m.fixed["low_exp_amp"] = True
        m.fixed["low_exp_decay"] = True
        m.values["low_exp_amp"] = 0.0
        m.values["low_exp_decay"] = 1.0
        m.fixed["low_lin_amp"] = True
        m.fixed["low_lin_slope"] = True
        m.values["low_lin_amp"] = 0.0
        m.values["low_lin_slope"] = 0.0

    # High tail test (independent of low side)
    print("Testing high exponential tail...")
    m.fixed["high_exp_amp"] = False
    m.fixed["high_exp_decay"] = False
    m.limits["high_exp_amp"] = (0, gaus_amp * 2)
    m.limits["high_exp_decay"] = (0.1, 50)
    m.values["high_exp_amp"] = gaus_amp * 0.15
    m.values["high_exp_decay"] = 1.0
    m.migrad()
    print(f"  High tail NLL = {m.fval:.2f} (delta = {m.fval - best_nll:.2f})")
    if m.fval < best_nll - 1:
        best_nll = m.fval
        print("  High tail ACCEPTED")
    else:
        m.fixed["high_exp_amp"] = True
        m.fixed["high_exp_decay"] = True
        m.values["high_exp_amp"] = 0.0
        m.values["high_exp_decay"] = 1.0
        print("  High tail REJECTED")

    # Final fit + errors
    print(f"Final fit with selected components...")
    m.migrad()
    m.hesse()
    print(f"Final NLL = {m.fval:.2f}")
    print(m)
    fig, ax = plt.subplots()
    ax.hist(energy, bins=200, range=(40, 100), histtype="step", density=True)
    x = np.linspace(40, 100, 5000)
    ax.plot(x, pdf(x, *[m.values[p] for p in m.parameters]))
    ax.set_xlabel("Energy [keV]")
    ax.set_yscale('log')
    fig.savefig("plots/fit.png")


def main():
    set_root_preferences()
    #for dataset in C.ALL_DATASETS:
    #    test_fit(dataset)
    test_fit(C.POSTREACTOR_AM241_20260115)


if __name__ == "__main__":
    main()
