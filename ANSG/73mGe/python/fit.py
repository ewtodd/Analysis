from analysis_utils import load_cpp_library
from analysis_utils.init import set_root_preferences
from analysis_utils.io import load_tree_data
from analysis_utils.fitting import single_peak_pdf, estimate_peak_params, fit_single_peak
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

    m = fit_single_peak(energy, 45, 75, 59.5)
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
