"""Tests for the kinetics fitting layer (no GUI)."""
import csv
from pathlib import Path

import numpy as np
import pytest

from fvm.fitting import PARAMETERS, BedData, FitCase, fit

REFERENCE = Path(__file__).resolve().parents[1] / (
    "docs/kesten_claude/vapor_reference.csv")
pytestmark = pytest.mark.skipif(not REFERENCE.exists(),
                                reason="Kesten reference data absent")


# -- data loading -----------------------------------------------------------
def test_loads_kesten_csv_with_unit_conversion():
    d = BedData.from_csv(REFERENCE)
    assert len(d) == 6
    assert "T" in d.channels and "X" in d.channels
    assert d.z[0] < d.z[-1], "positions must be sorted"
    assert d.z[-1] == pytest.approx(0.25 * 0.3048, rel=1e-6)   # z_ft -> m
    assert 400.0 < d.T.min() < 500.0                            # degR -> K


def test_frac3d_is_converted_to_the_two_step_convention():
    """Kesten's f = 0 is X = 0.25; loading must apply (3f+1)/4.

    Checked on the minimum rather than the first row, because BedData sorts by
    position and his stations are not listed in axial order.
    """
    d = BedData.from_csv(REFERENCE)
    assert d.X.min() == pytest.approx(0.25, abs=1e-6)
    assert d.X.max() == pytest.approx(0.734, abs=1e-3)


def test_csv_without_a_position_column_is_rejected(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("T_K,X\n900,0.5\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no position column"):
        BedData.from_csv(p)


def test_csv_with_nothing_to_fit_is_rejected(tmp_path):
    p = tmp_path / "bare.csv"
    p.write_text("z_mm,p_bar\n0,8.4\n10,8.3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="nothing to fit"):
        BedData.from_csv(p)


def test_generic_si_columns_load(tmp_path):
    p = tmp_path / "si.csv"
    p.write_text("z_mm,T_K,X\n0,500,0.10\n10,900,0.50\n20,950,0.70\n",
                 encoding="utf-8")
    d = BedData.from_csv(p)
    assert d.z[-1] == pytest.approx(0.02)
    assert d.T[-1] == pytest.approx(950.0)


# -- the window -------------------------------------------------------------
def test_window_selects_points():
    case = FitCase.kesten_demo()
    assert case.mask().sum() == 6
    case.z_window = (1e-3, None)
    assert case.mask().sum() == 2


def test_excluding_the_liquid_vapour_zone_transforms_the_residual():
    """The model does not represent that zone, so fitting against it is wrong.

    Including the first millimetre puts the temperature residual at ~59 K even
    at Kesten's own parameters; excluding it gives 1.5 K. That is the
    difference between the tool being usable and not.
    """
    case = FitCase.kesten_demo()
    kesten = dict(A_NH3=1.0e11, n_H2=1.0)
    wide, _ = case.errors(kesten)
    case.z_window = (1e-3, None)
    narrow, _ = case.errors(kesten)
    assert wide["T"] > 40.0
    assert narrow["T"] < 5.0


def test_model_reproduces_the_developed_bed(case=None):
    """Where the model applies, it should match Kesten to a couple of kelvin."""
    case = FitCase.kesten_demo()
    case.z_window = (1e-3, None)
    errs, _ = case.errors(dict(A_NH3=1.0e11, n_H2=1.0))
    assert errs["T"] < 5.0
    assert errs["X"] < 0.05


# -- objective --------------------------------------------------------------
def test_residuals_are_normalised_across_channels():
    """Temperature in kelvin must not swamp a dimensionless fraction."""
    case = FitCase.kesten_demo()
    case.z_window = (1e-3, None)
    r = case.residuals(dict(A_NH3=1.0e11, n_H2=1.0))
    assert np.all(np.abs(r) < 10.0)


def test_residuals_require_an_enabled_channel():
    case = FitCase.kesten_demo()
    with pytest.raises(ValueError, match="no channels"):
        case.residuals(dict(A_NH3=1e11), weights={"T": 0.0, "X": 0.0, "p": 0.0})


def test_kestens_values_beat_the_published_ones_on_his_own_data():
    """Consistent with the reference output coming from MAIN.f, not PARAM.f."""
    case = FitCase.kesten_demo()
    case.z_window = (1e-3, None)
    main_f, _ = case.errors(dict(A_NH3=1.0e11, n_H2=1.0))
    eq43, _ = case.errors(dict(A_NH3=0.3e11, n_H2=1.6))
    assert main_f["T"] < eq43["T"]


# -- fitting ----------------------------------------------------------------
@pytest.mark.slow
def test_autofit_improves_on_a_wrong_start():
    case = FitCase.kesten_demo()
    case.z_window = (1e-3, None)
    start = dict(A_NH3=2.0e10, n_H2=1.3)
    before = float(np.sqrt(np.mean(case.residuals(start) ** 2)))
    res = fit(case, fit_names=("A_NH3", "n_H2"), start=start, max_nfev=30)
    # The objective is the combined normalised residual over all weighted
    # channels, so assert on that rather than on temperature alone -- a fit
    # may trade a little temperature error for a better dissociation match.
    assert res["rms"] <= before + 1e-9
    assert res["errors"]["T"] < 10.0


@pytest.mark.slow
def test_fit_reports_sensitivity_for_each_parameter():
    case = FitCase.kesten_demo()
    case.z_window = (1e-3, None)
    res = fit(case, fit_names=("A_NH3", "n_H2"), max_nfev=12)
    assert set(res["sensitivity"]) == {"A_NH3", "n_H2"}
    assert all(v >= 0.0 for v in res["sensitivity"].values())


def test_unknown_parameter_is_rejected():
    case = FitCase.kesten_demo()
    with pytest.raises(ValueError, match="unknown parameter"):
        fit(case, fit_names=("nonsense",))


def test_parameter_table_is_self_consistent():
    for name, spec in PARAMETERS.items():
        lo, hi = spec["bounds"]
        assert lo < spec["default"] < hi, f"{name} default outside its bounds"
        if spec["log"]:
            assert lo > 0.0, f"{name} is logarithmic so bounds must be positive"
