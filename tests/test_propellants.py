"""Tests for tabulated propellants."""
import pytest

from fvm.mechanism import HydrazineShell405
from fvm.propellants import (PROPELLANT_DIR, TabulatedPropellant, available,
                             load, write_template)


def test_example_table_is_present():
    assert "hydrazine_X084" in available()


def test_tabulated_path_reproduces_the_reacting_mechanism():
    """The whole point of the shared interface: same answer, either route."""
    p = load("hydrazine_X084")
    tab = p.chamber_conditions(8.44e5)
    mech = HydrazineShell405().chamber_conditions(0.84)
    for key in ("T", "MW", "gamma", "cstar"):
        assert tab[key] == pytest.approx(mech[key], rel=1e-4), key


def test_perfect_gas_matches_the_table():
    p = load("hydrazine_X084")
    c = p.chamber_conditions(8.44e5)
    gas = p.perfect_gas(8.44e5)
    assert gas.gamma == pytest.approx(c["gamma"])
    assert gas.MW == pytest.approx(c["MW"])


def _table(tmp_path, body, name="demo"):
    path = tmp_path / f"{name}.csv"
    path.write_text(body, encoding="utf-8")
    return path


def test_single_row_is_pressure_independent(tmp_path):
    p = TabulatedPropellant.from_csv(
        _table(tmp_path, "p_bar,T_K,MW,gamma\n8.0,2100,23.8,1.24\n"))
    assert p.chamber_conditions(2e5)["T"] == pytest.approx(2100.0)
    assert p.chamber_conditions(50e5)["T"] == pytest.approx(2100.0)


def test_multiple_rows_interpolate(tmp_path):
    p = TabulatedPropellant.from_csv(_table(
        tmp_path, "p_bar,T_K,MW,gamma\n5,2000,24,1.25\n20,2100,23,1.23\n"))
    c = p.chamber_conditions(10e5)
    assert 2000.0 < c["T"] < 2100.0
    assert 23.0 < c["MW"] < 24.0


def test_extrapolation_is_refused_rather_than_guessed(tmp_path):
    p = TabulatedPropellant.from_csv(_table(
        tmp_path, "p_bar,T_K,MW,gamma\n5,2000,24,1.25\n20,2100,23,1.23\n"))
    with pytest.raises(ValueError, match="outside the tabulated range"):
        p.chamber_conditions(40e5)


def test_missing_columns_name_what_is_required(tmp_path):
    with pytest.raises(ValueError, match="missing required column"):
        TabulatedPropellant.from_csv(_table(tmp_path, "p_bar,T_K\n8,2100\n"))


def test_lowercase_headers_are_accepted(tmp_path):
    p = TabulatedPropellant.from_csv(
        _table(tmp_path, "p_bar,t_k,mw,gamma\n8.0,2100,23.8,1.24\n"))
    assert p.chamber_conditions()["T"] == pytest.approx(2100.0)


def test_comment_lines_are_ignored(tmp_path):
    p = TabulatedPropellant.from_csv(_table(
        tmp_path, "# from CEA run 2026-08-12\np_bar,T_K,MW,gamma\n8,2100,23.8,1.24\n"))
    assert p.chamber_conditions()["MW"] == pytest.approx(23.8)


def test_cstar_is_computed_when_absent_and_used_when_present(tmp_path):
    without = TabulatedPropellant.from_csv(_table(
        tmp_path, "p_bar,T_K,MW,gamma\n8,2100,23.8,1.24\n", "a"))
    with_cea = TabulatedPropellant.from_csv(_table(
        tmp_path, "p_bar,T_K,MW,gamma,cstar_m_s\n8,2100,23.8,1.24,1601\n", "b"))
    assert with_cea.chamber_conditions()["cstar"] == pytest.approx(1601.0)
    assert without.chamber_conditions()["cstar"] != pytest.approx(1601.0, rel=1e-6)
    assert 1000.0 < without.chamber_conditions()["cstar"] < 2000.0


def test_mole_fraction_columns_are_carried(tmp_path):
    p = TabulatedPropellant.from_csv(_table(
        tmp_path, "p_bar,T_K,MW,gamma,x_N2,x_H2O\n8,2100,23.8,1.24,0.32,0.41\n"))
    assert p.composition == {"n2": pytest.approx(0.32), "h2o": pytest.approx(0.41)}


def test_unknown_propellant_lists_what_is_available():
    with pytest.raises(FileNotFoundError, match="no propellant table"):
        load("unobtainium")


def test_template_round_trips(tmp_path):
    path = write_template("newprop", directory=tmp_path)
    assert path.exists()
    text = path.read_text(encoding="utf-8")
    for col in ("p_bar", "T_K", "MW", "gamma"):
        assert col in text
    with pytest.raises(FileExistsError):
        write_template("newprop", directory=tmp_path)


def test_template_is_blank_rather_than_prefilled(tmp_path):
    """It must not ship guessed numbers that could be mistaken for data."""
    text = write_template("x", directory=tmp_path).read_text(encoding="utf-8")
    data = [ln for ln in text.splitlines()
            if ln.strip() and not ln.startswith("#") and "p_bar" not in ln]
    assert all(set(ln) <= {","} for ln in data), "template must have no values"
