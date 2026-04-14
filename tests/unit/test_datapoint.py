import pytest
from molbench.molecule import Datapoint


def test_add_same_unit():
    result = Datapoint(1.0, "eV") + Datapoint(2.0, "eV")
    assert result.value == pytest.approx(3.0)
    assert result.unit == "eV"


def test_add_different_unit_raises():
    with pytest.raises(ValueError):
        Datapoint(1.0, "eV") + Datapoint(1.0, "au")


def test_sub():
    result = Datapoint(5.0, "au") - Datapoint(3.0, "au")
    assert result.value == pytest.approx(2.0)
    assert result.unit == "au"


def test_sub_different_unit_raises():
    with pytest.raises(ValueError):
        Datapoint(1.0, "eV") - Datapoint(1.0, "au")


def test_mul_scalar():
    result = Datapoint(2.0, "eV") * 3
    assert result.value == pytest.approx(6.0)
    assert result.unit == "eV"


def test_mul_float_scalar():
    result = Datapoint(2.0, "eV") * 0.5
    assert result.value == pytest.approx(1.0)


def test_mul_datapoint_raises():
    with pytest.raises(ValueError):
        Datapoint(2.0, "eV") * Datapoint(3.0, "eV")


def test_truediv():
    result = Datapoint(6.0, "au") / 2
    assert result.value == pytest.approx(3.0)
    assert result.unit == "au"


def test_truediv_datapoint_raises():
    with pytest.raises(ValueError):
        Datapoint(6.0, "au") / Datapoint(2.0, "au")


def test_floordiv():
    result = Datapoint(7.0, "au") // 2
    assert result.value == pytest.approx(3.0)
    assert result.unit == "au"


def test_floordiv_datapoint_raises():
    with pytest.raises(ValueError):
        Datapoint(6.0, "au") // Datapoint(2.0, "au")


def test_abs_negative():
    result = abs(Datapoint(-3.5, "au"))
    assert result.value == pytest.approx(3.5)
    assert result.unit == "au"


def test_abs_positive():
    result = abs(Datapoint(3.5, "au"))
    assert result.value == pytest.approx(3.5)


def test_eq_same():
    assert Datapoint(1.0, "eV") == Datapoint(1.0, "eV")


def test_eq_case_insensitive_unit():
    assert Datapoint(1.0, "eV") == Datapoint(1.0, "EV")


def test_neq_different_value():
    assert Datapoint(1.0, "eV") != Datapoint(2.0, "eV")


def test_neq_different_unit():
    assert Datapoint(1.0, "eV") != Datapoint(1.0, "au")


def test_neq_non_datapoint():
    assert Datapoint(1.0, "eV") != 1.0


def test_repr():
    dp = Datapoint(1.5, "au")
    assert repr(dp) == "1.5 au"
