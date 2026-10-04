import pytest
from phytclust.algo.bins import define_bins
from phytclust.exceptions import InvalidKError


class _PC:
    def __init__(self, n):
        self.num_terminals = n


def test_define_bins_basic_properties():
    pc = _PC(100)
    bins = define_bins(pc, num_bins=3, k_lo=1, k_hi=100)

    assert len(bins) >= 1
    assert bins[0][0] == 1
    assert bins[-1][1] == 100

    for (lo1, hi1), (lo2, hi2) in zip(bins, bins[1:]):
        assert lo1 <= hi1
        assert lo2 <= hi2
        assert hi1 + 1 == lo2

    assert all(1 <= lo <= hi <= 100 for lo, hi in bins)


def test_define_bins_rejects_invalid_range():
    pc = _PC(10)
    with pytest.raises(InvalidKError, match="k_hi must be > k_lo"):
        define_bins(pc, num_bins=3, k_lo=5, k_hi=5)


@pytest.mark.parametrize("parameter", ["num_bins", "k_lo", "k_hi"])
@pytest.mark.parametrize("value", [0, -1, 1.5, True, float("nan"), float("inf"), "3"])
def test_bin_parameters_require_positive_integers(parameter, value):
    from phytclust.exceptions import ConfigurationError
    arguments = {"num_bins": 3, "k_lo": 1, "k_hi": 10}
    arguments[parameter] = value
    error = ConfigurationError if parameter == "num_bins" else InvalidKError
    with pytest.raises(error, match=parameter):
        define_bins(_PC(10), **arguments)


def test_rounding_reduces_bin_count_without_gaps():
    bins = define_bins(_PC(3), num_bins=10)
    assert len(bins) < 10
    assert [k for lower, upper in bins for k in range(lower, upper + 1)] == [1, 2, 3]
    assert all(type(boundary) is int for pair in bins for boundary in pair)


def test_custom_range_is_covered_exactly():
    bins = define_bins(_PC(100), num_bins=4, k_lo=3, k_hi=20)
    assert [k for lower, upper in bins for k in range(lower, upper + 1)] == list(range(3, 21))
