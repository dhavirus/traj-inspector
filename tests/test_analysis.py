"""Numerical building blocks on hand-built inputs; expected values are worked out by hand in the comments."""
import math

import numpy as np
import pytest

from trajinspector import analysis as A


def test_count_transitions_matches_cln025_selftest():
    q = np.array([0.1, 0.5, 0.9, 0.85, 0.5, 0.2, 0.5, 0.9])
    t = A.count_transitions(q, 0.3, 0.8)
    assert (t["folds"], t["unfolds"]) == (2, 1)
    assert list(t["state"]) == [0, 0, 1, 1, 1, 0, 0, 1]
    assert t["folded_fraction"] == pytest.approx(0.5)


def test_series_that_never_leaves_the_middle_has_no_state():
    t = A.count_transitions(np.array([0.5, 0.79, 0.81, 0.79, 0.81]), 0.3, 0.8)
    assert t["folds"] == 0 and list(t["state"]) == [-1, -1, 1, 1, 1]


def test_min_image_cubic_box():
    # atom at x=1 and x=9 in a 10 A box: the image of the second sits at x=-1, 2 A from the first
    xyz = np.array([[1.0, 5, 5], [9.0, 5, 5]])
    assert A.min_image_distance_angstrom(xyz, np.array([10, 10, 10, 90, 90, 90.0])) == pytest.approx(2.0)


def test_min_image_triclinic_uses_the_diagonal_shift():
    # gamma=150: |a+b|^2 = 100+100+2*100*cos(150) = 26.8 -> 5.18 A, shorter than any axis (10 A)
    xyz = np.zeros((1, 3))
    d = A.min_image_distance_angstrom(xyz, np.array([10, 10, 10, 90, 90, 150.0]))
    assert d == pytest.approx(math.sqrt(200 + 200 * math.cos(math.radians(150))))


def test_q_values():
    keys = [(i, "CA") for i in range(6)]
    ref = np.array([[0, 0, 0], [30, 0, 0], [60, 0, 0], [90, 0, 0], [120, 0, 0], [4.0, 0, 0]])   # only residues 0 and 5 touch
    contacts, dropped = A.native_contacts(keys, ref, keys)
    assert dropped == 0 and len(contacts.first) == 1 and contacts.r0_angstrom[0] == pytest.approx(4.0)
    assert A.best_hummer_q(ref, contacts) == pytest.approx(1 / (1 + math.exp(5 * (4 - 1.8 * 4))))   # 0.982
    far = ref.copy()
    far[5, 0] = 1.8 * 4.0
    assert A.best_hummer_q(far, contacts) == pytest.approx(0.5)


def test_native_contacts_counts_unmatched_reference_atoms():
    ref_keys = [(0, "CA"), (5, "CA"), (5, "O1")]
    contacts, dropped = A.native_contacts(ref_keys, np.array([[0, 0, 0], [4.0, 0, 0], [4.0, 1, 0]]), [(0, "CA"), (5, "CA")])
    assert dropped == 1 and len(contacts.first) == 1


def test_extended_reference_is_a_hard_stop():
    keys = [(i, "CA") for i in range(6)]
    with pytest.raises(ValueError, match="no native contacts"):
        A.native_contacts(keys, np.array([[30.0 * i, 0, 0] for i in range(6)]), keys)


def test_block_ratio_flags_correlated_series():
    rng = np.random.default_rng(1)
    white = rng.normal(size=8000)
    walk = np.cumsum(rng.normal(size=8000))
    assert A.block_se_ratio(walk, 8) > 1.25   # random walk: SE still growing
    assert 0.5 < A.block_se_ratio(white, 8) < 1.5
    assert math.isnan(A.block_se_ratio(np.ones(100), 8))   # no variance


def test_slope_units():
    t = np.arange(10.0)
    assert A.slope(t, 3 * t + 2) == pytest.approx(3.0)
    assert math.isnan(A.slope(t[:2], t[:2]))


def test_ca_pairs_skip_neighbours_and_subsample():
    pairs, sub = A.ca_feature_pairs(5)
    assert not sub and len(pairs) == 6 and all(b - a >= 2 for a, b in pairs)   # (0,2)(0,3)(0,4)(1,3)(1,4)(2,4)
    pairs, sub = A.ca_feature_pairs(200)
    assert sub and pairs.max() == 199 and len(np.unique(pairs)) <= A.MAX_CA_FEATURE_ATOMS


def test_msm_summary_on_two_state_chain():
    rng = np.random.default_rng(0)
    state = np.cumsum(rng.random(6000) < 0.02) % 2          # switches on average every 50 frames
    feats = [np.column_stack([state, 1 - state, rng.normal(size=6000) * 0.1]) + rng.normal(size=(6000, 3)) * 0.05]
    out = A.msm_summary(feats, lag_frames=5, seed=0)
    assert out["vamp2"] > 1.5 and out["its_frames"][5] > 20 and out["ck_error"] < 0.1


def test_constant_series_has_no_inefficiency():
    out = A.stat_inefficiency(np.full(100, 0.7))
    assert math.isnan(out["g"]) and math.isnan(out["n_eff"])
