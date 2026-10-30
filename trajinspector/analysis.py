"""Numerical building blocks for gates 2-4 (pure numpy / pymbar / deeptime; no file access).

Ported from reference/cln025.py (best_hummer_q, count_transitions, stat_inefficiency, tica_msm), with MDAnalysis
coordinates in Angstrom instead of mdtraj nm.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.lib.mdamath import triclinic_vectors

# Best, Hummer & Eaton 2013 PNAS 110, 17874: native pair < 4.5 A, beta = 5 /A, lambda = 1.8 (heavy atoms, |i-j| > 3)
Q_CONTACT_CUTOFF_ANGSTROM = 4.5
Q_BETA_PER_ANGSTROM = 5.0
Q_LAMBDA = 1.8
Q_MIN_RESIDUE_SEPARATION = 3
MIN_IMAGE_SEARCH_ANGSTROM = 30.0   # shortcut: distances beyond this are reported as the cap; raise for protein-image checks above 3 nm
MAX_CA_FEATURE_ATOMS = 50          # shortcut: C-alpha subsampling keeps the pair count bounded; marked in the report
ITS_LAG_FACTORS = (0.25, 0.5, 1, 2, 4)
CK_MODELS = 4


@dataclass(frozen=True)
class NativeContacts:
    first: np.ndarray    # indices into the analysed selection
    second: np.ndarray
    r0_angstrom: np.ndarray


def atom_keys(residue_indices: np.ndarray, names: np.ndarray) -> list[tuple[int, str]]:
    """(residue ordinal inside the selection, atom name): the key used to match atoms between files."""
    ordinal = np.unique(residue_indices, return_inverse=True)[1]
    return list(zip(ordinal.tolist(), names.tolist()))


def native_contacts(ref_keys: list, ref_xyz: np.ndarray, target_keys: list) -> tuple[NativeContacts, int]:
    """Native heavy-atom contacts of the reference, expressed in target indices. Also returns how many reference
    atoms have no partner in the target (they are dropped)."""
    target_index = {k: i for i, k in enumerate(target_keys)}
    found = [(target_index[k], n) for n, k in enumerate(ref_keys) if k in target_index]
    target_idx = np.array([f[0] for f in found], dtype=int)
    ref_pos = ref_xyz[[f[1] for f in found]]
    ordinal = np.array([target_keys[i][0] for i in target_idx])
    pairs, dist = capped_distance(ref_pos, ref_pos, max_cutoff=Q_CONTACT_CUTOFF_ANGSTROM, return_distances=True)
    keep = (pairs[:, 0] < pairs[:, 1]) & (np.abs(ordinal[pairs[:, 0]] - ordinal[pairs[:, 1]]) > Q_MIN_RESIDUE_SEPARATION)
    if not keep.any():
        raise ValueError("the reference structure has no native contacts: pass the folded structure as reference")
    return (NativeContacts(target_idx[pairs[keep, 0]], target_idx[pairs[keep, 1]], dist[keep]),
            len(ref_keys) - len(found))


def best_hummer_q(xyz: np.ndarray, contacts: NativeContacts) -> float:
    """Fraction of native contacts in one frame (xyz in Angstrom, molecule already unwrapped)."""
    r = np.linalg.norm(xyz[contacts.first] - xyz[contacts.second], axis=1)
    return float((1.0 / (1.0 + np.exp(Q_BETA_PER_ANGSTROM * (r - Q_LAMBDA * contacts.r0_angstrom)))).mean())


def count_transitions(q: np.ndarray, unfolded: float, folded: float) -> dict:
    """Two-threshold transition count. Frames between the thresholds keep the previous state (-1 before the first)."""
    state, last, folds, unfolds = np.full(len(q), -1), -1, 0, 0
    for i, x in enumerate(q):
        new = 1 if x > folded else 0 if x < unfolded else last
        folds += last == 0 and new == 1
        unfolds += last == 1 and new == 0
        state[i] = last = new
    assigned = state >= 0
    return {"folds": int(folds), "unfolds": int(unfolds), "state": state,
            "folded_fraction": float((state[assigned] == 1).mean()) if assigned.any() else float("nan")}


def stat_inefficiency(x: np.ndarray) -> dict:
    """pymbar equilibration cut t0 (information only) and g, N_eff of the *whole* series.

    detect_equilibration maximises N_eff and assumes one relaxing basin. On a series that hops between two states
    it can cut away almost everything (a telegraph series of 4000 frames gave t0 = 3835), so g and N_eff are
    taken from all frames. A constant series has no defined g: both are NaN.
    """
    from pymbar import timeseries
    from pymbar.utils import ParameterError
    x = np.asarray(x, dtype=float)
    t0 = int(timeseries.detect_equilibration(x)[0])
    try:
        if np.ptp(x) == 0:   # float round-off in the mean would otherwise give a spurious g
            raise ParameterError("constant series")
        g = float(timeseries.statistical_inefficiency(x))
    except ParameterError:
        return {"t0_frame": t0, "g": float("nan"), "n_eff": float("nan")}
    return {"t0_frame": t0, "g": g, "n_eff": len(x) / g}


def block_se_ratio(x: np.ndarray, min_blocks: int) -> float:
    """SE with the largest block size that still gives min_blocks blocks / SE with half that block size.

    ~1 means the standard error has plateaued; larger means correlation is not yet averaged out.
    Returns NaN when the series is too short or has no variance.
    """
    def standard_error(size: int) -> float:
        n = len(x) // size
        means = x[: n * size].reshape(n, size).mean(axis=1)
        return float(means.std(ddof=1) / np.sqrt(n))
    size = len(x) // min_blocks
    if size < 2:
        return float("nan")
    small = standard_error(size // 2)
    return standard_error(size) / small if small > 0 else float("nan")


def min_image_distance_angstrom(xyz: np.ndarray, dimensions: np.ndarray) -> float:
    """Smallest distance between the selection and its own periodic image (26 lattice shifts, triclinic aware)."""
    box = triclinic_vectors(dimensions)
    best = MIN_IMAGE_SEARCH_ANGSTROM
    for shift in itertools.product((-1, 0, 1), repeat=3):
        if shift == (0, 0, 0):
            continue
        _, dist = capped_distance(xyz, xyz + np.array(shift) @ box, max_cutoff=best, return_distances=True)
        if len(dist):
            best = min(best, float(dist.min()))
    return best


def slope(t: np.ndarray, y: np.ndarray) -> float:
    """Least-squares slope dy/dt; NaN with fewer than 3 finite points."""
    ok = np.isfinite(t) & np.isfinite(y)
    return float(np.polyfit(t[ok], y[ok], 1)[0]) if ok.sum() > 2 else float("nan")


def ca_feature_pairs(n_ca: int) -> tuple[np.ndarray, bool]:
    """Index pairs (i, j >= i+2) over at most MAX_CA_FEATURE_ATOMS evenly spaced C-alpha atoms; flag if subsampled."""
    chosen = np.arange(n_ca) if n_ca <= MAX_CA_FEATURE_ATOMS else np.unique(np.linspace(0, n_ca - 1, MAX_CA_FEATURE_ATOMS).astype(int))
    pairs = [(a, b) for k, a in enumerate(chosen) for b in chosen[k + 2:]]
    return np.array(pairs, dtype=int).reshape(-1, 2), n_ca > MAX_CA_FEATURE_ATOMS


def msm_summary(features: list[np.ndarray], lag_frames: int, seed: int) -> dict:
    """VAMP-2 score, slowest implied timescale (frames) at lag * ITS_LAG_FACTORS, and the max CK error.

    Missing entries (None / absent) mean there was too little data at that lag.
    """
    from deeptime.clustering import KMeans
    from deeptime.decomposition import TICA, VAMP
    from deeptime.markov import TransitionCountEstimator
    from deeptime.markov.msm import MaximumLikelihoodMSM
    feats = [f.astype(np.float64) for f in features]
    vamp2 = float(VAMP(lagtime=lag_frames, dim=2).fit(feats).fetch_model().score(r=2))
    proj = [TICA(lagtime=lag_frames, dim=2).fit(feats).fetch_model().transform(f) for f in feats]
    n_clusters = int(min(100, max(2, sum(len(p) for p in proj) // 10)))
    kmeans = KMeans(n_clusters=n_clusters, fixed_seed=seed, max_iter=200).fit(np.concatenate(proj)).fetch_model()
    dtrajs = [kmeans.transform(p) for p in proj]

    def msm_at(lag: int):
        counts = TransitionCountEstimator(lag, "sliding").fit(dtrajs).fetch_model().submodel_largest()
        return MaximumLikelihoodMSM().fit(counts).fetch_model()

    its = {}
    for lag in sorted({max(1, round(f * lag_frames)) for f in ITS_LAG_FACTORS}):
        try:
            timescales = msm_at(lag).timescales(k=1)
        except ValueError:   # too little data at this lag
            timescales = []
        its[lag] = float(timescales[0]) if len(timescales) else None   # a one-state MSM has no timescale
    try:
        models = [msm_at(lag_frames * k) for k in range(1, CK_MODELS + 1)]
        test = models[0].ck_test(models, n_metastable_sets=2)
        ck_error = float(np.abs(test.predictions - test.estimates).max())
    except ValueError:
        ck_error = None
    return {"vamp2": vamp2, "its_frames": its, "ck_error": ck_error, "n_clusters": n_clusters}
