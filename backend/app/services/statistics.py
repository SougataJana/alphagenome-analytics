"""Deterministic, reproducible statistics for user-supplied score and gene sets."""

from datetime import datetime, timezone
from uuid import uuid4

import numpy as np
from scipy.stats import hypergeom


def _bh(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = np.argsort(p_values)
    adjusted = np.ones(n, dtype=float)
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * n / rank)
        adjusted[original_index] = running
    return adjusted.tolist()


def compare_groups(case_scores: list[float], control_scores: list[float], iterations: int, seed: int) -> dict:
    case = np.asarray(case_scores, dtype=float)
    control = np.asarray(control_scores, dtype=float)
    if not np.isfinite(case).all() or not np.isfinite(control).all():
        raise ValueError("Scores must be finite numeric values.")
    observed = float(case.mean() - control.mean())
    pooled_sd = float(np.sqrt(((case.size - 1) * case.var(ddof=1) + (control.size - 1) * control.var(ddof=1)) / (case.size + control.size - 2)))
    effect = observed / pooled_sd if pooled_sd else None
    rng = np.random.default_rng(seed)
    pooled = np.concatenate((case, control))
    exceedances = 0
    for _ in range(iterations):
        permuted = rng.permutation(pooled)
        delta = float(permuted[:case.size].mean() - permuted[case.size:].mean())
        exceedances += abs(delta) >= abs(observed)
    p_value = (exceedances + 1) / (iterations + 1)
    boot_rng = np.random.default_rng(seed + 1)
    boot = np.empty(iterations, dtype=float)
    for index in range(iterations):
        boot[index] = boot_rng.choice(case, size=case.size, replace=True).mean() - boot_rng.choice(control, size=control.size, replace=True).mean()
    return {
        "analysis_id": str(uuid4()),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "method": "two-sided permutation test of mean difference; percentile bootstrap confidence interval",
        "case_n": int(case.size), "control_n": int(control.size),
        "case_mean": float(case.mean()), "control_mean": float(control.mean()),
        "mean_difference": observed, "cohens_d": effect, "empirical_p_value": p_value,
        "confidence_interval_95": [float(v) for v in np.quantile(boot, [0.025, 0.975])],
        "iterations": iterations, "seed": seed,
        "correction": "plus-one empirical p-value correction; no multiple-test correction for this single comparison",
    }


def gene_enrichment(genes: list[str], background: list[str], gene_sets: dict[str, list[str]]) -> dict:
    selected = {str(g).strip().upper() for g in genes if str(g).strip()}
    universe = {str(g).strip().upper() for g in background if str(g).strip()}
    selected &= universe
    if not selected or len(universe) < len(selected):
        raise ValueError("The selected gene set must overlap the supplied background universe.")
    rows = []
    for name, members in gene_sets.items():
        pathway = {str(g).strip().upper() for g in members if str(g).strip()} & universe
        overlap = selected & pathway
        p = float(hypergeom.sf(len(overlap) - 1, len(universe), len(pathway), len(selected)))
        rows.append({"gene_set": name, "overlap_n": len(overlap), "set_n_in_background": len(pathway), "selected_n": len(selected), "background_n": len(universe), "genes": sorted(overlap), "p_value": p})
    adjusted = _bh([row["p_value"] for row in rows]) if rows else []
    for row, q_value in zip(rows, adjusted):
        row["q_value_bh"] = q_value
    rows.sort(key=lambda item: item["q_value_bh"])
    return {"analysis_id": str(uuid4()), "method": "one-sided hypergeometric over-representation test", "multiple_testing": "Benjamini-Hochberg FDR", "selected_genes_in_background": len(selected), "background_genes": len(universe), "results": rows}
