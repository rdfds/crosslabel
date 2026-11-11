"""Regenerate the matched-moment training, reference, and evaluation fields."""
import argparse
import json
from pathlib import Path

import numpy as np

from crosslabel import field_check
from crosslabel.field_screen import density_score, reference_tail
import simulation as prior
from marked_pairs import Library, summary
from spinna_design import generate as geometry_field, schedule as geometry_schedule


CONFIG = prior.NEW["matched_moments"]
NULL = prior.restricted(CONFIG["null"])
FAMILY_INDEX = 6


def train(bank):
    endpoints = []
    with np.load(bank, allow_pickle=False) as saved:
        for j, name in enumerate(("monomer", "dimer")):
            histograms = []
            for i in range(64):
                rng = np.random.default_rng(np.random.SeedSequence([730100, 1, FAMILY_INDEX, j, i]))
                xy, codes, _ = prior.generate(CONFIG, prior.restricted(float(j)), rng)
                histograms.append(summary(xy, codes, 4))
            endpoint = np.stack(histograms)
            np.testing.assert_array_equal(endpoint, saved[name])
            endpoints.append(endpoint)
    return Library.train(*endpoints)


def measure(seed, fractions, model):
    xy, codes, target = prior.generate(CONFIG, fractions, np.random.default_rng(seed))
    q, score = model.fit(summary(xy, codes, 4))
    interval = field_check(xy, codes, code_count=4, radius_nm=30., seed=seed)
    return float(q[0]), float(score[0]), interval, target, len(xy)


def check_density(path, matched_scores):
    saved = json.loads(path.read_text())
    fields = {}
    for plan in geometry_schedule():
        xy, _, _ = geometry_field(plan)
        fields[plan["key"]] = density_score(xy)["score"]
    reference = [fields[p["key"]] for p in geometry_schedule()
                 if p["family"] == "reference" and p["adequate"]]
    np.testing.assert_allclose(reference, saved["reference"], rtol=1e-9, atol=1e-12)
    fields.update(matched_scores)
    if len(saved["rows"]) != len(fields):
        raise ValueError("Density screen roster differs")
    for row in saved["rows"]:
        score = fields.pop(row["key"])
        np.testing.assert_allclose(score, row["density"]["score"], rtol=1e-9, atol=1e-12)
        if reference_tail(score, reference) != row["screen_p"]:
            raise ValueError(f"Density rank differs at {row['key']}")
    print(f"{len(saved['rows'])} density-screen scores and ranks verified.")


def replay(bank, rows):
    model = train(bank)
    reference = json.loads(rows.with_name("REFERENCE.json").read_text())
    lower_ref, pair_ref = [], []
    for seed in range(946001, 946101):
        _, score, field, _, _ = measure(seed, NULL, model)
        lower_ref.append(field["interval_95_model"][0])
        pair_ref.append(score)
    np.testing.assert_allclose(lower_ref, reference["matched_lower"], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(pair_ref, reference["marked_pair"], rtol=1e-9, atol=1e-12)
    outcomes = {True: [0, 0], False: [0, 0]}
    saved_rows = [json.loads(line) for line in rows.read_text().splitlines()]
    schedule = [(984000 + (0 if adequate else 500) + i, adequate)
                for i in range(100) for adequate in (True, False)]
    if [(r["seed"], r["adequate"]) for r in saved_rows] != schedule:
        raise ValueError("Expected the original 200 held-out field identities")
    density = {}
    for row in saved_rows:
        adequate = row["adequate"]
        q, score, field, target, records = measure(
            row["seed"], NULL if adequate else CONFIG["alt"], model,
        )
        interval = field["interval_95_model"]
        if (target, records, field["status"]) != (row["target"], row["records"], row["status"]):
            raise ValueError(f"Regenerated field differs at {row['seed']}")
        np.testing.assert_allclose(
            [q, score, field["estimate"], *interval],
            [row["pair_q"], row["pair_score"], row["point"], *row["ci"]],
            rtol=1e-9, atol=1e-12,
        )
        ranks = reference_tail(interval[0], lower_ref), reference_tail(score, pair_ref)
        if ranks != (row["nested_p"], row["pair_p"]):
            raise ValueError(f"Regenerated ranks differ at {row['seed']}")
        density[str(row["seed"])] = field["density_screen"]["score"]
        for i, p in enumerate(ranks):
            outcomes[adequate][i] += p <= .05
    print(f"128 training, 100 reference, and 200 held-out fields verified. "
          f"False alerts {outcomes[True]}; detections {outcomes[False]}.")
    check_density(rows.with_name("DENSITY.json"), density)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bank", type=Path)
    parser.add_argument("rows", type=Path)
    args = parser.parse_args()
    replay(args.bank, args.rows)
