"""Verify the record-disjoint witness certificate and its edit lower bound."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

RADIUS = 30.

def packing(xy, codes):
    tree = cKDTree(xy)
    available = np.ones(len(xy), dtype=bool)
    triples = []
    for anchor in range(len(xy)):
        if not available[anchor]:
            continue
        near = np.asarray(tree.query_ball_point(xy[anchor], RADIUS, return_sorted=True))
        near = near[available[near] & (codes[near] != codes[anchor])]
        labels = np.unique(codes[near])
        if len(labels) < 2:
            continue
        first = near[codes[near] == labels[0]][0]
        second = near[codes[near] == labels[1]][0]
        triple = [anchor, int(first), int(second)]
        available[triple] = False
        triples.append(triple)
    triples = np.asarray(triples, dtype=int).reshape(-1, 3)
    if len(np.unique(triples)) != triples.size:
        raise ValueError("Packing reused a record")
    if len(triples):
        marks = codes[triples]
        if np.any((marks[:,0]==marks[:,1]) | (marks[:,0]==marks[:,2]) | (marks[:,1]==marks[:,2])):
            raise ValueError("Packing contains a repeated code")
        if np.any(np.linalg.norm(xy[triples[:,1:]]-xy[triples[:,0,None]], axis=2)>RADIUS):
            raise ValueError("Packing contains an invalid anchor neighborhood")
    return triples

def verify(regenerate=False):
    folder = Path(__file__).resolve().parent / "data/witness_packing"
    with np.load(folder / "observations.npz") as data:
        xy, codes, groups = data["xy"], data["codes"], data["groups"]
    with np.load(folder / "certificate.npz") as data:
        certificate = data["triples"]
    ids = certificate[:, 1:]
    if len(np.unique(ids)) != ids.size:
        raise ValueError("Certificate reuses records")
    if np.any(groups[ids] != certificate[:, :1]):
        raise ValueError("Certificate crosses recorded groups")
    labels = codes[ids]
    if np.any((labels[:, 0] == labels[:, 1]) | (labels[:, 0] == labels[:, 2])
              | (labels[:, 1] == labels[:, 2])):
        raise ValueError("Certificate repeats codes")
    if np.any(np.linalg.norm(xy[ids[:, 1:]]-xy[ids[:, 0, None]], axis=2) > RADIUS):
        raise ValueError("Certificate exceeds the anchor radius")
    summary = json.loads((folder / "SUMMARY.json").read_text())
    if len(certificate) != summary["lower_bound_record_edits"] or len(xy) != summary["records"]:
        raise ValueError("Lower-bound count differs")
    if regenerate:
        rebuilt = []
        for group in np.unique(groups):
            indices = np.flatnonzero(groups == group)
            triples = packing(xy[indices], codes[indices])
            rebuilt.extend([[int(group), *map(int, triple)] for triple in indices[triples]])
        np.testing.assert_array_equal(rebuilt, certificate)
    print(f"{len(certificate)} disjoint triples verified in {len(xy)} recorded observations")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate", action="store_true")
    verify(parser.parse_args().regenerate)
