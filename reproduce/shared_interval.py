"""Apply the prespecified interval ablation to all saved follow-up fields."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np
from scipy.stats import binomtest

from dtc_interval import estimate, value_gradient

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
import robustness_design as benchmark




def save(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")
    temporary.replace(path)


def fields():
    saved = [json.loads(s) for s in
             (DATA/"robustness/EVALUATION.jsonl").read_text().splitlines()]
    by_seed = {r["seed"]: r for r in saved}
    for family, (setting, config) in enumerate(benchmark.SETTINGS.items()):
        for group, q in enumerate((*benchmark.GRID, *benchmark.OFF_GRID, None)):
            fractions = benchmark.ALTERNATIVE if q is None else benchmark.prior.restricted(q)
            for replicate in range(80 if q is None else 20):
                seed = 1200000+10000*family+100*group+replicate
                old = by_seed[seed]
                xy, codes, target = benchmark.prior.generate(
                    config, fractions, benchmark.rng(family, 2, group, replicate)
                )
                assert (len(xy), target, setting, q) == (
                    old["records"], old["target"], old["setting"], old["q"])
                yield xy, codes, old


def run(output):
    output.mkdir()
    (output/"arrays").mkdir()
    state = dict(status="running", completed=0, started_epoch=time.time())
    save(output/"STATE.json", state)
    def expired(*_):
        raise TimeoutError("Twenty-minute ablation limit")
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(1200)
    try:
        with (output/"collector.lock").open("x") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
            with (output/"FIELDS.jsonl").open("x", buffering=1) as stream:
                for xy, codes, old in fields():
                    result, arrays, stats = estimate(xy, codes, old["seed"])
                    point, _ = value_gradient(
                        stats["fine"].sum(0), stats["W"], stats["area"], stats["disc"])
                    np.testing.assert_allclose(
                        point, old["dtc_connected_per_anchor"],
                        rtol=1e-10, atol=2e-12)
                    interval = result.get("spatial_bootstrap_t_ci")
                    row = dict(seed=old["seed"], setting=old["setting"], q=old["q"],
                               target=old["target"], records=len(xy), result=result,
                               flag=None if interval is None else interval[0]>0,
                               covered=None if interval is None else interval[0]<=old["target"]<=interval[1],
                               crosslabel_interval=old["interval"])
                    np.savez_compressed(output/"arrays"/f'{old["seed"]}.npz',
                                        **arrays, **stats)
                    stream.write(json.dumps(row, allow_nan=False)+"\n")
                    state["completed"] += 1
                    save(output/"STATE.json", state)
                    if state["completed"] % 60 == 0:
                        print(f'{state["completed"]}/540 fields saved', flush=True)
        state["status"] = "complete"
    except BaseException as exc:
        state.update(status="failed_or_partial", error=repr(exc))
        raise
    finally:
        signal.alarm(0)
        state["ended_epoch"] = time.time()
        save(output/"STATE.json", state)


def summarize(folder):
    state = json.loads((folder/"STATE.json").read_text())
    if state["status"] != "complete" or state["completed"] != 540:
        raise ValueError("Complete 540-field comparison required")
    rows = [json.loads(s) for s in (folder/"FIELDS.jsonl").read_text().splitlines()]
    groups = []
    for setting in benchmark.SETTINGS:
        for adequate in (True, False):
            group = [r for r in rows if r["setting"]==setting and (r["q"] is not None)==adequate]
            support = [r for r in group if r["flag"] is not None]
            cross = [r["crosslabel_interval"][0]>0 for r in support]
            direct = [r["flag"] for r in support]
            a = sum(x and not y for x, y in zip(cross, direct))
            b = sum(y and not x for x, y in zip(cross, direct))
            groups.append(dict(setting=setting, adequate=adequate, planned=len(group),
                               supported=len(support), dtc_alerts=sum(direct),
                               crosslabel_alerts=sum(cross), crosslabel_only=a, dtc_only=b,
                               both=sum(x and y for x,y in zip(cross,direct)),
                               paired_p=binomtest(a,a+b).pvalue if a+b else 1.,
                               dtc_covered=sum(r["covered"] for r in support)))
    save(folder/"SUMMARY.json", dict(groups=groups, fields=len(rows)))
    print(json.dumps(groups, indent=2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--summarize", type=Path)
    args = parser.parse_args()
    if args.output:
        run(args.output)
    else:
        summarize(args.summarize)
