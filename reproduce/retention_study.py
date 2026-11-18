"""Reproduce the fixed Crosslabel retention and anchor-count experiment."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time

import numpy as np
import scipy
from scipy.stats import binomtest
from crosslabel.field import field_check
from generator import generate as simulate

RETENTION = (.15, .25, .50, .75, 1.)
COUNTS = (3500, 7000, 14000)
REPEATS = 80

def save(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def schedule():
    for replicate in range(1, REPEATS + 1):
        for pi, p in enumerate(RETENTION):
            for ni, count in enumerate(COUNTS):
                for alternative in (0, 1):
                    seed = 18000000 + 100000*pi + 10000*ni + 1000*alternative + replicate
                    yield dict(seed=seed, replicate=replicate, retention=p,
                               expected_records=count, alternative=bool(alternative),
                               entropy=[18112026, pi, ni, alternative, replicate])


def generate(plan):
    config = dict(shape="triangle", spacing=12., K=4, mask=False,
                  p=plan["retention"], density=plan["expected_records"]/100.)
    fractions = (.7, .2, .1) if plan["alternative"] else (.5, .5, 0.)
    return simulate(config, fractions, np.random.default_rng(plan["entropy"]))


def catalogue_hash(xy, codes):
    digest = hashlib.sha256(np.ascontiguousarray(xy, dtype="<f8").tobytes())
    digest.update(np.ascontiguousarray(codes, dtype="<i8").tobytes())
    return digest.hexdigest()


def run(output):
    output.mkdir()
    (output / "fields").mkdir()
    start = time.time()
    state = dict(status="running", started_epoch=start, attempted=0, completed=0)
    save(output / "STATE.json", state)
    hardware = dict(platform=platform.platform(), machine=platform.machine(),
                    processor=platform.processor(), python=sys.version,
                    numpy=np.__version__, scipy=scipy.__version__,
                    openblas_threads=os.environ.get("OPENBLAS_NUM_THREADS"),
                    omp_threads=os.environ.get("OMP_NUM_THREADS"))
    if sys.platform == "darwin":
        query = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                               text=True, capture_output=True)
        hardware["cpu"] = query.stdout.strip() if query.returncode == 0 else None
        hardware["cpu_query_error"] = query.stderr.strip() if query.returncode else None
    save(output / "ENVIRONMENT.json", hardware)

    def expired(*_):
        raise TimeoutError("Fixed 40-minute execution limit")

    signal.signal(signal.SIGALRM, expired)
    signal.alarm(2400)
    try:
        with (output / "collector.lock").open("x") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with (output / "FIELDS.jsonl").open("x", buffering=1) as rows:
                for plan in schedule():
                    if time.time() - start >= 2400:
                        expired()
                    if shutil.disk_usage(output).free < 2 << 30:
                        raise OSError("Two GiB write reserve reached")
                    state["attempted"] += 1
                    state["active_seed"] = plan["seed"]
                    xy, codes, target = generate(plan)
                    wall, cpu = time.perf_counter(), time.process_time()
                    result = field_check(xy, codes, code_count=4, radius_nm=30., seed=plan["seed"])
                    cpu_seconds = time.process_time() - cpu
                    wall_seconds = time.perf_counter() - wall
                    interval = result["interval_95_model"]
                    row = dict(plan, records=len(xy), target=target,
                               catalogue_sha256=catalogue_hash(xy, codes),
                               result=result, flag=result["positive_lower_limit"],
                               covered=None if interval is None else interval[0] <= target <= interval[1],
                               wall_seconds=wall_seconds, cpu_seconds=cpu_seconds)
                    np.savez_compressed(output / "fields" / f'{plan["seed"]}.npz',
                                        xy=xy, codes=codes)
                    rows.write(json.dumps(row, allow_nan=False) + "\n")
                    state["completed"] += 1
                    state["elapsed_seconds"] = time.time() - start
                    save(output / "STATE.json", state)
                    if state["completed"] % 30 == 0:
                        print(f'{state["completed"]}/2400 fields saved', flush=True)
        state["status"] = "complete"
    except BaseException as exc:
        state.update(status="failed_or_partial", error=repr(exc))
        raise
    finally:
        signal.alarm(0)
        state["ended_epoch"] = time.time()
        state["elapsed_seconds"] = time.time() - start
        save(output / "STATE.json", state)


def summarize(output, replay):
    state = json.loads((output / "STATE.json").read_text())
    if state["status"] != "complete" or state["completed"] != 2400:
        raise ValueError("Complete fixed study required")
    rows = [json.loads(line) for line in (output / "FIELDS.jsonl").read_text().splitlines()]
    planned = list(schedule())
    if len(rows) != len(planned):
        raise ValueError("Wrong field roster")
    replays = []
    for row, plan in zip(rows, planned):
        if any(row[k] != v for k, v in plan.items()):
            raise ValueError("Schedule mismatch")
        with np.load(output / "fields" / f'{plan["seed"]}.npz') as arrays:
            if catalogue_hash(arrays["xy"], arrays["codes"]) != row["catalogue_sha256"]:
                raise ValueError("Saved catalogue hash mismatch")
            if replay and plan["replicate"] == 1:
                xy, codes, target = generate(plan)
                if not np.array_equal(xy, arrays["xy"]) or not np.array_equal(codes, arrays["codes"]):
                    raise ValueError("Regenerated catalogue mismatch")
                if target != row["target"]:
                    raise ValueError("Regenerated target mismatch")
                result = field_check(xy, codes, code_count=4, radius_nm=30., seed=plan["seed"])
                if result != row["result"]:
                    raise ValueError("Replayed field interval mismatch")
                replays.append(plan["seed"])
    cells = []
    for p in RETENTION:
        for count in COUNTS:
            for alternative in (False, True):
                group = [r for r in rows if (r["retention"], r["expected_records"], r["alternative"])
                         == (p, count, alternative)]
                supported = [r for r in group if r["flag"] is not None]
                flags = sum(r["flag"] for r in supported)
                ci = binomtest(flags, len(supported)).proportion_ci(method="exact") if supported else (None, None)
                cells.append(dict(retention=p, expected_records=count, alternative=alternative,
                                  planned=REPEATS, completed=len(group), supported=len(supported),
                                  alerts=flags, ci95=list(ci),
                                  covered=sum(r["covered"] is True for r in group),
                                  mean_target=float(np.mean([r["target"] for r in group])),
                                  mean_records=float(np.mean([r["records"] for r in group]))))
    timings = []
    for count in COUNTS:
        group = [r for r in rows if r["expected_records"] == count]
        timings.append(dict(expected_records=count, calls=len(group),
                            records_quantiles=np.quantile([r["records"] for r in group], [.25,.5,.75]).tolist(),
                            wall_quantiles=np.quantile([r["wall_seconds"] for r in group], [.25,.5,.75]).tolist(),
                            cpu_quantiles=np.quantile([r["cpu_seconds"] for r in group], [.25,.5,.75]).tolist()))
    summary = dict(status="complete", cells=cells, timings=timings,
                   regenerated_and_replayed_seeds=replays, independent_fields=len(rows),
                   verified_catalogue_hashes=len(rows))
    save(output / "SUMMARY.json", summary)
    return summary

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--summarize", type=Path)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    if args.output:
        run(args.output)
    else:
        print(json.dumps(summarize(args.summarize, args.replay), indent=2))
