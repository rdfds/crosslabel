"""Compare field area at fixed retained source density."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import signal
import time

import numpy as np
from scipy.stats import binomtest

from crosslabel.field import field_check
from area_generator import generate
from collect_robustness import equivalent

HERE = Path(__file__).resolve().parent
BASE = HERE/"data/retention/FIELDS.jsonl"
RETENTION = (.25, .5, .75)
AREAS = (100., 200., 400.)


def save(path, data):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    temporary.replace(path)


def digest(xy, codes):
    result = hashlib.sha256(np.ascontiguousarray(xy, dtype="<f8").tobytes())
    result.update(np.ascontiguousarray(codes, dtype="<i8").tobytes())
    return result.hexdigest()




def schedule():
    for replicate in range(1, 81):
        for ai in (1, 2):
            for pi, p in enumerate(RETENTION):
                for alternative in (False, True):
                    yield dict(area=AREAS[ai], retention=p, replicate=replicate,
                               alternative=alternative,
                               entropy=[19112026, ai, pi, int(alternative), replicate],
                               seed=19000000+100000*ai+10000*pi+1000*int(alternative)+replicate)


def field(plan):
    config = dict(shape="triangle", spacing=12., K=4, mask=False,
                  density=70., p=plan["retention"], side_nm=1000*np.sqrt(plan["area"]))
    fractions = (.7, .2, .1) if plan["alternative"] else (.5, .5, 0.)
    return generate(config, fractions, np.random.default_rng(plan["entropy"]))


def baseline():
    rows = [json.loads(s) for s in BASE.read_text().splitlines()]
    rows = [r for r in rows if r["retention"] in RETENTION and r["expected_records"]==7000]
    if len(rows) != 480:
        raise ValueError("Wrong reused baseline roster")
    return rows


def engineering():
    for row in baseline():
        xy, codes, target = field(dict(row, area=100.))
        if digest(xy,codes) != row["catalogue_sha256"] or target!=row["target"]:
            raise ValueError("Baseline catalogue differs under generalized area generator")
    assert len(list(schedule())) == 960
    assert len({r["seed"] for r in schedule()}) == 960
    print("All 480 baseline catalogues and targets agree; 960 new seeds are distinct")


def run(output):
    output.mkdir()
    (output/"fields").mkdir()
    state = dict(status="running", started_epoch=time.time(), new_completed=0, reused=480)
    save(output/"STATE.json",state)
    def expired(*_):
        raise TimeoutError("Fixed 45-minute area-study limit")
    signal.signal(signal.SIGALRM,expired)
    signal.alarm(2700)
    try:
        with (output/"collector.lock").open("x") as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            with (output/"FIELDS.jsonl").open("x",buffering=1) as stream:
                for old in baseline():
                    stream.write(json.dumps(dict(old,area=100.,reused=True),allow_nan=False)+"\n")
                for plan in schedule():
                    if shutil.disk_usage(output).free < 2<<30:
                        raise OSError("Two GiB reserve reached")
                    xy,codes,target = field(plan)
                    result=field_check(xy,codes,code_count=4,radius_nm=30.,seed=plan["seed"])
                    interval=result["interval_95_model"]
                    row=dict(plan,records=len(xy),target=target,result=result,
                             catalogue_sha256=digest(xy,codes),reused=False,
                             flag=result["positive_lower_limit"],
                             covered=None if interval is None else interval[0]<=target<=interval[1])
                    np.savez_compressed(output/"fields"/f'{plan["seed"]}.npz',xy=xy,codes=codes)
                    stream.write(json.dumps(row,allow_nan=False)+"\n")
                    state["new_completed"]+=1
                    save(output/"STATE.json",state)
                    if state["new_completed"]%24==0:
                        print(f'{state["new_completed"]}/960 new fields saved',flush=True)
        state["status"]="complete"
    except BaseException as exc:
        state.update(status="failed_or_partial",error=repr(exc))
        raise
    finally:
        signal.alarm(0)
        state["ended_epoch"]=time.time()
        save(output/"STATE.json",state)


def summarize(folder,replay):
    state=json.loads((folder/"STATE.json").read_text())
    if state["status"]!="complete" or state["new_completed"]!=960:
        raise ValueError("Complete area study required")
    rows=[json.loads(s) for s in (folder/"FIELDS.jsonl").read_text().splitlines()]
    if len(rows)!=1440:
        raise ValueError("Wrong combined field count")
    replayed=[]
    for row in rows:
        if row["reused"]:
            continue
        with np.load(folder/"fields"/f'{row["seed"]}.npz') as arrays:
            if digest(arrays["xy"],arrays["codes"])!=row["catalogue_sha256"]:
                raise ValueError("Raw catalogue hash mismatch")
            if replay and row["replicate"]==1:
                xy,codes,target=field(row)
                if not np.array_equal(xy,arrays["xy"]) or not np.array_equal(codes,arrays["codes"]) or target!=row["target"]:
                    raise ValueError("Regenerated catalogue differs")
                result=field_check(xy,codes,code_count=4,radius_nm=30.,seed=row["seed"])
                if not equivalent(result,row["result"]):
                    raise ValueError("Complete field interval replay differs")
                replayed.append(row["seed"])
    groups=[]
    for p in RETENTION:
        for area in AREAS:
            for alternative in (False,True):
                group=[r for r in rows if (r["retention"],r["area"],r["alternative"])==(p,area,alternative)]
                assert len(group)==80
                support=[r for r in group if r["flag"] is not None]
                k=sum(r["flag"] for r in support)
                ci=list(binomtest(k,len(support)).proportion_ci(method="exact")) if support else [None,None]
                groups.append(dict(retention=p,area=area,alternative=alternative,planned=80,
                                   supported=len(support),alerts=k,ci95=ci,
                                   covered=sum(r["covered"] for r in support),
                                   mean_records=float(np.mean([r["records"] for r in group]))))
    save(folder/"SUMMARY.json",dict(groups=groups,new_fields=960,reused_fields=480,replayed_seeds=replayed))
    print(json.dumps(groups,indent=2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--engineering", action="store_true")
    mode.add_argument("--output", type=Path)
    mode.add_argument("--summarize", type=Path)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    if args.engineering:
        engineering()
    elif args.output:
        run(args.output)
    else:
        summarize(args.summarize, args.replay)
