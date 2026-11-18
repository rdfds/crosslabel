"""Regenerate every area-study catalogue and verify saved counts and intervals."""
import argparse
import json
from pathlib import Path

from scipy.stats import binomtest
import area_study as study
from verify_retention import equivalent


def main(replay=False):
    folder = Path(__file__).resolve().parent / 'data/area'
    rows = [json.loads(s) for s in (folder/'FIELDS.jsonl').read_text().splitlines()]
    summary = json.loads((folder/'SUMMARY.json').read_text())
    baseline = {r['seed']: r for r in study.baseline()}
    plans = {r['seed']: r for r in study.schedule()}
    if len(rows) != 1440 or len({r['seed'] for r in rows}) != 1440:
        raise ValueError('Incomplete or repeated study roster')
    replayed = []
    for row in rows:
        if row['reused']:
            old = baseline.pop(row['seed'])
            if any(row[k] != v for k, v in old.items()) or row['area'] != 100.:
                raise ValueError('Reused baseline differs')
        else:
            plan = plans.pop(row['seed'])
            if any(row[k] != v for k, v in plan.items()):
                raise ValueError('Changed area-study seed or setting')
        xy, codes, target = study.field(row)
        if study.digest(xy, codes) != row['catalogue_sha256'] or target != row['target']:
            raise ValueError('Regenerated catalogue or target differs')
        interval = row['result']['interval_95_model']
        if interval is None or row['flag'] != (interval[0] > 0):
            raise ValueError('Interval or decision differs')
        if row['covered'] != (interval[0] <= target <= interval[1]):
            raise ValueError('Target coverage differs')
        if replay and not row['reused'] and row['replicate'] == 1:
            result = study.field_check(xy, codes, code_count=4, radius_nm=30., seed=row['seed'])
            if not equivalent(result, row['result']):
                raise ValueError('Full interval replay differs')
            replayed.append(row['seed'])
    if baseline or plans:
        raise ValueError('Missing study fields')
    for cell in summary['groups']:
        group = [r for r in rows if (r['retention'],r['area'],r['alternative']) ==
                 (cell['retention'],cell['area'],cell['alternative'])]
        k = sum(r['flag'] for r in group)
        ci = list(binomtest(k,len(group)).proportion_ci(method='exact'))
        if len(group) != 80 or k != cell['alerts'] or not equivalent(ci, cell['ci95']):
            raise ValueError('Count or binomial interval differs')
        if sum(r['covered'] for r in group) != cell['covered']:
            raise ValueError('Coverage count differs')
    print(json.dumps(dict(catalogues_regenerated=len(rows), cells_checked=len(summary['groups']),
                          full_intervals_replayed=len(replayed), status='pass'), indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true')
    main(parser.parse_args().replay)
