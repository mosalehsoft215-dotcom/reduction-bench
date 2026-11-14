#!/usr/bin/env python3
"""Reproduce and record the 30-seed robustness check; see reseed_v3.py for one worker."""
import json,subprocess,sys,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
rows=[]
for seed in range(94001,94031):
    p=subprocess.run([sys.executable,str(ROOT/'authoring/evidence/reseed_v3.py'),str(seed)],capture_output=True,text=True,check=True)
    rows.append(json.loads(p.stdout.strip().splitlines()[-1]))
summary={'count':len(rows),'rmse_max_um':max(x['rmse'] for x in rows),'p95_max_um':max(x['p95'] for x in rows),'rmse_median_um':statistics.median(x['rmse'] for x in rows),'p95_median_um':statistics.median(x['p95'] for x in rows),'gauge_rms_max_um':max(x['gauge_rms'] for x in rows),'gauge_rms_median_um':statistics.median(x['gauge_rms'] for x in rows),'all_signs_correct':all(x['sign']==-1 for x in rows)}
ev=ROOT/'authoring/evidence'
(ev/'noise_realisations_raw.jsonl').write_text('\n'.join(json.dumps(x,separators=(',',':')) for x in rows)+'\n',encoding='utf-8')
(ev/'noise_realisations.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
