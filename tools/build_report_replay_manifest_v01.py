#!/usr/bin/env python3
"""Bind completed mode/contrast products to explicit directional DE inputs."""
import argparse,csv,json
from pathlib import Path

def build(results,de_manifest,output):
 root=Path(results)
 if not (root/'qualification.json').is_file():raise ValueError('Qualified input release receipt missing')
 receipt=json.loads((root/'qualification.json').read_text())
 if not receipt.get('status','').startswith('QUALIFIED'):raise ValueError('Input release is not qualified')
 with Path(de_manifest).open() as f:contrasts=list(csv.DictReader(f,delimiter='\t'))
 rows=[]
 for mode in sorted(p for p in root.iterdir() if p.is_dir()):
  for contrast in contrasts:
   table=mode/(contrast['source']+'_to_'+contrast['target'])/'tables'
   if not table.is_dir():continue
   ranking=table/'influence.txt';network=table/'influence_diffnetwork.tsv';de=Path(contrast['path'])
   if not all(p.is_file() for p in [ranking,network,de]):raise ValueError('Required report input missing')
   rows.append(dict(mode=mode.name,source=contrast['source'],target=contrast['target'],ranking=str(ranking),diffnetwork=str(network),de=str(de)))
 if not rows:raise ValueError('No qualified ranking products discovered')
 output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
 with output.open('w') as f:
  writer=csv.DictWriter(f,['mode','source','target','ranking','diffnetwork','de'],delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
 return len(rows)
if __name__=='__main__':
 p=argparse.ArgumentParser()
 for key in ['results','de-manifest','output']:p.add_argument('--'+key,required=True)
 print(build(**vars(p.parse_args())))
