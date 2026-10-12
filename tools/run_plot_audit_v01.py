#!/usr/bin/env python3
import argparse,hashlib,json,subprocess,csv
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--image',required=True);a=p.parse_args();root=Path(a.root)
results=[]
for plan in json.loads((root/'plan.json').read_text()):
    folder=(root/plan['script']).parent
    command=['apptainer','exec','--cleanenv','--no-home','--bind','/vast/projects/Immune_Aging,/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow',a.image,'Rscript','--vanilla',str(root/plan['script'])]
    with (folder/'run.log').open('w') as log:r=subprocess.run(command,cwd=folder,stdout=log,stderr=subprocess.STDOUT)
    pdf=root/plan['pdf'];valid=r.returncode==0 and pdf.is_file() and pdf.stat().st_size>1000 and pdf.read_bytes().startswith(b'%PDF')
    result=dict(plan,exit_code=r.returncode,status='PDF_CREATED' if valid else 'FAILED',pdf_bytes=pdf.stat().st_size if pdf.exists() else 0,pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest() if valid else None,command=command)
    results.append(result);(root/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(plan['id'],result['status'],flush=True)
with (root/'manifest.tsv').open('w') as f:
    fields=['id','tf','direction','section','source_notebook','source_lines','selection','layout','status','exit_code','pdf','script','pdf_sha256'];w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(results)
(root/'README.md').write_text('# Plot audit: Maf and Gfi1\n\nEvery active ggsave section in the archived notebooks was run independently. Young2Old has five plotting sections; Old2Young has one. Selectable layouts were enumerated. The first Young2Old section is shown both with its original top-10 selection and with explicit Maf-only selection.\n\nUse manifest.tsv to map each PDF to its notebook, exact source lines, selection and layout. Each folder contains executed.R and run.log; results.json records all replacements, inputs, commands and checksums. H13CORE DE tables supply all plots, as requested. Original filters and aesthetics remain section-specific. Undefined interactive colour limits are supplied explicitly and recorded; audit-only seed=1 is recorded. Failed variants remain labelled FAILED rather than replaced by another graph or relaxed filters. Production plotting/report policy is unchanged.\n\nThe package contains PDFs and executable source for inspection; its absolute HPC inputs are recorded for provenance, not shipped as biological data.\n')
