#!/usr/bin/env python3
"""Verify copied PDFs and build a local, traceable comparison index."""
import argparse,csv,hashlib,html,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',default='plot_audit');a=p.parse_args();root=Path(a.root);rows=[]
for folder in ['initial','section01_edge_attribute_fix','section01_star_api_fix']:
    for entry in json.loads((root/folder/'results.json').read_text()):
        row=dict(entry,audit_run=folder,pdf=folder+'/'+entry['pdf'],script=folder+'/'+entry['script'])
        if row['status']=='PDF_CREATED':
            pdf=root/row['pdf']
            if hashlib.sha256(pdf.read_bytes()).hexdigest()!=row['pdf_sha256']:raise ValueError('Copied PDF checksum differs')
            row['pdf_has_embedded_font_unicode_map']=b'/ToUnicode' in pdf.read_bytes()
        elif (root/row['pdf']).exists():
            partial=root/(row['pdf']+'.failed.partial');(root/row['pdf']).rename(partial);row['failed_partial_artifact']=str(partial.relative_to(root))
        rows.append(row)
(root/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
fields=['id','tf','direction','section','source_notebook','source_lines','selection','layout','audit_run','status','exit_code','pdf','script','pdf_sha256']
with (root/'manifest.tsv').open('w') as f:w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
success=[r for r in rows if r['status']=='PDF_CREATED']
parts=['<!doctype html><meta charset="utf-8"><title>Maf and Gfi1 plot audit</title><style>body{font:16px system-ui;margin:2rem auto;max-width:1300px;padding:0 1rem}section{margin:3rem 0;border-top:1px solid #bbb}object{width:100%;height:650px}code{overflow-wrap:anywhere}table{border-collapse:collapse}td,th{border:1px solid #bbb;padding:.4rem}nav a{display:block;margin:.4rem 0}.note{background:#fff4d2;padding:1rem}</style><h1>Maf and Gfi1 plot audit</h1>',f'<p>{len(success)} verified PDFs from {len(rows)} attempts. Each plot records the notebook section, source lines, layout, executed R script and input files. H13CORE DE tables were used throughout.</p>','<p class="note">The initial section-1 attempts failed because width_raw is missing; repaired variants restore it from the existing abs_weight attribute. Star-layout variants additionally remove an unsupported mode argument. Filters and edges are unchanged. The selected older two-step section uses weighted distances for its table and displays zero neighbourhood counts for Maf despite the visible graph; the later force-layout section explicitly uses unweighted hop distances. Long original titles may clip. These properties are retained for inspection, and production reporting rules are unchanged.</p>','<p><a href="manifest.tsv">Manifest TSV</a> · <a href="results.json">Complete provenance JSON</a> · <a href="README.md">Notes</a></p><nav>']
for i,r in enumerate(success):parts.append(f'<a href="#plot{i}">{html.escape(r["id"])}</a>')
parts.append('</nav>')
for i,r in enumerate(success):
    label=r['id'];source=r['source_notebook'];first,last=r['source_lines']
    parts.append(f'<section id="plot{i}"><h2>{html.escape(label)}</h2><p>{html.escape(source)}, lines {first}–{last}. Selection: {html.escape(r["selection"])}; layout: {html.escape(r["layout"])}; audit: {html.escape(r["audit_run"])}.</p><p><a href="{html.escape(r["pdf"])}" download>Download vector PDF</a> · <a href="{html.escape(r["script"])}">Executed R script</a> · <a href="{html.escape(str(Path(r["script"]).with_name("run.log")))}">Run log</a></p><object type="application/pdf" data="{html.escape(r["pdf"])}"><a href="{html.escape(r["pdf"])}">Open PDF</a></object></section>')
parts.append('<h2>Retained failed attempts</h2><table><tr><th>Attempt</th><th>Log</th></tr>')
for r in rows:
    if r['status']=='FAILED':parts.append(f'<tr><td>{html.escape(r["id"])}</td><td><a href="{html.escape(str(Path(r["script"]).with_name("run.log")))}">Failure log</a></td></tr>')
parts.append('</table>');(root/'index.html').write_text('\n'.join(parts))
(root/'README.md').write_text('# Maf / Gfi1 plot audit\n\nOpen index.html for a comparison index with PDF viewing windows and download links, or open PDFs directly. The folder is a portable package: transfer the entire plot_audit directory.\n\n27 completed vector PDFs cover all five active plotting sections in Young2Old for Maf and the one active section in Old2Young for Gfi1, including the selectable layouts. Section 1 includes its original top-10 selection and a separate Maf-only selection. Forty-three attempts are recorded: 16 failures are retained with logs and separately labelled partial artifacts.\n\nmanifest.tsv maps plots to exact source lines. Each plot folder includes executed.R and run.log. results.json records input paths, notebook/script/PDF checksums, transformations and commands. The original notebooks remain unchanged in historical/scripts/09_plots.\n\nSection-1 repair: alias width_raw to the graph’s existing abs_weight; star API repair: remove unsupported mode argument. No edge filters, weights or graph-selection thresholds were changed by these repairs. Explicit audit seed=1 and accepted colour scales for undefined interactive variables are recorded. All DE tables are H13CORE, as requested.\n\nInspection findings: the older selected two-step section reports zero Maf neighbourhood counts because its table uses weighted distances; the later section uses unweighted hop distances. Some long original titles clip at their original output width. These are preserved for review, not silently repaired in production.\n\nAll calculations ran in the pinned R container on HPC scratch. Only requested audit plots, executed source, logs and provenance were copied to the laptop. This folder is excluded from repository publication.\n')
print(json.dumps(dict(attempts=len(rows),verified_pdfs=len(success),failed_attempts=len(rows)-len(success),all_verified_pdfs_have_font_unicode_maps=all(r['pdf_has_embedded_font_unicode_map'] for r in success)),indent=2))
