#!/usr/bin/env python3
"""Validate portable report payloads and directional node values before publication."""
import argparse,base64,csv,hashlib,html,json,math,re,shutil
from pathlib import Path

def sha(path):
 digest=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):digest.update(block)
 return digest.hexdigest()
def tsv(path):
 with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def qualify(results,manifest,destination=None):
 root=Path(results);expected=tsv(manifest);records=[]
 for row in expected:
  contrast=row['source']+'_to_'+row['target'];base=root/row['mode']/contrast
  variants=sorted(p for p in base.iterdir() if p.is_dir() and (p/'plots').is_dir())
  if len(variants)!=2 or not (base/'main').is_dir():raise ValueError(f'Expected main and alternate report: {base}')
  de={r['gene_id']:r['log2FoldChange'] for r in tsv(row['de'])}
  ranking=tsv(row['ranking'])[:20]
  for variant in variants:
   plots=variant/'plots';settings=json.loads((plots/'plot_settings.json').read_text())
   statuses=tsv(plots/'plot_status.tsv')
   if [r['factor'] for r in statuses]!=[r['factor'] for r in ranking]:raise ValueError('Top20 order changed')
   if any(r['plot_status'] not in ['SELECTED_FORCE_DIRECTED_PLOT','CENTER_ONLY_FILTERED_EDGES'] for r in statuses):raise ValueError('Plot failure')
   reports=list((variant/'reports').glob('*.html'))
   if len(reports)!=1:raise ValueError('Expected one HTML report per contrast/filter variant')
   report=reports[0];body=report.read_text()
   match=re.search(r'<script id="pdf-data" type="application/json">(.*?)</script>',body,re.S)
   if match is None:raise ValueError('Embedded PDF payload missing')
   payload=json.loads(match.group(1));checks=json.loads(Path(str(report)+'.pdf_checksums.json').read_text())
   has_heatmap='expression' in payload
   count=22 if has_heatmap else 21
   expected_payload={'top10',*map(str,range(1,21))}
   if has_heatmap:expected_payload.add('expression')
   if set(payload)!=expected_payload or len(checks)!=count:raise ValueError('Unexpected report PDF payloads')
   if len(list(plots.glob('*.pdf')))!=count:raise ValueError('Unexpected number of source PDFs')
   for entry in payload.values():
    data=base64.b64decode(entry['bytes'],validate=True);pdf=plots/entry['name']
    if not data.startswith(b'%PDF-') or b'%%EOF' not in data[-1024:]:raise ValueError('Invalid PDF payload')
    if hashlib.sha256(data).hexdigest()!=checks[entry['name']] or sha(pdf)!=checks[entry['name']]:raise ValueError('Embedded PDF differs from native plot')
   if len(re.findall(r'<iframe\b',body))!=count or len(re.findall(r'data-id="top10" title=',body))!=1:raise ValueError('Unexpected PDF viewer count')
   if has_heatmap:
    expression=tsv(plots/'top20_tf_expression.tsv')
    heatmap=json.loads((plots/'heatmap_settings.json').read_text())
    if [r['factor'] for r in expression]!=[r['factor'] for r in ranking]:raise ValueError('Heatmap TF ranking mismatch')
    if heatmap['width_mm']!=90 or heatmap['font_family']!='Arial':raise ValueError('Heatmap size/font mismatch')
    for r in expression:
     if r['source']!=row['source'] or r['target']!=row['target']:raise ValueError('Heatmap contrast mismatch')
     if not math.isclose(float(r['log2FoldChange']),float(de[r['factor']]),rel_tol=1e-12,abs_tol=1e-12):raise ValueError('Heatmap directional logFC mismatch')
   if re.search(r'<(?:script|iframe)\b[^>]*src=["\']https?://',body):raise ValueError('Report requires external assets')
   if html.escape(row['source']) not in body or html.escape(row['target']) not in body:raise ValueError('Contrast labels missing')
   nodes=tsv(plots/'top10_nodes.tsv');gene_count=sum(n['type']=='Gene' for n in nodes)
   if gene_count==0:raise ValueError('No non-TF targets in overview')
   for n in nodes:
    supplied=n['log2FoldChange'];expected_fc=de.get(n['name'],'NA')
    missing=lambda x:x in ['', 'NA','NaN']
    if missing(expected_fc):
     if not missing(supplied):raise ValueError('Invented fold change')
    elif missing(supplied) or not math.isclose(float(supplied),float(expected_fc),rel_tol=1e-12,abs_tol=1e-12):raise ValueError(f'Directional logFC mismatch: {n["name"]}')
   records.append(dict(mode=row['mode'],source=row['source'],target=row['target'],variant=variant.name,settings=settings,
     report=str(report.relative_to(root)),report_sha256=sha(report),pdfs=count,overview_gene_nodes=gene_count,
     center_only_TFs=[s['factor'] for s in statuses if s['plot_status']=='CENTER_ONLY_FILTERED_EDGES'],
     input_sha256={k:sha(row[k]) for k in ['ranking','diffnetwork','de']}))
 receipt=dict(status='QUALIFIED_SELECTED_REPORTS',report_count=len(records),pdf_count=sum(r['pdfs'] for r in records),records=records)
 (root/'selected_report_qualification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 if destination:
  dest=Path(destination)
  if dest.exists():raise ValueError('Refuse to overwrite an existing report release')
  for record in records:
   src=root/record['report'];relative=src.parent.parent.relative_to(root)
   shutil.copytree(src.parent.parent,dest/relative)
  shutil.copy2(root/'selected_report_qualification.json',dest/'selected_report_qualification.json')
  shutil.copy2(manifest,dest/'report_input_manifest.tsv')
  for p in dest.rglob('*'):
   if p.is_file() and p.name not in ['report_input_manifest.tsv']:
    counterpart=root/p.relative_to(dest)
    if sha(p)!=sha(counterpart):raise ValueError('Publication copy differs')
 return receipt
if __name__=='__main__':
 parser=argparse.ArgumentParser()
 parser.add_argument('--results',required=True);parser.add_argument('--manifest',required=True);parser.add_argument('--destination')
 print(json.dumps(qualify(**vars(parser.parse_args())),indent=2))
