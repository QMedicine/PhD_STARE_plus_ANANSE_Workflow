#!/usr/bin/env python3
"""Create a local navigation page; each linked report is independently portable."""
import argparse,html,json
from pathlib import Path

def build(root):
 root=Path(root);receipt=json.loads((root/'selected_report_qualification.json').read_text())
 if receipt['status']!='QUALIFIED_SELECTED_REPORTS':raise ValueError('Unqualified report collection')
 rows=[]
 for record in receipt['records']:
  if not (root/record['report']).is_file():raise ValueError('Report missing')
  settings=record['settings']
  values=[record['mode'],record['source']+' → '+record['target'],record['variant'],
          'edge ≥ '+str(settings['min_edge_weight_plot'])+'; padj ≤ '+str(settings['padj_cutoff'])+'; |logFC| ≥ '+str(settings['abs_lfc_cutoff'])]
  rows.append('<tr>'+''.join('<td>'+html.escape(x)+'</td>' for x in values)+'<td><a href="'+html.escape(record['report'],quote=True)+'">Open self-contained report</a></td></tr>')
 body='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Selected plotting reports</title><style>body{font:16px system-ui;max-width:1400px;margin:2em auto;padding:1em}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:10px}p{line-height:1.5}</style><h1>Selected plotting reports</h1><p>Each HTML file embeds 20 focal PDFs and one top-10 overview PDF. Main and alternate variants share the same ranking scores. Known empty filtered graphs show centred warning panels; missing DE values appear grey. Plot settings, executed scripts, node/edge exports and PDFs accompany the reports.</p>'''
 body+='<p>Qualified reports: '+str(receipt['report_count'])+'; embedded PDFs: '+str(receipt['pdf_count'])+'.</p><table><tr><th>ATAC mode</th><th>Contrast</th><th>Variant</th><th>Display thresholds</th><th>Report</th></tr>'+''.join(rows)+'</table></html>'
 (root/'index.html').write_text(body)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root');build(p.parse_args().root)
