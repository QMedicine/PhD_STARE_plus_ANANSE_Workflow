#!/usr/bin/env python3
"""Self-contained HTML with embedded PDF bytes and native vector PDF viewers."""
import argparse,base64,csv,hashlib,html,json
from pathlib import Path

def methods_text(context):
    """Describe the supplied route and resolved parameters, including reused checkpoints."""
    esc=html.escape
    def value(key):return esc(str(context[key]))
    if context.get('union_windows'):
        peaks='The preserved named summit-window BED was reused as an explicit peak-universe checkpoint.'
    else:
        peaks=f'Summit windows were generated from 1-bp summits by adding {value("summit_flank")} bp on each side ({2*int(context["summit_flank"])+1} bp before boundary clipping), without merging or deduplication.'
    stage=context['stage']
    if stage in ['historical_replay','network_replay']:
        route='Archived affinities or networks and directional H13CORE DE tables were reused. Upstream reconstruction was not rerun in this replay.'
    elif stage=='atac_historical_hic_checkpoint':
        route='ATAC BAMs were recounted against the supplied peak universe; preserved corrected Hi-C H5, RNA TPM and H13CORE DE tables were reused.'
    else:
        route='RNA counts, TPM and DE tables were reconstructed, and ATAC BAMs were recounted. Preserved corrected Hi-C H5 matrices were reused.'
    if context.get('process_fastq'):route+=' RNA/ATAC samples entered through their declared FASTQ or BAM route.'
    preparation=f'{route} {peaks} Blacklisted alignment overlaps were removed and BAM overlaps counted with bedtools. ATAC activities used csaw normOffsets and calculateCPM(use.offsets=TRUE, use.norm.factors=FALSE, log=FALSE). State-specific contacts informed STARE 1.0.5 with {value("chosen_res")}-bp bins, a {value("stare_window")}-bp scoring window and {value("enhancer_window")}-bp enhancer window, gABC threshold {value("gabc_threshold")}, all_tss, -q True and -d True, the configured PSEM motifs and the archived modified shell script. <a href="https://doi.org/10.1093/bioinformatics/btad062">Hecker et al., 2023</a>.'
    expression=f'RNA expression used mean merged gene length per uppercase symbol for TPM; gene filters required sum of counts ≥{value("min_count")} and mean TPM ≥{value("min_tpm")}, and TF filters used ≥{value("min_count_tf")} and ≥{value("min_tpm_tf")}, respectively. The converter consumed unfiltered TPM. DESeq2 1.34.0 (R 4.1.3) used filtered raw counts, ratio size factors, a cell-means design (~0 + condition), Wald tests, parametric dispersion fitting and ashr shrinkage. Reconstructed reverse comparisons negate the primary log2 fold changes; checkpoint replays use the explicitly supplied directional DE tables.'
    if stage in ['historical_replay','network_replay','atac_historical_hic_checkpoint']:
        expression='RNA TPM and DE were supplied as checkpoints; their paths and the configured reconstruction settings are recorded below. These filters and DE fitting were not reapplied in this replay. The documented reconstruction method is: '+expression
    expression+=' Gene-filtered counts supply DESeq2. The separately computed TF-filter files are retained for compatibility; the custom converter uses unfiltered TPM and does not consume those TF-filter files.'
    if context.get('rna_counts_mode')=='annotation_symbol_sum':
        expression+=' This run uses the separately labelled annotation_symbol_sum method: counts are summed across all annotated IDs sharing an uppercase symbol, including distinct loci. A stable representative ID indexes each summed symbol. Missing/Ensembl-only names are excluded and excluded count totals are recorded. If no reference map or lengths are supplied, these are derived from featureCounts annotations and merged exon lengths. This method does not claim equivalence to the PhD retained-ID selection.'
    return preparation,expression

def build_report(source,target,mode,ranking,samplesheet,pdfs,versions,context,software_lock,output):
    esc=html.escape
    resolved=json.loads(Path(context).read_text())
    methods,expression=methods_text(resolved)
    with Path(ranking).open() as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)[:20]
    required={'factor','influence_score','target_score_scaled','G_score_scaled'}
    if not required.issubset(fields or []):raise ValueError('Ranking component columns missing')
    if len(rows)!=20:raise ValueError('Expected 20 ranked TFs for the report')
    with Path(samplesheet).open() as f:
        samples=list(csv.DictReader(f,delimiter='\t'))
    cards=[];payload={};checksums={}
    for i,row in enumerate(rows,1):
        tf=row['factor'];file=Path(pdfs)/f'ANANSE_2step_branch_{tf}_{source}_to_{target}.pdf'
        data=file.read_bytes()
        if not data.startswith(b'%PDF-'):raise ValueError(f'Invalid PDF {file}')
        payload[str(i)]={'bytes':base64.b64encode(data).decode(),'name':file.name}
        checksums[file.name]=hashlib.sha256(data).hexdigest()
        cards.append(f'<section><h2>{i}. {esc(tf)}</h2><button type="button" class="show" data-id="{i}">Open vector PDF viewer</button> <a class="download" data-id="{i}">Download PDF</a><iframe data-id="{i}" title="{esc(tf)} differential GRN" loading="lazy"></iframe></section>')
    table='<table><thead><tr>'+''.join('<th>'+esc(x or 'row')+'</th>' for x in fields)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(r.get(x) or '')+'</td>' for x in fields)+'</tr>' for r in rows)+'</tbody></table>'
    conversion='''STARE TF–gene affinities were averaged across replicates after intersecting genes and TF columns. Each TF affinity column was min–max scaled. The same RNA TPM vector, averaged across all RNA samples supplied to the converter, was used for both networks. Edge prob = scaled_affinity × (0.5 + 0.5 × scaled_TF_TPM) × (0.5 + 0.5 × scaled_target_TPM). Constant affinity columns were set to zero and zero-affinity edges excluded. GTF mapping converted targets to uppercase symbols, removed duplicate TF suffixes and retained the maximum probability for duplicate edges. These custom probabilities are regulatory scores, without a probability-calibration claim.'''
    conversion=conversion.replace('zero-affinity edges excluded',f'edges at or below scaled binding {esc(str(resolved.get("min_binding",0.0)))} excluded')
    ranking_method=f'''ANANSE 0.5.1 selected the top {esc(str(resolved['influence_edges']))} edges before subtraction and retained positive target-minus-source differences. Candidate TFs required positive target/source log2 fold change. G_score is absolute TF log2 fold change gated by adjusted p&lt;0.05. target_score sums the analogous DE gene scores weighted by the best cumulative differential-edge path probability divided by its path length, within the two-hop ego graph. Each component is dense-ranked then min–max scaled. The two scaled components contribute equally through their sum; influence_score is the dense rank of that sum, min–max scaled again. influence_score_raw is the sum of the unscaled components and is distinct from the ranking score. Direct target and plot neighbourhood counts annotate the networks and receive no additional weight in the final ranking. Reverse comparisons use target/source DE signs in the corresponding direction. <a href="https://doi.org/10.1093/nar/gkab598">Xu et al., 2021</a>. Rankings nominate TFs for investigation; they do not establish inhibition effects or experimentally verified rejuvenation.'''
    sample_text='<p>The declared cohort is shown in full because shared RNA expression and ATAC normalization can include samples beyond the displayed transition. Replay stages reuse the corresponding upstream checkpoints.</p><ul>'+''.join('<li>'+esc(' | '.join(r[k] for k in ['sample_name','replicate','species','assay','paired','cell_type','condition','cell_state']))+'</li>' for r in samples)+'</ul>'
    payload_json=json.dumps(payload).replace("<", chr(92)+"u003c")
    methods+='</p><p>'+expression
    methods+=' <a href="https://doi.org/10.1186/s13059-014-0550-8">Love et al., 2014</a>. ATAC offset normalization follows the csaw implementation, <a href="https://doi.org/10.1093/nar/gkv1191">Lun and Smyth, 2016</a>.'
    methods+='</p><p>Plots use the selected two-step branching design. Table neighbourhood counts use the unfiltered differential network; displayed edges and gene nodes use the original edge-weight and DE filters. TFs are squares, genes are circles and node fill indicates directional log2 fold change. Edge widths use the original direct-edge or two-step product proxy. The right panel displays scaled influence components and neighbourhood counts. A top-20 TF with no displayed edges causes strict plotting failure and withholds the complete report.'
    methods+='</p><h3>Resolved run settings and image lock</h3><pre>'+esc(json.dumps(dict(settings=resolved,software_lock=json.loads(Path(software_lock).read_text())),indent=2))+'</pre><p>'
    body=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>{esc(source)} → {esc(target)}</title><style>body{{font:16px system-ui;max-width:1300px;margin:2em auto;padding:1em;color:#17212b}}table{{border-collapse:collapse;font-size:12px;display:block;overflow:auto}}td,th{{padding:6px;border:1px solid #ccd2d8}}iframe{{width:100%;height:650px;border:1px solid #bbb;display:block;margin-top:1em}}section{{margin:3em 0}}button,a.download{{cursor:pointer;padding:8px}}p{{line-height:1.6}}pre{{white-space:pre-wrap}}</style><h1>{esc(source)} → {esc(target)}</h1><p>phd_stare_plus_ananse_workflow | ATAC mode: <strong>{esc(mode)}</strong> | The same directional DE table supplies ranking and plots.</p><h2>Input samples</h2>{sample_text}<h2>Methods</h2><p>{methods}</p><p>{conversion}</p><p>{ranking_method}</p><h3>Recorded software</h3><pre>{esc(Path(versions).read_text())}</pre><h2>Top 20 TFs and component scores</h2>{table}<h2>Differential GRN plots</h2><p>The PDFs below are stored inside this HTML. The browser PDF viewer preserves vector graphics and selectable text, with zoom controls. Download is available even when the browser has no PDF viewer.</p>{''.join(cards)}<script id="pdf-data" type="application/json">{payload_json}</script><script>const payload=JSON.parse(document.getElementById('pdf-data').textContent), urls={{}};function url(id){{if(!urls[id]){{const data=Uint8Array.from(atob(payload[id].bytes), c=>c.charCodeAt(0));urls[id]=URL.createObjectURL(new Blob([data],{{type:'application/pdf'}}));}}return urls[id];}}document.querySelectorAll('a.download').forEach(a=>{{a.href=url(a.dataset.id);a.download=payload[a.dataset.id].name;}});document.querySelectorAll('button.show').forEach(b=>b.addEventListener('click',()=>{{document.querySelector('iframe[data-id="'+b.dataset.id+'"]').src=url(b.dataset.id)+'#toolbar=1';}}));const observer=new IntersectionObserver(entries=>entries.forEach(e=>{{if(e.isIntersecting){{e.target.src=url(e.target.dataset.id)+'#toolbar=1';observer.unobserve(e.target);}}}}),{{rootMargin:'200px'}});document.querySelectorAll('iframe').forEach(f=>observer.observe(f));</script></html>'''
    Path(output).write_text(body)
    Path(str(output)+'.pdf_checksums.json').write_text(json.dumps(checksums,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['source','target','mode','ranking','samplesheet','pdfs','versions','context','software-lock','output']:p.add_argument('--'+key,required=True)
    build_report(**vars(p.parse_args()))
