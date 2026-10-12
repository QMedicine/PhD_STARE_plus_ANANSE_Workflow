#!/usr/bin/env python3
"""Render every plotting section as a traceable, independently executed audit."""
import argparse,hashlib,json,re
from pathlib import Path

def make(source,out,repair_edge_attribute=False,repair_star_mode=False):
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    plans=[]
    an=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_100kb_w5000000_m5000000_H13_motifsV2_loess_offsets_CPM')
    rn=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results')
    for direction,tf in [('Young2Old','MAF'),('Old2Young','GFI1')]:
        notebook=source/'historical/scripts/09_plots'/('GRN_plotting_'+direction+'.Rmd')
        lines=notebook.read_text().splitlines();chunks=[];start=None
        for i,line in enumerate(lines,1):
            if line.startswith('```{r'):start=i
            elif line.startswith('```') and start is not None:
                code='\n'.join(lines[start:i-1])+'\n'
                if re.search(r'^ggsave\(',code,re.M):chunks.append((start+1,i-1,code))
                start=None
        stem='STARE_YoungHSC_to_OldHSC' if direction=='Young2Old' else 'STARE_OldHSC_to_YoungHSC'
        de=rn/('HSC_50kbp_window_H13CORE_motifs' if direction=='Young2Old' else 'HSC_50kbp_window_H13CORE_motifs_Young_over_Old')/'01_counts/05_deseq2'/('A_Young_HSC:B_Old_HSC' if direction=='Young2Old' else 'A_Old_HSC:B_Young_HSC')
        de=de/(de.name+'.deseq2.results.tsv')
        for section,(first,last,base_code) in enumerate(chunks,1):
            global_plot='n_tfs             <-' in base_code
            if repair_edge_attribute and not global_plot:continue
            force_plot=bool(re.search(r'^plot_layout\s*<-',base_code,re.M))
            variants=[('as_written_top10',layout) for layout in ['nicely','fr','kk','circle','random','sugiyama','star']]+[('focal_TF_only',layout) for layout in ['nicely','fr','kk','circle','random','sugiyama','star']] if global_plot else [('focal_TF',layout) for layout in (['nicely','fr','kk','circle','random'] if force_plot else ['notebook_default'])]
            for selection,layout in variants:
                if repair_star_mode and layout!='star':continue
                audit_id=f'{tf}_{direction}_section{section:02d}_L{first}-{last}_{selection}_{layout}'
                if repair_edge_attribute:audit_id+='_edge_attribute_fix'
                if repair_star_mode:audit_id+='_star_api_fix'
                folder=out/audit_id;folder.mkdir(exist_ok=True);code=base_code;changes=[]
                assignments=dict(influence_file=str(an/('TF_Influence_'+direction)/(stem+'_TF_influence.txt')),diffnet_file=str(an/('TF_Influence_'+direction)/(stem+'_TF_influence_diffnetwork.tsv')),deseq_file=str(de),out_dir=str(folder))
                if not global_plot:assignments['tf_of_interest']=tf
                if force_plot:assignments['plot_layout']=layout
                if global_plot:assignments['layout_name']=layout
                for key,value in assignments.items():
                    code,n=re.subn(r'^'+key+r'\s*<-.*$',lambda m:key+' <- '+json.dumps(value),code,count=1,flags=re.M)
                    if n!=1 and key!='out_dir':raise ValueError('Missing assignment '+key)
                    if n:changes.append(dict(assignment=key,value=value))
                if selection=='focal_TF_only':
                    code,n=re.subn(r'^n_tfs\s*<-.*$','n_tfs <- 1',code,count=1,flags=re.M)
                    code,n2=re.subn(r'(top_tfs_tbl\s*<-\s*inf\s*%>%)',r'\1\n  filter(factor == "'+tf+'") %>%',code,count=1)
                    if n!=1 or n2!=1:raise ValueError('Missing focal selection sites')
                    changes.append(dict(assignment='top_tfs_tbl',value='Restrict focal-only variant to '+tf+'; n_tfs=1'))
                pdf=folder/(audit_id+'.pdf')
                code,n=re.subn(r'^out_pdf\s*<-.*$',lambda m:'out_pdf <- '+json.dumps(str(pdf)),code,count=1,flags=re.M)
                if n!=1:raise ValueError('Missing PDF path')
                prefix='set.seed(1) # Explicit audit-only seed for comparable layouts\n'
                if repair_edge_attribute:
                    before='# Layout selection\n'
                    after='E(g)$width_raw <- E(g)$abs_weight # Audit repair: restore the raw edge attribute referenced by the original scale\n\n'+before
                    if code.count(before)!=1:raise ValueError('Unexpected edge-attribute repair site')
                    code=code.replace(before,after)
                    changes.append(dict(before=before,after=after,reason='Original graph stores abs_weight but geom_edge_link references missing width_raw; filters, weights and selection are unchanged'))
                if repair_star_mode:
                    before='layout_as_star(g, center = center_idx, mode = "out")';after='layout_as_star(g, center = center_idx)'
                    if code.count(before)!=1:raise ValueError('Unexpected star API site')
                    code=code.replace(before,after).replace('layout_as_star(g, mode = "out")','layout_as_star(g)')
                    changes.append(dict(before=before,after=after,reason='Remove unsupported mode argument from igraph star layout; preserve star center and graph edges'))
                if 'fc_limits' in code and not re.search(r'^fc_limits\s*<-',code,re.M):
                    prefix+='fc_limits <- c(-2,4)\nfc_breaks <- c(-2,-1,0,2,4)\n'
                    changes.append(dict(assignment='missing interactive fc_limits/fc_breaks',value='Explicit accepted Young-to-Old scale [-2,4], ticks [-2,-1,0,2,4]'))
                script=folder/'executed.R';script.write_text(prefix+code)
                plans.append(dict(id=audit_id,tf=tf,direction=direction,section=section,source_notebook=str(notebook.relative_to(source)),source_notebook_sha256=hashlib.sha256(notebook.read_bytes()).hexdigest(),source_lines=[first,last],selection=selection,layout=layout,script=str(script.relative_to(out)),script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),pdf=str(pdf.relative_to(out)),changes=changes,seed=1,input_files=assignments))
    (out/'plan.json').write_text(json.dumps(plans,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True);p.add_argument('--repair-edge-attribute',action='store_true');p.add_argument('--repair-star-mode',action='store_true');a=p.parse_args();make(a.source,a.out,a.repair_edge_attribute,a.repair_star_mode)
