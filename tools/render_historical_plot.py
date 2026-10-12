#!/usr/bin/env python3
"""Render the archived R chunk with explicit run paths; preserve plotting logic."""
import argparse,re
from pathlib import Path

def render(template,influence,diffnetwork,deseq,tf,outdir,transition):
    code=Path(template).read_text()
    if '# - Uses FILTERED diffnetwork for PLOT (min_edge_weight)' not in code:
        raise ValueError('Expected the generic user-selected branching template')
    values=dict(influence_file=influence,diffnet_file=diffnetwork,deseq_file=deseq,tf_of_interest=tf,out_dir=outdir)
    for key,value in values.items():
        if any(c in value for c in ['"','\\','\n']):raise ValueError('Unsupported character in plot argument')
        code,n=re.subn(r'^'+key+r'\s*<-\s*.*$',key+' <- "'+value+'"',code,count=1,flags=re.M)
        if n!=1:raise ValueError('Missing assignment '+key)
    code,n=re.subn(r'^out_pdf <- .*$',f'out_pdf <- file.path(out_dir, paste0("ANANSE_2step_branch_", tf_of_interest, "_{transition}.pdf"))',code,count=1,flags=re.M)
    if n!=1:raise ValueError('Missing PDF assignment')
    source,target=transition.split('_to_',1)
    import json
    return 'source_label <- '+json.dumps(source)+'\ntarget_label <- '+json.dumps(target)+'\n'+code+'\n'

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['template','influence','diffnetwork','deseq','tf','outdir','transition','output']:p.add_argument('--'+key,required=True)
    args=vars(p.parse_args());output=args.pop('output');Path(output).write_text(render(**args))
