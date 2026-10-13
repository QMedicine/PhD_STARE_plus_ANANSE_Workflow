#!/usr/bin/env python3
"""Inventory every archived/current script and locate label-association review sites."""
import hashlib,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
EXT={'.sh','.bash','.py','.R','.r','.Rmd','.rmd','.nf','.config','.swift','.def'}
PATTERNS={
 'file_order':r'\bsort\b|\bsorted\s*\(|list\.files\s*\(|\.glob\s*\(|\bfind\s+|\bls\s+',
 'positional_bind':r'\bpaste\b|\bcbind\s*\(|bind_cols\s*\(|\bzip\s*\(|\.iloc\b|withIndex\s*\(',
 'label_assignment':r'colnames\s*\(|setnames\s*\(|rownames\s*\(|names\s*\([^)]*\)\s*(?:<-|=)|\bLABELS\s*=|\bBAMS\s*=|sample_names|sample_labels|condition.*(?:factor|rep\s*\()',
 'contrast_direction':r'contrast\s*(?:=|<-)|source\s*(?:=|<-)|target\s*(?:=|<-)|log2FoldChange.*(?:-\s*|\*\s*-1)|Young2Old|Old2Young',
 'id_rewrite':r'renameSeqlevels|seqlevels.*<-|\.str\.(?:replace|upper)|sub\s*\(|gsub\s*\(|\.rename\s*\(|basename\s*\(',
 'keyed_association':r'sample_name|cell_state|\.to_dict\s*\(|\.set_index\s*\(|\.reindex\s*\(|\.loc\b|\bmatch\s*\(|\.findAll\b|\.index\s*=|available\s*=|metadata\s*=',
}

def scan():
 paths=sorted(p for name in ['historical','bin','modules','workflows','tools','tests','templates','containers'] for p in (ROOT/name).rglob('*') if p.is_file() and p.suffix in EXT)
 paths+= [ROOT/'main.nf',ROOT/'report_replay.nf',ROOT/'nextflow.config']
 inventory=[];sites=[];groups={}
 for path in sorted(set(paths)):
  raw=path.read_bytes();sha=hashlib.sha256(raw).hexdigest();relative=str(path.relative_to(ROOT));lines=raw.decode(errors='replace').splitlines()
  scope='historical' if relative.startswith('historical/') else 'current'
  flags=[]
  for i,line in enumerate(lines,1):
   types=[key for key,pattern in PATTERNS.items() if re.search(pattern,line)]
   if types:
    flags.extend(types);sites.append(dict(path=relative,line=i,categories=types,text=line[:600]))
  role=('production' if relative.startswith(('bin/','modules/','workflows/')) or relative in ['main.nf','report_replay.nf','nextflow.config']
        else 'archived_dependency' if relative.startswith(('historical/software/','historical/upstream_diffhic_standard/'))
        else 'historical_stage_or_alternative' if scope=='historical'
        else 'container_recipe' if path.suffix=='.def'
        else 'test_or_development_tool')
  inventory.append(dict(path=relative,scope=scope,role=role,sha256=sha,lines=len(lines),review_categories=sorted(set(flags)),duplicate_group=sha,static_screen='COMPLETE',runtime_qualification='SEE_FINDINGS_AND_NATIVE_RECEIPTS_NOT_INFERRED_FROM_SCREEN'))
  groups.setdefault(sha,[]).append(relative)
 out=ROOT/'evidence/audit/label_identity_audit_v01';out.mkdir(parents=True,exist_ok=True)
 (out/'script_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
 (out/'review_sites.json').write_text(json.dumps(sites,indent=2)+'\n')
 (out/'duplicate_groups.json').write_text(json.dumps({k:v for k,v in groups.items() if len(v)>1},indent=2)+'\n')
 print(json.dumps(dict(scripts=len(inventory),unique_contents=len(groups),historical=sum(x['scope']=='historical' for x in inventory),current=sum(x['scope']=='current' for x in inventory),review_sites=len(sites)),indent=2))
if __name__=='__main__':scan()
