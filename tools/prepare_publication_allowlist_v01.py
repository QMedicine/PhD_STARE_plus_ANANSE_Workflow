#!/usr/bin/env python3
"""Audit an explicit code/provenance staging set; never stage data by directory sweep."""
import json,re,subprocess
from pathlib import Path
files=(subprocess.check_output(['git','ls-files','-z'])+subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'])).decode().split('\0')
roots=['bin/','containers/','docs/','historical/','modules/','templates/','tests/','tools/','workflows/']
standalone={'.gitattributes','.gitignore','AGENTS.md','README.md','main.nf','nextflow.config','nextflow_schema.json'}
allowed=[p for p in files if p and (p in standalone or any(p.startswith(root) for root in roots)) and not p.endswith(('.html','.pdf'))]
receipts=['evidence/audit/rna_locus_groups_v01/','evidence/audit/rna_numerical_regression_v01/','evidence/audit/native_contracts_v01/','evidence/audit/native_contracts_v02/','evidence/audit/selected_plot_v05/','evidence/audit/deseq_cellmeans_v02/','evidence/audit/tpm_boundary_v01/','evidence/audit/influence_regression_v01/','evidence/audit/rebuilt_affinity_regression_v01/','evidence/audit/saved_union_affinity_regression_v01/','evidence/audit/rna_fastq_regression_v01/','evidence/audit/rna_fastq_regression_v02/','evidence/audit/paired_geometry_ATAC/','evidence/audit/paired_geometry_RNA/','evidence/audit/qualified_bundle_v01/','evidence/audit/rna_modes_validation_v01/','evidence/audit/estimated_ice_validation_v01/','evidence/audit/verified_ATAC_input_v01/','evidence/audit/fastq_byte_diagnostic_v01/','evidence/audit/failed_fastq_identity_v01/','evidence/audit/source_snapshot_v12/','evidence/audit/atac_fastq_regression_v01/','evidence/audit/atac_header_diagnostic_v01/','evidence/audit/atac_fastq_regression_v02/','evidence/audit/hic_presplit_fixture_v01/','evidence/audit/atac_alignment_diff_v01/','evidence/audit/atac_header_seed_probe_v01/','evidence/audit/atac_fastq_regression_v03/','evidence/audit/diffhic_count_fixture_v01/','evidence/audit/native_contracts_v21/','evidence/audit/native_contracts_v28/','evidence/audit/fixed_contacts_fixture_v01/','evidence/audit/connector_fixture_v01/','evidence/audit/plot_audit_v01/']
for p in files:
 if p.startswith('evidence/audit/') and any(p.startswith(x) for x in receipts):allowed.append(p)
allowed.extend(['evidence/publication_allowlist_v01.txt','evidence/publication_audit_v01.json','evidence/native_validation_summary_v01.json','evidence/audit/curated_assets_observed_v01.jsonl','evidence/audit/discovery_manifest_v01.json'])
secret_patterns=[re.compile(x.encode()) for x in [r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',r'gh[pousr]_[A-Za-z0-9]{30,}',r'github_pat_[A-Za-z0-9_]{30,}',r'(?i)authorization:\s*bearer\s+[A-Za-z0-9._-]{20,}']]
flagged=[p for p in allowed if Path(p).is_file() and any(pattern.search(Path(p).read_bytes()) for pattern in secret_patterns)]
if flagged:raise SystemExit('Potential credentials in files: '+str(flagged))
Path('evidence/publication_allowlist_v01.txt').write_text('\n'.join(sorted(set(allowed)))+'\n')
Path('evidence/publication_audit_v01.json').write_text(json.dumps(dict(allowed_files=len(set(allowed)),allowed_bytes=sum(Path(p).stat().st_size for p in set(allowed) if Path(p).is_file()),secret_pattern_matches=len(flagged),status='PREPARED_NOT_STAGED',excluded_formats=['reads','BAM','contact matrices','images','large scientific outputs','PDF/HTML reports']),indent=2)+'\n')
print(Path('evidence/publication_audit_v01.json').read_text())
