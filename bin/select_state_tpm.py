#!/usr/bin/env python3
"""Select RNA TPM by exact state/sample identity; retain explicit historical replay."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

GENE_COLUMNS = ('gene','gene_id','Gene','Geneid','symbol','gene_name','ENSEMBL','EnsemblID','Ensembl ID')

def read_table(path):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        if not fields or len(set(fields)) != len(fields):
            raise ValueError('Empty or duplicate table headings: ' + str(path))
        rows = list(reader)
    if not rows or any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError('Empty table or malformed row: ' + str(path))
    return fields, rows

def state_vector(tpm, samplesheet, state, mode='state_specific'):
    fields, rows = read_table(tpm)
    gene = next((name for name in GENE_COLUMNS if name in fields), fields[0])
    _, samples = read_table(samplesheet)
    if not {'sample_name','assay','cell_state'}.issubset(samples[0]):
        raise ValueError('Samplesheet requires sample_name, assay and cell_state')
    names = [row['sample_name'] for row in samples]
    if len(set(names)) != len(names):
        raise ValueError('Duplicate samplesheet sample_name')
    state_samples = [row['sample_name'] for row in samples if row['assay']=='RNA' and row['cell_state']==state]
    if not state_samples:
        raise ValueError('No RNA samples for cell_state ' + state)
    if mode == 'state_specific':
        selected = state_samples
    elif mode == 'historical_pooled':
        selected = [name for name in fields if name != gene]
    else:
        raise ValueError('Unsupported TPM expression mode')
    if any(name not in fields for name in state_samples):
        raise ValueError('RNA sample missing from TPM columns: ' + ', '.join(name for name in state_samples if name not in fields))
    seen = set()
    vector = {}
    for row in rows:
        identifier = row[gene]
        if not identifier or identifier in seen:
            raise ValueError('Empty or duplicate TPM gene identifier: ' + identifier)
        seen.add(identifier)
        values = [float(row[name]) for name in selected]
        if not all(math.isfinite(value) and value >= 0 for value in values):
            raise ValueError('Non-finite or negative TPM: ' + identifier)
        vector[identifier] = sum(values)/len(values)
    return gene, selected, rows, vector

def select(tpm, samplesheet, state, output, metadata, mode='state_specific'):
    gene, selected, rows, vector = state_vector(tpm, samplesheet, state, mode)
    with Path(output).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, [gene]+selected, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows({name:row[name] for name in [gene]+selected} for row in rows)
    record = dict(tpm_expression_mode=mode, cell_state=state, selected_rna_samples=selected,
                  n_genes=len(vector), mean_operation='arithmetic within selected sample columns',
                  input_tpm_sha256=hashlib.sha256(Path(tpm).read_bytes()).hexdigest(),
                  samplesheet_sha256=hashlib.sha256(Path(samplesheet).read_bytes()).hexdigest(),
                  selected_tpm_sha256=hashlib.sha256(Path(output).read_bytes()).hexdigest())
    Path(metadata).write_text(json.dumps(record, indent=2)+'\n')
    return record

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for key in ['tpm','samplesheet','state','output','metadata']:
        parser.add_argument('--'+key, required=True)
    parser.add_argument('--mode', choices=['state_specific','historical_pooled'], default='state_specific')
    print(json.dumps(select(**vars(parser.parse_args()))))
