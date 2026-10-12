#!/usr/bin/env python3
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

def compare(expected,observed):
    def read(p):return pd.read_csv(p,sep=',' if str(p).endswith('.csv') else '\t',index_col=0)
    a,b=read(expected),read(observed)
    if a.index.has_duplicates or b.index.has_duplicates:raise ValueError('Duplicated table index')
    same_keys=set(a.index)==set(b.index);same_columns=set(a.columns)==set(b.columns)
    result=dict(expected=str(expected),observed=str(observed),expected_sha256=hashlib.sha256(Path(expected).read_bytes()).hexdigest(),observed_sha256=hashlib.sha256(Path(observed).read_bytes()).hexdigest(),same_keys=same_keys,same_columns=same_columns,expected_shape=list(a.shape),observed_shape=list(b.shape))
    if same_keys and same_columns:
        b=b.loc[a.index,a.columns]
        av=a.apply(pd.to_numeric,errors='raise').to_numpy(dtype=float);bv=b.apply(pd.to_numeric,errors='raise').to_numpy(dtype=float)
        equal=(av==bv)|(np.isnan(av)&np.isnan(bv));finite=np.isfinite(av)&np.isfinite(bv);delta=np.abs(av[finite]-bv[finite])
        result.update(exact_numeric=bool(equal.all()),n_different_values=int((~equal).sum()),max_absolute_error=float(delta.max(initial=0)),max_relative_error=float((delta/np.maximum(np.abs(av[finite]),1e-15)).max(initial=0)),within_1e_12=bool(np.allclose(av,bv,rtol=1e-12,atol=1e-12,equal_nan=True)))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--expected',required=True);p.add_argument('--observed',required=True);p.add_argument('--output',required=True);a=p.parse_args();Path(a.output).write_text(json.dumps(compare(a.expected,a.observed),indent=2)+'\n')
