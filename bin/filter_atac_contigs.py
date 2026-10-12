#!/usr/bin/env python3
"""Preserve the legacy MAPQ/unmapped/contig filter with explicit contig parameters."""
import argparse,subprocess

def filter_bam(bam,output,mapq,exclude,threads,mode,samtools='samtools'):
    # Apply genomic contig exclusions by exact RNAME identity, retaining valid SAM headers.
    excluded=set(filter(None,exclude.split(',')))
    reader=subprocess.Popen([samtools,'view','-@',str(threads),'-h','-q',str(mapq),'-F','4',bam],stdout=subprocess.PIPE,text=True)
    writer=subprocess.Popen([samtools,'view','-@',str(threads),'-b','-o',output,'-'],stdin=subprocess.PIPE,text=True)
    try:
        for line in reader.stdout:
            keep=not any(contig in line for contig in excluded) if mode=='legacy_line_substring' else (line.startswith('@') or line.split('\t',3)[2] not in excluded)
            if keep:writer.stdin.write(line)
        writer.stdin.close();reader.stdout.close()
        if reader.wait()!=0 or writer.wait()!=0:raise RuntimeError('samtools BAM filter failed')
    finally:
        for proc in [reader,writer]:
            if proc.poll() is None:proc.terminate()
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bam',required=True);p.add_argument('--output',required=True);p.add_argument('--mapq',type=int,default=10);p.add_argument('--exclude',default='');p.add_argument('--threads',type=int,default=1);p.add_argument('--mode',choices=['legacy_line_substring','rname'],default='legacy_line_substring');p.add_argument('--samtools',default='samtools');a=p.parse_args();filter_bam(a.bam,a.output,a.mapq,a.exclude,a.threads,a.mode,a.samtools)
