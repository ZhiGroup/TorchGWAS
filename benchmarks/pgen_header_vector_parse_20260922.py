"""Ordered real-PGEN comparison of blockwise and global index assembly."""
import argparse
import json
from pathlib import Path
import struct
import time

import numpy as np

from torchgwas.pgen_reader import (PGEN_MAGIC, VBLOCK_SIZE, PgenFormatError,
    PgenHeader, read_header, vrtype_index_layout, _decode_lengths)


def vector_header(path):
    with open(path,'rb') as stream:
        head=stream.read(12)
        if head[:2]!=PGEN_MAGIC:raise PgenFormatError('Invalid PGEN magic')
        mode=head[2]
        if mode not in (0x10,0x11):raise PgenFormatError('Unsupported mode')
        m,n=struct.unpack('<II',head[3:11])
        ctrl=head[11];bits,width=vrtype_index_layout(ctrl)
        blocks=(m+VBLOCK_SIZE-1)//VBLOCK_SIZE
        raw=stream.read(8*blocks)
        if len(raw)!=8*blocks:raise PgenFormatError('Truncated vblock offsets')
        block_offsets=np.frombuffer(raw,dtype='<u8').copy()
        types=[];lengths=[]
        for index in range(blocks):
            count=min(VBLOCK_SIZE,m-index*VBLOCK_SIZE)
            raw=stream.read((count*bits+7)//8)
            packed=np.frombuffer(raw,dtype=np.uint8)
            if bits==4:
                unpacked=np.empty(len(packed)*2,dtype=np.uint8)
                unpacked[0::2]=packed&15;unpacked[1::2]=packed>>4
                types.append(unpacked[:count])
            else:
                types.append(packed)
            lengths.append(_decode_lengths(stream.read(count*width),count,width))
    vrtypes=np.concatenate(types)
    record_lengths=np.concatenate(lengths)
    record_offsets=np.empty(m,dtype=np.uint64)
    record_offsets[0]=block_offsets[0]
    np.cumsum(record_lengths[:-1],dtype=np.uint64,out=record_offsets[1:])
    record_offsets[1:]+=block_offsets[0]
    starts=np.arange(blocks,dtype=np.int64)*VBLOCK_SIZE
    if not np.array_equal(record_offsets[starts],block_offsets):
        raise PgenFormatError('Variant-block length/index mismatch')
    return PgenHeader(storage_mode=mode,variant_ct=m,sample_ct=n,header_ctrl=ctrl,
        vrtype_bits=bits,length_bytes=width,vblock_offsets=block_offsets,
        vrtypes=vrtypes,record_offsets=record_offsets,record_lengths=record_lengths)


def timed(call):
    w=time.perf_counter();c=time.process_time()
    value=call()
    return value,dict(wall_seconds=time.perf_counter()-w,cpu_seconds=time.process_time()-c)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--path',required=True);p.add_argument('--report',required=True)
    args=p.parse_args();path=Path(args.path).resolve(strict=True)
    observations=[]
    for name,call in [('original',lambda:read_header(path)),
                      ('vector',lambda:vector_header(path)),
                      ('vector',lambda:vector_header(path)),
                      ('original',lambda:read_header(path))]:
        header,timing=timed(call)
        observations.append(dict(name=name,**timing,header=header))
    reference=observations[0]['header']
    for record in observations[1:]:
        other=record['header']
        for field in ('storage_mode','variant_ct','sample_ct','header_ctrl',
                      'vrtype_bits','length_bytes'):
            if getattr(reference,field)!=getattr(other,field):
                raise RuntimeError('Header field differs: '+field)
        for field in ('vblock_offsets','vrtypes','record_offsets','record_lengths'):
            if not np.array_equal(getattr(reference,field),getattr(other,field)):
                raise RuntimeError('Header vector differs: '+field)
    report=dict(path=str(path),samples=reference.sample_ct,markers=reference.variant_ct,
        observations=[{k:v for k,v in row.items() if k!='header'} for row in observations],
        exact_vectors=True,scope='Ordered single-process parser comparison under shared load; no GWAS speedup claim.')
    dest=Path(args.report);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2))
    print(json.dumps(report))


if __name__=='__main__':main()
