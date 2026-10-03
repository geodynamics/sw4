#!/usr/bin/env python3
"""Audit executed interface propagation inputs; case names do not prove coverage."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
from review_regressions import interface_selection


def scientific_identity(work):
    lines=[]
    for line in (work/'case.in').read_text().splitlines():
        tokens=line.split('#',1)[0].split()
        if not tokens or tokens[0]=='fileio':continue
        fields={}
        for token in tokens[1:]:
            key,sep,value=token.partition('=')
            if not sep:raise ValueError(f'Unexpected fixture input: {line}')
            if tokens[0]=='rechdf5' and key=='outfile':continue
            if tokens[0]=='rechdf5' and key=='infile':
                # Compare station scientific content, not HDF5 container bytes.
                with h5py.File(work/value) as stream:
                    stations={name:{field:stream[name][field][()].tolist()
                                    for field in stream[name]} for name in stream}
                value=stations
            else:
                try:value=float(value)
                except ValueError:pass
            fields[key]=value
        lines.append((tokens[0],fields))
    canonical=json.dumps(lines,sort_keys=True,separators=(',',':'),allow_nan=False)
    return hashlib.sha256(canonical.encode()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-dir',type=Path,action='append',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    inventory=[];identities={};positions=set()
    for root in args.evidence_dir:
        results=json.loads((root/'comparison.json').read_text())
        for name in results:
            if not name.startswith('interface-'):continue
            work=root/name/'cpu'
            selection=interface_selection(name,work/'run.log')
            identity=scientific_identity(work)
            if identity in identities:
                raise ValueError(f'Duplicate scientific setup: {identities[identity]} and {root/name}')
            identities[identity]=str(root/name)
            gpu=root/name/'gpu'
            if not gpu.is_dir():raise ValueError(f'Missing GPU comparison: {gpu}')
            if scientific_identity(gpu)!=identity or interface_selection(name,gpu/'run.log')!=selection:
                raise ValueError(f'CPU/GPU scientific inputs or grids differ: {name}')
            kind='force' if name.startswith('interface-force-') else 'moment'
            positions.add((kind,selection['depth']))
            inventory.append({'name':name,'root':str(root),'scientific_sha256':identity,**selection})
    if len(inventory)!=22:raise ValueError(f'Expected 22 unique interface setups, found {len(inventory)}')
    for kind in ('force','moment'):
        for depth in (2000.-1e-7,2000.,2000.+1e-7):
            if (kind,depth) not in positions:raise ValueError(f'Missing {kind} at depth {depth:.17g}')
        selected=[item for item in inventory if item['name'].startswith(f'interface-{kind}-')]
        stencils={(item['grid'],item['kc'] if item['grid']==0 else item['kc']-item['Nz'])
                  for item in selected}
        if not {(0,1),(0,2),(0,3),(1,-3),(1,-2),(1,-1)}<=stencils:
            raise ValueError(f'Missing special interface stencils for {kind}: {stencils}')
    args.output.write_text(json.dumps({'unique_interface_setups':22,'cases':inventory},indent=2)+'\n')
    print('PASS: 22 unique scientific setups; exact interface and both epsilon sides present')


if __name__=='__main__':main()
