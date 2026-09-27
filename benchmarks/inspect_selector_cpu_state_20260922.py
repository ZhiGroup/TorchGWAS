"""Read-only CPU/NUMA availability observations; no affinity or system changes."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import time


def snapshot():
    cpu={}
    for line in Path('/proc/stat').read_text().splitlines():
        parts=line.split()
        if parts and parts[0].startswith('cpu') and parts[0][3:].isdigit():cpu[int(parts[0][3:])]=[int(x) for x in parts[1:]]
    vm={key:int(value) for key,value in (line.split() for line in Path('/proc/vmstat').read_text().splitlines())
        if key.startswith(('numa_','pgmigrate','pgfault','pgmajfault'))}
    return dict(utc=datetime.now(timezone.utc).isoformat(),cpu=cpu,vm=vm)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    topology={}
    for number in sorted(os.sched_getaffinity(0)):
        base=Path('/sys/devices/system/cpu')/('cpu'+str(number))
        fields={key:(base/'topology'/key).read_text().strip() for key in ['physical_package_id','core_id','thread_siblings_list']}
        fields['nodes']=[p.name for p in base.glob('node[0-9]*')]
        topology[number]=fields
    rows=[snapshot()]
    for _ in range(7):time.sleep(1.);rows.append(snapshot())
    intervals=[]
    for previous,current in zip(rows,rows[1:]):
        busy={}
        for number in topology:
            delta=[b-a for a,b in zip(previous['cpu'][number],current['cpu'][number])]
            total=sum(delta[:8]);idle=delta[3]+delta[4]
            busy[number]=(total-idle)/total if total else None
        intervals.append(dict(busy_fraction=busy,vm_delta={key:current['vm'][key]-previous['vm'][key] for key in previous['vm']}))
    report=dict(topology=topology,snapshots=rows,intervals=intervals,loadavg=Path('/proc/loadavg').read_text().strip(),
        numa_balancing=Path('/proc/sys/kernel/numa_balancing').read_text().strip(),scope=__doc__)
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    for number in range(12,20):
        print(json.dumps(dict(cpu=number,topology=topology.get(number),busy=[row['busy_fraction'].get(number) for row in intervals])),flush=True)
    print(json.dumps(dict(numa_balancing=report['numa_balancing'],loadavg=report['loadavg'],vm_delta=intervals[-1]['vm_delta'])),flush=True)


if __name__=='__main__':main()
