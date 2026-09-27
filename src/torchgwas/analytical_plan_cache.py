"""Bounded reuse of deterministic analytical plans, never scan observations."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
import time

FORMAT = 'torchgwas-analytical-plan-cache-v1'
MAX_ENTRY_BYTES = 64 << 20
MAX_ENTRIES = 32
_NAME = re.compile(r'[0-9a-f]{64}\.json\Z')


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def input_identity(path):
    path=Path(path).resolve(strict=True);stat=path.stat()
    return dict(path=str(path),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,
                ctime_ns=stat.st_ctime_ns,device=stat.st_dev,inode=stat.st_ino)


def input_is_stable(identity):
    """Do not trust timestamp keys within a filesystem timestamp tick.

    A one-second grace interval covers coarse subsecond clocks without an
    input-size-dependent hash pass. Future timestamps are ineligible too.
    Revalidate identity after planning/loading before publishing or execution.
    """
    return time.time_ns()-max(identity['mtime_ns'],identity['ctime_ns'])>=1_000_000_000


class AnalyticalPlanCache:
    """Optional cache of caller-validated plans with atomic, bounded JSON entries.

    This is a trusted local artifact, not an authentication boundary. A checksum
    detects corruption. Ordinary input edits invalidate through filesystem
    identity, including ctime; it is not a full-content input hash. No lock is
    held during planning, so concurrent misses may compute the same plan twice.
    """
    def __init__(self,directory,request):
        self.directory=Path(directory).expanduser()/'torchgwas-analytical-plans-v1'
        self.key=hashlib.sha256(canonical(dict(format=FORMAT,request=request))).hexdigest()
        self.path=self.directory/(self.key+'.json')
        self.state='miss'
        self.write_status='not_attempted'

    def load(self):
        try:
            # fstat and a bounded read also protect against a file growing after
            # stat. Caches from other versions or incomplete writes are misses.
            with self.path.open('rb') as stream:
                if os.fstat(stream.fileno()).st_size>MAX_ENTRY_BYTES:
                    self.state='oversized';return None
                raw=stream.read(MAX_ENTRY_BYTES+1)
            if len(raw)>MAX_ENTRY_BYTES:
                self.state='oversized';return None
            record=json.loads(raw)
            if (set(record)!={'format','key','plan','sha256'} or record['format']!=FORMAT
                or record['key']!=self.key or not isinstance(record['plan'],dict)
                or record['sha256']!=hashlib.sha256(canonical(record['plan'])).hexdigest()):
                self.state='invalid';return None
            plan=record['plan']
            if (not isinstance(plan.get('selected'),dict)
                or not isinstance(plan.get('search_space'),dict)
                or 'input_file_identity' not in plan['search_space']):
                self.state='invalid';return None
            try:os.utime(self.path,None)
            except OSError:pass
            self.state='hit'
            return plan
        except FileNotFoundError:
            self.state='miss'
        except (OSError,ValueError,TypeError,RecursionError):
            self.state='invalid'
        return None

    def store(self,plan):
        temporary=None
        try:
            payload=canonical(plan)
            record=canonical(dict(format=FORMAT,key=self.key,plan=plan,
                                  sha256=hashlib.sha256(payload).hexdigest()))
            if len(record)>MAX_ENTRY_BYTES:
                self.write_status='oversized';return
            self.directory.mkdir(parents=True,exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=self.directory,prefix='.pending-',suffix='.json',delete=False) as stream:
                temporary=Path(stream.name);stream.write(record);stream.flush();os.fsync(stream.fileno())
            os.replace(temporary,self.path);temporary=None
            # Only recognized entries in this version's dedicated directory are
            # evicted. Racing readers handle removal as a cache miss.
            entries=[]
            for path in self.directory.iterdir():
                if path.is_file() and _NAME.fullmatch(path.name):
                    try:entries.append((path.stat().st_mtime_ns,path))
                    except FileNotFoundError:pass
            excess=max(0,len(entries)-MAX_ENTRIES)
            for _,path in sorted(entries,key=lambda row:(row[0],row[1].name)):
                if not excess:break
                if path==self.path:continue
                try:path.unlink(missing_ok=True);excess-=1
                except OSError:pass
            self.write_status='stored'
        except (OSError,ValueError,TypeError,RecursionError):
            # Cache availability never changes the selected analytical answer.
            self.write_status='unavailable'
        finally:
            if temporary is not None:
                try:temporary.unlink(missing_ok=True)
                except OSError:pass

    def audit(self):
        return dict(status=self.state,write_status=self.write_status,key=self.key,
                    directory=str(self.directory),format=FORMAT,
                    max_entry_bytes=MAX_ENTRY_BYTES,max_entries=MAX_ENTRIES,
                    scope='Analytical plan reuse only; no association timings. Source/calibration, input identity, workload, search and output settings form the key. Live context and memory are checked again before execution.')
