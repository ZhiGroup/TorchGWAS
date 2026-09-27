"""Summarize py-spy raw (collapsed) profiles by thread role.

Usage: python pyspy_summary_20260924.py PROFILE.txt [--rate 250] [--wall SECONDS] [--top 12]

Roles are inferred from frames: main (run_linear_gwas), shard producer
(run_shard), tile worker (sumstats_tiled worker / _blocked), decode reader
(fill / produce_direct), hub dispatcher (_dispatch), other. For a --gil
profile, samples are GIL-holding samples, so samples / (wall * rate) is the
fraction of wall time some thread held the GIL.
"""
import argparse
import re
from collections import Counter, defaultdict

ROLES = [('shard', ('run_shard',)), ('tile', ('_blocked', 'worker (torchgwas/sumstats_tiled')),
         ('reader', ('fill (torchgwas/streaming', 'produce_direct', 'produce (torchgwas/streaming')),
         ('hub', ('_dispatch',)), ('main', ('run_linear_gwas', '<module>'))]


def role(frames):
    text = ';'.join(frames)
    for name, keys in ROLES:
        if any(key in text for key in keys):
            return name
    return 'other'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('profile')
    parser.add_argument('--rate', type=float, default=250)
    parser.add_argument('--wall', type=float)
    parser.add_argument('--top', type=int, default=12)
    args = parser.parse_args()
    by_role = Counter()
    inclusive = defaultdict(Counter)
    leaf = defaultdict(Counter)
    threads = defaultdict(set)
    total = 0
    for line in open(args.profile):
        line = line.rstrip()
        if not line or ' ' not in line:
            continue
        stack, count = line.rsplit(' ', 1)
        count = int(count)
        frames = stack.split(';')
        thread = next((f for f in frames if f.startswith('thread')), 'thread ?')
        frames = [f for f in frames if not f.startswith(('process', 'thread'))]
        name = role(frames)
        by_role[name] += count
        threads[name].add(thread)
        total += count
        for frame in set(re.sub(r':\d+\)', ')', f) for f in frames):
            inclusive[name][frame] += count
        if frames:
            leaf[name][re.sub(r':\d+\)', ')', frames[-1])] += count
    print(f'{args.profile}: {total} samples ({total/args.rate:.1f} s at {args.rate:g} Hz)')
    if args.wall:
        print(f'  samples / (wall x rate) = {total/(args.wall*args.rate):.2f}')
    for name, count in by_role.most_common():
        print(f'\n[{name}] {count} samples = {count/args.rate:.1f} thread-seconds over {len(threads[name])} threads')
        print('  inclusive:')
        for frame, value in inclusive[name].most_common(args.top):
            print(f'    {value/args.rate:7.2f} s  {frame[:110]}')
        print('  leaf (self):')
        for frame, value in leaf[name].most_common(6):
            print(f'    {value/args.rate:7.2f} s  {frame[:110]}')


if __name__ == '__main__':
    main()
