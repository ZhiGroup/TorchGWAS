"""Per-thread self and inclusive time from a py-spy speedscope profile.

    python benchmarks/speedscope_threads.py profile.json [--min-samples 200] [--top 10] [--match torchgwas]
"""
import argparse
import collections
import json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('profile')
    parser.add_argument('--min-samples', type=float, default=200)
    parser.add_argument('--top', type=int, default=10)
    parser.add_argument('--match', action='append', default=None,
                        help='substring of a frame file to list inclusively (repeatable)')
    parser.add_argument('--thread', default=None, help='substring of the thread name to show')
    args = parser.parse_args()
    matches = args.match or ['torchgwas']
    data = json.load(open(args.profile))
    frames = data['shared']['frames']

    def label(index):
        frame = frames[index]
        return f"{frame['name']} ({frame.get('file', '').split('/')[-1]}:{frame.get('line', '')})"

    for profile in data['profiles']:
        samples, weights = profile['samples'], profile['weights']
        total = sum(weights)
        if total < args.min_samples or (args.thread and args.thread not in profile['name']):
            continue
        own, inclusive = collections.Counter(), collections.Counter()
        for stack, weight in zip(samples, weights):
            if stack:
                own[stack[-1]] += weight
            for index in set(stack):
                inclusive[index] += weight
        print(f"===== {profile['name']}  samples {total:g}")
        print('  self:')
        for index, count in own.most_common(args.top):
            print(f'   {count:8g} {100 * count / total:5.1f}%  {label(index)}')
        print('  inclusive:')
        shown = 0
        for index, count in inclusive.most_common():
            if any(match in frames[index].get('file', '') for match in matches):
                print(f'   {count:8g} {100 * count / total:5.1f}%  {label(index)}')
                shown += 1
                if shown >= args.top:
                    break


if __name__ == '__main__':
    main()
