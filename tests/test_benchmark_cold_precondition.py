"""A timed scan can only start after every private input has zero resident pages."""
import ast
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock


class ColdPreconditionTest(unittest.TestCase):
    def make_checker(self, directory, observations):
        root = Path(directory)
        prefix = root/'fixture'
        for name in ('fixture.pgen', 'fixture.pvar', 'fixture.psam', 'phenotype.npy', 'covariates.npy'):
            (root/name).write_bytes(b'fixture')
        source = Path(__file__).resolve().parents[1]/'benchmarks/direct_calculator_physical_multi.py'
        if not source.exists():
            self.skipTest(f'benchmarks/{source.name} is not in this repository')
        tree =ast.parse(source.read_text())
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'cold_inputs')
        operating_system = SimpleNamespace(fsync=Mock(), posix_fadvise=Mock(), POSIX_FADV_DONTNEED=4)
        timer = SimpleNamespace(sleep=Mock())
        namespace = dict(os=operating_system, time=timer, json=json, out=root,
                         resident_pages=Mock(side_effect=observations))
        exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), 'exec'), namespace)
        return namespace['cold_inputs'], prefix, operating_system, timer

    def test_immediately_cold_needs_no_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            check, prefix, operating_system, timer = self.make_checker(directory, [(7, 10), (0, 10)]*5)
            rows = check(prefix.parent, prefix)
            self.assertEqual(len(rows), 5)
            self.assertTrue(all(row['resident_at_launch'] == 0 for row in rows))
            self.assertTrue(all(len(row['eviction_attempts']) == 1 for row in rows))
            self.assertEqual(operating_system.posix_fadvise.call_count, 5)
            timer.sleep.assert_not_called()

    def test_transient_pages_are_recorded_before_zero_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            check, prefix, operating_system, timer = self.make_checker(directory, [(7, 10), (2, 10), (0, 10)]+[(0, 10)]*8)
            rows = check(prefix.parent, prefix)
            self.assertEqual([attempt['resident_pages'] for attempt in rows[0]['eviction_attempts']], [2, 0])
            self.assertTrue(all(row['resident_at_launch'] == 0 for row in rows))
            self.assertEqual(operating_system.posix_fadvise.call_count, 6)
            timer.sleep.assert_called_once_with(0.05)

    def test_persistent_pages_stop_scan_and_preserve_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            check, prefix, operating_system, timer = self.make_checker(directory, [(7, 10)]+[(2, 10)]*5)
            with self.assertRaisesRegex(RuntimeError, 'run not started'):
                check(prefix.parent, prefix)
            evidence = json.loads((prefix.parent/'cold_precondition_failure.json').read_text())
            self.assertEqual(evidence['files'][0]['resident_at_launch'], 2)
            self.assertEqual(len(evidence['files'][0]['eviction_attempts']), 5)
            self.assertEqual(operating_system.posix_fadvise.call_count, 5)
            self.assertEqual(timer.sleep.call_count, 4)