import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import unittest

EXECUTABLE = Path(sys.argv.pop(1)).resolve() if len(sys.argv) > 1 else Path('./main').resolve()
MASS = 6.646473667973e-27


def configuration(particles=4097, bins=500, threads=4, steps='3 4', resume='no',
                  dump='yes', mass=MASS, dt=2.215491295991e-7, sigma=0.001):
    return f'''# Regression fixture
Particles: {particles}
Bins: {bins}
Time step: {dt}
Mass: {mass}
Threads: {threads}
Batch steps:
 {steps}
# Steps are incremental.
Resume: {resume} input.dmp
Save: {dump} output.dmp
Wall uncertainty: {sigma}
'''


def checkpoint(positions, momenta, step=0):
    return struct.pack('=i', step) + struct.pack(f'={len(positions)}d', *positions) + struct.pack(f'={len(momenta)}d', *momenta)


def read_checkpoint(data, count):
    step = struct.unpack_from('=i', data)[0]
    positions = struct.unpack_from(f'={count}d', data, 4)
    momenta = struct.unpack_from(f'={count}d', data, 4 + 8 * count)
    return step, positions, momenta


class SimulationTest(unittest.TestCase):
    def run_case(self, config=None, data=None, seed='42', setup=None, succeeds=True):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'datos.in').write_text(config or configuration(), encoding='utf-8')
            if data is not None:
                (root / 'input.dmp').write_bytes(data)
            if setup:
                setup(root)
            env = dict(os.environ, KINETICBOX_SEED=seed, OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE')
            env.pop('OMP_THREAD_LIMIT', None)
            result = subprocess.run([str(EXECUTABLE)], cwd=root, env=env, capture_output=True, timeout=30)
            self.assertEqual(result.returncode, 0 if succeeds else 1, result.stderr.decode(errors='replace'))
            self.assertNotIn(b'Sanitizer', result.stderr)
            self.assertNotIn(b'runtime error:', result.stderr)
            files = {path.name: path.read_bytes() for path in root.iterdir() if path.is_file()}
            return result, files

    def assert_histogram(self, data, particles, bins):
        lines = data.decode().splitlines()
        self.assertEqual(len(lines), 2 * bins + 6)
        rows = [line.split() for line in lines[1:]]
        self.assertEqual(sum(int(row[1]) for row in rows), particles)
        self.assertEqual(sum(int(row[3]) for row in rows), particles)
        self.assertTrue(all(math.isfinite(float(value)) for row in rows for value in row))
        values = re.findall(r'=\s*([-+\w.eE]+)', lines[0])
        self.assertEqual(len(values), 7)
        self.assertTrue(all(math.isfinite(float(value)) for value in values))
        return rows

    def test_fresh_odd_population_and_every_batch_count(self):
        _, files = self.run_case()
        for filename in ['X0000000.dat', 'X0000003.dat', 'X0000007.dat']:
            self.assert_histogram(files[filename], 4097, 500)
        step, positions, momenta = read_checkpoint(files['output.dmp'], 4097)
        self.assertEqual(step, 7)
        self.assertTrue(all(math.isfinite(value) for value in positions + momenta))

    def test_fixed_seed_repeats_at_one_and_four_threads(self):
        for threads in [1, 4]:
            config = configuration(threads=threads)
            _, first = self.run_case(config)
            _, second = self.run_case(config)
            self.assertEqual(first['output.dmp'], second['output.dmp'])
            self.assertEqual(first['X0000007.dat'], second['X0000007.dat'])

    def test_requested_thread_count_overrides_environment(self):
        result, four = self.run_case(configuration(threads=4))
        _, one = self.run_case(configuration(threads=1))
        self.assertIn(b'requested threads = 4', result.stdout)
        self.assertNotEqual(four['output.dmp'], one['output.dmp'])

    def test_alternative_bin_counts(self):
        for bins in [2, 16, 250]:
            _, files = self.run_case(configuration(bins=bins))
            rows = self.assert_histogram(files['X0000007.dat'], 4097, bins)
            self.assertAlmostEqual(float(rows[0][0]), -0.5 - 1.0 / bins, places=5)

    def test_initial_position_diagnostic_uses_shifted_bins(self):
        for bins in [16, 500]:
            _, files = self.run_case(configuration(bins=bins))
            text = files['X0000000.dat'].decode()
            counts = [int(row.split()[1]) for row in text.splitlines()[1:]]
            expected_full = 4097 / bins
            statistic = 0.0
            for index in range(bins + 2, 2 * bins + 3):
                expected = expected_full / 2 if index in [bins + 2, 2 * bins + 2] else expected_full
                statistic += (counts[index] - expected) ** 2 / expected
            statistic /= bins + 1
            observed = float(re.search(r'chi2x =\s*([\d.eE+-]+)', text).group(1))
            self.assertAlmostEqual(observed, statistic, places=6)

    def test_resume_legacy_file_and_unambiguous_large_filenames(self):
        data = checkpoint([0.0, 0.1], [0.0, 0.0], 100000000)
        _, files = self.run_case(configuration(particles=2, resume='sí', steps='1 1'), data)
        self.assertIn('X100000001.dat', files)
        self.assertIn('X100000002.dat', files)
        self.assertEqual(read_checkpoint(files['output.dmp'], 2)[0], 100000002)

    def test_ballistic_motion_without_collision(self):
        data = checkpoint([0.0, 0.1, -0.1], [1e-24, -1e-24, 0.0], 2)
        _, files = self.run_case(configuration(particles=3, resume='yes', steps='10', sigma=0), data)
        step, positions, momenta = read_checkpoint(files['output.dmp'], 3)
        self.assertEqual(step, 12)
        self.assertAlmostEqual(positions[0], 10e-24 * 2.215491295991e-7 / MASS, places=14)
        self.assertAlmostEqual(positions[1], 0.1 - positions[0], places=14)
        self.assertEqual(momenta, (1e-24, -1e-24, 0.0))

    def test_reflection_parity_for_both_directions_and_multiple_crossings(self):
        for distance in [0.75, 1.75, 2.75]:
            data = checkpoint([0.0, 0.0], [3e-23, -3e-23])
            _, files = self.run_case(configuration(particles=2, resume='yes', steps='1', sigma=0,
                                                    dt=distance * MASS / 3e-23), data)
            _, positions, momenta = read_checkpoint(files['output.dmp'], 2)
            crossings = math.trunc(distance + 0.5)
            parity = -1 if crossings % 2 else 1
            self.assertAlmostEqual(positions[0], parity * (distance - crossings), places=13)
            self.assertAlmostEqual(positions[1], -positions[0], places=13)
            self.assertTrue(math.isclose(momenta[0], parity * 3e-23, rel_tol=1e-14))
            self.assertTrue(math.isclose(momenta[1], -momenta[0], rel_tol=1e-14))

    def test_invalid_configuration_is_rejected(self):
        cases = [dict(particles=0), dict(particles=-1), dict(bins=1), dict(bins=501),
                 dict(threads=0), dict(mass=0), dict(dt=-1), dict(sigma=-0.1),
                 dict(mass='nan'), dict(dt='inf'), dict(particles='999999999999999999999'),
                 dict(steps='0'), dict(steps='1oops'), dict(steps=''), dict(resume='maybe')]
        for change in cases:
            with self.subTest(change=change):
                result, _ = self.run_case(configuration(**change), succeeds=False)
                self.assertIn(b'Error:', result.stderr)

    def test_more_than_fifty_batches_is_rejected(self):
        result, _ = self.run_case(configuration(steps=' '.join(['1'] * 51)), succeeds=False)
        self.assertIn(b'50 batches', result.stderr)

    def test_long_or_missing_fields_are_rejected(self):
        for text in [configuration().replace('input.dmp', 'a' * 500),
                     configuration().replace('Resume: no input.dmp', 'Resume: no'),
                     configuration() + 'unexpected\n', configuration().split('Save:')[0],
                     configuration().replace('Particles: 4097', 'Particles: ' + '1' * 5000)]:
            self.run_case(text, succeeds=False)

    def test_invalid_seed_is_rejected(self):
        for seed in ['-1', ' -1', '+1', ' 1', '4294967296', '', 'abc']:
            self.run_case(seed=seed, succeeds=False)

    def test_bad_checkpoint_lengths_and_step_counts_are_rejected(self):
        good = checkpoint([0.0, 0.1], [1e-24, -1e-24])
        for data in [b'', good[:-1], good + b'extra', checkpoint([0.0, 0.1], [0.0, 0.0], -1)]:
            result, _ = self.run_case(configuration(particles=2, resume='yes'), data, succeeds=False)
            self.assertIn(b'checkpoint', result.stderr)

    def test_invalid_particle_states_are_rejected(self):
        for position, momentum in [(float('nan'), 0), (0, float('inf')), (100, 0), (0, 1e-20)]:
            result, _ = self.run_case(configuration(particles=1, resume='yes'),
                                      checkpoint([position], [momentum]), succeeds=False)
            self.assertIn(b'outside histogram support', result.stderr)

    def test_total_step_overflow_is_rejected_before_saving(self):
        result, files = self.run_case(configuration(particles=1, resume='yes', steps='1 1'),
                                      checkpoint([0], [0], 2147483646), succeeds=False)
        self.assertIn(b'INT_MAX', result.stderr)
        self.assertNotIn('output.dmp', files)

    def test_nonfinite_trajectory_is_rejected_before_saving(self):
        result, files = self.run_case(configuration(particles=1, resume='yes', steps='1', dt=1e300, mass=1e-300),
                                      checkpoint([0], [3e-23]), succeeds=False)
        self.assertIn(b'non-finite trajectory', result.stderr)
        self.assertNotIn('output.dmp', files)

    def test_negative_energy_is_reported_without_clamping(self):
        config = configuration(particles=1, resume='yes', steps='1', dt=0.1 * MASS / 1e-30, sigma=0)
        result, files = self.run_case(config, checkpoint([0.49], [1e-30]), seed='1', succeeds=False)
        self.assertIn(b'negative momentum squared', result.stderr)
        self.assertNotIn('output.dmp', files)

    def test_histogram_open_error_is_reported(self):
        result, _ = self.run_case(setup=lambda root: (root / 'X0000000.dat').mkdir(), succeeds=False)
        self.assertIn(b'histogram output', result.stderr)

    def test_checkpoint_replace_error_is_reported_and_temp_removed(self):
        result, files = self.run_case(setup=lambda root: (root / 'output.dmp').mkdir(), succeeds=False)
        self.assertIn(b'checkpoint', result.stderr)
        self.assertFalse(any('.tmp.' in name for name in files))

    def test_disabled_checkpoint_does_not_write(self):
        _, files = self.run_case(configuration(dump='no'))
        self.assertNotIn('output.dmp', files)


if __name__ == '__main__':
    unittest.main(verbosity=2)