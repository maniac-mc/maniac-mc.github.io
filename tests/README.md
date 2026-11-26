# MANIAC tests

## Layout

- `units/` - compiled Fortran unit tests (`test_*.f90`), linked against the
  maniac objects in `build/`. No external data needed.
- `readers/` - parsing tests: `prescan/` and `scan/` are compiled Fortran
  tests; `data/` and `inputs/` run the full binary on good/bad files.
- `writers/` - end-to-end tests that run maniac + LAMMPS and check the
  written topology (requires LAMMPS, set `LMP` if not in `PATH`).
- `integration/` - end-to-end runs of the `maniac` binary (energy vs
  references, GCMC moves, energy drift).
- `cli/` - command-line argument handling.

## Running

From the repo root, after `./build.sh`:

```bash
tests/run-all.sh            # everything
cd tests/units && ./run-test.sh   # one suite
```

`run-all.sh` exits non-zero if any suite fails. Writer tests are skipped
with a warning if no LAMMPS binary is found.

## Writing a compiled unit test

Drop a `test_*.f90` file containing a `program` into `tests/units/` (or any
suite using `lib.sh`). The runner compiles it with the maniac modules and
links it against all objects except `main.o`. Fail with `stop 1`
(preferably `error stop 1`); exiting 0 means PASS. Seed the RNG with
`call seed_rng(<fixed>)` (from `random_utils`) so the test is reproducible.

All compiled-test runners share `tests/lib.sh`; do not copy-paste the
compile/link logic into new `run-test.sh` files.

## House rules

- Tests must be able to fail: a runner script must exit non-zero when a
  test fails.
- Keep test artifacts (`build/`, `log.*`, `fort.*`, `outputs/`) out of the
  git tree; runners clean up after themselves and `.gitignore` covers them.
