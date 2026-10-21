# Repository Guidelines

## Project Structure & Module Organization

Goma is a finite element solver written primarily in C11 and C++17, with some
Fortran. Solver implementations live in `src/`, with headers in `include/`.
Preserve established module prefixes such as `mm_`, `rf_`, and `dp_` when
extending related functionality. `tests/` contains unit tests; `guts/` is the
regression-test submodule. `cmake/` holds dependency-discovery modules, `tpls/`
contains third-party library build tools and bundled dependencies, and
`scripts/` provides maintenance utilities. Sphinx documentation and figures
live in `docs/`.

## Build, Test, and Development Commands

Install and configure dependencies using `BUILD.md` and `tpls/README.md` first.
When using the TPL installer, source its generated `config.sh` before
configuring Goma.

- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON`: configure the solver and Catch2 tests.
- `cmake --build build -j 2`: compile the configured targets.
- `./build/goma -v`: verify the executable and inspect enabled features.
- `ctest --test-dir build --output-on-failure -j 2`: run unit tests.
- `CLANG_FORMAT=clang-format-21 scripts/run-clang-format.sh`: format C/C++ sources, headers, and tests in place; inspect the resulting diff.
- `uv sync --locked`, then `uv run make -C docs github SPHINXOPTS="-W --keep-going -n"`: install documentation dependencies and build with strict checks.

## Documentation

Extensive sphinx documentation is available in `docs/` directory and should be 
updated when adding new equations, boundary conditions, and other user options.

## Coding Style & Naming Conventions

Follow `.clang-format`: two-space indentation, no tabs, attached braces, and a
100-column limit. Use clang-format 21 to match CI. Match neighboring naming
conventions, including snake_case C functions and module-prefixed filenames.
Keep changes focused and preserve existing license notices.

## Testing Guidelines

Unit tests use Catch2 v3 with descriptive `TEST_CASE` names and tags such as
`[gds][gds_vector]`. Place tests under the corresponding subsystem in `tests/`
and register new sources in `tests/CMakeLists.txt`. Add regression coverage for
numerical behavior changes. Initialize GUTS with `git submodule update --init
guts`, then run `python3 scripts/guts.py -g ../build/goma -v --json
results.json` from `guts/`. No numeric coverage threshold is configured.

`guts/` is an extremely expensive regression suite it is better to do a quick check
with a small problem  `python3 scripts/guts.py -g ../build/goma -v --json
results.json -p p_C_T_U` before running the whole suite.

## Commit & Pull Request Guidelines

Recent commits use short, action-oriented subjects such as “Cleanup AutoDiff
files”; merged subjects commonly include a PR number. Use a concise imperative
subject. Describe the problem, resulting behavior, affected solver options, and
validation in each PR. Link relevant issues and include numerical comparisons
or documentation screenshots when they help reviewers assess the change. Ensure
applicable build, test, formatting, and documentation checks pass.
