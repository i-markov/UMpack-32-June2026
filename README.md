![C++](https://img.shields.io/badge/C%2B%2B-20-blue)
![GCC](https://img.shields.io/badge/GCC-13%2B-orange)
![Build](https://img.shields.io/badge/build-passing-brightgreen)
![Status](https://img.shields.io/badge/status-stable-brightgreen)

# UMpack

Contact email: `igor.markov1@gmail.com`

UMpack is a collection of layout and utility packages used by tools
for hypergraph partitioning, placement and floorplanning.  This checkout
contains 32 top-level package directories with buildable libraries and
tools; the top-level Perl configuration script is the intended way to
configure and build the full set.

## Package Inventory

Top-level package directories:

`ABKCommon`, `AnalytPl`, `Capo`, `ClusteredHGraph`, `Combi`,
`CongestionMaps`, `Constraints`, `Ctainers`, `FMPart`, `FilledHier`,
`FixOr`, `GenHier`, `GeomTrees`, `Geoms`, `HGraph`, `HGraphWDims`, `MLPart`,
`ParquetDBFromRBP`, `ParquetFP`, `PartEvals`, `PartLegality`,
`Partitioners`, `Partitioning`, `PlaceEvals`, `Placement`, `RBPlace`,
`RowIroning`, `ShellPart`, `SmallPart`, `SmallPlacement`, `SmallPlacers`, `Stats`.

## Requirements

The build is driven by the checked-in Makefiles and shell wrappers.  A
typical Linux build needs:

- Perl 5
- GNU make, either as `make` or `gmake`
- `g++`
- `/bin/sh`
- standard regression utilities used by the local scripts, including
  `diff`, `wc`, `egrep`, `sed`, `cat`, `rm` and `ls`

Regression tests also need the UMpack libraries themselves.  Let
`perl config` build the component libraries first; it creates `lib/`
and symlinks each package's static and shared libraries there.  If a
test executable cannot load a shared library, set:

```sh
export LD_LIBRARY_PATH="$PWD/lib:$LD_LIBRARY_PATH"
```

## Configure

Run the configuration script from the repository root:

```sh
perl config
```

The script updates package Makefiles and `make-gnu`/`makeOpt-gnu`
wrappers for the local system, creates the top-level `lib/` symlinks,
checks the compiler, and then offers to build and test every component.
On Linux it can also use multiple compile jobs when more than one
processor is detected.

## Build

For a full build, use the configuration script and answer yes when it
asks whether to compile UMpack components:

```sh
perl config
```

When optimization is enabled, `config` switches to the `makeOpt-gnu`
wrappers.  Without optimization it uses `make-gnu`.  Individual packages
can still be rebuilt from their own directories with the same wrappers,
for example:

```sh
cd HGraph && ./makeOpt-gnu
```

## Test

The recommended path is to let `perl config` run tests after the full
library build.  Internally it invokes each package's `test` target in
dependency order, which builds the package's test executables and then
runs its local regression script.

Each component directory has a `regression` script.  These scripts run
the local test executables, write fresh output files such as `new.out`,
filter volatile lines such as timing, platform and seed information, and
compare the result with the checked-in expected output.  A successful
test reports an empty diff summary, usually as zero lines, words and
bytes in `diffs.notime`.

## Random Seeds

Several packages use the ABKCommon seed handler to make randomized code
repeatable.  The mechanism is simple:

- If `seeds.in` is present in the working directory, it overrides the
  normal nondeterministic seed source.  The first line is an in-use flag
  and the second line is the initial seed value.  The common pattern is
  `0` on the first line and a chosen seed on the second line.
- If `seeds.in` is absent, the first random object in the process picks a
  nondeterministic seed from the environment and records it in
  `seeds.out`.
- Once the initial seed is fixed, every randomized object created in that
  process consumes the same deterministic seed stream, so repeated runs
  produce the same randomized decisions and the same output.
- To replay a previous run, capture the seed file it used and feed that
  back as `seeds.in` before rerunning the test.
- On a normal exit, the seed handler rewinds `seeds.out` and marks it as
  no longer in use.  The file is not appended to at shutdown; it is the
  run-time log of the seeds that were actually consumed.
- If a process crashes before shutdown, the existing `seeds.out` may still
  contain the seed stream from that run.  Copying it to `seeds.in` can
  help reproduce the crash, because the next run will replay the same
  initial seed and the same deterministic seed sequence.

This is why many regression scripts create or copy `seeds.in` before
running randomized tests, and why `seeds.out` is treated as a transient
run artifact.

The detailed implementation commentary is in [`ABKCommon/abkseed.h`](./ABKCommon/abkseed.h)
and [`ABKCommon/seed.cxx`](./ABKCommon/seed.cxx).

To rerun one package after it has been built:

```sh
cd HGraph
./regression
```

If a package was cleaned first, rebuild its library and test executables
before running `./regression`, or run the package wrapper's `test` target:

```sh
cd HGraph
./makeOpt-gnu test
```

## CapoTest1 Outputs

`CapoTest1` is the flexible Bookshelf driver for Capo.  In its normal mode it
reads the placement described by `-f design.aux`, runs Capo, and can write
optional outputs:

- `-savePl file.pl` writes the final Bookshelf placement.
- `-plot basename` writes a gnuplot visualization with cells and nets.
- Other `-plot*` options select narrower visualizations, such as nodes only,
  rows, sites, or node names.

`CapoTest1` also has a one-dimensional ordering mode.  When `-saveOrder` is
present, the driver reads the netlist from `-f`, forces every node to unit size,
builds a single-row placement with one slot per movable node, runs Capo, and
writes the movable node names in final left-to-right order:

```sh
cd Capo
./CapoTest1.exe -f TESTS/saurabh1.aux -saveOrder ordering.out -savePl ordering.pl -plot ordering
```

The same run can save the final Bookshelf placement with `-savePl` and a
gnuplot visualization with any existing `-plot*` option.  In this mode,
`-saveOrder` is the option that selects the synthetic one-row layout; no
separate mode flag is needed.

## 2026 Modernization

Significant updates in 2026 funded by IonQ, Inc:

1. Build system modernization for GCC 13+ compatibility
2. Removal of deprecated SGI STL extensions and custom allocators
3. Migration to standard C++ containers and algorithms
4. Comprehensive test infrastructure improvements
5. Enhanced documentation (CLAUDE.md, FIXES_AND_LESSONS.md)

This work was performed by Igor L. Markov (IonQ contractor) and
Daiwei Chen (IonQ employee).

For detailed engineering notes on the modernization effort,
see [FIXES_AND_LESSONS.md](./FIXES_AND_LESSONS.md).

## Notes

- It is recommended that you let the installation program build all the
  libraries for you and then let it test them.
- No compiler errors or warnings should appear.
- Regression output should not differ from the precomputed output files
  after the scripts remove machine-dependent lines.
