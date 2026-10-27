# Capo

Capo is the top-down partitioning and placement controller in UMpack.
It drives recursive partitioning, block construction, and the search
logic used by the Capo placer.

## Key Files

- `capoPlacer.h`
- `capoPlacer.cxx`
- `capoMainLoop.cxx`
- `partProbForCapo.h`
- `partProbCtor.cxx`
- `DOCS/`

## Build

Optimized builds use `./makeOpt-gnu`.

Use the package wrapper from this directory:

```sh
./makeOpt-gnu
```

For a debug build, use `./make-gnu`.

## Test

Run the package regression script:

```sh
./regression
```

For placement-oriented runs, these test executables are the most useful entry
points:

- `CapoTest0`: Bookshelf-based smoke test with fallback to a bundled `.aux`
  file when `-f` is omitted.
- `CapoTest1`: Flexible Bookshelf placement run.  It reads the design from
  `-f design.aux`, runs Capo, and can save the final placement with
  `-savePl file.pl` or a gnuplot visualization with `-plot*` options.
- `CapoTest2`: DB/DEF-based placement run that exercises the database import
  path.
- `CapoTest4`: DB-based placement run that can also check consistency and
  save a DEF.
- `CapoTest6`: Placement recovery utility that normalizes a `.pl` file and
  writes `out-recovered.pl`.

## CapoTest1 One-Row Ordering

`CapoTest1` can also derive a linear node ordering from a netlist.  Passing
`-saveOrder file` selects a synthetic one-row mode:

- the netlist is read from `-f design.aux`, which may be an `HGraph` aux
  containing only `.nodes` and `.nets`;
- every node is forced to unit size;
- Capo builds a single-row placement with one slot per movable node;
- the ordering file contains one movable node name per line in final
  left-to-right row order.

Example:

```sh
./CapoTest1.exe -f TESTS/saurabh1-hg.aux -saveOrder ordering.out -savePl ordering.pl -plot ordering
```

The same run can still write a Bookshelf `.pl` file with `-savePl` and a
gnuplot script with any existing `-plot*` option.  The package regression
checks this mode by comparing `capo1.order` with `expected1.order`.

Congestion shifting (`-congestionShifting`) uses routing demand divided by
routing supply. If either side of a cut has no routing supply, Capo uses
uniform whitespace weights instead. This includes the synthetic one-row
layout, which has no routing-resource data. The regression script checks
that a fully occupied 12-cell chain remains legal in both one-row and
multirow modes with congestion shifting enabled.

## Notes

- Legacy design notes and examples live under `DOCS/`.
- `README` in this directory is an older package note preserved for context.
