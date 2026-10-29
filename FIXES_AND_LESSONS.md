# Engineering Lessons and Cleanup Record

This file has two purposes:

- **Reusable lessons** documents practices to apply when maintaining or
  modernizing the remaining packages.
- **Completed cleanup** records repository-wide and package-specific work that
  has already been performed.

## Reusable lessons

### Build and linking

- Use the GNU toolchain consistently in active package Makefiles:
  - `CC = g++ $(OPTFLAGS) $(CCFLAGS)`
  - `LD = g++ -shared`
  - `LDFINAL = g++ -Wl,--no-as-needed`
  - `CCFLAGS = -fPIC`
  - `SODIRKEY = -Wl,-rpath,`
  - `SONAMEKEY = -Wl,-soname,`
- Keep that GNU spelling in Makefiles on every platform. On macOS,
  `gnu-make-wrapper.sh` routes C++ compiles and links through `cxx-wrapper.sh`,
  which calls Apple Clang and translates the GNU-only options (`-shared`,
  `-soname`, `--start-group`, `--no-as-needed`). Fix new platform differences
  in that wrapper, not by adding per-package Makefile branches. Do not add
  GCC-only flags such as `-ffloat-store` unconditionally.
- Put `-D` string macros that contain spaces in octal (`\040`), not hex
  (`\x20`). Hex escapes consume every following hex digit, so `"\x20Darwin"`
  parses as `\x20D` and Clang rejects it; octal escapes stop after three
  digits.
- Keep inline comments off linker-flag definitions such as `SONAMEKEY`. A
  trailing space in `-Wl,-soname,` becomes an empty SONAME argument under GNU
  `ld`.
- Use `$(CURDIR)` for package-local library paths and output names. `$(PWD)` can
  still refer to the parent shell directory under `make -C`, leaving clean
  rules and shared-library paths stale.
- Prefer pattern rules when prerequisites matter. If one-off explicit compile
  rules remain, verify that they use the same toolchain variables as the
  generic rules; stale variables such as `$(CC2)` may break only clean builds.
- When a shared build wrapper accumulates package-specific branches, collapse
  them aggressively once multiple paths now share the same invocation. Keep
  only the explicit exceptions that still need different compiler variables.
- Remove package-local wrappers with hardcoded paths when shared, self-locating
  wrappers provide the same behavior.
- When a top-level config step creates a symlink farm such as `lib/`, make the
  setup idempotent. Detect stale or missing links after a checkout move or
  rename, prompt before rebuilding only when the tree is actually stale, and
  refresh the individual links on reruns instead of leaving old targets in
  place.
- When a GNU link fails after a toolchain migration, inspect library order and
  transitive shared-library dependencies before changing source code. Move a
  provider library later in the link line when appropriate.
- If local libraries have cyclic dependencies, wrap the dependency cluster in
  `-Wl,--start-group ... -Wl,--end-group` instead of repeatedly reordering
  individual `-l` entries.
- After changing a shared header or compatibility shim, rebuild every dependent
  library. Leaf executables can otherwise retain stale behavior from old shared
  objects.
- In packages that link tests against local shared libraries, header-only API
  changes may require an explicit `make lib` even when the test executable
  itself recompiles cleanly. Otherwise link failures or stale runtime behavior
  can come from an out-of-date package library rather than the edited source.
- If a regression starts crashing after build-system cleanup, perform a clean
  library rebuild before debugging source behavior. Package drivers do not
  always rebuild their local shared libraries.
- When driving the top-level `perl config` workflow, it is safe to accept
  defaults through the `regression-only version` prompt. Let the package build
  phase finish even if compiler errors appear, so you capture the full list of
  failing packages from one run.
- In the same `perl config` workflow, answer yes to
  `Would you like to test UMpack components now?` only if all packages built
  cleanly. If some builds failed, answer no and skip tests unless test
  execution is specifically needed to debug an already-known test failure.
- Make `test` or `regression` targets depend on `lib` where practical. If that
  is not possible, regression scripts should clearly report that `make lib`
  may be required after `make clean`.
- Listing `lib` beside the test executables (`test: lib $(TEST_EXES)`) does
  not order them under `make -jN`: tests can link while the package library is
  still being rebuilt. List the executables in `TEST_EXES` and add
  `$(TEST_EXES): lib`. With the old rule, a clean `makeOpt-gnu -j 14 test`
  failed in 15 of 19 packages on macOS.
- Remove inactive compiler, linker, suffix-rule, and test variants together.
  Leaving static, Purify, Solaris, or pseudo-build fragments behind creates
  configuration paths that drift.

### Headers and compatibility

- Include the system headers that declare legacy APIs such as `strcasecmp()`;
  do not rely on incidental transitive includes.
- Removing a monolithic compatibility header often exposes hidden transitive
  dependencies such as `RandomRawUnsigned`, `CHAR_MAX`, or `sqrt()`. Expect to
  add explicit includes for RNG, limits, and math facilities as part of the
  same cleanup rather than treating those follow-up fixes as unrelated bugs.
- When removing an umbrella header from a public interface, make exposed
  dependencies explicit there too: standard containers, inline assertion
  macros, and forward declarations for stream types used in method signatures.
- When code prints standard containers through project stream helpers, include
  `ABKCommon/abkio.h` explicitly. Dropping the umbrella header without replacing
  that dependency can break otherwise safe header trims in tests and library
  implementations alike.
- Do not reintroduce the removed `uofm::*` allocator-backed containers in new
  code or examples. Prefer standard library containers and streams, plus focused
  project helpers where legacy behavior still requires them.
- Be careful with legacy headers that do not protect their declarations against
  multiple inclusion through different paths. `paramproc.h`, for example, can
  collide if included explicitly while another remaining umbrella include still
  pulls in `abkcommon.h`; remove the transitive umbrella provider first or keep
  the explicit parser include at a boundary that is not duplicated.
- If a translation unit still needs a header that transitively includes
  `abkcommon.h`, you can often remove its own direct `abkcommon.h` include
  without yet adding the corresponding narrow header explicitly. That keeps the
  dependency surface shrinking while avoiding duplicate inclusion bugs from
  legacy unguarded headers such as `paramproc.h`.
- libc++ (Apple Clang) is stricter than libstdc++ about transitive includes
  and GCC-internal names. Include headers such as `<sstream>` explicitly,
  qualify `std::ios_base` and `std::setiosflags`, and do not use libstdc++
  typedefs such as `__istream_type`.
- Keep compatibility behavior in project-local shims rather than adding
  forbidden specializations to `std`.
- Delete compatibility alias headers once live includes have moved to standard
  headers. Leaving empty alias headers behind invites new code to depend on
  removed migration surfaces.
- When replacing legacy hash-table aliases, convert call sites to
  `std::unordered_map` and `std::unordered_set` directly, but keep explicit
  project hashers for legacy content-key semantics such as C strings.
- For hash tables keyed by `char*` or `const char*`, hash and compare string
  contents rather than pointer addresses.
- For maps keyed by `std::string`, use the matching string hash. Hashing
  `value.c_str()` hashes an address and can make equivalent keys miss during
  lookup.
- Replace deprecated algorithms such as `random_shuffle` with their modern
  equivalents when doing so preserves behavior.
- When replacing a shared typedef such as `bit_vector`, audit public signatures
  as well as member fields. A caller using the typedef may stop converting to a
  callee that still accepts the old allocator-specific container type.
- Legacy parser support is split across headers: command-line parameter classes
  live in `ABKCommon/paramproc.h`, while stream parsing manipulators such as
  `eathash()`, `needword()`, `needcaseword()`, and `skiptoeol` live in
  `ABKCommon/abkio.h`. Removing `abkcommon.h` often exposes both dependencies
  separately.
- Remove dead template helpers only after checking for unqualified call sites.
  Member functions named `clone()` or `next()` are not evidence that a free
  helper is still used.

### Testing and regressions

- Do not assume `make test` executes tests. In this tree it often only builds
  test executables; inspect the package for a regression driver or required
  direct invocation.
- When adding a test executable, update both the Makefile test target and the
  regression driver so the test is actually run.
- Regression scripts must check each test's exit status and the result of
  `diff`. Scripts that print `Ok` unconditionally or finish with an unrelated
  successful command can hide failures.
- Use `/bin/sh` for regression scripts unless another shell is genuinely
  required. When converting legacy shell scripts, preserve stderr capture if
  warnings may be part of the golden output.
- Do not use `echo -n` in regression scripts. macOS `/bin/sh` prints a literal
  `-n` and a newline. Use `printf '%s' "..."`, which behaves the same under
  every POSIX `sh`.
- When a tool's output file name changes (for example, gnuplot scripts now
  written as `.gpl` instead of `.plt`), update the regression script and
  golden output in the same change.
- Remove stale `TEST_IN`, `TEST_OUT`, and `TEST_DIFF` variables, including their
  clean-rule references, after moving regression logic into a script.
- Add assertion-based unit tests for semantics and invariants that are not
  visible in golden-output regressions.
- When replacing legacy `abkfatal` or `abkassert` checks with exceptions, add
  negative-path tests that exercise mismatched sizes, empty inputs, and
  out-of-range queries explicitly. Golden-output regressions rarely cover
  release-mode validation behavior.
- After changing a public header that affects template aliases or typedefs,
  rebuild stale regression executables as well as libraries. Otherwise a test
  binary can keep an old mangled signature and report a runtime symbol lookup
  failure even though the library was rebuilt correctly.
- For randomized APIs, test invariants and deterministic edge cases rather than
  exact sampled sequences.
- For floating-point results, use tolerances or narrowly scoped normalization
  when differences are numerical noise. Do not mask algorithmic changes.
- Normalize diagnostic source line numbers in golden-output regressions when
  the line number itself is not under test. For warnings formatted as
  `Warning in file:line`, strip only the numeric suffix so real warning text
  still participates in the diff.
- Normalize environment-sensitive runtime diagnostics in golden-output
  regressions when the diagnostic itself is not under test. Memory-thrashing
  and page-fault notices can vary with host pressure and should be filtered
  alongside timing and platform banners.
- The optimized default path in `perl config` enables `memory-saving options`
  unless answered otherwise, and that temporarily removes the
  `SIGNAL_DIRECTIONS` block from `ABKCommon/config.h`. This can change
  regression output, so restore `ABKCommon/config.h`, rebuild affected
  packages, and avoid committing that temporary config change unless the user
  explicitly wants the memory-saving state preserved.
- Call the intended overload explicitly when testing legacy APIs with ambiguous
  overload sets.
- Exercise copy construction using both file-backed and in-memory objects.
  File readers may populate optional metadata that direct constructors omit.
- Refresh a golden output only after confirming that the changed result is
  intentional, deterministic, and stable. If the sole mismatch is a trailing
  blank line, update the golden file rather than the production code path.
- For parser helper modernization, preserve line-number accounting, comment
  skipping, EOL checks, and exact failure points. Broad `operator>>` manipulator
  chains should be migrated with targeted parser tests rather than by changing
  only the central helper implementation.
- If a legacy regression target has a broken wrapper path, use direct compile
  checks for representative translation units, but record that the full
  regression remains unverified until the wrapper is repaired.
- Regression drivers should invoke the executable that is actually built by the
  package Makefile and fail early with a clear message if it is missing.
- Test and utility translation units often relied on `abkcommon.h` for parser,
  timer, and banner helpers. After trimming umbrella includes, add explicit
  `ABKCommon/paramproc.h`, `ABKCommon/abktimer.h`, `ABKCommon/abkio.h`, or
  `ABKCommon/infolines.h` includes instead of restoring a broad dependency.

### Parsing, strings, and data structures

- Prefer bounded formatting for fixed buffers and direct character assignment
  for one-byte delimiters. Replace fixed parser token buffers with `std::string`
  or line-oriented parsing when accepted syntax remains unchanged.
- Legacy `vector<bool>` expressions that relied on side-effect-free indexing can
  turn into `[[nodiscard]]` warnings under newer libstdc++ releases. Treat a
  bare `_flags[idx];` expression as suspicious and replace it with the intended
  assignment or remove it if it was dead code.
- If a valid root or parent lookup fails, inspect assumptions inherited from
  legacy key semantics and self-parented root handling.
- Do not mutate a tree while recursively traversing it for annotation. Compute
  annotations first, then perform structural cleanup in a separate phase.
- A `const&` member initialized from a temporary (`_ref(Type(0))`) dangles as
  soon as the constructor returns. Bind it to a function-local static or to an
  object that outlives the instance.
- Backward iteration over a `std::set` must not decrement `begin()`. Loops of
  the form `for (--g; g != s.begin(); --g)` are wrong when `g` is already the
  first element; use `while (g != s.begin()) { --g; ... }`.
- Guard ratios against zero denominators before normalization. Demand divided
  by zero routing supply yields infinity and later NaN, which silently breaks
  comparisons such as a cutline search. Fall back to a neutral weight instead.

### Source and documentation maintenance

- Keep modernization edits local and mechanical unless existing algorithmic
  behavior is demonstrably broken.
- Remove commented-out code and stale comments about unsupported compilers,
  shells, and build variants.
- Keep the top-level package inventory synchronized with the actual directory
  structure.
- Each package README should briefly describe the package, important files,
  build procedure, and test entry point.
- Keep historical package names in archived context, but use current names in
  user-facing documentation.
- Shared wrapper logic must locate itself relative to its own script rather
  than relying on an absolute checkout path.

## Completed cleanup

### Repository-wide work

- Converted active package Makefiles from Solaris-oriented settings to GNU
  `g++` builds and removed `CC`, `-G`, and `-R` assumptions.
- Replaced obsolete suffix rules where GNU make requires pattern rules with
  explicit prerequisites.
- Added linker groups where old builds depended on Solaris library-order
  behavior.
- Added missing system headers and replaced SGI STL assumptions with local
  compatibility shims where still required.
- Corrected content hashing for C strings and project string types.
- Removed ABKCommon's legacy hash-map and hash-set compatibility aliases in
  favor of direct standard unordered containers, while retaining explicit hash
  helpers for legacy string-keyed containers.
- Replaced the GCC-era `abklimits.h` compatibility shim with `<limits>` and
  converted `InvalidUnsigned` to `std::numeric_limits<unsigned>::max()`.
- Removed unused ABKCommon template helpers (`square`, `deleteAllPointersIn`,
  free `next`, and free `clone`) after replacing the active `square()` call
  sites with direct multiplication.
- Converted the shared `bit_vector` typedef to `std::vector<bool>` and adjusted
  affected APIs that still accepted allocator-specific vector types.
- Replaced selected `random_shuffle`, `sprintf`, `strcpy`, `strcat`, and fixed
  parser-buffer uses with safer modern equivalents.
- Removed inactive code, stale compatibility comments, and obsolete build
  variants.
- Collapsed repeated `make-gnu` and `makeOpt-gnu` flag handling into a shared,
  self-locating wrapper while preserving the package-facing command names.
- Made the top-level `config` script refresh `lib/` symlinks safely after a
  checkout rename by checking for stale or unexpected entries before asking to
  rebuild the directory.
- Added native macOS builds with Apple Clang: `config` detects Darwin and
  offers parallel builds, `cxx-wrapper.sh` translates GNU linker options, and
  the Linux GCC command line is unchanged. Related portability fixes: octal
  escapes in `ABKSYSTEM`, explicit libc++ includes and `std::` qualification,
  and `printf` instead of `echo -n` in all regression scripts. All
  regressions exit 0 on macOS and Linux; ClusteredHGraph, GenHier, and FMPart
  still show known platform-dependent diffs on macOS.
- Made test executables depend on `lib` in every package Makefile so that
  parallel `make -jN test` no longer links tests against a library that is
  still being rebuilt.

### Package-specific record

#### Cross-package entries

- `Capo` and `GeomTrees`: removed C++17 `register` warnings from FastSteiner,
  dropped misleading TU-local `inline` annotations that triggered optimizer
  warnings, replaced fixed-buffer `sprintf` calls in `capoPlacer.cxx` with
  bounded `snprintf`, and restored the intended `_placed[b1Cell] = true`
  update in `swapCells()`.
- `Capo` and `MLPart`: normalized warning source line numbers in regression
  comparisons so harmless source-line drift no longer creates golden-output
  diffs while preserving the warning text itself.
- `ClusteredHGraph`, `FMPart`, and `MLPart`: restored green regressions after
  GNU-only Makefile updates and clean dependent rebuilds.
- `Combi`, `Geoms`, `Stats`, and `Constraints`: added focused assertion-based
  tests for permutation, geometry, statistics, and constraint behavior,
  including deterministic edge cases, invariants, and legacy overload paths.
- `Combi` and `Stats`: after the validation cleanup, downstream verification
  stayed green through `Partitioning/regression` and a full `Capo` rebuild plus
  regression run, confirming that the refactoring did not break immediate
  consumers.
- `PartLegality` and `Partitioners`: refreshed deterministic golden outputs
  after confirming the corrected baselines were stable.

#### Package entries

- `Capo`: congestion shifting (`-congestionShifting`) now uses uniform
  whitespace weights when either side of a cut has no routing supply, such as
  the synthetic one-row layout. Before the fix, zero supply produced NaN
  weights. The new `regression-congestion` script checks that a fully occupied
  12-cell chain stays legal in one-row and multirow modes.
- `FMPart`: fixed `BucketArrayHL::getHighest()` decrementing
  `_weightList.begin()` and inserting the zero-gain bucket into itself when no
  negative-gain buckets exist. Added `FMPartTest9` to the test build and
  regression.
- `RBPlace`: `RBPCoreRow` default and `y`-only constructors bound
  `_placement` to a temporary; they now bind it to a shared static empty
  placement. The regression now reads gnuplot output from `out.gpl`.

- `Combi`: replaced the remaining local `abkFactorial()` uses in gray
  permutation generation with a package-local helper and made the permutation
  bit-vector API explicit `std::vector<bool>`, then verified the change with
  `CombiTest7` and the package regression.
- `Combi`: replaced the core `mapping`, `permut`, `grayPermut`, and `grayPart`
  validation dependencies on `abkcommon`/`abkassert` with local exception-based
  checks, added negative-path semantic coverage to `CombiTest7`, and made the
  RNG and limit dependencies explicit where the old umbrella header had been
  providing them transitively.
- `FilledHier`: removed recursive tree mutation from the area pass. Rebuilding
  `HGraph` after the hash-shim correction also removed the stale crash.
- `FMPart`: removed redundant direct `paramproc.h` includes from test drivers
  that already receive the parser declarations through remaining umbrella
  includes, and made `std::cout`/`std::endl` usage explicit in auxiliary tests.
- `GenHier`: fixed root and parent lookup, corrected hierarchy-key hashing, and
  refreshed the deterministic expected output after cluster ordering changed.
  Later header cleanup also confirmed that `genHier.h` only needs assertions,
  hash support, limits, and standard containers; parser and metadata helpers
  belong in `genHier.cxx`, and `GenHierTest0.cxx` must include
  `ABKCommon/paramproc.h` directly instead of relying on the public header to
  pull command-line parsing in transitively.
- `GeomTrees`: added focused tests for Manhattan MSTs.
- `GeomTrees`, `HGraphWDims`, `CongestionMaps`, `SmallPlacement`, and
  `SmallPlacers`: removed direct `abkcommon.h` dependencies where the code only
  needed focused parser, timer, assertion, stream, allocator, or metadata
  headers, then revalidated with package regressions and downstream rebuilds
  through `RBPlace`, `Capo`, and `ParquetFP`.
- When replacing shared container types, start with private members and helper
  signatures. That keeps the change mechanically safe, and a clean forced
  library build plus downstream rebuild will show whether any public API work is
  still justified.
- `RBPlace` and `Capo`: replaced remaining public-header `abkcommon.h`
  dependencies with focused standard and `ABKCommon` includes, then removed
  the remaining direct `abkcommon.h` usage throughout `Capo` tests and
  implementation files while preserving green `makeOpt-gnu -B test` and
  `regression` runs after each tranche. During that work, `paramproc.h` again
  proved that direct parser-header additions must stop at boundaries that do
  not already receive it through a surviving umbrella include. Finishing the
  remaining `RBPlace` source files also confirmed that public `bit_vector`
  members require `ABKCommon/abktypes.h` in the owning header, and that the
  placement/AUX readers need `ABKCommon/abkio.h` explicitly because parser
  manipulators are not declared by `paramproc.h`.
- `ParquetFP`, `ParquetDBFromRBP`, and `AnalytPl`: replaced direct
  `abkcommon.h` usage with focused verbosity, parser, timer, assertion, and type
  headers; fixed `ParquetFP`'s hidden reliance on `FPcommon.h` as a transit
  source for allocator-related types; normalized environment-sensitive memory
  diagnostics in `ParquetDBFromRBP/regression`;
  and revalidated with package regressions plus a downstream `Capo` rebuild and
  regression.
- `HGraph`: added focused tests for in-memory hypergraphs and fixed copying
  unnamed in-memory hypergraphs.
- `ParquetFP`: removed obsolete static and Purify Makefile paths and the
  hardcoded local `makeOpt` wrapper; replaced local C-style name construction;
  converted the regression script to POSIX `sh`; and restored a green
  regression after the GNU-only Makefile update and clean rebuild.
  The later clang-format pass also exposed a broader lesson: some headers in
  this codebase depend on direct base includes rather than transitive ones, so
  automated formatting or include reordering should always be followed by a
  package rebuild to catch type-visibility regressions early.
- `PartEvals`: its regression became green after rebuilding `HGraph`.
- `Ctainers`/`BitBoard`: converting an inline shared utility to `std::vector`
  requires rebuilding every library that consumes the affected object layouts,
  not just the directly edited package.
  Stale `HGraphWDims` and `PlaceEvals` shared libraries manifested as
  downstream `RBPlace`/`Capo` segfaults even though the source fixes were
  correct. Rebuild dependency chains before diagnosing these as logic bugs.
- Regression scripts run under `/bin/sh` should not rely on brace expansion.
  Use explicit paths such as `seeds.in seeds.out`; otherwise stale `seeds.out`
  files can create noisy seed warnings or hide abnormal prior termination.
- `Partitioning`: a clean rebuild of `libPartitioning.so` cleared its regression
  crash. Added tests for partitioning buffers and restored the orphaned
  `PartitioningTest12` to the package test build.
- `Placement`: normalized the printed bounding box in `PlacementTest1.cxx`,
  aligned its expected output, and added tests for geometry and transforms.
- `SmallPart`: restored a green regression through the clean rebuild path.
- `Capo`, `FMPart`, `MLPart`, `Partitioners`, `PlaceEvals`, and `SmallPart`:
  removing umbrella includes exposed many hidden test-driver dependencies on
  `ABKCommon/paramproc.h`, `ABKCommon/abkio.h`, `ABKCommon/infolines.h`, and
  `ABKCommon/abktypes.h`. The practical lesson is that after an include
  cleanup, package libraries are not enough: rebuild the package tests and then
  rerun the regressions, because the first failures often come from smoke
  tests, not from the libraries themselves.
- `Stats`: replaced several `abkcommon.h` validation dependencies in
  `trivStats`, `linRegre`, `rancor`, `freqDistr`, `loadedDie`, and
  `keyCounter` with standard exceptions and limits, fixed the zero-threshold
  `KeyCounter` constructor bug, fixed `TrivialStats` reading `data[0]` before
  handling empty floating-point inputs, and extended `StatsTest21` with
  negative-path validation checks.
- `Stats`: continued the same cleanup into `expMins` and `multiRegre`,
  replacing `UINT_MAX` and fatal precondition checks with standard limits and
  exception-based validation, while extending `StatsTest21` to cover expected
  minima and multiple-regression positive and negative paths.
