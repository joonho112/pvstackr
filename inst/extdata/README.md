# pvstackr extdata

These files are small synthetic data made for the examples and tests of
pvstackr.

- `pisa_tiny.csv` has the column layout of a PISA file (12 students, two
  plausible values, four replicate weights) but contains no real PISA records.
- `pisa_tiny_manifest.dcf` records how the data and the example fit were made,
  with checksums.
- `examples/pisa_tiny_stack_direct.rds` is an example `stack_direct` fit of
  these data. It was made without a sampler, with fitting functions that return
  draws built around the target, so it shows what a fit looks like, not how
  well the method works.

The files are kept under version control and are rebuilt from the package's
development sources.

They serve the examples and tests only. They are not suitable for real
inference, coverage claims or performance benchmarking.
