# Reconstruction record

This is a new reconstruction and rerun of work reported earlier in the same
conversation. The earlier local repository, bundle, patch, result files and
test logs were unavailable when publication resumed. The remote repository
did not contain the reported head `9c939d4` or the research branch. We do not
claim to have recovered or independently authenticated those missing objects.

The prior conversation records protocol commits `e75dac6`, `ed41b82` and
`0434f63`, a primary failed gate, an exploratory donor-only sensitivity, a
sequencing-header audit, and numerical results. Those records inform the
reconstruction; they are not externally timestamped preregistration evidence.
The old hashes are historical identifiers, not working provenance links.

The specification below is committed before this rerun. The earlier results
are already known, so this run is a reproduction with known targets, not a
new blinded test. Reconstructed code, newly retrieved source bytes, new test
logs and newly calculated results supply the evidence for this contribution.
No previous test count is inherited as a current test result. A mismatch must
be disclosed and investigated without changing the design to force agreement.

## Known comparison values (not fabricated outputs)

- HCTbc primary fold gains: 0.14424247849217992, 0.04254064448634634,
  0.9093840407316872, 0.14649162119995474 nats per read.
- HCTbc equal-culture mean: 0.31066469622754206; fixed-predictor adverse
  1% outcome-mass lower mean: 0.3032144246585811.
- HCTbc donor-only `1e-5` threshold: all four gains zero, all donor exponents
  zero; donor mass fractions removed in culture order 0.00457817568539598,
  0.003465069460411585, 0.0033411302679190458, 0.0019795896037803073.
- HCTbc outcome-only `1e-5` threshold mean: 0.3114309239165016.
- SW6bc primary gain zero; strict `1e-3` floor unavailable.
- Exploratory temporal diagnostic: three of four positive home means;
  equal-home mean 1.4803276373883425; nine of 24 evaluation halves have one
  occupied receiver barcode. These are conditional selected-half scores,
  not whole-library scores or independent biological replicates.

Source workbook expected SHA256:
`5c61f93a15134b795e6655ecab18f2d30e7d52daed33168ca0a60406303ee100`.
The prior report's source and cache fingerprints are checked against newly
retrieved bytes; a retained self-declared metadata label is insufficient.

Publication of the completed branch and a draft PR was explicitly authorized
by the user. This does not authorize merging or promoting findings to a
reviewed or clinical status.
