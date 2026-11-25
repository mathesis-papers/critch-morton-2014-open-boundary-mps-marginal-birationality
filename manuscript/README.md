# Proof of Critch–Morton Conjecture 8

Posing locator: Andrew Critch and Jason Morton, *Algebraic Geometry of Matrix
Product States*, SIGMA 10 (2014), 095, Section 2, p. 6, Conjecture 8.

For the affine closures of translation-invariant open-boundary binary bond-two
matrix product states, every consecutive three-site amplitude-sum map is
birational at every length N ≥ 4. The note gives an explicit nonempty rational
chart, both inverse identities and prime vanishing ideals. Written corollaries
give rationality in dimension seven, arbitrary fixed nonzero site-dependent
coherent contractions onto every consecutive block of length at least three,
minimality of three retained amplitude indices, and exact membership equations
on the reconstruction chart. The N−2 pulled-back quartics are necessary but
insufficient for global membership. These are amplitude statements; the
reduced-density Conjectures 9 and 10 are not covered.

## Files and reproduction

- `note.tex`: complete statements, proofs, coverage table and six references.
- `metadata.json`: scientific description, creators, license and declaration link.
- `verify.py` and `checks.json`: quartic coefficient identity, integer witness and
  exact rational reconstruction examples.
- `verify_extensions.py` and `extension-checks.json`: ordered coherent contractions,
  longer blocks, normalized seven-parameter coordinates and chart equations.
- `source-manifest.json`: exact sizes and SHA256 hashes of the preceding files.

Build with Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, with a TeX distribution containing the packages used by `note.tex`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded in the main source. Both checkers use Python 3's
standard library, exact integers and `fractions.Fraction`:

```sh
python3 verify.py
python3 verify_extensions.py
```

`checks.json` gives 120 source-inverse cases, 120 target-inverse cases and 47,040
coordinate equalities, for six parameter choices and all triple windows at
lengths 4–8. `extension-checks.json` gives 42 weighted source cases, 42 weighted
target cases, 60 longer-block cases, 42 chart cases and 42 seven-parameter cases,
with 14,160 coordinate equalities, at lengths 4–7. These are finite evaluations;
the all-length statements and corollaries have the written proofs in the note.

## Declaration coverage

The declaration `D5.S3.Quantum.TensorNetworks.OpenBoundaryMPSMarginalBirationality.result`
in the [pinned Lean source](https://github.com/the-omega-institute/trureturing/blob/b95cdbb8cccca425964251ba38cc0b0e517e478a/D5/S3/Quantum/TensorNetworks/OpenBoundaryMPSMarginalBirationality.lean)
has type `claim`: the two prime vanishing ideals and `IsBirationalVia`, with
polynomial rational coordinates, a nonempty chart and both inverse identities
for every N ≥ 4 and every consecutive triple. The four corollaries are written
proofs rather than additional declarations in that file. The quartic coefficient
comparison is exhaustive; rational evaluations have the finite scope above.

## Creators and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Creator order: Ma, then Zhang. Licensed under Creative Commons Attribution 4.0
International (CC BY 4.0, `cc-by-4.0`).
