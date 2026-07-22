# HBN Feasibility Check For WP07

Date: 2026-07-13.

Question: Can Healthy Brain Network (HBN) support a p-factor psychopathology example where apparent extra latent dimensions may reflect nonlinear residual dependence?

## Finding

HBN is a strong substantive candidate, but the immediately open HBN-EEG/OpenNeuro releases are not sufficient for the desired analysis.

## What Is Open Without DUA

The HBN-EEG OpenNeuro/GitHub releases expose `participants.tsv` files with:

- participant ID,
- release number,
- sex,
- age,
- handedness,
- commercial-use flag,
- full-phenotype flag,
- derived `p_factor`,
- derived `attention`,
- derived `internalizing`,
- derived `externalizing`,
- EEG task availability/quality fields.

Across the currently checked releases (`ds005505`, `ds005506`, `ds005507`, `ds005508`, `ds005509`, `ds005510`, `ds005511`, `ds005512`, `ds005514`, `ds005515`, `ds005516`), `participants.tsv` covers 3155 participants. The file trees checked for `ds005507` and `ds005516` contain only `participants.tsv` and `participants.json` among phenotype-like files; no raw CBCL item or syndrome-scale tables are present there.

These derived scores are useful for downstream analyses, but they are not enough to test whether a one-factor p model, a two-factor internalizing/externalizing model, or a bifactor model is needed. The relevant factor model has already been applied before release.

## What Requires HBN DUA / LORIS Access

The HBN phenotypic portal states that full phenotypic data require completing the HBN Data Usage Agreement and then accessing the LORIS database. The portal also explains that users can query and download selected instruments/fields as CSV after access is approved.

The Release 11 data dictionary zip is publicly downloadable and confirms that raw CBCL instruments are documented:

- `CBCL.xlsx`: Child Behavior Checklist age 6-18, with raw item variables `CBCL_01`, `CBCL_02`, etc., scored 0-2.
- `CBCL_Pre.xlsx`: preschool CBCL, with raw item variables `CBCL_Pre_01`, `CBCL_Pre_02`, etc., scored 0-2.

Therefore, HBN should support the desired p-factor residual-diagnostic analysis if we obtain LORIS/DUA access to raw CBCL data.

## Suitability For The Manuscript

Scientific suitability: high.

- HBN is transdiagnostic pediatric mental-health data.
- The relevant symptom instrument is CBCL, which is central in youth p-factor/internalizing/externalizing work.
- HBN is more open/reproducible than ABCD if DUA access can be obtained and documented.

Immediate reproducibility: limited.

- Without DUA/LORIS access, we cannot fit the one-factor, two-factor, or bifactor models ourselves from raw indicators.
- The derived OpenNeuro p/internalizing/externalizing scores are circular for the manuscript question because they presuppose the bifactor scoring model.

## Recommended Next Step

If the user is willing to use DUA-gated data, request/obtain HBN LORIS access and download `CBCL` for age 6-18 participants. Then run:

1. one-factor p model,
2. two-factor internalizing/externalizing model,
3. bifactor p + internalizing/externalizing model,
4. residual diagnostics on the one-factor model,
5. nonlinear follow-up for localized residual dependencies,
6. comparison of whether nonlinear terms reduce the apparent need for multiple latent dimensions.

If DUA-gated data are not acceptable for the paper, HBN should not replace the fully reproducible Wage example.
