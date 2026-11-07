# Referee requests and revision scope for discussion

25 October 2026. Author-internal triage, not an agreed execution plan or a
record of completed revisions. Based on the referee comment boxes preserved
in `papers/II/referee-response-v1/response-draft.tex`, distinguished from our
proposed replies and author-side checks.

The referee describes most points as suggestions intended to strengthen the
paper. Every point still needs a response. The wording does not require us
to adopt every suggested method or every additional test in our draft reply.

## What the referee asks

| ID | Explicit request or concern | Status of the suggested method | Proportionate response to discuss |
|---|---|---|---|
| R1 | State the physical scales associated with the ratios. Explain that they measure mass concentration rather than classical efficiencies | “Mass transfer efficiency” is offered as terminology to consider, not a required replacement | Define the scales and instantaneous ratios consistently; choose terminology that does not imply a measured transfer rate |
| R2a | Motivate the study through unusually extensive coverage at high resolution; acknowledge that data combination is established | Narrative revision | Reframe the Introduction around the mosaic and science questions |
| R2b | Explain whether the scientific results would change without TP, given its small apparent effect on compact cores | A counterfactual scientific question, not an explicit instruction to rerun the whole paper | Identify which results depend on TP and quantify relevant differences using existing paired products where suitable |
| R3 | Give quantitative definitions, boundaries and uncertainties for the young/intermediate/evolved regions; justify their interpretation as an evolutionary sequence | No particular boundary-perturbation algorithm is prescribed | Show the actual boundaries and activity diagnostics, distinguish measurement apertures, and address uncertainty in the classification/selection |
| R4a | Quantify consistency in the overlapping uv range rather than just showing coverage | Cited studies/tools are examples | Provide a quantitative, common-response consistency diagnostic and explain what it measures |
| R4b | Compare almica with conventional 12m+7m imaging on a common PB footprint; clarify the apparent suggestion that CASA cannot do this; consider showing the combined PSF | Comparison and PSF are framed as suggestions; they still warrant a substantive response | Inspect existing CASA products first; compare suitable products and show the dirty PSF if available; revise claims about general method superiority |
| R4c | Quantify “marginal” TP effects in an appendix; compare J-comb with a conventional method such as feathering | Feathering is suggested; **additional TP joint deconvolution is explicitly optional** | Give paired-source statistics and a same-input feather comparison if feasible; reserve more detailed tests for results that need them |
| R5 | Explain PB treatment of 12m7m and 12m7mTP, especially how a supposedly uncorrected combined image treats unattenuated TP. Show PB=0.2 on Figs. 2 and 3 | Direct explanation and figure requests; no prescribed rerun | Establish the actual image-response and photometry conventions, correct the account, and add the selection contours. Change processing only if needed |
| R6 | Assess whether different 7m/high-resolution sensitivities and footprints bias the mass ratios or account for fewer cores near map edges | Injection/recovery and parent matching are **not explicitly prescribed** | Establish matched selection and size-dependent sensitivity; assess the effect on the reported ratios. Decide how much recovery testing is needed to support the retained claim |
| R7 | Justify neglecting compact free–free in the interferometric data; assess possible bias in core masses/ratios | In-band spectral index **or** matched-resolution radio comparison are suggested examples, not a demand for both | Use a defensible constraint or bounded contamination scenario; state unresolved limitations and their effect on the conclusions |
| R8a | Explain why sub-beam fitted sizes are not direct size measurements; give the deconvolution method and treatment of unresolved sources | Direct clarification request | Describe the implemented method accurately and how size limits enter subsequent analysis; correct implementation if the audit requires it |
| R8b | Define opacity, gas-to-dust ratio, brightness temperature and optical-depth correction in the mass calculation | Direct clarification request | Make conventions and quantities explicit and consistent |
| R8c | Explain why retain a maximum optical-depth correction of 1.084 amid much larger mass uncertainties; give brightness temperatures and consider intrinsically thick cores | Direct explanation and calculation request | Report beam-averaged brightness temperatures, justify the correction, and discuss filling-factor limitations |
| R8d | Address the major systematic in applying 18″ environmental temperatures to subarcsecond cores; justify or revise the 10–20% uncertainty | No particular scenario grid or new temperature observation is demanded | Demonstrate mass/ratio sensitivity to plausible temperature assumptions and narrow conclusions where necessary |
| R9a | Explain the distinction between the two intercore surface-density estimators directly in the text and caption | No explicit covariance/refitting request | Define exactly how each is measured and illustrate their spatial scales |
| R9b | Assess whether multiple velocity components affect the CO moment-based quantities | No specific decomposition method is demanded | Examine representative spectra/velocity selection and show or discuss the impact at the level needed for the kinematic claims |
| R10 | Add units and improve the lower panel of Fig. C1. Make the SIMBA comparison quantitative | A pixel scatter plot **or another quantitative comparison** is requested; Python is a suggestion | Produce a clear bandpass panel and a quantitative comparison with compatible frequency, units, beam and footprint |
| Presentation | Supply marked changes and a latexdiff PDF | Explicit referee request; editor deliverables are tracked separately | Prepare these after the clean manuscript is finalized, following the frozen-variant workflow |

## Additional work introduced by our draft

The following may be scientifically useful, or become necessary to defend a
particular retained claim, but must not be attributed to an explicit referee
instruction. They have not all been agreed as revision commitments.

- Extensive injection/recovery grids in flux, size, background, crowding and
  mosaic position, with completeness corrections for summed mass.
- A new parent-associated core/fragment analysis. Defining the current
  regional ratio remains necessary; switching estimands would be a separate
  scientific choice.
- A specific suite of regional-boundary perturbations. The referee asks for
  quantitative definitions and uncertainty, without prescribing this test.
- Full covariance, spatial-independence and likelihood reanalysis of the
  environmental correlations, plus a universal-factor interpretation audit.
- A prescribed set of source-level flags, contamination exclusions,
  temperature grids, or a full component-decomposed velocity product.
- All of the scale/offset fits, aperture statistics and residual diagnostics
  proposed for SIMBA, rather than a focused quantitative comparison.
- Complete historical recovery of every processing script. It helps verify
  R5, but the referee asks for a correct, supported account of the processing
  and response; the exact original script is not itself a requested deliverable.

Numerical consistency and avoiding a double PB/opacity correction remain
author responsibilities. Calling these checks additional does not justify
retaining an established error. Scope should follow the claims we retain
and the actual evidence, rather than treating every draft “we will” as a
fixed obligation.

## Provenance availability and current R5 evidence

The user reports that Sihan performed J-comb and the original script is not
available; getsf preparation was done on a remote server that is presently
inaccessible. Do not repeatedly request those scripts as though they were
locally available. This limits historical verification of R5 but does not
prevent discussion or independent work on the other comments.

The local audit establishes that the current getsf inputs match the
unmultiplied restored 12m7m image and the dust-only combined image, and that
the maintained catalogue retains the raw getsf flux columns. It also
distinguishes the selection gain from the gain used in the delivered
`pbant` image. It does not establish the residual response of the original
MIRIAD mosaic or certify the full J-comb chain. A response must distinguish
these observations from reconstruction or assumptions about missing steps.

Subsequent author confirmation establishes the reported processing order.
PB-corrected 12m7m and TP are combined first, both extraction images are
attenuated with the same combined 12m7m PB, and extracted source fluxes are
divided by that gain. The R5 reply, clean-v2 text and Fig. 2/3 boundaries
have now been updated accordingly. See the later entry in `progress.md`.
This confirmation does not constitute a newly recovered script or a measured
noise-uniformity test.

## Proposed discussion order

1. Agree the scope of the central claims, especially what the regional mass
   ratios are intended to establish. This determines how much R6–R8 testing
   is necessary.
2. For R5, distinguish a verified description correction from a processing
   correction. Discuss whether available products can establish the response
   or whether a controlled reconstruction is needed. Do not assume a full
   imaging rerun is required.
3. Choose focused evidence for R4, R6, R7 and R8d; consider available data and
   explicitly bounded conclusions. Inspect existing comparison products first.
4. Proceed with independent definitions, text and presentation items after
   agreeing their wording. Revise the draft's future-tense promises to match
   the chosen work, and finalize the response only after that work is done.

No new imaging, extraction, mass recalculation or manuscript edit is
authorized by this discussion note itself.
