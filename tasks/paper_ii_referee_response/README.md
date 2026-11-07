# Paper II referee response

This coordination task owns the revision of LANCET Paper II, AA62184-26,
and the point-by-point response to the first referee report. It connects
referee diagnostics to the existing LANCET analysis tasks; it does not create
a second scientific pipeline.

The task configuration and subtask status are in [task.toml](task.toml).
Record completed work and evidence in [progress.md](progress.md).

## Manuscript workflow

The only working manuscript source is:

`papers/II/submission_a&a_v2/lancet-II-v2.tex`

Edit this file directly. Wrap new or replaced scientific text in
`\revise{...}` so that it appears in red. The macro changes the color of the
actual text, without adding a bracketed author note or making it bold.
Use separate revision spans for separate paragraphs; preserve labels,
citations, equations and LaTeX syntax. Describe deletions in the change log;
the final latexdiff will show removed text.

The clean source has `\markrevisionstrue` during preparation. Once the
manuscript is finalized, `\markrevisionsfalse` produces ordinary submission
text from the same source. This is a finalization step, not an instruction
to remove the current red markings now.

These two reserved files are frozen during preparation:

- `papers/II/submission_a&a_v2/lancet-II-v2-bold.tex`
- `papers/II/submission_a&a_v2/lancet-II-v2-diff.tex`

Do not edit, regenerate, compile as the working target, or synchronize these
files while revising the clean source. Generate the requested bold and
latexdiff versions only after the user says the clean manuscript is final.
Use the submitted v1 manuscript as the referee comparison baseline.
Do not automatically run latexdiff after an edit.

The existing reply draft is
`papers/II/referee-response-v1/response-draft.tex`. Its folder name identifies
the first response draft; it is not the active manuscript version. It was
prepared from the submitted v1 and must be checked against the current v2
before its numbers or proposed changes are reused. The formatting reference
remains `papers/II/response-template.tex`.

## Building the response in VS Code

When `papers/II/` is the open VS Code folder, its local
`.vscode/settings.json` selects the installed Tectonic compiler for LaTeX
Workshop. This avoids the missing `latexmk` executable. The recipe uses an
absolute compiler path, writes the PDF beside the selected root source, and
produces SyncTeX information for source/PDF navigation.

Use **Build LaTeX project** or **Build with recipe: Tectonic**. The workspace
setting applies to the `papers/II/` folder; opening a different workspace
does not automatically inherit this nested settings file.

Open `papers/II/` using **File > Open Folder** before building. A file-only
window whose Explorer says **No Folder Opened** does not load this folder's
workspace settings and can fall back to `latexmk`. The correct window shows
`II` as the Explorer root. If inspecting settings manually, use the workspace
settings rather than the read-only `defaultSettings.jsonc` reference.

The manuscript folder also has its own LaTeX Workshop settings for a
multi-root workspace that opens `submission_a&a_v2` and `referee-response-v1`
separately. Its default recipe runs PDFLaTeX through the local TeX Live
`latexmk` supplied with the Realtime TeX extension. The manuscript recipe
uses an absolute executable path and places output beside the source.

## Two kinds of scientific work

**Standalone diagnostics** answer a specific referee question without changing
the adopted analysis. Examples include a representative joint-imaging
comparison, feathering comparison, boundary perturbation, temperature
scenario, or a small source-injection pilot. Write their products under
`results/paper_ii_referee_response/diagnostics/<subtask>/<run-id>/`. Record
inputs, settings, commands, masks, units and outcome. A diagnostic is evidence
only for the conditions it actually tests.

**Maintained updates** belong in the existing owning task's `task.toml`, its
runner in `scripts/`, and reusable functions in `src/lancet/`. Examples
include a corrected PB convention, catalogue selection, deconvolution/limit
treatment, contamination handling, uncertainty model, or publication figure.
Verify the affected behavior, rerun the affected stages and downstream
consumers, and record the actual resulting products. Do not leave an adopted
correction only in a temporary script, a figure copied into the paper, or a
second response-specific implementation.

Subtasks marked `diagnostic_then_update` start with an isolated diagnostic;
if the result changes the adopted analysis, update the listed owner. A null
result can instead support a bounded explanation or limitation in the reply.
No full imaging run or mass-production rerun is implied by initializing this
coordination task.

## Referee work packages

| ID | Report item | Work | Owning tasks / integration target |
|---|---|---|---|
| R1 | 1 | Physical scales, snapshot mass-fraction terminology, consistent numerical accounting | `core_formation`; clean manuscript |
| R2a | 2a | Scientific motivation and mosaic coverage/resolution | Clean manuscript |
| R2b | 2b | Separate ACA and TP contributions; state result provenance | `source_catalogue`, `synthesize_continuum`, `column_density` |
| R3 | 3 | Region definitions, L/M apertures, boundary sensitivity | `source_catalogue`, `column_density`, `core_formation` |
| R4a | 4a | Quantitative overlap consistency under a matched sky response | `combination_diagnostic` |
| R4b | 4b | Conventional joint imaging versus almica; dirty PSF | `combination_diagnostic`; isolated continuum-imaging diagnostic |
| R4c | 4c | J-comb versus feathering; source/background differences | `combination_diagnostic`, `source_catalogue`, `synthesize_continuum` |
| R5 | 5 | PB construction, noise conventions, PB=0.2 figure overlays | `prepare_continuum_images`, `plot_allscale_continuum`, `column_density`, `source_catalogue` |
| R6 | 6 | Common footprint, parent associations, source recovery, mass completeness | `source_catalogue`, `core_formation`; isolated injection/recovery pilot |
| R7 | 7 | Compact free–free constraints and affected mass sums | `construct_freefree_model`, `synthesize_continuum`, `physical_properties`, `core_formation` |
| R8 | 8a–d | Sizes/limits, opacity, brightness temperature, temperature systematics | `source_catalogue`, `physical_properties`, `background_effects`, downstream consumers |
| R9 | 9a–b | Intercore estimators, spatial/error covariance, velocity components | `environment_measurements`, `correlation`, `c18o_moments` |
| R10 | 10 | SIMBA bandpass presentation and quantitative comparison | `synthesize_continuum` |
| E1 | Introductory/editorial requests | Final reply, clean/marked/diff PDFs, references and language | Finalization only; reserved manuscript variants stay frozen |

`tasks/c18o_phangs_combination/` is an existing spectral-line pilot. It is not
a substitute for the continuum imaging comparison in R4b and must not be
rerun or repurposed merely because both involve data combination.

## R4 combination-diagnostic ownership

At the author's request, R4 uv coverage, visibility diagnostics and subsequent
combination comparisons extend the maintained `combination_diagnostic` task.
Use the local `combination_diagnostic` symlink and the corresponding result
symlink to reach its configuration and evidence. Keep the scientific code and
products there; this coordination task records status and manuscript/reply
links. Figure 1 placement is a later decision.

## Per-subtask completion record

For each item, record:

1. The precise question and whether the work is a diagnostic, maintained
   update, or manuscript clarification.
2. The input versions, configuration, common footprint and physical units.
3. The command and the actual evidence, including unsuccessful or null tests.
4. The owning task/script changes and affected downstream products, if any.
5. The clean-manuscript location and red revision span.
6. The matching reply text, with actual outcomes and verified locations.

Use `pending`, `in_progress`, `needs_evidence`, `ready_for_review`, and
`complete` as status values. Completed code does not automatically mean that
the scientific referee item is complete. Keep proposed tests in future tense
until performed; distinguish catalogue consistency from injection recovery,
environmental temperature from core temperature, and beam-averaged radio
constraints from core-scale constraints.

## Dependencies and manuscript figures

The comment-by-comment dependency review and proposed starting sequence are
in [dependencies.md](dependencies.md). Start with the R5 PB/selection audit
and the R3 boundary definitions; distinguish work that can start immediately
from dependencies needed to finalize a scientific result.

Follow the dependencies of the existing task configurations. In particular,
changes to the source catalogue can affect physical properties, environmental
measurements, core mass fractions and correlations. Changes to the opacity,
temperature or contamination convention must be propagated consistently to
all relevant mass/surface-density products. Shared uncertainties should not
be treated as independent merely because quantities appear in different
tables.

Keep diagnostic figures in their diagnostic directory and maintained figures
under the owner's `results/` directory. Copy a verified adopted figure into
the v2 `figures/` directory deliberately and log its source. Do not bulk-copy
all results or overwrite an existing paper figure with an exploratory one.

Read-only scientific inputs remain read-only. This task follows the repository
output policy: new referee diagnostics belong under `results/`, not `data/`.

## Finalization checklist

- [ ] Scientific work and the clean manuscript are agreed final.
- [ ] All response claims and numerical values match the adopted task outputs.
- [ ] All referee items have recorded evidence or an explicit limitation.
- [ ] Final figure/table/section references are verified.
- [ ] Save the red working source/marked PDF and produce the clean submission
  PDF from the same manuscript source.
- [ ] Generate and check the reserved bold version from the finalized source.
- [ ] Generate and check latexdiff against the submitted v1 source.
- [ ] Compile and visually check the clean, marked and diff deliverables.
- [ ] Finalize the response and cover letter; check current submission deadline.

Initialization of this task does not satisfy or trigger these finalization
steps.
