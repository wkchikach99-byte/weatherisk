# v24 submission to-do (audit of 2026-11-12)

This audit covers four areas: Thorsten's comments C1 to C18, the current Guide for Authors and 2026 template, the numbers and figures, and the repository and release.
- Nothing in the paper or the code was changed for it.
- Line numbers refer to `softwarex_paper_v24.tex` at commit 53b6db3.
- "Verified" means checked against the code, the template or the output files.

## Order of work

1. Paper fixes (section B), after Widad's go-ahead.
2. Decisions on what becomes public (A1).
3. Repo preparation (section C).
4. Merge, push, release, DOI (A2, A3).
5. Paste the DOI, rebuild the PDF, upload (A6).

Zenodo archives every tracked file of the release, and an archive cannot be withdrawn. A1 therefore comes before the tag.

## A. Only Widad can do these

| # | Item | Why |
|---|---|---|
| A1 | Decide what is public before tagging v0.2.0. The candidates are `r_code/` (J. Contzen's code; it has no licence from him), `r_code/main.pdf` (his pre-publication manuscript), `judgment_calls_for_review.md`, the paper drafts and cover letter in `docs/software_x/final/`, and the personal paths and emails in `hpc/*.slurm` and `hpc/setup_albedo.sh`. | The MIT LICENSE names only the three authors. v24 no longer says Contzen published the code, so nothing in the paper backs hosting it. The link at l.511 points to the pinned commit 4cf18b2 and keeps working if `r_code/` is removed from the tip later. |
| A2 | Merge `remove-r-runtime` into `main` (a fast-forward is possible) and push. Do **not** squash. | `main` still contains R and Rust, a broken `pyproject.toml` and no pure-Python HPC script. The paper cites commit 7561a6e (l.688), which a squash would orphan. CI has never run on GitHub, so check that it is green after the push. It has not been verified on Ubuntu with older NumPy and SciPy. |
| A3 | Check that the Zenodo–GitHub integration is on, publish GitHub Release v0.2.0 and note the new DOI. | The DOI in the paper (zenodo.20108563, l.121) is the **v0.1.1** archive, which is the R-dependent code. Verified on zenodo.org. |
| A4 | C13 wording. Choose option A (Thorsten's master text, verbatim) or option B (the same text, but "simulation study" instead of "empirical illustration"). Ideally ask Thorsten. | v24 dropped his sentence "served its intended purpose very well". Whether Contzen ran `r_code/` on climate data for his Figure 9 cannot be checked from the repo. Both options are spelled out in section K of `judgment_calls_for_review.md` and in the audit report. C4 is fine as it is; restoring his full clause would overstate the validation. |
| A5 | Confirm the facts that Claude would write in: the postal addresses (B3), whether Dickhaus has funding to name (B4), and the AI declaration ("the authors" or "the first author", and which tools; l.817). | These are facts only Widad knows. |
| A6 | Editorial Manager: complete the declarations tool and upload its .docx, enter the funders, and upload the PDF, the source zip, Figures 1–3 as separate files and `highlights.txt`. | The guide says the declarations tool must always be completed. |
| A7 | `LICENSE.txt` and `src/` layout: rename or copy, or ask softwarex@elsevier.com. | The guide asks for both. The repo has `LICENSE` and a flat `weatherisk/` package. Claude can prepare either change locally. |

## B. Paper fixes Claude can make on request

Every item has been verified against the template or the code. **TW** marks items that touch supervisor-marked wording. Those items keep the wording and only add to it or correct a fact.

| # | Line(s) | Problem | Fix |
|---|---|---|---|
| B1 | 102–149 | **Desk-rejection risk.** The metadata table is the old 9-row version, with a "Reproducible Capsule" row and different labels. The 2026 template has 8 rows, and its instructions say to "leave the left column as it is". | Use the template headings `Required Metadata` / `Current code version` and its 8 left-column labels verbatim. C2 becomes the GitHub link for v0.2.0. The Zenodo DOI moves to a software reference. |
| B2 | 599 | **TW.** AWI is used in the body (l.590, 595, 599, 699) but expanded only in the abstract and the Acknowledgements. Thorsten's C3 body definition was lost in the cut. | Restore his phrase "Alfred Wegener Institute (AWI)" at l.599 (+3 words). |
| B3 | 70–72 | Affiliations have no postal addresses. The guide requires the full address. | Am Handelshafen 12, 27570 Bremerhaven; Bibliothekstraße 1, 28359 Bremen (A5). |
| B4 | 833–836 | The funding statement has no sentence saying the funders had no role. | Add the standard Elsevier sentence. |
| B5 | 384–401 | Figure 1 module labels are wrong (verified): data loading is `cmip6_data`, not `io`; stage 1a is `stl`; local MLE is in `density`, not `estimation`; stage 6 is `estimation`, not `pipeline`. l.426 claims the figure "names the module of each stage". | Correct the four labels. |
| B6 | 87–90; highlight 5 | The abstract and highlight 5 give only the run without detrending, which matches the published 104. The paper itself says Contzen et al. used STL detrending, which gives 142 and 14. | Add "with STL detrending, 142 and 14" (about 8 words). Define EDC in the highlight. |
| B7 | 570, 590, 594, 696 and others | **TW.** "Figure 9" and "their Figure 3" refer to Contzen's figures, but the paper has its own Figures 1–3. The Figure 3 caption reads "The Figure 9 case study…". | At the first mention and in the caption, write "Figure 9 of Contzen et al. [4]". Keep the "case study" wording. |
| B8 | 299, 552–558 | Example 1 never states its settings. Section 2 gives the CMIP6 defaults (ε = 5, smoothing 2), but `paper_stripes` uses ε = 4 and smoothing 4 (verified, `parameters.py:84–85`). | Add one sentence on the preset settings, and write "default in `CMIP6Config`" at l.299. |
| B9 | 431–438 | `weatherisk cmip6` has no `--no-detrend` option (verified), so the headline run cannot be reproduced through the CLI. `validate` is only a sanity check. The heading names `scalable`, which nothing imports (verified). The checkpoint claim covers only the CMIP6 pipeline. | Reword four short phrases. |
| B10 | 760–762 | "The full case study ran in 26 to 36 hours" covers only stages 3–6, resumed after stage 1 (l.634–636). | "Stages 3–6 of the case study took 26 to 36 hours on 64 cores." |
| B11 | 288, 599–601 | ν is used before it is defined (it is first defined in Table 3). The text says the SLURM script is "(Listing 1)", but the listing shows the Python calls it runs. | Define ν at first use, and write "which runs the commands of Listing 1". |
| B12 | `cover_letter.tex` 75–79 | The cover letter mentions "metadata tables" (plural) and "a code appendix", but v24 has one table and no appendix. | Replace with the v24 contents. |
| B13 | 18–39, 97, 115, 121, 208–215, 817 | The source zip goes to the editors. It still has the "THORSTEN-18 TRACKING" header, the `\rev*` macros, the RELEASE BLOCKER comments, and `\revt` wrapped around new text in C13 that is not Thorsten's. | Make a clean copy for upload. Keep the tracked v24 for the supervisors. |
| B14 | references | Eight bibitems have no DOI, although DOIs exist (checked against Crossref). The software references (pyextremes, SpatialExtremes) lack versions. | Add them. References are not counted in the word limit. |

After B1 to B14 the count rises from 3675 to about 3750 words, under the 4000 limit.

## C. Repo preparation Claude can do locally (Widad pushes)

| # | Item |
|---|---|
| C1 | Bump the version to 0.2.0 in `pyproject.toml:7`, `weatherisk/__init__.py:3` and `CITATION.cff`. Fix the CITATION date (2026-05-07; Zenodo says 2026-05-10). Rename the CHANGELOG `[Unreleased]` section to `[0.2.0]` and remove its TODO comment. |
| C2 | README: the preset table says `stripes` has resolution 10; the code says 51 (verified). The HPC workflow hard-codes `/albedo/work/user/wikchi001/...` and `--account=paleodyn.paleodyn` (verified). The Example 1 command used for the paper is missing. |
| C3 | After A3, paste the new DOI into the paper (l.121), `pyproject.toml`, the README badge and `CITATION.cff`, then rebuild the PDF. |
| C4 | Before the merge, rerun ruff and the tests with R truly hidden, and count the tests in a clean export. |

## Checked and fine

- **Thorsten's comments.** C1, C5–C12, C14–C17 are kept verbatim. C9 and C18 keep their intent. No deleted claim has come back ("asymptotic", "regardless", "centuries", "stationary signal").
- **Length and structure.** 3675 words, 3 figures, 5 template sections. The end statements are in the guide's order. All references are cited and all `\ref`s resolve. The competing-interests text is the standard wording. danek2020 is a correct [dataset] entry.
- **Numbers.** Checked against `pipeline_results.npz` and the frozen R fixtures:
  - cluster counts 104/23 and 142/14;
  - Table 2 deviations;
  - run times 36:02:48 / 25:58:03;
  - memory figures;
  - Example 1 (k = 7/5, MAE 0.37/0.86);
  - test tolerances.
- **Tests.** With R truly removed from PATH: 235 tracked tests pass in 137 s.
  - The command used before did not hide R. On this Mac `Rscript` is in `/usr/local/bin`.
  - The paper's claim at l.504 holds.
- **Figures.** Figure 3's panels match its caption. Both PNGs meet the guide's resolution minimum. The Figure 2 panel titles and ticks are very small at print width. Enlarging them is optional.
