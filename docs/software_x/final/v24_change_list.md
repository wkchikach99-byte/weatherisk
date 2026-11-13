# v24 change list, item by item (2026-11-13)

This list collects every point from the four audits of 2026-11-12 that concerns the paper, the highlights or the cover letter. Repository and release items are listed only where the paper depends on them; the rest are in `submission_todo_v24.md`.

- Line numbers refer to `softwarex_paper_v24.tex` at commit 53b6db3.
- Nothing has been changed yet.

**Markers**
- **[R]** could lead to a desk return or rejection.
- **[F]** is a factual error a reviewer can check.
- **[C]** is clarity or consistency.
- **[T]** is a small point.
- **TW** means the passage carries supervisor wording (`\revt` is Thorsten; `\rev`, `\revb` and `\revg` are earlier supervisor edits). The fix keeps their words and only adds to them or corrects a fact.
- **W** means the item needs a fact or decision from Widad first.

## 1. Front matter

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 1 | [R] W | 70–72 | `Bremerhaven, Germany` and `Bremen, Germany` only. The guide asks for the full postal address. | `Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research, Am Handelshafen 12, 27570 Bremerhaven, Germany` and `University of Bremen, Bibliothekstraße 1, 28359 Bremen, Germany`. Add departments if you want. |
| 2 | [F] | 84–90 | The abstract gives only the run without detrending, which matches the published 104. The paper itself (l.662–667) says Contzen et al. used STL detrending, which gives 142 and 14. | Append to the last abstract sentence: "; with the detrending reported in that study, the counts are 142 and 14." (+12 words) |
| 3 | [T] | 75–91 | The abstract has 121 words. The template asks for about 100; the guide allows up to 250. | No change needed. If you want it shorter, item 2 can be offset by shortening the first sentence. |
| 4 | [T] TW | 87 | "Alfred Wegener Institute (AWI)-ESM-1-1-LR" reads awkwardly. | No change, because this is Thorsten's C3 wording. |

## 2. Metadata table (the main desk-rejection risk)

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 5 | [R] | 102 | `\section*{Metadata}` | `\section*{Required Metadata}` followed by `\section*{Current code version}`, as in the 2026 template. |
| 6 | [R] | 109–149 | 9 rows with old labels: C3 "Reproducible Capsule", C8 "Link to developer documentation", C9 email. The template says "leave the left column as it is". | Exactly 8 rows, with the left column copied from the template: <br>C1 Current code version: v0.2.0 <br>C2 Permanent GitHub link to code/repository used for this code version: `https://github.com/wkchikach99-byte/weatherisk/tree/v0.2.0` <br>C3 Legal Code License: MIT License <br>C4 Code versioning system used: git <br>C5 Software code languages, tools, and services used: Python; SLURM (optional, for HPC runs) <br>C6 Compilation requirements, operating environments \& dependencies: "Python ≥3.10 on Linux or macOS;" followed by the current C7 text <br>C7 If available Link to developer documentation/manual: `https://github.com/wkchikach99-byte/weatherisk#readme` <br>C8 Support email for questions: widad.kchikach@awi.de |
| 7 | [R] | 121 | The Zenodo DOI sits in the old "Capsule" row, and that DOI (20108563) is the v0.1.1 archive, which is the R-dependent code. | Remove it from the table. Add a software reference for weatherisk v0.2.0 with the new DOI (item 49) and cite it in the first sentence of Section 2. |
| 8 | [R] W | 115, 121 | RELEASE BLOCKER comments. The paper says v0.2.0, but the code is at 0.1.1. | After the release, delete both comments and paste the new DOI (repo steps in `submission_todo_v24.md`, A2–A3 and C1–C3). |
| 9 | [T] | C7 text | "R is not required. The validation fixtures were generated once in R…" | No change. It is correct and matches `pyproject.toml`. |

## 3. Motivation

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 10 | [C] TW W | 208–215 | C13: the replacement text dropped Thorsten's sentence "That implementation served its intended purpose very well". | Choose option A (Thorsten's master text, verbatim) or option B (the same text, but "simulation study" instead of "empirical illustration"). Both are spelled out in `judgment_calls_for_review.md` §K. Ideally ask Thorsten. |
| 11 | [C] | 208–215 | `\revt{...}` is wrapped around new sentences that are not Thorsten's. | Remove `\revt` from the text that is not his. This matters only if you build a tracked PDF for him. |
| 12 | [T] TW | 202–212 | Thorsten's "fits max-stable processes to given climate datasets" (C11) sits next to "does not contain the preprocessing of climate-model output". | This is settled by item 10. If you take option A or B, both statements are true. |
| 13 | [T] | 228 | scikit-extremes is cited only as a URL footnote. | Make it a reference with an access date (item 51). |

## 4. Software description

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 14 | [C] | 288 | ν is used before it is defined. Its first definition is in Table 3. | "raised to the power~$\nu$, the degrees of freedom of the extremal-$t$ model, as in Opitz~\cite{opitz2013}." |
| 15 | [F] | 299 | "(default $\varepsilon=5$)" is a CMIP6 default. Example 1 uses ε = 4. | "(default $\varepsilon=5$ in \texttt{CMIP6Config})" |
| 16 | [F] TW | 306 | "(default 2.0 in grid-point units)" is also a CMIP6 default. Example 1 uses 4. | "(default 2.0 grid points in \texttt{CMIP6Config})" |
| 17 | [F] TW | 333 | "(pairs up to 4.0ε apart, three starting points per cluster)" is the CMIP6 path. The synthetic path uses pairs up to 4 grid spacings (verified, `estimation.py:188`). | "(in the CMIP6 pipeline, pairs up to $4.0\,\varepsilon$ apart and three starting points per cluster)" |
| 18 | [T] | 319–322 | $D_2$ is "1 − overlap". The code returns 100 × (1 − overlap). | No change. Average linkage and the quantile threshold do not depend on the scale. |
| 19 | [F] | 386–401 | Figure 1 module labels are wrong (verified against the imports). | Data loading: `cmip6_data` (not `cmip6_pipeline + io`). Stage 1a: `stl`. Stage 3: `density + backend` (not `density + estimation`). Stage 6: `estimation + backend` (not `pipeline`). Stages 1b, 1c, 4 and 5 stay as they are. |
| 20 | [T] | 404–406 | The labels "Step 1a / 1b / 1c" on the right repeat the box text, and the paper elsewhere says "stage". | Delete the three labels. |
| 21 | [T] TW | 411 | "step~1a can be switched off" | "stage~1a can be switched off" |
| 22 | [C] | 426 | "Figure 1 names the module of each stage." | No change once item 19 is done. Today the sentence is false. |
| 23 | [C] | 431 | The heading "(\texttt{scalable}, \texttt{cli})" names `scalable`, which no pipeline or script imports (verified). | "\paragraph{Scalability and execution (\texttt{cli})}" |
| 24 | [F] | 433 | "\texttt{validate} for the synthetic designs". `validate` is only a quick stationary sanity check; the synthetic designs run through `run_nonstationary_pipeline`. | "\texttt{validate} for a quick synthetic sanity check" |
| 25 | [F] | 434 | "\texttt{cmip6} for the Figure~9 case study". `weatherisk cmip6` has no `--no-detrend` option (verified), so the run that gives 104 cannot be made through the CLI. | "\texttt{cmip6} for the case study of Section~\ref{sec:ex-fig9} with STL detrending (the run without detrending uses \texttt{scripts/reproduce\_fig9.py --no-detrend})" |
| 26 | [T] | 435 | NOAA is not expanded. | "the US National Oceanic and Atmospheric Administration (NOAA) Climate Prediction Center (CPC)" (+6 words), or drop "NOAA". |
| 27 | [F] | 436–438 | "The pipeline writes a checkpoint after each step, and after each cluster within the final step". Only the CMIP6 pipeline writes checkpoints. | "The CMIP6 pipeline writes a checkpoint after each stage, and after each cluster within the in-cluster re-estimation" |
| 28 | [T] | 245, 739 | CLI and API are not defined. | Optional: write "command-line interface (CLI)" at l.245. Leave API as it is. |

## 5. Validation (§2.3)

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 29 | [T] | 504 | "With R removed from the search path, all 235 tests pass in about 2 min". | No change. Rechecked with R truly hidden: 235 passed in 137 s. |
| 30 | [T] TW | 511–512 | `\revt` is wrapped around the pinned URL. Thorsten asked for a URL; the pinned commit is ours. | Remove `\revt` here when cleaning the source (item 54). The text stays. |
| 31 | [T] | 532 | Ellipse overlap: only the relative bound is given. The test also has an absolute bound of 1e-14. | Optional: "relative $\leq 5\times10^{-15}$, absolute $\leq 10^{-14}$". |
| 32 | [T] | Table 1 | Overfull box of 2.6 pt. | Invisible in print. No change. |

## 6. Example 1

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 33 | [F] | 556–558 | The run settings are not stated, and the Section 2 defaults do not apply to this example. | After "took 8.4~hours of wall time", add: "The preset sets $\varepsilon$ and the smoothing radius to 4 grid points, and the in-cluster re-estimation uses pairs up to 4 grid points apart." (+25 words) |
| 34 | [C] | 570 | "(their Figure~3)" next to our own Figure 3. | "(Fig.~3 in~\cite{contzen2025})" |
| 35 | [T] | Figure 2 | Panel titles and ticks are very small at print width. | Optional: regenerate with a larger font. |

## 7. Example 2 ("Figure 9" means Contzen's figure; the paper has only Figures 1–3)

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 36 | [C] TW | 590 | Title: "Revisiting the Figure~9 case study" | "Revisiting the published case study" (the "--- AWI-ESM-1-1-LR global precipitation" part stays) |
| 37 | [T] TW | 594 | "revisit the Figure~9 analysis of Contzen et al." | No change. It clearly refers to their figure. |
| 38 | [C] TW | 599 | Thorsten's C3 body definition was lost in the cut. AWI is used in the body four times without being defined. | "of the \revt{Alfred Wegener Institute (AWI)} \emph{albedo} HPC cluster" (+3 words, his wording) |
| 39 | [C] | 599–601 | "…with the SLURM script \texttt{hpc/run\_fig9\_pyonly.slurm}, once with STL detrending and once without (Listing~\ref{lst:fig9})". The listing shows the Python calls, not the SLURM script. | "…with the SLURM script \texttt{hpc/run\_fig9\_pyonly.slurm}, which runs the two commands of Listing~\ref{lst:fig9}." |
| 40 | [C] TW | 616 | Table caption: "for the Figure~9 case-study analysis" | "for the case study of Example~2" |
| 41 | [C] TW | 621 | Listing caption: "Running the Figure~9 case study on HPC" | "Running the case study of Example~2 on HPC" |
| 42 | [C] | 658 | Table caption: "Cluster counts for the Figure~9 case study" | "Cluster counts for the case study of Example~2" |
| 43 | [T] | 687–689 | "the only later change, the R ranking of tied annual maxima in the madogram". Two later commits touch the code; only one changes the computation. | "the only later change to the computation, R's rank rule for tied annual maxima in the madogram," |
| 44 | [C] | 696 | Figure 3 caption: "The Figure~9 case study of Contzen et al. recomputed…" | "The case study of Contzen et al.~\cite{contzen2025} (their Figure~9) recomputed…" |
| 45 | [R] W | 688 | Cites commit `7561a6e`, which exists only locally. | Merge without squashing so the commit survives on GitHub. |
| 46 | [T] | Figure 3 | Coastlines only, no borders. | No map disclaimer needed. |

## 8. Impact and Conclusions

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 47 | [T] TW | 728 | "has $n^2$ entries". l.446 says $n(n-1)/2$. Both are true: the full matrix has n² entries and n(n−1)/2 distinct ones. | Optional: "has $n(n-1)/2$ distinct entries". |
| 48 | [F] TW | 760–762 | "The full Figure~9 case study … ran in 26 to 36~hours on 64~cores." This covers stages 3–6 only. | "Stages 3--6 of the case study on the AWI-ESM-1-1-LR T63 grid took 26 to 36~hours on 64~cores." |

## 9. References

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 49 | [R] W | new | There is no reference to weatherisk itself. The DOI moves here from the metadata table (item 7). | `[software] W.~Kchikach, T.~Dickhaus, G.~Lohmann, weatherisk: spatial clustering of climate extremes via max-stable processes, version 0.2.0, Zenodo (2026), \url{https://doi.org/10.5281/zenodo.<new>}`. Use the repo's title. |
| 50 | [C] | 866–874 | pyextremes and SpatialExtremes have no version and no access date. The guide wants the URL and the last-accessed date for web references. | pyextremes: "version 2.x [software], … [accessed DD November 2026]". Use the version we compared against. SpatialExtremes: "R package version 2.1-0 [software], CRAN (2022), \url{https://doi.org/10.32614/CRAN.package.SpatialExtremes}". |
| 51 | [T] | 228 | scikit-extremes is a footnote. | Make it a bibitem with "[accessed DD November 2026]" and cite it. This also removes the footnote from the word count. |
| 52 | [C] | refs | Eight references have DOIs that are not given (checked on Crossref). | davison2012 10.1214/11-STS376; padoan2010 10.1198/jasa.2009.tm08577; dehaan2006 10.1007/0-387-34471-3; schlather2002 10.1023/A:1020977924878; opitz2013 10.1016/j.jmva.2013.08.008; cooley2006 10.1007/0-387-36062-X_17; saunders2021 10.1007/s10687-020-00395-y; bauer2021 10.1038/s41558-021-00986-y |
| 53 | [T] | 859–864, 878–889 | contzen2025 has no issue number. The journal abbreviations are not in LTWA form. | "Extremes 28\,(4) (2025) 713--737"; "J.\ Off.\ Stat."; "J.\ Multivar.\ Anal." |
| 53a | [T] W | 912 | danek2020 gives "version 20200212". | Confirm that this is the version of the files you used on albedo. DataCite lists 20250527 as the metadata version. |

## 10. End statements

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 54 | [C] W | 815–822 | AI declaration: "used Claude (Anthropic) in order to assist". The guide's text is "used [TOOL] to [REASON]". The AUTHOR CHECK comment is still there. | "used Claude (Anthropic) to assist with…". Confirm "the authors" or "the first author", and delete the comment at l.817. |
| 55 | [R] W | 833–836 | The funding paragraph is not in the Elsevier form and has no sentence on the funders' role. The guide says to state that the funders had no involvement. | "This work was supported by the Helmholtz Association [grant number HIDSS-0005, MarDATA Helmholtz School for Marine Data Science] and the European Research Council [ERC Synergy Grant i2B, grant number 101118519]. The funders had no role in study design; in the collection, analysis and interpretation of data; in the writing of the report; or in the decision to submit the article for publication." Confirm whether Dickhaus has funding to add. |
| 56 | [C] | 829 | "on which the Figure~9 case-study pipeline was executed" | "on which the case study of Example~2 was run" |
| 57 | [C] W | 830–831 | "We thank J. Contzen for the original R implementation". Nothing in the paper says he agreed to his code being in our public repo. | Once he agrees, add "and for permission to include it in the repository". If he does not, see `submission_todo_v24.md` A1. |
| 58 | [T] | 806–812 | Data availability: the Figure 9 cluster outputs are not deposited. | Optional: deposit them on Zenodo as a [dataset] and cite them. |

## 11. Length (after the changes)

| # | Mark | Now | Change to |
|---|---|---|---|
| 59 | [R] | The count is 3675 words. That count includes the abstract, the text and the captions. It does not include the end statements (235 words) or the two result tables and the parameter table (112 words). Under the strictest reading of the guide, the total is already about 4022. Items 2, 14, 25, 26, 33, 38 and 55 add about 110 words. | When updating, offset the additions with cuts of about 150 words in non-supervisor text. The longest candidates are the STL-variants passage (l.669–680) and the computational-cost detail (l.634–645). The target is under 4000 even under the strict reading. |

## 12. Source file for upload

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 60 | [C] | 18–39 | "THORSTEN-18 TRACKING VERSION" header and the `\rev*` macro definitions. | Upload a clean copy: macros inlined, header removed. Keep the tracked v24 for the supervisors. |
| 61 | [T] | 97, 115, 121, 817 | Internal comments. | Remove them in the clean copy. |
| 62 | [T] | preamble | The template sets `\setlength{\parindent}{0pt}`. | Add it. |

## 13. Highlights (`highlights.txt`)

| # | Mark | Now | Change to |
|---|---|---|---|
| 63 | [C] | Bullet 5: "Recovers the published 104 EDC clusters for CMIP6 precipitation, no detrending". EDC is not defined in the highlights, and the detrended result is left out (as in item 2). | "Recovers the published 104 clusters without detrending; STL detrending gives 142" (80 characters) |

## 14. Cover letter (`cover_letter.tex`)

| # | Mark | Line | Now | Change to |
|---|---|---|---|---|
| 64 | [C] | 75–79 | "the mandatory metadata tables" and "a code appendix with the code for both examples". v24 has one table and no appendix. | "The manuscript contains the mandatory code metadata table, a description of the six-stage pipeline, two worked examples (a synthetic validation and the AWI-ESM-1-1-LR case study) with their result figures, and a short command listing; further code examples are in the repository README. The source code is publicly available on GitHub and permanently archived on Zenodo." |
| 65 | [T] | 10–11, 51–53 | Short address and the old DOI. | Use the same address as item 1 and the new DOI. |
| 66 | [T] | 51, 74, 82 | Internal comments. | Delete them. |

## 15. Upload package (Editorial Manager)

| # | Mark | Item |
|---|---|---|
| 67 | [R] W | Complete the declarations tool and upload its .docx. The guide says this must always be done. |
| 68 | [C] | Upload Figures 1–3 as separate files. Figure 1 is TikZ, so export it to `Figure_1.pdf`. |
| 69 | [C] | Upload: the manuscript PDF; a source .zip (clean .tex and both PNGs; no .bib is needed, since the references are inline); `highlights.txt`; the cover letter. Also enter the funders and the GitHub URL. |

## Decisions only Widad (or the supervisors) can make

- Item 10, C13 wording. Ideally ask Thorsten.
- Item 57 and Contzen's permission for `r_code/`.
- Items 1, 53a, 54 and 55: the facts.
- 13 green sentences that were reworded in v23 (listed in `judgment_calls_for_review.md` §K). Ask Gerrit whether to restore them. None of them is wrong.
- C4: no change recommended. Restoring Thorsten's full clause would overstate what the validation covers.
- Not a paper change, for your note to Thorsten: the run that matches 104 fits stationary GEVs to 156 years without detrending. That is the situation behind his C9 and C14. The paper makes no validity claim for it.
