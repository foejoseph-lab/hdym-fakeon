# Brief 14 (DRAFT by Claude, 31 Oct 2026 night, for the author to correct): rung 1 closed four ways

**This brief stands alone.** Briefs 10–13 are not in git and are not re-supplied by the author. Every rule a
session is expected to follow is written out below. Rules marked **[R]** are reconstructed from how the repo's
notes and code cite earlier briefs (e.g. "brief 12 §2.3"), not from the briefs themselves; the author corrects
them here if the reconstruction is off. Earlier citations of "brief 12 §x" in the notes are left as they are;
this file is what they now resolve to (§2 maps them). From this brief on, briefs live in `notes/briefs/` and
each one restates everything.

## 1. The goal (unchanged, and first)

> A research paper, backed by simulation, on a part of the study of the spin-2 graviton in the presence of a
> fakeon that the community at large has not yet addressed.

(Reworded from brief 13 for grammar only; the author to confirm.) "Unaddressed" means as of 2026 and is
rechecked at each step that leans on it. Comparison theory agravity (Salvio & Strumia 1403.4226). Quantum
treatment; classical work is the ħ → 0 check. Perturbative route. The paper is the focus, not a deadline;
whether one results is open (author, 31 Oct), "a mountain of salt" on the word paper.

## 2. Standing rules

### 2a. Method and order

1. **Historical order, technique of the day first, then each later technique, until the methods run out.** At
   each paper on the path: do the computation with the technique current at its date (computer-aided), then
   with what was used next, and so on toward 2026; expect rough agreement; where they disagree that is a
   finding. Pre-1977 results that a paper rests on are presumed fact unless historically known wrong. The sweep
   is a vacuum: forward pass, back pass for dead vines, forward for alternative paths that rejoin (walked one
   level back, not to the invention of the wheel), repeat. "Vacuum forward then back, then accept what a steam
   vac gives or rip up as a final" (author).
2. **Stumble legitimately before switching.** Diagrams / method P before anything else; the era's common
   technique only after a real stumble, with the stopping criterion written down before the attempt starts.
   Where rules 1 and 2 conflict the ladder table decides and says so in its status column.
3. **Method P and method M [R, cited as brief 12 §2.1].** Method P = Feynman parameters with the average
   continuation; method M = momentum-space residues, the check and balance. Every number that matters has
   both, or says why not.
4. **Yang–Mills and purely classical work enter only with the author's written yes at the time [R, §2.2]**, and
   only if a paper on the path requires them. Finite QED (Lee–Wick 1970), Einstein–Dirac, Einstein–Yang–Mills
   are cited as lines, not computed, unless he says otherwise. hdym and hdym-su3 are the Yang–Mills side and
   stay off the ladder.
5. **Each paper is summarised to the author before anything is built [R, §2.3]**: the setup note, the
   computation pinned down so two people could do it independently, what counts as agreement.
6. **Standing checks vs. tools [R, §2.4].** A consistency check that cannot fail by construction (Schwarz
   reflection of a transcribed expression, etc.) is a tool to call when something looks odd, not a standing
   check in the run script.
7. **A folder appears when there is code, where the author says [R, §2.5]**; folders are named by concept,
   never by year, author or method; techniques accumulate inside; a new folder only for a new concept.
   Current: `fakeon/toy/`, `bubble/`, `matter_loop/`.
8. **One run, then talk [R, §2.6].** No batch of variants before the first result has been discussed.
9. **Review against intent first [R, §2.7], then the independent review.** The independent review is a
   product/bias review run as a **fresh agent, before a patch is sent**: did the framing lead the computation,
   were expectations written after the answers, does the narrative exceed the logs. Findings are folded in
   and the patch re-made. This is a visible gate in §3 of every brief. (Tonight's patch went out once before
   its review; it was replaced.)
10. **Check, don't assume** — about files on his machine, about what ran, about what a result means. "Assume
    is an ass of you and me."
11. **Predictions written in code headers stay in the record when they fail** (`results/superseded/`). A failed
    prediction relabelled as "expected" without saying so is exactly what the bias review is for. A predicted
    number whose sign depends on a convention is written with the convention, or scored as a slip.
12. **Single-method or single-graph numbers are not checks.** Only the full sum means anything.
13. **Controls are logged like results.** A run made to show what goes wrong without an ingredient gets its own
    log and a paragraph in the results note.
14. **The fakeon enters only at the finite parts.** Divergences and beta functions are prescription-independent
    (1803.07777: the renormalization "coincides with the one of the Euclidean version of the theory").

### 2b. Running and hardware

15. **Smoke test here, full run there (author, 31 Oct night).** Claude writes the script with its predictions in
    the header, runs it once in the sandbox at the smallest setting that exercises every code path and confirms
    a sane log. The full run is the author's: exact command line(s) for `C:\GitHub\hdym-fakeon`, predicted wall
    clock, `--jobs 22` where it applies, output straight into `results\`. The author commits the logs; Claude
    reads them from `origin/main` next session. Nothing is written up as a result from a smoke test.
16. **Time limit:** any run under 20 minutes predicted wall clock; every log ends
    `# END <name> elapsed=<s>s status=...`. Hardware: 22 cores, RTX 4080; sympy is single-threaded and rational;
    GPU only if numerics move to float arrays. The sandbox has 2 cores.
17. **Conventions (from hdym):** outputs written straight into `results/` by the code; signed and complex
    values kept; the ghost treated exactly wherever the theory allows and a naive operator labelled as such;
    every number in a README traces to a committed log, script or figure named beside it; superseded code kept
    in `results/superseded/`.

### 2c. Sources and access

18. **Paywalls:** not buying, not borrowing, not scanning. Forge on within open sources (arXiv, INSPIRE,
    numdam, author pages); mark what is unseen; re-score in one line if a text ever arrives. "Any true
    head-scratch later I'll pay the $40" is the author's call alone. Links Claude asks the author to paste are
    arXiv, INSPIRE, numdam or an author's own page, nothing else; anything else gets "is that open?" first.
    DeWitt 1965 (~$600) is the standing example of what this forbids; cited as unseen.
19. **How a page becomes readable in the sandbox:** the author pastes the URL into the chat. It costs nothing.
    Open-access PDFs he downloads stay on his side, out of git — cite the URL.
20. **Reference marks:** [v] verified against a listing this session; [m] from memory, textbook-stable; [m?]
    from memory, check before leaning on it. Verified 31 Oct night: 1803.07777 = Anselmi & Piva, JHEP 05
    (2018) 027. Still [m?]: Fradkin & Vilkovisky PRD 8, 4241 (1973); Modesto–Shapiro 1512.07600 = PLB 755, 279.

### 2d. Documents

21. **Briefs are in git from now on** (`notes/briefs/briefNN.md`), each self-contained; a new brief
    supersedes the previous one in full. Drafts by Claude are marked DRAFT until the author corrects them.
22. **Instrument checks [R, brief 12 §4]:** `bubble/` at p⃗ = 0 is an instrument check, done and kept; it
    confirms the average continuation is under control on the simplest integral. Not a result.
23. **Housekeeping [R, brief 12 §5]:** naming (method P/M, "principal-value (FW)"), logs regenerated after
    renames, histories kept at the foot of notes.
24. **Working conventions [R, brief 10 §5]:** literal paths, one command per line, no wildcards or environment
    variables (the author is a month into git); patches `git am`-ed from `C:\Users\Owner\Downloads`.

## 3. Repositories and the patch gate

| repo | state 31 Oct night |
|---|---|
| github.com/foejoseph-lab/hdym-fakeon | `main` 90ea4f4 on origin. **Patch 0001 (revised after the fresh-agent review) sent tonight, to be `git am`-ed**: rung 1 (b2), LS rerun read, 1803.07777 read, sweep row, this brief. **Then the author regenerates the two `rung1_dewitt_*.log` files with the `matter_loop/README.md` commands and commits them** (first use of rule 15). |
| github.com/foejoseph-lab/hdym | `main` 3a34e0f; `bianchi9/results/ls_rerun.log` read tonight (reconciliation note, new section). |
| github.com/foejoseph-lab/hdym-su3 | `main` 9e99e19, untouched; off the ladder (rule 4). |

Patch flow: made against the current `main` → fresh-agent review (rule 9) → findings folded in → sent →
`git am` by the author → he pushes promptly; check `git log origin/main` at session start.

## 4. The ladder (chronology note §5), status

| rung | concept folder | status |
|---|---|---|
| 0 | `fakeon/toy/` — Pais–Uhlenbeck 1950, the fakeon projection on one mode | done 30 Oct |
| 1 | `matter_loop/` — one-loop matter on a classical background (Utiyama–DeWitt 1962) | **done 31 Oct, four techniques**: (a) 1962 cutoff, one placement; (b1) proper time on the graphs; (b2) DeWitt 1965 recursion on generic backgrounds; (c) dim reg. Log order identical four ways = 't Hooft–Veltman 3.40; (b2) adds the Euler coefficient −1/360. Open: 1/p⁴ resummation printed and fitted; unshifted-cutoff placement; spinor/vector matter |
| 2 | `bubble/` — the one-loop bubble under ghost prescriptions (Lee–Wick 1969 / Nakanishi 1971) | p⃗ = 0 done 30 Oct; **p⃗ ≠ 0 (Nakanishi 1971 Lorentz-invariance test) and the Minkowski real-axis prescription (Aglietti–Anselmi 1612.06510) next** |
| 3 | (new folder when it starts) — one-loop graviton on a background ('t Hooft–Veltman 1974) | not started; PDF pp. 80–94 unread; Einstein–Maxwell (PRD 10, 401) as the algorithm check |
| — | **pause** before Tomboulis 1977 / Stelle 1977 | review against intent with the author; fix the next rungs |

Later: Stelle 1977 proper; the 1978–85 beta-function history (Julve–Tonin diagrams first, by rule 2); Lee–Wick →
fakeon at the finite parts; the h/χ/φ multiplet against agravity; the unaddressed C ≠ 0 piece.

## 5. First tasks next session

1. **Rung 2 proper in `bubble/`**: the bubble at p⃗ ≠ 0 under the Lee–Wick contour with real spatial momenta
   (Nakanishi 1971), CLOP, the 1703.04584 region-wise continuation, the average continuation; then the
   Minkowski real-axis case of 1612.06510 as the fifth prescription. Summaries of Nakanishi PRD 3, 811 and
   3, 3235 and Lee–Wick PRD 3, 1046 to the author first (open sources only; if paywalled, from abstracts and
   from 1703.04584's account of them, marked so). Smoke test here, full run there.
2. **Read 1803.07777 §3 in full** for the discrepancies with agravity's eqs. (54a) and (55), and check whether
   Salvio 1804.09944 (our anchor) absorbed them. Agravity is the comparison theory; not optional. (Author pastes
   `https://arxiv.org/pdf/1803.07777v2` when asked.)
3. Rung 1 leftovers if time allows: the finite p⁴ log p² part of Π and the 1/p⁴ resummation; the
   unshifted-cutoff placement.
4. INSPIRE citation walk from 1901.09273 and 1909.01169, 2024 onward — the sweep's flagged gap. Author pastes
   the two INSPIRE search URLs when asked, or does it by hand and pastes the titles.

## 6. Open with the author

- The [R] rules above: confirm or correct.
- `notes/` vs hdym's folder-README + `MATH.md` + `rabbit_holes/`. Not forced yet.
- Whether spinor and vector (abelian) matter are added to `matter_loop/` before rung 3.
- The three remaining suggested actions of the reconciliation note (eq. (11) check against `L`; the paper
  clause; the 1903.01044 sign match).
- The wording of §1.

## 7. Things the author said that a future Claude should keep hearing

- "If this brief reads too restrictive it's since I went too far restricting it. Just point out if anything
  doesn't make sense to do and why." — Push back in writing; don't comply silently, don't go forth alone.
- "Legit stumble before initial avoidance."
- "Flogging a dead horse, but not without reason."
- "Assume is an ass of you and me" — don't be Mulder.
- The vacuum metaphor for the sweep; the plasma/QGP metaphor for the fakeon kink (regionwise analytic — a
  phase boundary, not a branch point of any one function).
- "Everything is a vacuum forward then back, then accept what a steam vac gives or rip up as a final."
- "If Riem² needs anything then ok" — tooling is the author's to provide when asked, after the in-house route
  has been tried (tonight sympy did (b2) in 6.5 min; no new tool).
- He is learning git as this goes; literal paths, one command per line, no wildcards or env variables.
- He is old school, watches the sessions closely, and notices when one gets sharp. Point out, don't snap.
- "Why if not there": a brief must not call on priors that are not in front of the reader.
