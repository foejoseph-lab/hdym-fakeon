# hdym-fakeon — the spin-2 graviton in the presence of a fakeon

A research paper, backed by simulation, on a part of the study of the spin-2 graviton with a fakeon that the
community has not yet addressed. The comparison theory is agravity (Salvio & Strumia 2014).

## The rule of the project (author, 31 Oct 2026)

Every result on the path to the goal, from Utiyama–DeWitt 1962 up to at least agravity, is redone with the technique of its day and then with each later technique **until the
techniques stale out**, with the intent of extracting usable and missed mathematical techniques. Agreement
between techniques is the pass condition, not the product; the product is what each technique sees that the
others do not, recorded in `notes/TECHNIQUES.md` at the end of every rung. A rung is not left with a known open
computation while the next one starts. Results before 1962 that the path rests on are presumed fact unless
history marked them wrong; Pais–Uhlenbeck 1950 is such an input. **The route is quantum.** Classical work
(`fakeon/toy/`, any Bianchi IX projection) is ħ → 0 check material for the quantum result, done only when that
result needs it, and is never a step of the plan.

## The question

Quadratic gravity, S = ∫√−g [R/2 + R²/(6f₀²) − C²/(4M²)], carries a massive spin-2 mode of mass M with the
wrong-sign kinetic term. Agravity keeps it as a genuine degree of freedom and argues its runaway is metastable.
Anselmi's fakeon prescription removes it from the spectrum by a quantization rule whose classical limit is a
nonlocal, lower-order equation of motion: the fakeon's own equation is solved with the half-retarded,
half-advanced Green function and substituted back (the "projection"). That projection has been carried out only
on flat space, on FLRW (where the Weyl tensor vanishes and the spin-2 mode drops out), and for linear tensor
perturbations about de Sitter. On any background with C ≠ 0 the spin-2 fakeon has never been projected.

The Bianchi IX minisuperspace — the Mixmaster billiard — is the simplest background with C ≠ 0, chaotic, and
with the unprojected ghost's behaviour already mapped (hdym `bianchi9/`: adiabatic below H ≈ 0.2M, a power-law
repeller above, the ghost ends the chaos). The paper asks, in the quantum theory, what the spin-2 fakeon does on
a background with C ≠ 0 and how that compares with agravity's unprojected ghost; the classical projection on
the same flow is the ħ → 0 check of that answer, not the answer.

## Where things stand

| step | status | where |
|---|---|---|
| Literature: is it unaddressed, what are the exact equations | done 30 Oct 2026 | `notes/2026-10-30_literature_sweep.md` |
| Reconciliation of `bianchi9/` with Lehners & Stelle 1909.01169 | mapped; LS rerun done on the author's side 31 Oct (hdym `bianchi9/results/ls_rerun.*`) and read 31 Oct evening (note §"Read") | `notes/2026-10-30_lehners_stelle_reconciliation.md` |
| ħ → 0 check material (not the route): the classical fakeon projection on a one-mode toy | done 30 Oct 2026; kept for the ħ → 0 check of the quantum result, used only when that needs it | `fakeon/toy/`, `notes/2026-10-30_step1_toy.md` |
| Chronology 1939–2014 and the ladder | written 31 Oct, for the author to correct | `notes/2026-10-31_chronology_1939_2014.md` |
| Rung 1: one-loop matter on a classical background (Utiyama–DeWitt 1962), four techniques | done 31 Oct 2026 (DeWitt's 1965 covariant recursion added in the evening) | `matter_loop/`, `notes/2026-10-31_rung1_utiyama_dewitt_{setup,results}.md`, `notes/2026-10-31_rung1_b2_dewitt_recursion_setup.md` |
| Rung 2: the one-loop bubble under the ghost prescriptions (Lee–Wick 1969 / Nakanishi 1971) | p⃗ = 0 case done 30 Oct (the instrument check); p⃗ ≠ 0 and the Minkowski real-axis case next | `bubble/`, `notes/2026-10-30_rung1_fakeon_scalar_bubble.md` |
| Rung 3: one-loop graviton on a background ('t Hooft–Veltman 1974) | — | — |
| Pause before 1977; then Stelle and the ladder's later rungs (chronology note §5) | — | — |

## Layout

```
README.md            this file
notes/               dated working notes and write-ups (lab-notebook style lives here, not in READMEs)
fakeon/toy/          one-degree-of-freedom fakeon projection: code, logs, a short README (rung 0)
bubble/              the one-loop bubble with the fakeon, Lee–Wick and principal-value prescriptions, two methods,
                     certified against the published closed forms (rung 2 at p⃗ = 0; the instrument check)
matter_loop/         one-loop quantized matter on a classical background: the 1962 computation redone with the
                     techniques of 1962, 1951/65 and 1972–74 (rung 1)
```

Folders are named by **concept**, not by year, author or method, and techniques accumulate inside them (the
author, 31 Oct 2026); a new folder appears only when a new concept does. The ladder itself is historical and is
in `notes/2026-10-31_chronology_1939_2014.md` §5: the toy (1950) → matter on a background (1962) → the bubble
under prescriptions (1969–71) → the graviton loop (1974) → pause → Stelle 1977 and onward to agravity. Each rung is
redone with the technique of its day, then with each later technique; the goal — a quantum treatment of the spin-2
mode with a fakeon against agravity — is unchanged; the classical projection work is ħ → 0 check material, not the route.

Each folder follows hdym's pattern: a short README (what it asks, how to run, what is in `results/`), code at the
top of the folder, outputs written straight into `results/`.

## Relation to the other repositories

Nothing here is a copy of hdym or hdym-su3. Material from them is referenced by path, never duplicated:

- hdym `bianchi9/` is the Bianchi IX machinery (Cartan derivation, generated right-hand sides, symplectic form,
  Kasner tools) and the unprojected-ghost results this paper compares against. Code here that needs it imports
  from the hdym checkout at `HDYM_ROOT` (default `C:\GitHub\hdym`), read-only; anything that has to diverge is
  derived afresh here and says what it was derived from.
- hdym `islands/code/` is the certification pattern (untrusted generator, trusted dyadic checker, Lean 4 as the
  design target) to follow when a result here is certified.
- hdym `paper/main.tex`, item (v) of its outlook, points to this repository; the Lee–Wick Yang–Mills paper stays
  what it is.

## Conventions

As in hdym: logs end `# END <name> elapsed=<s>s status=...`; outputs go straight into `results\` from the code;
signed and complex values are kept; the ghost is treated exactly wherever the theory allows and a naive operator
is labelled as such; any run is kept under 20 minutes predicted wall clock, and every log carries its elapsed
seconds so later changes can be judged against them. Physics and judgement calls are the author's.
