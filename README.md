# hdym-fakeon — the spin-2 graviton in the presence of a fakeon

A research paper, backed by simulation, on a part of the study of the spin-2 graviton with a fakeon that the
community has not yet addressed. The comparison theory is agravity (Salvio & Strumia 2014).

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
repeller above, the ghost ends the chaos). The paper asks what the fakeon-projected spin-2 mode does on the same
flow, where the projection stops being usable, and how both compare with agravity's unprojected ghost.

## Where things stand

| step | status | where |
|---|---|---|
| Literature: is it unaddressed, what are the exact equations | done 30 Oct 2026 | `notes/2026-10-30_literature_sweep.md` |
| Reconciliation of `bianchi9/` with Lehners & Stelle 1909.01169 | mapped; rerun pending | `notes/2026-10-30_lehners_stelle_reconciliation.md` |
| Path B step 1: the projection on a one-mode toy | done 30 Oct 2026 | `fakeon/toy/`, `notes/2026-10-30_step1_toy.md` |
| Step 2: the projection on Bianchi I + one Taub wall, constraint reduced alongside | deferred: the ħ → 0 check of the quantum result, run when that needs it | — |
| Step 3: full Bianchi IX, add R² for agravity, compare on the same backgrounds | — | — |
| Route 1 (perturbative, fakeon against agravity at one loop), rung 1: the scalar bubble two ways against the published numbers | done 30 Oct 2026 | `loops/scalar/`, `notes/2026-10-30_rung1_fakeon_scalar_bubble.md` |
| Route 1, rung 2: fakeon QED, then Yang–Mills self-energies and cutting equations | next | — |
| Route 1, rungs 3–4: the h/χ/φ multiplet against agravity, then the unaddressed piece | — | — |

## Layout

```
README.md            this file
notes/               dated working notes and write-ups (lab-notebook style lives here, not in READMEs)
fakeon/toy/          one-degree-of-freedom fakeon projection: code, logs, a short README
loops/scalar/        rung 1 of the perturbative route: the scalar bubble with the fakeon, Lee-Wick and
                     Feynman-Wheeler prescriptions, two methods, certified against the published closed forms
```

The goal is a quantum treatment of the spin-2 mode with a fakeon against agravity; the classical projection work
(`fakeon/toy/`, the Bianchi steps) is its ħ → 0 check. The perturbative route is a ladder — scalar, QED,
Yang–Mills, the gravitational multiplet — each rung reproducing published results before anything new; later
rungs get sibling folders `loops/qed/`, `loops/ym/`, `loops/gravity/`.

Later folders follow hdym's pattern: one folder per sub-question, each with a short README (what it asks, how to
run, what is in `results/`), code at the top of the folder, outputs written straight into `results/`.

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
