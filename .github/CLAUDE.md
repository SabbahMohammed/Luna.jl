# LupoAirOsc.jl — working context

Place at repo root of `LupoAirOsc.jl`. Branch: `adding-Luna-rerun`. Companion
Luna fork: `SabbahMohammed/Luna.jl`, branch `dissociation-support`, pinned into
LupoAirOsc via a local `Pkg.develop(path=...)` (check `Manifest.toml`'s
`[[deps.Luna]] path = ...` line — if it's missing or points elsewhere, the
pin has drifted and needs re-`develop`ing before anything below is trustworthy).

**If you are a new session picking this up: read this whole file before
touching code.** It exists so you don't have to re-derive what's already been
found. See "Keeping this file current" at the bottom — you are expected to
extend it, not just read it.

## What this simulates

Companion code to Sabbah, Brahms, Belli & Travers et al., *"Periodic ultrafast
soliton dynamics induced by the molecular dissociation inside an air-filled
hollow-core optical fiber."*

A 30 fs / 800 nm / 1 kHz Ti:Sapph pump self-compresses in a single-ring AR-HCF
(30 µm core, `a = 15e-6`, ~210 nm walls) down to ~1.7 fs and ~10¹⁴ W/cm². At
that intensity O₂ and N₂ dissociate **via the near-IR pump field itself**
(multiphoton/tunnelling, ADK — not linear UV photolysis; the paper explicitly
drops the λ<240 nm term as negligible because that light scatters out through
the fibre resonance). Freed O atoms build O₃ through the Chapman cycle. O₃'s
Hartley band (~255 nm) sits on the UV resonant dispersive wave (RDW), so
accumulating ozone both absorbs the RDW and detunes its phase matching for the
*next* pulse.

**The central scientific question** (unchanged since the project started, and
**not yet worked on this session** — see "Current stage" below): the paper's
own model only reproduces the **He-O₂** case (pure oxygen chemistry gives a
monotonic RDW red-shift settling to equilibrium, no periodicity). The
oscillation only appears experimentally with **air** (N₂ present); the paper
calls full N-chemistry "exceedingly difficult" and explicitly leaves the
mechanism unresolved. `ReactionDiffusion.jl` carries a 21-reaction Chapman +
O(¹D) + NOx network over 10 species, so reproducing the oscillation is the
actual goal — everything so far has been getting the simulation to a state
where that's even attemptable.

Two reference targets from the paper (Sec. III E, Fig. 5-7):
- **He-O₂**: 79%-21%, 12 bar, 2.5 µJ, 22 cm fibre. Ozone reaches
  **9.3×10²⁴ m⁻³ after 5 s** (~3% of local gas density).
- **Air / N₂-O₂**: 79%-21%, 6 bar, 2.2 µJ. This is the case with the
  unexplained ~100 s oscillation period.

## Current stage

All work in this session has been on the **He-O₂ validation case** —
deliberately, since it's the case the paper's own model reproduces, so it's
the right thing to get trustworthy before attempting air. The air/oscillation
question itself has not been touched.

**What's now believed solid:**
- Environment: both repos build and run cleanly from a fresh clone (given the
  dev-path pin above is intact).
- Self-compression + RDW generation reproduces cleanly in bare Luna
  propagation (N₂/O₂ *and* He/O₂), confirmed against a pristine
  pre-session Luna baseline.
- O₃'s dispersion (real + imaginary refractive index) is genuinely wired into
  the propagation, not just present in `PhysData`, **and is now applied along
  z rather than frozen at the fibre entrance** — see the linop fix below,
  which is the most important thing in this file. Confirmed by matched
  propagations with/without ozone (`examples/o3_dispersion_check.jl`) and by
  the localized-vs-uniform comparison (`examples/o3_localization_check.jl`).
- σ(Hartley peak) matches the literature cross-section to 0.2%; Φ(O¹D)
  branching at the peak matches literature (~0.9).
- Atom conservation in the O₂/O₃ photolysis bookkeeping (`apply_lunaion`) is
  now verified tight (was silently destroying O atoms every pulse before).

**Load-bearing finding, not yet fully chased down:** the propagation is
*extremely* sensitive to O₃ density near some threshold. A 2% O₃ fill doesn't
fade the RDW in place — it **relocates the entire band** (~250 nm →
~310-320 nm), a phase-matching (real-index) effect, not just Hartley
absorption. This means small errors in O₃ bookkeeping can flip qualitative
behaviour, not just shift numbers slightly. Treat any O₃-density result near
this range with real suspicion until cross-checked.

**Cautionary tale, worth reading before trusting any O₃ result:** for most of
this session the answer to "why doesn't the real run's RDW shift, when a
uniform fill at lower density shifts it dramatically?" was pursued as
*physics* — localization, O₂ depletion, spectral downsampling, Monitor
timing — and a confident, wrong conclusion ("localization is the whole
story") was reached and written down. The actual cause was the const-linop
bug below: an infrastructure defect that only manifests for a NON-UNIFORM
fill, which is exactly why every uniform control case looked fine and
appeared to corroborate the wrong story. The lesson to carry: when a
comparison isolates "localized vs uniform" and the localized case is the odd
one out, suspect the code path that only the localized case exercises before
inventing physics to explain it.

**Fixed this session, most recent first (see git log for full detail —
this is the "why", not a duplicate of commit messages):**
- **The linear operator was frozen at z=0** (`PropAir.jl`, Luna
  `Capillary.jl`). `setup_propair` used `LinearOps.make_const_linop`, which
  evaluates `Modes.β`/`Modes.α` once with no `z` kwarg and whose `βfun!`
  ignores z — so every propagation ran the whole 22 cm with the dispersion
  AND the loss of the mixture at the fibre entrance. Uniform fills are exact
  under that, which is why the test cases never caught it. Chemistry-produced
  ozone is a narrow bump at z ~ 10 cm whose fraction at z=0 is ~1e-34, so its
  refractive index — the real part that moves RDW phase-matching and the
  imaginary part that IS the Hartley band — was discarded in full; only the
  nonlinear polarisation, via `densityfun(z)`, ever saw it. Measured gap at
  the ozone peak: 8233 rad/m. Now `make_linop`. Switching required three
  fixes in Luna's `Capillary.jl`: cached `neff_wg`/`neff` methods for
  `loss=:gas` (which had none, so `make_linop` was a hard MethodError), a
  per-z cache in `gas_mixture`'s `coren` (it rebuilt the whole Sellmeier
  closure per (ω,z) pair), and clamping the pressure splines' z into [0, L]
  (they extrapolate cubically, and the adaptive stepper probes past the fibre
  end, producing negative partial pressures that make CoolProp throw).
  **Everything previously concluded about localized-vs-uniform ozone is void.**
  With the fix the real localized profile moves the RDW 250 → 312 nm; it had
  appeared not to move it at all. Costs ~4x per solve.
- `RunMonitor` now also writes the LAST pulse's full `lunaoutput` (every z,
  not just `zslices`) to `<path>_last.jld2` on the same throttled cadence as
  its lightweight history, overwriting each time —
  `Monitor.load_propagation(path)` reconstructs a real, working
  `Luna.Output.MemoryOutput` from it, from any process, while the run is
  still going. Before this there was no way to get a true
  `Luna.Plotting.prop_2D` for a long run except from a still-alive process's
  in-memory state; the lightweight history was never enough (no full field).
  `examples/plot_last_propagation.jl` is the standalone script for it.
  **A `he_o2.jl` run started before this landed won't have `_last.jld2`** —
  only runs started after this commit do.
- `Monitor` (the live-run plotter) recorded the O₃(z) profile *after* a
  pulse's chemistry was applied, but the spectrum it paired that with was
  computed *before* — a one-pulse-stale mismatch. Given the sensitivity
  above, this alone could produce a misleading monitor plot. Fixed; **not yet
  re-validated with a fresh long run** (the fix landed after the last
  `he_o2.jl` run that produced `examples/he_o2_monitor.jld2`).
- O₃ no longer gets a `PlasmaCumtrapz` (ionisation/free-electron) response.
  Modelling choice: ionised ozone is assumed to dissociate before a
  free-electron response would matter, so what would be "ionisation" is
  costed entirely as dissociation (`DissCumtrapz` — energy loss only, no
  free-electron phase) at ionisation's own energy. `PhysData.jl`'s
  `:O3_diss` is deliberately set equal to `:O3`'s real ionisation potential
  (12.53 eV) to match — **this is intentionally different from the paper's
  own 9.5 eV ADK-fit value**, don't "fix" it back without checking here first.
  O₂ and N₂ are unaffected (still get both a real ionisation response and
  their own separate dissociation channel).
- `DissCumtrapz` (dissociation response, used for O₂/N₂/O₃) no longer carries
  a spurious free-electron phase term copied from `PlasmaCumtrapz` —
  dissociation produces neutral fragments, no plasma current.
- The three `DissCumtrapz` responses were structurally unreachable
  (appended past the end of the density vector `NonlinearRHS.Et_to_Pt!`
  zips responses against) — dissociation drew *zero* energy from the field
  until this was fixed. Now nested inside each gas's own response tuple.
- `neff_gasonly` (the `loss=:gas` waveguide-index path `PropAir.jl` needs)
  was dropping the Marcatili waveguide dispersion term entirely, not just
  its wall-loss part — this silently killed self-compression/RDW generation
  for every run using it (i.e. every LupoAirOsc run). Fixed.
- Two different, incompatible species orderings existed in LupoAirOsc
  (`ReactionState`'s field order vs. the ODE/PDE unknown-vector order) and
  were being hand-indexed inconsistently. `src/Species.jl` is now the single
  source of truth (`Species.SOLVE`, `Species.STATE`, `Species.IDX`); index by
  name, never by integer literal.
- `apply_lunaion`'s O₂/O₃ photolysis bookkeeping had a sign/ordering bug
  destroying O atoms every pulse (O₃ + hν → O₂ + O releases an O₂; the code
  subtracted it instead, from an already-decremented O₃).
- Luna's `Dissfrac_*` stat is the fraction dissociated **this pulse**, not a
  running total — the old delta-vs-previous-pulse pattern went to zero once
  the field settled, silently switching the O source off after ~2 pulses.
  New `diss_mode=:per_pulse` fixes this. **Default is still `:delta`**
  (backward compat) — `examples/he_o2.jl` already opts into `:per_pulse`;
  deciding whether to flip the global default (or delete `:delta` entirely)
  is still open.
- Generalised the chemistry (`State.jl`/`ReactionDiffusion.jl`/`Rates.jl`)
  for an inert bath gas (He) alongside air's reactive N₂ —
  `GasFill`/`air_fill`/`he_o2_fill`. `RateParams` takes an explicit
  third-body density + `M_efficiency`, `DiffusionParams` a `D_scale`.
  **Both are still placeholder 1.0** (N₂-equivalent) for He — real values
  are needed before He-O₂ numbers mean anything quantitatively; ozone
  production rate depends on `M_efficiency` directly.
- Added `UPPETrigger`: adaptive UPPE re-solve (re-use the previous Luna
  solve while composition drift stays under a tolerance, default 0.1%)
  replacing a hand-set `skip` pulse-count constant. Chemistry/diffusion
  still advance every 1 ms regardless. `examples/uppe_trigger_validation.jl`
  is the framework to demonstrate this converges to strict per-pulse UPPE
  as tolerance tightens — **written but not yet run to completion**
  (strict-UPPE reference runs are expensive; needs someone to actually run it
  and record the convergence numbers here).
- Added `Monitor.jl`: file-backed live plotting for long runs (JLD2 history +
  periodic PNG), so a run can be watched from another process without
  blocking. `examples/watch_run.jl` tails it.

**Concretely next, in rough priority order:**
1. There was a `he_o2.jl` run in progress (user-run, outside any assistant
   session) as of this writing, started **before** the Monitor-timing fix,
   the O₃-ionisation-removal change, and `save_full`/`load_propagation` all
   landed — its `examples/he_o2_monitor.jld2`/`.png` reflect the *old*
   physics and won't have a `_last.jld2`. Don't treat that run's numbers as
   validating the current code. Once it's done (or a fresh one is started),
   check whether the RDW-relocation behaviour now shows up in the monitor
   plot (or doesn't, which would itself be informative given how the
   standalone localized-profile test behaved — see
   `examples/o3_localization_check.jl`), and use
   `examples/plot_last_propagation.jl` to get a real `prop_2D` from it.
2. Get real `M_efficiency`/`D_scale` values for He (literature third-body
   efficiency and diffusivity), or at least bound how much they matter.
3. Run `uppe_trigger_validation.jl` and record the actual convergence
   numbers here.
4. Decide `diss_mode` default; consider deleting `:delta` if `:per_pulse` is
   adopted everywhere.
5. Only once He-O₂ is trusted quantitatively: move to air/N₂-O₂ and the
   oscillation-period question.

## Architecture

```
Pump pulse (800 nm, 30 fs)
  └─ PropAir.propair()          Luna UPPE solve — rebuilt from scratch each call
       └─ Kerr + Plasma (O2,N2) + DissCumtrapz (O2,N2,O3)
            └─ State.apply_lunaion()   deltas vs prev_*frac → concentrations
                 └─ ReactionDiffusion  21 reactions, VoronoiFVM 1-D along z
                      └─ new gas fractions → next propair() call, 1 ms later
```

The 1 kHz loop runs at full per-pulse resolution
(`ReactionDiffusion.jl`'s `tstops = range(..., step=1e-3)`); what's adaptive
is *whether* each callback triggers a fresh UPPE solve (`UPPETrigger`), not
the chemistry/diffusion timestep, which never coarsens.

## File map

| File | Role |
|---|---|
| `src/Species.jl` | Single source of truth for species ordering (`SOLVE` vs `STATE`) — index by name via `IDX`, never by literal. |
| `src/PropAir.jl` | Luna setup + run. `setup_propair` / `run_propair` split. |
| `src/State.jl` | `ReactionState`, `GasFill`/`air_fill`/`he_o2_fill`, `UPPETrigger`, `runnew` main loop, `effective_densities` (dispersion decomposition, below), back-coupling. |
| `src/ReactionDiffusion.jl` | VoronoiFVM system, 21-reaction `reaction!`, `flux!`, Rodas5. |
| `src/Rates.jl` | k1–k21 from JPL Publication 19-5, Troe falloff; diffusion coefficients. |
| `src/PhotoChem.jl` | `PhotoReactor` — O₃/NO₂/NO₃/N₂O photolysis from local spectrum. |
| `src/Monitor.jl` | File-backed live plotting for long runs. |
| `src/Plotting.jl` | `plotfullsol` etc., PyPlot only. |
| `src/Old.jl`, `src/SimpleReact.jl` | superseded |
| `examples/he_o2.jl` | The paper's He-O₂ validation case. |
| `examples/o3_dispersion_check.jl` | Two single propagations (0%/2% O₃, atom-conserving), confirms ozone dispersion actually perturbs propagation. |
| `examples/o3_localization_check.jl` | Diagnostic: real (localized) vs. uniform O₃ profile comparison, using a real run's saved history. |
| `examples/uppe_trigger_validation.jl` | Strict-per-pulse vs. adaptive-trigger convergence check (not yet run). |
| `examples/watch_run.jl` | Watches a `Monitor` history from another process. |

## Dispersion: how the 10-species chemistry maps onto tracked gases

`State.effective_densities` collapses the 10-species chemistry onto the
gases Luna actually has dispersion data for, **conserving atoms**:

```
N2_eff = N2 + 0.5*(N + NO + NO2 + NO3) + N2O
O2_eff = O2 + 0.5*(O + O1D) + 0.5*NO + 1.0*NO2 + 1.5*NO3 + 0.5*N2O
O3_eff = O3
```

(e.g. NO₃ is 1 N + 3 O → 0.5 N₂ + 1.5 O₂.) For an inert bath (He) there's no
N-chemistry; the bath density is just held constant. This is implemented, not
proposed — if you're tempted to add an `:N2O`-as-dispersion-proxy shortcut,
don't; that was a real bug here before (fed NO₂ concentrations to Luna under
the N₂O key) and this is what replaced it.

## Conventions

- Concentrations in **cm⁻³**, time in **s**, inside the chemistry modules;
  Luna/PhysData work in **m⁻³**. Conversion sites are the usual bug source —
  check units explicitly when a density looks off by 10⁶.
- `clamp!(up, 0, 2*rs.mtot)` in `State.jl` is a crude physical bound (a
  stiffness signal, not a fix) — must key off *total* density, not N₂
  specifically, since `mN2 == 0` for a He-O₂ fill.
- Don't "fix" `apply_lunaion`'s delta-vs-`prev_*frac` pattern for O₂/N₂ — it
  correctly avoids double-counting Luna's *cumulative* dissociation output
  for those two. (This is *not* the same bug as the O₃ `diss_mode` issue
  above, which is about `Dissfrac_*` not being cumulative *across pulses* —
  different axis, don't conflate them.)
- When in doubt about whether a physical quantity matches the paper, extract
  the actual paper PDF text and check — several bugs this session were
  caught exactly that way (dissociation energies, cross-section values,
  fibre geometry).

## Keeping this file current

**If you make a change here, update this file before you finish the
session** — not just the git commit message. Future sessions (any model,
not just Claude) read this file first and treat it as ground truth for
"what's already known/fixed/open." A commit history is a record of what
happened; this file is supposed to be a record of what's *true right now* —
keep "Current stage" honest: move finished items out of "next", add new
findings (especially surprising or counter-intuitive ones, like the RDW
relocation above), and flag anything you leave half-done as clearly open
rather than letting it silently disappear. If you're not sure something's
worth recording, err toward recording it — the cost of a stale line here is
much lower than the cost of the next session re-discovering the same bug.
