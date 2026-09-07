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
  Re-checked end-to-end on a live run including the boundary: over 0 → 0.332 s
  the total O budget (`O + O¹D + 3·O₃ + NO + 2·NO₂ + 3·NO₃ + N₂O + 2·O₂`,
  summed over all 201 z-nodes) drifts by **+2.0e-6 relative**, i.e. +0.002% of
  the ozone made. The residual is a *gain*, not a loss, and is the O₂ Dirichlet
  reservoir refilling depletion near the ends — not atoms leaking out. Ozone's
  own boundary loss is ~5e-5 of the ozone made (end-node O₃ is ~2e-14 cm⁻³
  against a 1.8e19 peak). **Whatever is wrong with the rates, it is not
  atom bookkeeping.**
- **The ozone → dispersion → RDW feedback loop is live and observable** — for
  the first time, since the const-linop bug below had severed it. A run on the
  fixed code moves the RDW peak 251.1 nm → 306.9 nm within 6 ms and then creeps
  to 309.3 nm by 68 ms, landing where the standalone localized-profile test
  predicted (312 nm). An oscillation *requires* this loop, so no run predating
  the fix could have produced one even in principle.

**THE current blocker — the chemistry runs ~200x too fast.** The paper reaches
~3% O₃ after **5 s**; the current code reaches 3.36% after **13 ms**, and
settles toward ~6.5% against the paper's ~3%. Measured from a live run's saved
stats, per pulse at the compression point (z ≈ 10.5 cm):

| stat | peak value |
|---|---|
| `Dissfrac_O2` | **0.278 % per pulse** |
| `Dissfrac_O3` | 22.1 % per pulse |
| `ionfrac_O2` | 2.89 % per pulse |
| peak intensity | 6.29e13 W/cm² |

The paper's ozone curve implies a per-pulse O₂ dissociation of ~1.4e-5, so the
O-atom source is **~200x too large**. That single number accounts for
essentially the whole discrepancy. Note also that O₂ *ionisation* is 10x larger
than O₂ dissociation but feeds nothing chemically (plasma only), and that O₃
photolysis at 22%/pulse is **not** a net odd-oxygen sink — it cycles O₃ → O and
the chemistry converts it straight back.

Two candidate knobs, deliberately *not* both to be tuned (they are degenerate
against this one observable, so fitting both leaves you unable to say which was
wrong):
- **`:O2_diss` barrier, currently 15.5 eV** (superexcited state, Song et al. —
  NOT the 5.12 eV ground-state bond energy; that was tried and is wrong, see
  the comment in `PhysData.jl`). ADK is exponential in Ip: measured at the run's
  actual peak field (21.8 GV/m), 17 eV gives 0.049x, 18 eV gives 0.0061x, and
  **~18.4 eV supplies the whole 200x**. Three eV of barrier covers it.
- **`diss_yield` (φ), added this session**, default 1.0 so nothing changes
  unless passed. The branching ratio of the strong-field channel: the model as
  written assumes every superexcited molecule dissociates into two O atoms and
  both survive to make ozone. Structural expectation, **not yet measured**:
  φ scales the source linearly but odd oxygen Ox = O + O₃ is destroyed by
  reactions quadratic in it, so equilibrium O₃ ~ √φ and time-to-equilibrium
  ~ 1/√φ. If that holds, the 200x in the source buys only ~14x in ozone level,
  landing at ~0.46% (below the paper's 3%) while stretching the timescale to
  only ~0.7 s (short of 5 s) — i.e. **φ alone could not match both level and
  timescale**, which would point back at the barrier energy.
  `examples/diss_yield_scan.jl` measures the actual exponent (frozen field,
  `uppe_tol=Inf`, so it isolates the chemistry's response to φ);
  **written, handed to the user to run, result not yet recorded here.**

**The other blocker — cost.** At last measurement the UPPE trigger was
re-solving on essentially every pulse (58 solves / 57 pulses = 102%, drift
0.0057 against `uppe_tol=5e-3`), giving **3.84 hours of wall clock per
simulated second**. The paper's 5 s target is ~19 h; a single ~100 s
oscillation period would be ~16 days. The drift metric itself is well designed
(relative χ shift at 255/300/800 nm — the right quantity); it fires constantly
because O₃ genuinely changes that fast. **Fix the chemistry pace first and
re-measure — do not tune `uppe_tol` blind.** (Minor: `examples/he_o2.jl`
sets `uppe_tol = 5e-3` but its comment says 0.1%; it is 0.5%.)

**Oscillation mechanism hunt (started 2026-09-07) — what the data say, what
the model says, and where it stands.** Sources: the paper draft
(`~/Downloads/Ozone.pdf`), the 22.5 cm and 27 cm pressure/energy grids, the
thesis Ch. 7 (Dropbox `PhD_Thesis/reviews/MSabbah_PhD_Thesis_John_ozone.pdf`,
text dump in the session scratchpad), and the side-scattering images in
Dropbox `PhD/Ozone paper/images/side/`.
- **Geometry: the air/N₂-O₂ experiments used a 22.5 cm fibre** (thesis 7.1;
  the 27 cm grid is a second length). `PropAir`'s 15 cm default is NOT the
  experiment. 6 bar, 2.2 µJ, 79/21, dry synthetic N₂/O₂, sealed static cell.
- **What "disappearance" is.** Thesis p.101, confirmed by the side-scattering
  images (UV at 230–295 nm visible from the side while the output is dark):
  the RDW is *still generated*; it is *absorbed by ozone downstream* of the
  compression point (z_c ≈ 7.3 cm, so ~15 cm of ozone to cross). Luna agrees:
  at 22.5 cm, **0.03% uniform O₃ already kills the 268 nm band** — pure
  Beer–Lambert, σ(255 nm)·n·L ≈ 2.6e4 per unit fraction. The observable is
  therefore **transmission through the downstream ozone column**, not a
  phase-matching switch. (At 7–8 bar the RDW sits at 300–350 nm, on the
  Hartley tail, so it never vanishes — it shifts and then drifts *back*.)
- **The slow, N₂-only process is ozone *decaying*.** In every air panel the
  band returns to its original wavelength after ~50–130 s; in He-O₂ it never
  does (700 s). Ozone in the dark at 300 K has no volumetric sink except
  chemistry, and diffusion to the ends is ~4000 s. So something present only
  with N₂ removes ozone from the downstream column on ~100 s. Oscillation
  (period ~45 s at 5 bar, ~100 s at 6 bar for 22.5 cm; ~100–200 s for 27 cm)
  appears only in a *window* of pulse energy (2.2 µJ at 22.5 cm; 1.6/2.2 µJ
  at 27 cm); outside it the band returns and stays. Thesis: "no obvious
  relationship" with parameters — treat the period scalings as soft.
- **Reduced 0-D model (`examples/oscillation/reduced_nox_model.jl`) with the
  code's own network finds NO oscillation and NO slow ozone decay** in any
  regime tried (φ 1→0.005, τ_diff 3.8 s→∞, N source ×1→every N₂⁺). NOx is
  capped at ~2e14 cm⁻³ (~1 ppm) because N + NO₂ → N₂O + O (k17) beats N + O₂
  → NO + O (k14, ground-state 8.5e-17) once NO₂ reaches ~1e14, and N₂O is a
  dead end. That is 4–5 orders too weak to touch ozone. **Missing physics
  identified (none implemented in the full model yet):** (1) plasma
  recombination atom sources — O₂⁺ + e → O + O and N₂⁺ + e → N + N; the
  code discards all ionisation products as "plasma only", yet ionfrac_N2 at
  z_c is 1.0%/pulse, 3000× the ADK N₂ dissociation; (2) N(²D): a large
  fraction of those N atoms is excited, and N(²D) + O₂ → NO + O at ~5e-12
  bypasses the N₂O drain entirely (a `k14mult` knob now exists in the reduced
  model as a proxy); (3) N + N + M / N + O + M recombination, without which
  the ion source is absurd (N₂O reached 18% of the gas); (4) NO₃ + NO₂ + M →
  N₂O₅ (NOx sequestration) and its thermal/UV release.
- **Working hypothesis (not yet demonstrated): NOx titration of the
  downstream ozone column.** Filament literature (Petit 2010, Camino 2015)
  has NOx at tens of ppm alongside O₃ at hundreds of ppm — i.e. NOx
  *comparable to the ~100 ppm column that blacks out the RDW*. NO + O₃ →
  NO₂ + O₂ then NO₂ + O₃ → NO₃ needs no O atoms, so downstream of the zone
  each N removes ~2 O₃ stoichiometrically. Accumulating NOx over ~100 s
  titrates the column → band returns (He-O₂: no titrant → never returns).
  Oscillation would be the relaxation limit of that titration when NOx
  production ≈ ozone production (the energy window). Untested.
- **Side-scattering, calibrated spectrometer (user confirmed), lines never
  identified:** in `images/side/side.png` a very bright narrow line at 232 nm
  and a weaker one at 245 nm switch on at t ≈ 35 s and persist. These sit on
  the **NO γ-band v′=1 heads (1,2) 233.0 and (1,3) 244.0 nm** to ~1 nm; the
  strongest γ bands, (0,0) 226.9 and (1,0) 214.9, terminate on v″=0 and are
  self-absorbed by ground-state NO, which would explain their absence. If
  this assignment holds it is *direct* evidence of NO, appearing on a ~35 s
  timescale. Check: weak (1,4) 255.4 and (0,3) 258.7 nm (a faint line near
  256 nm is visible). The alternative reading is non-guided RDW light
  (λ < 240 nm scatters at the fibre resonance). Unresolved.
- `examples/oscillation/nox_titration_model.jl`: 3-box (upstream / zone /
  downstream) model with a Beer–Lambert readout and the four missing pieces
  of physics above as parameters (φ, η, f2D). Built to test the titration
  hypothesis. **First scan (φ=0.14, 300 s):** with the code's own network
  (η=0, f2D=0) the downstream column only grows — output −130…−160 dB, band
  never returns (He-O₂-like). With half the N atoms as N(²D) (η=0, f2D=0.5)
  NOx titrates the downstream column: output **−51 dB @50 s → −17 dB @150 s
  → −13 dB @300 s** — the experiment's dark-then-return on the right
  timescale. Net NOx supply = (2·f2D−1)·N-source, a knife-edge between
  N(²D)+O₂→NO and N(⁴S)+NO→N₂+O. Larger sources (f2D=1, or any η ≥ 1e-3)
  return in 3–30 s and NO runs away to 10–30% of the gas: **once ozone is
  gone this network has no NOx sink**, so the titrant is never consumed and
  nothing can oscillate. **Missing: the wall.** In a 30 µm bore the
  diffusion-limited loss is 5.78 D/r² ≈ 5e4 /s; NO₃/N₂O₅ (γ~1e-3 on silica)
  die on the wall in ~20 µs, making NO + O₃ → NO₂ → (+O₃) NO₃ → wall a
  *stoichiometric*, self-consuming titration — the closing step a
  relaxation cycle needs — and ~95% of N(⁴S) recombines on the wall
  (dissolving the knife-edge). γ(O₃) ≲ 1e-9 is bounded by He-O₂ keeping its
  ozone for 700 s. Wall terms are now a parameter (`wall`) of the model;
  **Wall scan (φ=0.14, η=0, 300 s; k_wall O 2.1e4, N 2.2e4, NO₃ 1.1e4,
  N₂O₅ 8.1e3, NO₂ 12 s⁻¹):** f2D=0.5 → **band returns at 102 s**
  (−40.6 dB @50 s, −10.6 @100 s, −3.8 @150 s, −2.1 @300 s) with NOx bounded
  at ppm levels (NO 3e14, NO₂ 1e15, N₂O₅ ~1e10) — the first physically
  sensible reproduction of the dark-then-return, and absent entirely with the
  code's own network. f2D=0.75 / 1 return at 48 / 29 s but NO runs away to
  0.3–4% (NO has no sink). **No oscillation in any run** — after the return
  the band stays on (matches the non-oscillating grid panels, not the 2.2 µJ
  ones). Structural reason: one slow variable (NOx) driving a fast one (O₃)
  cannot cycle. Candidates for the second slow variable / bistability, in
  order: (i) N₂O — already climbing past 1% in every run, and O(¹D)+N₂O→2NO
  is UV-gated, so its release switches with the RDW; (ii) the RDW wavelength
  crossing the Hartley band as ozone changes (the Beer–Lambert readout uses
  σ(255 nm) only, and box B is attenuated across the whole spectrum instead
  of only in the Hartley band — a known crudeness); (iii) spatial: in the
  no-wall run the downstream column collapsed abruptly at ~145 s when the
  NOx front arrived — a propagating front the 3-box model cannot resolve.
  Figures: `examples/oscillation/nox_*.png`. Each 300 s run costs ~20 min.
- **The 3-box model's "source only in the zone" was too kind.** The first 1-D
  frozen-optics run (`air_1d_frozen.jl`, 2 s) shows ozone made along the
  *whole* fibre — recompression peaks at 7, 10, 14 and 19 cm, and ~1e15 cm⁻³
  even at the entrance within 0.3 s — so the downstream column is 1e17–1e18
  everywhere by 2 s (output −1500 dB). From the frozen solve's stats along z:
  the ADK N₂ channel (21 eV) is 1e-5–1e-2 of the O source everywhere — never
  a titrant. The *ionisation* channel (PPT) gives N/O ≈ 0.5 at the peak and
  ≈ 0.08 off-peak for any η (the ratio ionfrac_N2·0.79 / ionfrac_O2·0.21 is
  η-independent), and ionfrac_O2 off-peak is ~900× Dissfrac_O2 (8e-4/pulse
  at z = 1 cm), so with η > 0 both ozone and NOx are made along the whole
  fibre. **η (recombination-to-atoms yield) is the lever, and it acts on
  both species.** Physics note: at 6 bar three-body attachment e + O₂ + M
  (~1e10 s⁻¹) competes with dissociative recombination, so η is genuinely
  uncertain, likely 0.1–0.5.
- **Where a limit cycle could come from, now that the sinks are in:** with
  ozone high, NO₂ is drained by NO₂ + O₃ → NO₃ → wall at k₂₀·[O₃] ≈ 30 s⁻¹, so
  NOx cannot accumulate; as ozone falls NO₂ rises ∝ 1/[O₃] and the catalytic
  O + NO₂ → NO + O₂ / NO + O₃ cycle strengthens — positive feedback, i.e. a
  switch, matching the abrupt "return". In the ozone-free state NOx then
  dies only by NO₂ wall loss (γ = 1e-6 → 12 s⁻¹, 80 ms — too fast; γ ~ 1e-8
  → ~10 s) after which ozone rebuilds. **The period would be set by the
  NO₂/NO wall uptake coefficients**, the least-known numbers in the model —
  a testable, physically meaningful outcome if the 1-D run bears it out.
- Timing of the 1-D frozen-optics run: nbundle=10, 201 nodes, reltol 1e-5 /
  abstol 1e-12 → 180 s wall per simulated second (1.8 s per 10 ms kick).
  Tolerances and node count are now exposed (`rd_reltol`, `rd_abstol`,
  `npoints`); background commands must use absolute paths (the harness resets
  the shell cwd to the parent directory).
- Also found and fixed on the way: O(¹D)+O₃ and O(¹D)+N₂O branching double
  counted in `Rates.jl` (each channel carried the total rate).

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
- **1-D model extended for the NOx mechanism (2026-09-07), all opt-in — defaults
  reproduce the old behaviour exactly.** `Species`: 11th species `:N2O5`
  (appended to `SOLVE`; in `STATE` after `:N2O`), so `NSPEC` is now 11 and
  `ReactionState` has an `N2O5` field. `Rates`: k22 N+N+M, k23 N+O+M, k24
  NO₂+NO₃+M→N₂O₅ (Troe, JPL 19-5), k25 N₂O₅+M→NO₂+NO₃ (= k24/K_eq, 0.05 s⁻¹);
  `D_N2O5`; new `WallParams(a, dp; scale, γ)` — per-species first-order wall
  loss `min(γ·v̄/2r, 5.78·D/r²)` in `SOLVE` order, defaults γ = 1e-3 for O, N,
  O(¹D), NO₃, N₂O₅ and 1e-6 for NO₂, **0 for O₃ and NO** (He-O₂ keeps its
  ozone for 700 s ⇒ γ(O₃) ≲ 1e-9). `ReactionDiffusion.reaction!` carries the
  four reactions and the wall term (atoms return as O₂/N₂, NOy sticks);
  `ReactionDiffusionSolver(...; wall=nothing)`. `State.apply_lunaion(...;
  ion_yield=0, f2D=0, npulses=1)`: O₂⁺+e→O+O and N₂⁺+e→N+N from Luna's
  `ionfrac_*` stats scaled by `ion_yield`; a fraction `f2D` of fresh N atoms
  is applied as N(²D)+O₂→NO+O instantaneously. `apply_photochem(...;
  attenuate=false, npulses=1)`: Beer–Lambert attenuation of the spectrum at
  each z by the O₃ column upstream of it, using the O₃ cross-section on the
  ω grid (new `PhotoReactor.σO3` field). `runnew(...; ion_yield, f2D,
  attenuate, nbundle=1, rd_dtmax=1e-4, rd_dt_after=1e-14)`: `nbundle` applies
  n pulses as one kick (fractions 1−(1−f)ⁿ, photon dose ×n, tstops every
  n ms — the paper's own "5 chemistry steps per solve" done honestly); the
  two `rd_*` knobs expose the RD solver's `dtmax` and its post-kick restart
  dt, which were hard-coded to 1e-4 and 1e-14 and cost ~1 s wall per pulse
  (2211 unknowns). `setup_run(...; wall=nothing)`. The RD runner now takes
  `N2O5` positionally after `N2O`; its two per-pulse `println`s are behind
  `verbose`. `examples/oscillation/air_1d_frozen.jl` is the frozen-optics
  driver with the transmission readout (T = exp(−σ₂₅₅∫_{z_c}^{L} n_O₃ dz)
  from the Monitor O₃(z) history). **Runs from a script must set
  `ENV["MPLBACKEND"]="Agg"` before loading the package** — `Monitor.render`
  from inside the solver callback crashes the Tk backend.
- **Bug: the N₂ source was silently off in `:per_pulse` mode.** `apply_lunaion`
  still used the delta-vs-`prev_N2frac` form for N₂ after `:per_pulse` was
  introduced for O₂/O₃, so with a reused field (frozen optics, or any pulse
  between UPPE re-solves) N₂ dissociated only on the first pulse. Now follows
  `diss_mode` like O₂. (The "Conventions" note below about not touching the
  delta pattern predates `diss_mode` and applies only to `:delta` mode.)
- Added **`diss_yield`** (φ) to `apply_lunaion`/`runnew` — see the blocker
  section above. Default 1.0, so no behaviour change unless passed. Applies to
  the O₂/N₂ strong-field channel only, NOT to O₃ photolysis (whose own implicit
  yield of 1.0 is a separate, untouched assumption). Committed together with
  `examples/diss_yield_scan.jl`.
- **Corrected the boundary-condition comment** in `ReactionDiffusion.jl`, which
  claimed "zero-flux (Neumann) boundaries, meaning no species can leave or
  enter". The live code is the exact opposite: **Dirichlet at both ends** — O₂
  and N₂ held at bulk fill (an infinite reservoir), every reactive species
  including O₃ held at **zero** (a perfect sink). The stale comment described a
  commented-out block the Dirichlet loop had replaced. Harmless for He-O₂
  (D_O3 ≈ 0.014 cm²/s at 12 bar → ~3.8 mm diffusion length over 5 s, against an
  ozone peak ~10 cm from either end; half-fibre diffusion time ~4000 s), but
  **not safe to assume for the air case**, where ~100 s is far closer to the
  diffusion timescale.
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
  appeared not to move it at all.
- **Made that affordable.** The naive switch cost ~3x per solve, because every
  frequency point at every z step re-entered the Sellmeier stack. Since
  `n(λ,z) = sqrt(1 + Σᵢ γᵢ(λ)·ρᵢ(z))` has all its wavelength dependence in the
  z-independent `γᵢ`, `gas_mixture` now returns a tabulatable `GasMixtureIndex`
  struct, `neff_β_grid` precomputes `γᵢ` on the frequency grid, and `ρᵢ(z)`
  comes from a `densityspline` rather than CoolProp (whose root-find was 60% of
  the whole linop at ~76 µs a call). Net: ~12 s per propagation vs ~9 s for the
  old incorrect const linop — correctness costs ~35%, not ~200%. Tabulated
  `neff` matches the generic path to 2.2e-16.
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
1. **Record the `diss_yield_scan.jl` result here** — the fitted exponent
   (peak O₃ ~ φ^α) and whether any φ reaches the paper's 3% at ~5 s. That
   decides between the two knobs above. If α ≈ 0.5, φ is the wrong lever and
   the barrier energy is the one to move.
2. **Get the chemistry pace right** (one knob, chosen on physical grounds),
   then re-run `he_o2.jl` and check the RDW red-shift is *monotonic* on the
   way to equilibrium, as the paper reports. Damped ringing would mean the
   feedback gain is still too high.
3. **Re-measure the UPPE resolve rate** after (2), not before. If it is still
   ~100%, that is the moment to think about `uppe_tol` or the drift metric.
   Any long run is impractical until this is fixed.
4. Get real `M_efficiency`/`D_scale` values for He (literature third-body
   efficiency ~0.6-0.7 and diffusivity). Only ~1.5x, so it does not fix (2),
   but it should be set correctly before quantitative comparison.
5. Run `uppe_trigger_validation.jl` and record the actual convergence
   numbers here.
6. Decide `diss_mode` default; consider deleting `:delta` if `:per_pulse` is
   adopted everywhere.
7. Only once He-O₂ is trusted quantitatively: write `examples/air.jl`
   (6 bar, 2.2 µJ, 79/21 N₂-O₂ — `air_fill` exists, no script does yet) and
   go after the oscillation. **It is an air-only phenomenon**: He-O₂ has no
   NOx branch and the paper reports it monotonic, so it cannot be reproduced
   in He-O₂ by construction. Revisit the Dirichlet boundary condition there.

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
