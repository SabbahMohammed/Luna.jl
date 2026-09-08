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
  The RD tolerances are still **hard-coded** in the runner (`Rodas5`, reltol
  1e-5, abstol 1e-12, `isoutofdomain` rejecting u < −1e-15); only `dtmax`
  and `dt_after` are exposed (`rd_dtmax`, `rd_dt_after`), and the node count
  via `setup_run(npoints=)`. Background commands must use absolute paths (the
  harness resets the shell cwd to the parent directory).
- **Result of the first 100 s frozen-optics 1-D runs (2026-09-07, 6 bar air,
  2.2 µJ, 22.5 cm, φ=0.14, f2D=0.5, wall×1, nbundle=20, 201 nodes, the
  runner's reltol 1e-5/abstol 1e-12, dtmax 1 ms, dt_after 1e-8; 4.3 h wall
  each = 150–155 s/s):**
  `examples/oscillation/air1d_T100_f2D0.5_eta{0.03,0.1}_wall1.0_phi0.14_nb20*`.
  Neither run ever transmits: the readout (z_c = 7.3 cm) is −942 dB (η=0.03)
  and −372 dB (η=0.1) at 100 s. What was learned:
  (i) **The rise is ~100× too fast.** The downstream column passes τ=1
  (1e17 cm⁻²) before the first save at 0.4 s and sits at 2–6e19 cm⁻² by 1 s;
  the experiment keeps the band visible for ~50 s. Every later timescale
  inherits this, so the absolute O-atom yield (φ, and η's O₂⁺ share) must be
  calibrated before any period means anything — the He-O₂ shift-rate data is
  the right constraint (no NOx there).
  (ii) **Titration works, ∝ η, but only where the pulse is weak.** At η=0.1 the
  15–21 cm lobe falls 19× (5.7e19 → 3e18 cm⁻²) in 100 s; at η=0.03 only 3×.
  But the entrance lobe (0–5 cm, 7.5e18 cm⁻³) and the last centimetre before
  the exit (final recompression, 5e18 cm⁻³ = 4.6e18 cm⁻², −229 dB on its own)
  are untouched. Reason: ADK dissociation (15.5 eV) is steeper in intensity
  than PPT ionisation, so at the intensity peaks O/N ≫ 1 and ozone wins; in
  the weaker lobes the ion channel gives N/O ≈ 0.1–0.5 and NOx wins. So with
  frozen optics the band can never return — the exit spike alone blocks it —
  and its fate depends on the optics *moving* (entrance ozone changes the
  dispersion, which relocates/weakens the downstream recompression). The next
  run must be unfrozen (`uppe_tol` ~10 % drift).
  (iii) **NOx never accumulates**: NO₂ and NO peak within 1 s (Σz ≈ 3e18 and
  1e18 cm⁻³ at η=0.03) and decay; N₂O keeps growing and overtakes total O₃ by
  ~20 s — N+NO₂→N₂O is still the dominant nitrogen sink, exactly as in the 0-D
  model. No switch, just steady titration by a small standing NOx level.
  Julia buffers stdout to a file, so a background run's log lags the run by
  many minutes — judge progress from the Monitor PNG mtime, not the log.
- **The propagation already absorbs the RDW — `attenuate` is a frozen-optics
  correction only.** Luna's ozone index is complex
  (`PhysData.ref_index_o3_analytical_reference`: two Gaussians in the imaginary
  part, Hartley 255 nm and Chappuis 600 nm), it reaches the linop through
  `sellmeier_gas(:O3)` → `γ_ozone_analytical` → the mixture index, and since the
  const-linop fix that index is rebuilt per z. Converting κ(λ) to a
  cross-section, σ = 2ωκ/c per unit density:

  | λ (nm) | 230 | 255 | 268 | 280 | 300 | 320 |
  |---|---|---|---|---|---|---|
  | σ (1e-18 cm²) | 4.95 | 11.7 | 8.59 | 4.07 | 0.438 | 0.0138 |

  which is the literature Hartley band (1.15e-17 cm² at the peak) to 2 %. So
  with live optics the RDW is absorbed inside the UPPE solve, and
  `apply_photochem(...; attenuate=true)` would double-count — it is correct
  ONLY for the frozen driver, whose single ozone-free solve would otherwise
  apply unabsorbed UV to the photochemistry.
- **Mistake in the frozen driver's readout (results above still stand):** it used
  σ(255 nm) = 1.15e-17 cm² for a band at 268–273 nm, where σ is 8.6e-18 — a 34 %
  overstatement of the absorption. Both runs were opaque by orders of magnitude,
  so no conclusion changes, but a wavelength-resolved readout is the right form.
  `examples/oscillation/air_1d_live.jl` avoids the issue entirely by reading
  `Monitor.uvenergy`, the UV energy that actually leaves the fibre.
- **The 1-D solver is ~6x faster (2026-09-07).** Profiling 10 kicks showed ~80 %
  of the time in VoronoiFVM's `eval_rhs!`, which re-assembles the full Jacobian
  with dual numbers on every call (0.89 ms for 11 × 201 unknowns) and then keeps
  only the residual; a Rosenbrock step makes ~9 such calls.
  `ReactionDiffusion.FastRHS` evaluates the same discretisation directly in
  ~4 µs (agreement 1e-17 relative, and it is now a package test). Jacobian,
  sparsity and mass matrix still come from VoronoiFVM. Wall time per 10 kicks:
  24.7 s baseline → 6.7 s (same tolerances, results to 3e-10) → 3.9 s at
  reltol 1e-4 / abstol 1e-10. A 100 s run is ~40 min instead of 4.3 h.
  New kwargs `rd_alg`, `rd_reltol`, `rd_abstol`, `rd_fast_rhs`; all defaults
  reproduce the old behaviour. **BDF integrators (FBDF, QNDF, KenCarp47) all die
  immediately with `DtLessThanMin`** on this penalty-Dirichlet mass-matrix form
  (the 1e30 penalty rows), so it is Rosenbrock only until the BCs change.
- **First live-optics 1-D run (3 s, φ=0.14, η=0.1, f2D=0.5, wall×1, tol=0.1,
  201 nodes; `air1dlive_T3_*`), and a new mechanism candidate.** Everything
  happens inside 0.1 s at this φ, so treat the timescale as meaningless until
  φ is calibrated — but the *spatial* behaviour is new and did not appear in any
  frozen run:
  (i) The ozone profile is no longer two lobes. It is near-uniform along the
  whole fibre (~2e18 cm⁻³) with **narrow, deep holes** where ozone falls below
  1 cm⁻³ over ~5 mm at z ≈ 7.5 cm, a second opening at z ≈ 11 cm by 2.6 s.
  **These are at least partly the pulse-bundling bug found straight afterwards
  (see below) and must be re-run before being believed** — that bug zeroed
  downstream ozone every kick wherever the per-pulse photolysis fraction times
  nbundle exceeded 1, which is exactly the high-intensity region where the holes
  appear. A real burn-out effect should still exist (at the compression point
  O₃ photolysis plus O + O₃ → 2O₂ can outrun O + O₂ + M), but its depth and
  position here are not trustworthy.
  (ii) **The holes migrate downstream** as ozone accumulates upstream and changes
  the dispersion that sets where the pulse recompresses. A transparent hole is
  transparent at 255 nm, so this is a spatial relaxation-oscillator candidate
  that no 0-D or 3-box model can express and that frozen optics structurally
  cannot show: the compression point walks along the fibre, opening a
  transparent channel; when it reaches far enough downstream the RDW gets out
  again. This is the earlier candidate (iii) (a propagating front) appearing on
  its own. Subject to the same caveat as (i) — re-run before relying on it.
  (iii) The RDW **relocates rather than fading**: centroid 273 → 327 nm within
  0.12 s, output down 20 dB. Consistent with the standing "load-bearing finding"
  that ~2 % O₃ moves the band to 310–320 nm.
  Cost: 207 s wall per simulated second, worse than frozen, because ozone at
  φ=0.14 drives 38 UPPE re-solves in 3 s. Calibrating φ down should fix the cost
  and the physics together, since the re-solve rate is driven by ozone drift.
- **φ is NOT the reason the model is too fast — my earlier diagnosis was wrong.**
  `examples/oscillation/he_o2_phi_calibration.jl` scans φ against the paper's own
  He-O₂ number (Sec. III E / Fig. 6b, quoted in `examples/he_o2.jl`): 9.3e24 m⁻³
  ozone, ~3 % of the local density, after ~5 s at 12 bar, 2.5 µJ, 22 cm.
  Frozen optics, 5 s each, ~80 s wall per φ:

  | φ | 0.14 | 0.05 | 0.015 | 0.005 | 0.0015 | 0.0005 |
  |---|---|---|---|---|---|---|
  | peak O₃ (1e23 m⁻³) | 29.0 | 15.5 | 6.40 | 2.53 | 0.89 | 0.38 |
  | % of n | 0.98 | 0.53 | 0.22 | 0.086 | 0.030 | 0.013 |
  | t₉₀ (s) | 0.04 | 0.04 | 0.02 | 0.02 | 0.02 | 0.02 |

  Two results, both against my earlier claim. **φ = 0.14 UNDER-produces ozone by
  3.2×**, not over-produces: the level scales as φ^0.78 and hitting the target
  would need φ ≈ 0.53. And **t₉₀ is 0.02–0.04 s at every φ over a 280× range** —
  the equilibration time does not depend on the source at all, so no φ can buy
  the paper's 5 s rise. (t₉₀ is resolution-limited at one 20 ms kick, so the true
  value may be faster still.)
- **Why the model equilibrates in ~8 pulses instead of ~5000, and why the scan
  above cannot answer the calibration question anyway.** Measuring the
  single-pulse ozone photolysis fraction along the He-O₂ fibre: 2.7e-4 at the
  entrance, 5e-4 at 6.6 cm, 0.024 at 11 cm, then **0.114–0.122 from 13 cm to the
  exit**. A 30 µm core puts ~1e16 photons/cm² on a 1e-17 cm² cross-section, so
  the RDW photolyses ~12 % of the ozone it passes, per pulse — 1/e in 8 pulses.
  With the optics FROZEN the RDW is nailed inside the Hartley band forever, so
  ozone hits a photolysis-limited equilibrium almost immediately. That is an
  artefact of my scan design, not of the model: photolysis is not a net odd-oxygen
  sink (O₃ + hν → O + O₂, then O + O₂ + M → O₃ in ~100 ns), so what should set
  the real 5 s timescale is ozone pushing the RDW *out* of its own absorption
  band — σ falls 100× between 255 nm and 320 nm — which frozen optics forbids by
  construction. **Any φ calibration must therefore be run with live optics.**
- **Bug fixed: pulse bundling in `apply_photochem` scaled the dose instead of
  saturating.** It multiplied the spectral energy density by `nbundle`, but
  photolysis removes a fraction of what is present, so n pulses leave (1−f)ⁿ,
  not 1−nf. With f = 0.122 and nbundle = 20 the old form asked to destroy 244 %
  of the ozone: O₃ went negative and was zeroed by `runnew`'s clamp, every kick,
  over the downstream half of the fibre. Correct bundling leaves 0.88²⁰ = 7.8 %.
  `PhotoReactor` now takes `npulses` and applies 1−(1−f)ⁿ per species, splitting
  between channels by their single-pulse yields (exact — the quantum yields are
  intensity-independent), f capped at 1; NO₂, NO₃ and N₂O get the same
  treatment. **This invalidates the ozone profiles of every nbundle=20 run so
  far** — both 100 s frozen runs and the 3 s live run. Their NOx conclusions
  (titration ∝ η, no NOx accumulation, N + NO₂ → N₂O drain) are qualitative and
  probably survive; the ozone profiles and all timings do not. `nbundle=1` runs
  were never affected. Now covered by a package test.
- **Bug: `runnew` returned a post-photolysis transient, not the between-pulses
  state — every ozone number ever read off a returned `rs` is ~115x low.** The
  RD runner built `tstops = range(0, tmax, step=pulse_dt)` and fired the
  chemistry callback at every one of them, including the last, which is `tmax`
  itself whenever tmax is a multiple of pulse_dt (essentially always). The
  solver saves a point on BOTH sides of a callback, so the last saved point —
  the state written back into the caller's arrays — is the instant AFTER a
  photolysis kick and BEFORE the ~100 ns O + O₂ + M → O₃ that puts the ozone
  back (photolysis is not a net odd-oxygen sink). The saved pairs on He-O₂,
  φ=0.05, nbundle=20 make it plain:

  | t (s) | 0.02 | 0.10 | 0.18 | 0.20 |
  |---|---|---|---|---|
  | pre-kick O₃ (1e18 cm⁻³) | 1.325 | 6.127 | 10.152 | 11.035 |
  | post-kick O₃ (1e18 cm⁻³) | 0.0096 | 0.0481 | 0.0864 | 0.0960 |

  Identical for `fast_rhs` true and false, so it long predates the speed-up.
  Kicks within one pulse interval of `tmax` are now excluded, leaving that much
  field-free relaxation, and `save_end=true` is explicit. **The monitor
  histories were always right** (they snapshot pre-kick), which is why the air
  runs' O₃(z) plots looked sane while the calibration's `maximum(rs.O3)` did
  not. Note also what the table shows about nbundle=20: the ozone is crashing
  by 99 % and rebuilding every 20 ms, which is the bundling-validity problem
  below, not a transient worth modelling.
- **Bug in the same three lines: `condition` indexed `tstops[ni]` unguarded**, so
  any `tmax` that is not an exact multiple of `pulse_dt` threw `BoundsError`
  once the solver stepped past the final tstop. Now bounds-checked. (This is why
  no run had ever been given a non-commensurate tmax.)
- **`nbundle` must satisfy f·nbundle ≪ 1, where f is the per-pulse photolysis
  fraction (0.12 downstream).** The saturation fix makes the arithmetic right,
  but bundling 20 pulses still applies 1−0.88²⁰ = 92 % of the ozone's photolysis
  in a single instant and then lets it recover over 20 ms. The real system never
  loses 92 % of its ozone at once. So `nbundle` ≲ 2 for these cases, not 20 —
  which costs 10–20x, roughly offset by the 6.3x solver speed-up. Every
  nbundle=20 result so far is affected.
- **Live-optics φ calibration (8 s, tol=0.1, nbundle=20, He-O₂ 12 bar / 2.5 µJ /
  22 cm), reading the MONITOR — the returned-state bug above makes the printed
  "peak O3" column of these runs wrong by 115x:** at φ=0.05 the monitor's peak
  ozone plateaus at 15.9e24 m⁻³, already 1.7x ABOVE the paper's 9.3e24 target
  (so **φ is well below 0.05**, the opposite of the frozen scan's φ ≈ 0.53), and
  the rise is 50 % at 0.2 s, 90 % at 0.64 s, 99 % at 1.2 s against the paper's
  ~5 s. The RDW centroid jumps 255 → 312 nm within the first snapshot and then
  sits there, so the band does leave the Hartley region, but in ~0.1 s rather
  than over seconds. Only 13–16 UPPE re-solves in 8 s. Superseded by the
  nbundle=2 re-run.
- **Why a lower φ should fix the rise time as well as the level, and how to
  read the answer.** The production per pulse is 2φ·f_ADK·[O₂]; at the
  compression peak f_ADK ≈ 9e-3, so φ=0.05 makes 5.6e16 cm⁻³ per pulse against
  an equilibrium of ~1.6e19 — 280 pulses, 0.28 s. φ=0.005 gives 2800 pulses,
  2.8 s. So the rise time scales as (level)/φ and a single φ can in principle
  match BOTH halves of the paper's target (9.3e24 m⁻³ AND ~5 s). If no φ matches
  both, the sink side of the network is wrong, not the source — that is the
  discriminating test this calibration is for.
- Sanity check on the photolysis fraction, since everything above turns on it:
  100 nJ of RDW at 260 nm is 1.3e11 photons, over a 7.07e-6 cm² core that is
  1.9e16 cm⁻², times σ = 1e-17 cm² gives f = 0.19. The measured 0.12 corresponds
  to ~65 nJ, i.e. 2.6 % conversion of a 2.5 µJ pump. So the model's per-pulse
  photolysis is physically right; it is a genuine consequence of a 30 µm core,
  not a coding error.
- **φ ≈ 0.05 is the calibration (2026-09-07), on the LEVEL.** With all three
  bugs fixed and nbundle=2, live optics, 8 s, He-O₂ 12 bar / 2.5 µJ / 22 cm
  (~800 s wall each, 9–16 UPPE solves):

  | φ | 0.05 | 0.015 | 0.005 |
  |---|---|---|---|
  | peak O₃ (1e24 m⁻³) | 8.78 | 5.69 | 3.66 |
  | % of local density | 2.97 | 1.93 | 1.24 |
  | fraction of the 9.3e24 target | 0.945 | 0.612 | 0.394 |
  | t₉₀ (s) | 0.25 | 0.45 | 0.71 |

  **φ = 0.05 hits the paper's ozone dead on**: 8.78e24 against 9.3e24, and
  2.97 % of the local density against the paper's "about 3 %". Use φ = 0.05.
  Scalings: level ∝ φ^0.38, t₉₀ ∝ φ^-0.45 — note the SIGN, lower φ is slower,
  the opposite of what I predicted from the production-rate argument above (that
  argument ignored that the equilibrium level falls with φ too).
- **Open question for the user: is the paper's "~5 s" a rise time or a
  measurement checkpoint?** The model reaches the right level in 0.25 s and then
  plateaus. If ~5 s is when the paper *sampled* the ozone, there is no
  discrepancy at all and the calibration is complete. If it is the observed rise
  time, then no single φ can match both halves (t₉₀ = 5 s would need φ ≈ 5e-5,
  giving a level 15x below target), which would put the error in the SINK side
  of the network, not the source. `examples/he_o2.jl`'s own header calls it "the
  paper's 9.3e24 m⁻³ checkpoint" and sets `tevolve = 1`, which leans towards
  checkpoint — but this is worth one look at Sec. III E / Fig. 6b before
  building on it.
- **First trustworthy air runs (2026-09-08): the model relaxes monotonically and
  cannot oscillate.** 150 s, live optics, φ=0.05 (calibrated), nbundle=2, 201
  nodes, f2D=0.5, wall×1, tol=0.1, at η = 0.1 and 0.3. ~6.5 h wall each
  (150–162 s/s). First runs in which optics, ozone level and pulse bookkeeping
  are all correct together. Result: **the RDW relocates 273 → 330 nm within ~1 s
  and stays**, output −17 to −22 dB, for the whole 150 s. No return, no cycling.

  | | η = 0.1 | | η = 0.3 | |
  |---|---|---|---|---|
  | t (s) | front (cm) | exit plug (cm) | front (cm) | exit plug (cm) |
  | 13 | 5.96 | 9.56 | 5.29 | 9.56 |
  | 52 | 4.28 | 6.53 | 2.59 | 5.40 |
  | 108 | 3.04 | 5.51 | 2.02 | 4.39 |
  | 150 | 2.93 | 4.61 | 2.14 | 3.71 |

  ("front" = how far from the entrance ozone stays above 1e17 cm⁻³; "exit plug" =
  the same measured back from the exit. Boundary-pinned nodes excluded.)
  The steady state is an ozone plug at the entrance (0–3 cm) and another at the
  exit (last 4–5 cm) with a depleted middle. **Those two plugs are what keep the
  band relocated**, and neither clears.
- **The depletion front is a NOx DIFFUSION front, and it stalls.** It decelerates
  (0.015 cm/s at η=0.1 over 50–150 s, 0.003 cm/s at η=0.3, both slowing) and
  parks at 2–3 cm — which is just sqrt(D·t) for NO₂ at 6 bar (D ≈ 0.017 cm²/s,
  sqrt(0.017×150) = 1.6 cm). NOx is made only where the pulse has compressed, and
  it reaches upstream only by diffusion, so the front creeps as sqrt(t) and never
  arrives. Higher η parks it further upstream, sooner. This is monotone
  relaxation, not a limit cycle.
- **Methodological blocker found in the same runs: `uppe_tol` quantises the
  observable, and at 0.1 the optics are effectively frozen.** 59 of the 63 UPPE
  solves at η=0.1 happen in the first 13 s; only 4 more in the remaining 137 s.
  The composition-drift metric is dominated by the fast initial ozone build, and
  once ozone plateaus the drift sits at 0.02–0.13 while the spatial profile
  still reorganises by several cm. So the UV output is a piecewise-constant
  staircase with steps tens of seconds long — **it is structurally incapable of
  showing a ~100 s oscillation**. If the mechanism is optical (the compression
  point switching between two positions as ozone slowly changes — and the
  experiment's narrow 2.2 µJ window does suggest bistability rather than
  chemistry), these runs could not have revealed it. Re-running at tol=0.01,
  101 nodes. A drift metric sensitive to spatial redistribution, not just
  composition, is probably the right longer-term fix.
- **tol=0.01 is not viable**: 217 UPPE solves for 1.4 simulated seconds,
  extrapolating to ~77 h for 150 s. Killed. The drift metric is dominated by the
  fast initial ozone build, so tightening it globally just re-solves the opening
  transient to death. Needs a metric that is loose early and tight once ozone
  plateaus, or one keyed to spatial redistribution.

**Experimental constraints from the user's OneNote (2026-09-08), and what they do
to the model.** Two pages: "Probe with oscillation" and "Nox investigation from",
both 19 Sep 2025.

- **The model overproduces ozone by ~10x, and TWO independent measurements say
  so.** (i) The pump spectrum at 800 nm is constant through the oscillation, so
  the soliton dynamics are unchanged. But ozone's index contribution is only
  weakly dispersive — at the model's 1.3 % ozone, Δn at 800 nm is 3.86e-5, which
  is **14.6 % of (n_air − 1)** — so the model's ozone would rewrite the IR
  dynamics completely. Holding the IR index perturbation to a few percent caps
  ozone at 0.1–0.3 %. (ii) "No visible sign for absorption around 650 nm": the
  Chappuis band (σ = 3.67e-21 cm² at 650 nm, in Luna's own index) gives 0.4–1.0 dB
  at the model's column, which should be plainly visible. Below 0.1 dB caps the
  column at 6.3e18 cm⁻², a mean of 0.19 %. Both bounds agree; the model runs at
  0.78–1.3 %. NOTE: still to confirm whether the 650 nm page is transmission or
  side-scattering — the dB figures assume transmission through the full 22.5 cm.
- **Even at the bounded ozone the band still goes dark**: 0.15 % gives tens of
  optical depths at 268 nm. What it does NOT do is drag phase-matching to 330 nm,
  which is exactly the failure mode of the runs above.
- **Panel (a) of "Probe with oscillation" is a pump-PROBE absorption measurement,
  not the RDW** (probe at the same rep rate, ~100 ns after the pump; the void at
  ~400 nm is the anti-resonant fibre's non-guided resonance band). This is the
  most diagnostic data available, because O₃ and NO₂ look nothing alike there:

  | λ (nm) | 280 | 300 | 320 | 340 | 360 | 450 |
  |---|---|---|---|---|---|---|
  | O₃ (dB at model column) | 382 | 42 | 2.9 | 0.12 | 0.01 | 0.02 |
  | NO₂ (dB at model column) | 0.016 | 0.039 | 0.074 | 0.11 | 0.15 | 0.13 |

  Ozone prints a knife-edge (×300 between 300 and 340 nm; model edge at ~327 nm),
  NO₂ a flat broad hump that is invisible at model levels. So: the edge position
  gives the column; the optical depth at 315–325 nm is linear in it and should
  breathe if ozone is what oscillates; any probe light below 300 nm is a third
  independent bound. A broad 350–500 nm hump instead would mean NO₂ far above
  anything the network makes. Figure: `examples/oscillation/probe_prediction.png`.
- **The ~100 ns probe delay is a direct measurement of the per-pulse photolysis
  fraction.** O + O₂ + M → O₃ at 6 bar air runs at 2.87e6 s⁻¹, i.e. 1/e in 348 ns,
  so at 100 ns only **25 % has reformed** — the probe samples the gas with ~75 %
  of the photolysed odd oxygen still atomic. The ozone deficit it sees is
  therefore ≈ 0.75·f, and f = 0.12 predicts a ~9 % dip in Huggins-band absorption
  between probe-delay and between-pulse conditions. That is a direct test of the
  single number driving the whole model.
- **Why tripling η did nothing, and what the real knob is.** Ground-state atomic
  N is a NOx SINK here, not a source. Pseudo-first-order fates of one N atom at
  the model's own densities (NO 3.6e15, NO₂ 2.8e15, O₃ 1.2e18 cm⁻³):

  | channel | rate (s⁻¹) | share |
  |---|---|---|
  | N + NO → N₂ + O (destroys an NO) | 1.09e5 | 75 % |
  | N + NO₂ → N₂O + O (buries the N) | 3.35e4 | 23 % |
  | N + O₂ → NO + O | 2.68e3 | 2 % |
  | N + O₃ → NO + O₂ | 2.34e2 | — |

  So a ground-state N is ~37x more likely to destroy an NO than to make one. Net
  NOx per N atom: **+0.135 at f2D = 0.5, +0.481 at 0.7, +0.827 at 0.9, +1.0 at
  1.0.** f2D is therefore ~7x more powerful than η over that range, and η = 0.1 vs
  0.3 gave the same ozone column (2.65e19 vs 2.52e19 cm⁻²) precisely because the
  extra nitrogen cancelled itself. **f2D was fixed at 0.5 arbitrarily; literature
  for N₂⁺ + e⁻ dissociative recombination (Peterson 1998 branchings) supports
  ~0.6–0.65.** Running f2D = 0.7 and 0.9 at η = 0.3 next.
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
