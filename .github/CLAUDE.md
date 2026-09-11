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

**Where things stand (2026-09-11).**

- The campaign is at the start of **Phase 0** of `LupoAirOsc/docs/OSCILLATION_PLAN.md`: nothing in
  `targets/` beyond the two raw extracts, no register, no reduced-model skeleton yet.
- Calibrated so far: the `:O2_diss` ADK barrier, **21.0 eV**, from the early-time slope of the
  255 nm probe transmission (2 bar, 33 cm, 1.3 µJ in fibre); φ = 1 follows from it. γ(O) ≲ 1e-5 is
  a bound from the same curve's late-time curvature, obtained with the wall off and
  `o3_from_adk=false`. The barrier was fitted to a by-eye digitisation and is refitted in Phase 3.
- Not calibrated: `:O3_diss` (12.53 eV; the ~14–14.5 eV figure is an extrapolation and the fit is
  gated on Phase 1), `:N2_diss` (21 eV, unconstrained), every wall coefficient other than γ(O),
  and the baseline sink settings themselves (undecided; the register will hold them).
- **Nothing is reproduced.** Air does not oscillate in the model. He-O₂ was matched only to the
  paper's *simulated* ozone density, at 15.5 eV with φ = 0.05, and has not been run since
  2026-09-04; its measured target is the RDW wavelength trajectory (plan section 2).
- The code defaults are stale relative to all of the above: `Rates.WallParams` γ(O) = 1e-3,
  `air_1d_live.jl` wall on and φ = 0.14, `he_o2.jl` 22 cm instead of 22.5 cm.
- What the model does show: the ozone → dispersion → RDW loop is live; the 6 bar run makes
  1.88e18 cm⁻² of ozone in 10 s, peaked at the compression point, and the RDW falls only 2.2 dB;
  the driver spectrum is identical in both RDW states (0.00 dB), so the switch is absorption of an
  already-generated RDW.

### Log of this section, 2026-09-02 to 2026-09-08 — kept as a record, superseded in places

Everything from here to the next `##` heading was written before the 21.0 eV calibration and is
kept as the record of how the model got here. Where it disagrees with the summary above or with
the dated sections below, it is superseded. In particular: the 15.5 eV barrier and the "200× too
fast" blocker were resolved by the 21.0 eV calibration; φ = 0.05 was fitted model-against-model
and is superseded by φ = 1 at 21.0 eV; the He-O₂ ozone density target was the paper's simulation.

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
both 19 Sep 2025. **Read the "THE PROBE IS AN NO2 MONITOR" block below first — it
is the single most load-bearing fact in this file.**

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

---

## THE ANTI-PHASE DATASET, MEASURED FROM THE h5 (2026-09-08)

`double_spectrum_timeseries_61deg_200mbarAr_Hisol_4uJ_6barAr_pcf.h5`, read directly.
**These supersede the provisional numbers read off the plot.** Extracted traces are saved
at `LupoAirOsc/targets/raw/antiphase_61deg_200mbarAr.tsv` (25000 rows: t, RDW band, driver
band, probe band) so nobody re-reads 400 MB. Plots: `examples/antiphase61_windows.png`,
`examples/antiphase61_spectra.png`.

### The time base is in the file — it was never uncalibrated

`timestamps` row 1 is a **Unix epoch in seconds**; row 2 the same in µs. First sample
1559238010.169 = 2019-05-30 18:40:10, last 1559238967.098 = 18:56:07, which matches the
`params/timestamp` string exactly. Both spectrometers start within **16 ms** of each other.

> **956.93 s over 25000 slices → 0.0383 s per slice, 26.12 Hz.**

Not 40000 s, and not the ~600 s estimated. The notebook comment "this does not represent
the actual time" is wrong: it does. Use row 1.

### THE SPECTROMETER NAMING IS SWAPPED — identify channels, never trust the names

Confirmed by the user: **the acquisition software labels the two spectrometers the wrong
way round.** In this file the group called `hisol` holds the **pump**, and `pcf` holds the
**probe** — the opposite of the names. Treat it as systematic across this rig's
`double_spectrum_timeseries_*` files, and *still verify per file*.

Two independent tests, both cheap. Run the first; fall back to the second.

**1. Onset order (primary).** The probe is unblocked first and the pump only afterwards, so
the run captures the full dynamics from the pump's first shot. Sum each channel over
200-1100 nm per frame, find the first frame above `dark + 0.2 × (working level − dark)`:
**the channel that comes on FIRST is the probe.**

| file | `pcf` onset | `hisol` onset | verdict |
|---|---|---|---|
| `61deg_200mbarAr` | **12.09 s** | 17.09 s | 5.0 s gap — decisive, `pcf` = probe |
| `65deg_100mbarAr` | 12.48 s | 12.45 s | 30 ms — **inconclusive**, use test 2 |

So the test works only when there is a real gap. Require ≳1 s before trusting it.

**2. Spectral signature (fallback and cross-check).**

| | spectrum |
|---|---|
| **pump** | broad 600-1100 nm supercontinuum **plus** a deep-UV band near **280 nm** that switches |
| **probe** | an isolated band near **370 nm** and little else, plus residual 800 nm |

On the 61deg file both tests agree with each other and both contradict the group names,
which is why the naming is now known to be swapped. The probe band at ~370 nm also
independently corroborates the earlier "THE PROBE IS AN NO2 MONITOR" section — that was
derived from a completely separate argument.

### THE DRIVER SPECTRUM IS IDENTICAL IN BOTH STATES: 0.00 dB

Median over the RDW-high slices versus the RDW-low slices, after amplitude locking:

| band | HIGH | LOW | difference |
|---|---|---|---|
| RDW 220-320 nm (`hisol`) | 2.026e4 | 3269 | **+7.92 dB** |
| driver 600-1100 nm (`hisol`) | 7.710e5 | 7.716e5 | **−0.00 dB** |
| probe 345-440 nm (`pcf`) | 2.084e5 | 3.276e5 | **−1.96 dB** |

**This is the single strongest constraint in the dataset.** The driver supercontinuum is
unchanged to two decimal places between the RDW-on and RDW-off states. The soliton
dynamics, the self-compression, and therefore the RDW *generation* are identical in both
phases. So **the RDW is generated every shot and then removed** — the switch is
post-generation attenuation, not a change in generation or phase matching.

That kills a whole class of candidate mechanisms (anything acting through the compression
point, the pulse energy, or the dispersion at 800 nm) and points squarely at Hartley-band
absorption of an already-generated RDW. It also explains why the ozone back-action scan
found the compression point immobile: it should be.

Consistency check to do: 7.92 dB of in-band absorption over the post-compression length,
at σ(O₃) ≈ 5.5e-18 cm² near 282 nm, needs a column of ~4.1e17 cm⁻² downstream. Compare the
model's 8.6e14 cm⁻² naive steady state at 6 bar — **~500× short**, which is the gap to close.
**[Superseded 2026-09-10.]** That 8.6e14 assumed the pump destroys 52 % of the ozone every
pulse, which the 255 nm match later showed is far too strong. With that fit's limiting sinks
the 6 bar run reached 1.88e18 cm⁻² in 10 s (see "The 6 bar RDW trace" below) — enough ozone in
total, but peaked at the compression point where the RDW is generated. The gap is location,
not amount.

### The oscillation, measured

- **Switch depth +7.92 dB** in the 220-320 nm sum; ~20 dB at the band peak (800 → 7 counts),
  since the band sum carries background. Quote the in-band figure.
- **Duty cycle 41 %** high, after locking. This is a *shape* target and constrains the ratio
  of the two slow branch traversal times.
- **Anti-phase confirmed quantitatively**: correlation of log(RDW) against log(probe) over
  the locked region is **−0.772**.
- **Period 147, 156, 159, 160 s** (Schmitt-triggered rising edges at 302, 449, 605, 764,
  924 s), lengthening ~9 % across four cycles. The first 285 s interval is the locking
  transient — the first two maxima never reach full amplitude.
- The waveform is **square**: flat top, flat bottom, fast transitions. Relaxation-like, not Hopf-like, **in this run** — *[qualified 2026-09-10: one run at one set of parameters; across the thesis parameter grid depth and duty cycle vary between cells, so the oscillation type is open — see "THE GOAL IS A MECHANISM"]*.
  Amplitude locks after 2-3 cycles while the period keeps drifting.

Note the period here (~150 s) is not the ~100 s quoted from the 6 bar spectrograms. Since
period is not a target this changes nothing, but do not carry one dataset's period to
another.

### Detrend the probe MULTIPLICATIVELY

The probe's slow decay is mostly the probe source decaying before the fibre, not chemistry.
Its absolute modulation shrinks with the mean while the *ratio* is roughly preserved, so
divide by the envelope rather than subtracting it; an additive detrend leaves a spurious
amplitude decay for the model to explain with chemistry it does not have.

## THE PERIOD IS NOT A TARGET — and its fragility is a clue (2026-09-08)

**[Partly superseded 2026-09-10.]** Still true that the period is not a target. But the
"fragile period, robust amplitude ⇒ relaxation" argument below no longer holds: the parameter
grid (see "THE GOAL IS A MECHANISM", 2026-09-10) shows depth and duty cycle vary across
parameters too, so the oscillation type is open.

**Do not fit the oscillation period.** It is not reproducible even between nominally
identical experiments; a slight change in coupled energy or fill pressure moves it. Fitting
to it is fitting to noise. The target is that **the oscillation exists**, robustly, with the
right switching depth and the anti-phase relation between pump and probe.

**The irreproducibility is itself a measurement, and it redirects the search.**

| | period | amplitude |
|---|---|---|
| Hopf bifurcation | robust — set by the linear eigenfrequency | fragile — grows as √(distance from threshold) |
| Relaxation / SNIC | fragile — set by slow passage, diverges near threshold | robust — set by the fast branches |

The experiment reports **fragile period, robust amplitude**: the relaxation column. So the
object to look for is **slow–fast structure with hysteresis between two branches**, not a
Hopf bifurcation of a smooth limit cycle, and the right analysis is nullclines and
slow-manifold geometry rather than linear stability about a fixed point.

This also explains the earlier negative result: `reduced_switch_model.jl` scanned twelve
parameter sets looking for a limit cycle, found every trajectory relaxing monotonically to
a steady state, and concluded no oscillation. It was looking for the wrong object. A
relaxation oscillator needs an N-shaped (or otherwise folded) nullcline; the finding that
the destruction rate rises monotonically — f ∝ O₃^−0.55, so f·O₃ ∝ O₃^0.45 — is precisely
the statement that **the fold is missing**. Finding what folds that nullcline is the
problem, and it is a sharper question than "does it oscillate".

## SECOND DATA LOCATION: the anti-phase pump-probe experiment (2026-09-08)

```
~/Library/CloudStorage/Dropbox-Heriot-WattUniversityTeam/RES_EPS_Lupo/Projects/Ozone/
  phd/spec/air/dispersive_wave_amplification_and_deamplification/probe/30-05-2019/
  another set of data/
```

`double_spectrum_timeseries_*.h5` — two spectrometers recording simultaneously (`hisol`
and `pcf`), which is the source of the two anti-phased traces. Pre-rendered PNG crops name
both their wavelength and time windows, so they index the h5 cheaply.

Configuration, confirmed by the user 2026-09-08:

| item | answer |
|---|---|
| the `Ar` in the filenames | argon generates the **probe**; the sample fibre is **air** |
| `hisol` | the **probe** |
| `pcf` | the **pump** |
| time axis | **uncalibrated and wrong** — 40000 s should read ~600 s |

**The probe's slow decay is mostly NOT chemistry.** It is largely the probe source decaying
before the light enters the fibre. Some part may be NOx, but the two are not separable a
priori, so **the probe trace must be detrended before it becomes a target** and only its
oscillatory component is a chemistry observable. Fitting the raw probe envelope would be
fitting a property of the probe generation stage — the same class of silent error as
fitting a pre-coupling pulse energy.

This supersedes the earlier note that "the slow decay in the probe is due to the source
decaying away", which could be read as the ozone/NOx source. It is the optical source.

## JULIA SESSIONS: AgentREPL AND MCPRepl (2026-09-10)

A fresh Julia process here spends **~42 s** before doing any work (`using LupoAirOsc` 34.2 s,
first propagation compile 6.4 s); the same short propagation in a warm session takes 0.5 s.
Two MCP servers are registered with Claude Code at user scope:

- **`agentrepl`** — AgentREPL.jl 0.7.1, dev checkout at `~/.julia/dev/AgentREPL`. For
  sub-agents. Stdio, no network port; each `session` is an isolated Julia worker with Revise
  auto-loaded. One session per agent, `activate` LupoAirOsc, call `revise` after edits.
- **`mcprepl`** — MCPRepl.jl 0.1.0, dev checkout at `~/.julia/dev/MCPRepl`. For sharing your
  own REPL: `using Revise, MCPRepl; MCPRepl.start!()`. It serves code execution on a
  localhost-only port (3000 for the first REPL) with **no authentication**, so any local
  process can run code while it is up; stop it when done. Its adapter also has a `spawn_repl`
  tool that starts such a REPL in tmux without you. tmux is not installed, so that fails
  today, but agents must never call `mcprepl` tools.

Revise 3.17 is in the global environment. Long runs (the 6 bar chemistry) stay detached
scripts writing to disk, never inside a session. Packages LupoAirOsc does not depend on
(HDF5, StructuralIdentifiability, GlobalSensitivity) go in a separate environment. Full setup:
`LupoAirOsc/docs/OSCILLATION_PLAN.md` section 8.3.

## USER DECISIONS, SECOND BATCH, AND WHAT THE DOCUMENTS SAY ABOUT He-O₂ (2026-09-11)

Answers to the plan review of 2026-09-11, plus a reading of the paper and thesis. All recorded in
`OSCILLATION_PLAN.md` (sections 0, 2, 3, 6, 7, 8.2, 8.5).

- **The 255 nm run's 1.3 µJ is the in-fibre energy** (user). Thesis 7.2.2 and paper II C say
  "∼1.3 µJ coupled pump pulses".
- **Probe pulse energy: a few nJ** (user); thesis Figs. 7.1 and 7.4 say a few tens of nJ. Either way
  no strong-field effect: the probe reads the chemistry and never drives it.
- **The waveplate angle in the anti-phase filenames belongs to the probe-generation setup**, a
  separate apparatus, and says nothing about the pump (user). The angles in the pump-only files
  (87.75–91.15deg, tracking the 3/4/5 µJ labels) are presumably the pump attenuator that thesis 7.1
  describes. No calibration of that angle exists (user): the in-fibre energy is what the user states
  or what the paper or thesis writes, nothing else.
- **Baseline sink settings are not decided** (user: the agents check; only the experimental results
  are truth, every model parameter may move). The code defaults are stale: `Rates.WallParams` still
  has γ(O) = 1e-3, `air_1d_live.jl` defaults to wall×1 and φ = 0.14, and the 255 nm match was run
  with the wall off and `o3_from_adk=false`. The register records the baseline; no run relies on a
  script default, and the agent preamble now says so.
- **He-O₂ was never validated against a measurement.** The 9.3e24 m⁻³ / 3 % / 5 s ozone figures
  (paper III B and Fig. 6b; thesis 7.3.3 and Fig. 7.9b) are explicitly simulation output. What was
  measured (paper II E, Figs. 5 and 7a; thesis 7.2.4, Figs. 7.6 and 7.10a): 2.5 µJ in the fibre (user; the
  thesis writes "launched"), 12 bar 79 %–21 % He-O₂, **22.5 cm** — `examples/he_o2.jl` uses 22 cm,
  taken from the paper's Fig. 8 simulation example, and must be corrected. The RDW appears weakly
  near 255–260 nm, vanishes, reappears near 275 nm, red-shifts to ~325 nm within ~10 s and holds
  there for the rest of the run (the thesis figure spans 50 s, the text says 700 s), with no
  periodicity. The He-O₂ outputs in `examples/` date from 2026-09-04, before the 21.0 eV
  calibration; the case has not been re-run since. He-O₂'s regression target is the measured
  wavelength trajectory, after a re-run at 21.0 eV, φ = 1, 22.5 cm.
- **The grid data is raw h5, not images**: `Ozone paper/images/grid/` holds 31 files, three per
  pressure per fibre length labelled 3, 4, 5 µJ, plus a 7-file 6 bar energy scan at 27 cm. Thesis
  7.4 confirms the fill is air and does not state the energy convention; the figures' 1.6/2.2/2.8 µJ
  are in-fibre and the files' 3/4/5 µJ labels map onto them in that order (user). Thesis 7.4 reads the panels slightly differently from the plan
  (27 cm: 6 bar at 1.6, 2.2 and 2.8 µJ, not 7 bar; 22.5 cm: nothing at 7 or 8 bar) and says "no
  obvious conclusion can be drawn" about the window; A1 settles it from the files.
- **Tier 1 counts in the reduced model** (user), once it passes the reduced-model convergence gate
  (tables, interpolation, zones, ODE tolerance, bundling interval); Phase 5 confirms in the full
  model.
- **A4 writes the reduced-model skeleton in Phase 0** (user: "the best agent"), so A2 has something
  to analyse in Phase 1; zone structure is A4's call, starting from the 3-box
  `examples/oscillation/nox_titration_model.jl`.
- **Only the orchestrator commits, and only when satisfied** (user). Agents never run git commands
  that change history.
- **The remaining review gaps were resolved by the user's defaults (2026-09-11)**, all in the
  plan: the `:O3_diss` fit is gated on A2's Phase 1 verdict; `:O2_diss` is refitted in Phase 3
  against A1's digitised curve; the reduced-versus-full check points are four steady cells with
  stated tolerances; A3 also tabulates per-zone photolysis rates; the mean probe-band absorption
  is judging-only unless `hisol_stability2` justifies the detrend; A1 extracts any NO signature
  for `:N2_diss`, else a 0.1×–10× bracket; the Phase 0 gate has a fixed feature list; file
  formats are fixed in plan section 8.6; answers return through `ANSWERED_` files and
  `SendMessage`; A5 rejects through `audit/REJECT_*`; a permission allow list sits in
  `.claude/settings.json` in both `LupoAirOsc/` and `~/.julia/dev/`; Phase 5 order is 6 bar 2.2,
  6 bar 2.8, 5 bar 2.2, 8 bar 2.2 with at most four Julia processes; this file's "Current stage"
  was rewritten and the dispersion section corrected to 11 species with the N₂O₅ terms the code
  already has; the plan's section 4 was folded into the section 8.2 briefs.

## AGENT MODELS, EFFORT AND USAGE LIMITS (2026-09-11)

User decision: the RDW-oscillation campaign runs **all on Fable**. The orchestrator, A4 DYNAMICS
and A5 AUDIT run at **max** effort; A1 ORACLE, A2 REGISTER and A3 OPTICS at **xhigh**. The five
sub-agent definitions (frontmatter `model: fable` plus `effort`) are
`LupoAirOsc/.claude/agents/rdw-a*.md`, copied identically to user level in
`~/.claude/agents/rdw/` so a chat opened anywhere finds them (project agents are found only
between a chat's working directory and its repository root).
(An earlier version of the plan had an Opus orchestrator and justified centralising
interpretation by "the sub-agents are Fable", as if Fable were the weaker model. It is the
stronger one; the rule stands because only the orchestrator sees every agent's work.)

Traps, checked against the Claude Code docs for v2.1.267:

- Spawn agents by type with **no `model` override** — a per-call model beats the frontmatter.
- **`max` effort lasts one session**: run `/effort max` in every orchestrator session.
  `effortLevel` in settings.json accepts only up to `xhigh`, so the `max` in
  `~/.claude/settings.json` does not provide it.
- **`CLAUDE_CODE_EFFORT_LEVEL` overrides agent frontmatter** — never set it.

**At a usage limit: wait for the reset, then continue** (user, 2026-09-11). Claude Code
(v2.1.234 or later) does this itself in an open interactive session signed in to claude.ai
(`autoContinueAtUsageLimit`, on by default). Four caveats:
- it does not start the wait on its own for a limit more than 24 h from resetting (a weekly
  limit) — pick "Wait here, then continue automatically" in `/rate-limit-options`;
- it re-arms at most twice in a row;
- it needs Enter after more than ~30 min of sleep;
- it still stops on permission prompts.

Detached Julia runs keep going through the wait. After the continuation, the orchestrator
resumes interrupted agents with `SendMessage` and never restarts finished or still-running
work. Full procedure: `OSCILLATION_PLAN.md` section 8.4.

## WHERE THE DOCUMENTS AND DATA ARE (2026-09-08)

Everything is under
`~/Library/CloudStorage/Dropbox-Heriot-WattUniversityTeam/RES_EPS_Lupo/Projects/Ozone/`:

```
Ozone.pdf                                  the paper, 12 pp
MSabbah_PhD_Thesis_Final_submission.pdf    the thesis, 144 pp; ozone work is CHAPTER 7

Ozone paper/Figures/                       18 h5, 20 ipynb, 12 pdf in nine topic folders.
                                           Every folder has its own plot.ipynb, the
                                           provenance chain from h5 to published panel.

phd/spec/air/dispersive_wave_amplification_and_deamplification/probe/30-05-2019/
  another set of data/                     the anti-phase pump-probe experiment

data_for_nature/linear/                    Serdyuchenko/Gorshelev ozone cross-sections and
                                           the real/imaginary index parts -- the provenance
                                           of Luna's ozone refractive index
```

Plus `LupoAirOsc/targets/raw/`, already-parsed traces — do **not** re-read the 400 MB h5.

**Use exactly those two documents.** Older paper drafts and other thesis versions exist
elsewhere in the tree; ignore them. Most of the experimental record is in the thesis chapter
and in `Figures/`.

**The paper's figures mix measurement with the paper's own model output**, and captions do
not always distinguish them. A simulation result is not ground truth — it is a prior model's
answer, and that model's ozone is now known to be ~9200× high. Two numbers already bit us
this way: the 9.3e24 m⁻³ He-O₂ ozone figure, and the ~27th-order energy scaling, which came
from the paper's simulated 1.2/1.3/1.4 µJ curves and not from data. Record for every target
whether it is measured or simulated; when a figure does not say, ask.

## THE PUMP-OFF EXPERIMENT: THE OSCILLATION'S CLOCK RUNS ON PUMPED TIME (2026-09-11)

**File:** `Projects/Ozone/phd/spec/air fiber 10/27-01/6bar_air_3uJ_225mm_fiber10_40000frames.h5`
(2020-01-27, 6 bar air, 22.5 cm, fibre 10, pump only; extracted to
`LupoAirOsc/targets/raw/pumpblock_27-01-2020_6bar_air.tsv`). Identified from the folder's notebook
and saved plots, which match the user's figures. Layout `/raw/spectrum`, `/raw/λ`, text timestamps;
62.65 Hz. **Use 2.5 µJ in simulation** for the `3uJ` label (user, 2026-09-11).

The pump was blocked twice mid-oscillation. RDW / pump band (cancels coupling drift, which is only
2–5 % here), relative to the first peak:

| stretch | pump | RDW |
|---|---|---|
| fresh gas | on 2 → 315 s | dark 20–140 s (0.03); above 25 % from 193 s; peak at 273 s; 0.64 and falling at the block |
| gap 1 | off ~114 s | — |
| | on 429 → 443 s | spike 1.04 within ~2 s, down to 0.29 within ~5 s |
| gap 2 | off ~110 s | — |
| | on 553 → 667 s | back at once at 0.40, no spike; peak at 632 s, 1.17× the first |

A 2 s burst at 396 s at 18 % of normal pump power produced no RDW.

- **The phase-setting state survives ~110 s without the pump** (~190 s onset from fresh gas versus
  an immediate return after each gap).
- **The phase advances with pumped time:** peaks 359 s apart on the clock but 133 s apart in
  full-power pumped time (41 + 13 + 79 s), about one period of comparable 6 bar runs. Only one peak
  precedes the first block, so this run's own period is not measured.
- **A fast quencher relaxes in the dark and rebuilds within seconds** (the spike after gap 1), while
  the slow memory does not. Local ozone at the compression point fits the fast role. The missing
  spike after gap 2 is unexplained.

**Mechanism filter:** the slow variable must change while pumped and barely change in 110 s of
darkness. N₂O₅ re-partitioning, axial diffusion and end exchange keep evolving in the dark, so each
must be shown not to move the phase; a long-lived inventory made or destroyed only by the pump
(total NOx, N₂O) passes. Reproducing this run is part of Tier 2.

## SECOND PLAN REVIEW: DECISIONS AND APPARATUS FACTS (2026-09-10)

**Apparatus facts from the user — inputs, never fit parameters:**
- **Fill pressure is stable but declines slowly**, e.g. 6 → 5.92 bar over a run (~1.3 %).
- **Pulse-energy RIN is below 1.5 %.** Whether the mean energy also drifts slowly over a run is
  not known.

**Decisions (user):**
- **Success has two tiers.** Tier 1, the goal: a physically constrained, numerically converged
  mechanism that makes the RDW oscillate. Tier 2, confidence: it also reproduces the 22.5 cm
  grid's qualitative map, the pump–probe anti-phase relation including which leads, and a
  response to the apparatus drift and noise consistent with the experiment. Failing Tier 2
  does not undo Tier 1.
- **The 27 cm grid is not modelled.** Use its dynamics for insight; simulate it only to
  understand something specific, never as a target or validation.
- **Calibration and judging data stay separate.** Parameters are calibrated against observables
  that do not require reproducing the oscillation; the oscillation's own features only judge a
  mechanism.
- The identifiability rule is refined from "one observable per parameter" to **"never estimate
  parameters the data cannot tell apart"**.

**Two bounds on where run-to-run variability can come from:**
- Intrinsic chemical noise is negligible: ~1e12 ozone molecules in the 0.16 µL core gives ~1e-6
  relative fluctuations.
- Energy RIN barely affects the ozone source: at up to ~5.5th-order energy scaling, 1.5 % noise
  raises the mean source by ~0.3 % and averages out over the thousands of pulses in any chemical
  timescale.

So variability more likely reflects sensitivity near the edges of the oscillating window, or slow
drift. The slow pressure decline is the concrete candidate to test, including whether it can
stretch the period within a run (the anti-phase intervals grow from ~2000 to ~4300 frames).

**Gates added to the plan:** numerical convergence (tolerances, at least 201 z points, `nbundle`,
UPPE re-solve cadence — the artefacts this project has actually hit); a feedback-loop gate before
any mechanism search, which also tests which species carries the loop; sensitivity to the
apparatus drift and noise; and a `conflicts/` ledger for disagreeing evidence.

## THE GOAL IS A MECHANISM; THE PARAMETER GRID MAKES EVERY NUMBER A GUIDE (2026-09-10)

**The goal is a physical mechanism that makes the RDW oscillate — regardless of timescale,
duty cycle or peak-to-valley depth** (user). No single experimental number is a target; the
numbers only say whether a mechanism is in roughly the right regime.

Why, from thesis Figures 7.11 and 7.12 (printed pp. 97–98, PDF pp. 111–112; PDF page = printed
page + 14): a grid of 5–8 bar × 1.6/2.2/2.8 µJ for 22.5 cm and 27 cm fibres, 700 s per cell.
Read off the figures, to be re-derived from `Ozone paper/images/grid`:

- **Oscillation is confined to a narrow window** — at 22.5 cm, 6 bar/2.2 µJ (strong),
  5 bar/2.2 µJ (shallow), possibly 7 bar/2.8 µJ; at 27 cm, 6 bar/1.6 and 2.2 µJ and
  7 bar/2.2 µJ. Almost every other cell settles to a steady RDW.
- **The window moves with fibre length**: lower energy at 6 bar, and 7 bar, at 27 cm.
- **Period (~70 to >200 s), depth (shallow beading to deep switching), onset time**, and
  whether the modulation is also spectral (the 6 bar cells swing toward shorter wavelengths)
  all vary from cell to cell.

So the best test of a mechanism is a **map**: steady for most parameters, oscillating in a
limited window that moves the right way with fibre length. The gas fill and energy convention
for the grid are not in the captions — confirm from the thesis main text or ask.

This supersedes the "robust amplitude ⇒ relaxation oscillator" inference in the 2026-09-08
sections on the period and the waveform: amplitude is not robust across parameters, so the
oscillation type is open.

**Thesis interpretation worth knowing** (printed p. 99, discussion — an interpretation, not a
measurement): nitrogen's dissociation energy is higher than oxygen's, so O atoms are produced
faster, ozone forms first, and the DW shift happens before O atoms are lost to reactions with
nitrogen. Consistent with the red flag that `:N2_diss` should sit above `:O2_diss`.

**Access.** On 2026-09-10 macOS blocked every Dropbox folder for Visual Studio Code
("Operation not permitted"), so `images/grid` and the h5 data were unreadable from Claude
Code. Grant VS Code access under System Settings → Privacy & Security and restart it. An
identical copy of the thesis (same size and page count) is at
`~/Zotero/storage/REG4PSGD/MSabbah_PhD_Thesis_Final_submission.pdf`.

## CORRECTIONS FROM THE PLAN REVIEW (2026-09-10)

- **Fibre ends are right as modelled.** The gas cells are sealed, and the user confirmed that
  holding both ends at fresh-fill composition — the model's Dirichlet boundary for every
  species — is correct. It is not a detail: over a ~1000 s run, exchange with the ends reaches
  7–10 cm into the fibre, comparable to the ~7 cm between the input end and the compression
  point, so the ends are a real ozone and NOx sink and reduced models must keep an
  end-exchange term.
- **NOx is not simply catalytic.** The network already forms N₂O₅ (k24/k25) and loses N₂O₅ and
  NO₃ to the wall with γ = 1e-3 each. The sink exists; its rate is the question, and 1e-3 is a
  guess that looks high against the γ(O) ≲ 1e-5 the 255 nm fit found on the same silica.
- **The ozone gap is location, not amount** — see the superseded note in the anti-phase section.
- **The model is 1-D along z** (201 points); radial mixing (~20 µs) is treated as instant and
  walls enter as first-order losses. Temperature is fixed at 298 K; gas heating is bounded at
  ≲0.03 K even if all 2.2 mW were absorbed, against ~0.4 K for a 1 % change in O + O₃.

## FIXED EXPERIMENTAL PARAMETERS — never fit these (2026-09-08)

These are measured properties of the apparatus, identical across every fibre and
every pump pulse in this work. They are **not** free parameters, and a fit that
moves them is a fit that has gone wrong:

| quantity | value | note |
|---|---|---|
| core **diameter** | 30 µm | so `a = 15e-6` (radius) in every script |
| pump pulse duration | ~30 fs | all pump pulses, all pressures |
| repetition rate | 1 kHz | pump and probe |
| probe delay | +100 ns after the pump | same rep rate |
| centre wavelength | 800 nm | |

Fibre length is 22.5 cm for the 6 bar air / N₂-O₂ oscillation experiments and
**33 cm** for the 2 bar 255 nm probe-transmission experiment. Pressure is stated per
experiment.

**Pulse energy needs care: the number in a filename is not always the energy in the
fibre.** The convention differs per dataset folder, and getting it wrong is silent —
the ADK source runs at ~5th order in energy at 6 bar, so a factor 1.8 in energy is a
factor ~20 in ozone, and the result looks plausible rather than erroneous.

| dataset | labelled | **use in simulation** |
|---|---|---|
| `N2-O2/6bar_4uJ_89.5deg_225mm_fiber10_raw_60000frames.h5` | 4 µJ | **2.2 µJ** |
| `air_results/air_6bar_2.2uJ_22.5cm.h5` | 2.2 µJ | **2.2 µJ** |
| `air fiber 10/27-01/6bar_air_3uJ_225mm_fiber10_40000frames.h5` (pump-off run) | 3 µJ | **2.5 µJ** |

The 4 → 2.2 µJ mapping absorbs coupling loss and other transmission factors that are
deliberately out of scope for the model. Apply it; do not try to derive it. Every other
folder's convention is unknown and must be asked about, not inferred.

Those two files are byte-identical in size (488,542,540), a week apart in July 2020, same
6 bar / 22.5 cm / 60000 frames — so at the same in-fibre energy they differ in **one
variable only, synthetic N₂/O₂ versus real air**, isolating argon, water and CO₂.

Data lives at `~/Library/CloudStorage/Dropbox-Heriot-WattUniversityTeam/RES_EPS_Lupo/
Projects/Ozone/Ozone paper/Figures`; every folder there carries its own `plot.ipynb`,
which is the provenance chain from h5 to published panel.

If the model disagrees with data, the thing to move is a *model* parameter — the
effective ADK barrier, a rate coefficient, a branching ratio — not the apparatus.

## THE ADK O₂ BARRIER IS CALIBRATED: 21.0 eV, and φ = 1 (2026-09-08)

**This supersedes every `diss_yield`/φ value used anywhere earlier in this repo.**

### The measurement

The 255 nm probe-transmission experiment (2 bar air, 33 cm, 1.3 µJ, 30 fs, 1 kHz)
is the cleanest calibration in the whole dataset, and it had not been used. It is a
direct Beer-Lambert readout of the ozone **column**:

- 255 nm is the Hartley peak, σ₂₅₅ = 1.15e-17 cm².
- **No RDW is generated at these parameters** (user, confirmed by the model: 1e-5 of
  the output energy below 400 nm here versus 0.12-0.19 at 6 bar / 2.2 µJ). No RDW
  means no UV photolysis and no ozone feedback on the optics, so the column is set by
  the ADK O₂ source and the Chapman sinks *alone*. Every other case in this repo
  tangles those together.

Digitised from the red "Exp 1.3 µJ" trace, with column = −ln(I/I₀)/σ₂₅₅:

| t (s) | I/I₀ | O₃ column (cm⁻²) |
|---|---|---|
| 0.5 | 0.78 | 2.16e16 |
| 1.0 | 0.58 | 4.74e16 |
| 2.0 | 0.31 | 1.02e17 |
| 4.0 | 0.125 | 1.81e17 |
| 6.0 | 0.055 | 2.52e17 |
| 8.0 | 0.02 | 3.40e17 |

The rise is close to linear, slightly sub-linear: **~4.3e13 cm⁻² per pulse** at 1 kHz.
(Linear in column ⇒ exponential in I/I₀; the late-time trace sits above that, so
there is a real sink, not just a source.)

### The discrepancy, and what it was

At the fixed geometry with the old 15.5 eV barrier the model produced
**3.97e17 cm⁻² per pulse — 9200× too much**. That is the true size of the ozone
overproduction, and it is far bigger than the ~10× previously inferred from the
indirect arguments (Chappuis null, constant 800 nm spectrum, sweep rate). It is why
every air run in this repo has needed `diss_yield` between 1e-4 and 0.05 to look sane.

The overproduction is **entirely in the barrier**, not in the geometry or the pulse.
Scan at 2 bar / 33 cm / 1.3 µJ / 30 fs / a = 15 µm (`examples/oscillation/o2diss_barrier_fit.jl`,
one propagation per energy, all trial barriers evaluated on the same field via the new
`PropAir` `diss_scan` kwarg):

| barrier (eV) | source (cm⁻²/pulse) | φ needed |
|---|---|---|
| 15.5 (old) | 3.97e17 | 1.08e-4 |
| 18.0 | 6.93e15 | 6.2e-3 |
| 20.0 | 2.34e14 | 0.183 |
| 20.8 | 5.84e13 | 0.736 |
| **21.0** | **4.11e13** | **1.045** |
| 21.4 | 2.03e13 | 2.12 |
| 23.0 | 1.15e12 | 37.3 |

**`:O2_diss` is now 21.0 eV in `PhysData.jl`** (old value commented in place), and
`diss_yield` should be **1.0** — the branching-ratio fudge disappears entirely.

Why the barrier and not φ: φ is a *branching ratio*, and a branching ratio of 1e-4 for
a superexcited state is not physical. The ADK "ionisation potential" for a dissociation
channel is by contrast an openly **effective** parameter — ADK is a hydrogenic
tunnelling formula pressed into service for molecular dissociation — so it is the
honest place to absorb a factor of 1e4.

### A bound worth keeping: the source cannot be steeper than ~10th order

The same figure carries the paper's own simulations at 1.2 / 1.3 / 1.4 µJ
(I/I₀ ≈ 0.87 / 0.30 / 0.02 at t = 2 s), i.e. columns in the ratio 1 : 8.8 : >29, which
is d ln(col)/d ln(E) ≈ **27**. The model gives 6.7 at 15.5 eV and 12.0 at 21.0 eV, and
**no barrier can reach 27**:

> For any tunnelling rate R = A·exp(−B/E) with A ~ 1e16-1e17 s⁻¹, the intensity order at
> an operating point is ½·ln(A/R). Producing the measured source needs a dissociated
> fraction ~2e-6 in 30 fs, i.e. R ~ 1e8 s⁻¹, so ln(A/R) ≈ 18-21 and the order is
> **9-10.5 — independent of the barrier**, which only sets *where* in z that rate is
> reached. Measured compression sensitivity d ln(I_peak)/d ln(E) is 0.9-1.4 at these
> parameters, so the total energy order is bounded at roughly 9-15.

So the remaining factor ~2 in energy sensitivity is **not** in the source term. It has
to come from somewhere else — most plausibly a threshold, e.g. the self-compression
point crossing the fibre exit. At 2 bar / 33 cm / 30 fs the model puts the compression
point at 21.9 cm (peak 1.18e14 W/cm²), moving 27.4 → 18.5 cm as the energy goes
1.0 → 1.6 µJ. Worth checking against where the experiment thinks it compresses; do not
"fix" it by moving `a` or `τ`.

### What this invalidates

Every air result in `examples/oscillation/` predating this used 15.5 eV together with a
φ chosen to compensate, and the compensation was applied as a flat multiplier on the
source. That is *not* equivalent to raising the barrier, because the barrier also
changes the **spatial profile** of the source (it concentrates it into the small region
where the intensity is highest) and its **energy scaling** (order 6.7 → 12.0). The
N₂ source moves with it, so the NOx supply — the thing the whole oscillation argument
turns on — needs re-measuring, not rescaling.

Measured at the 6 bar oscillation parameters (22.5 cm, 30 fs, `sixbar_recalibrated.jl`),
with the source columns in cm⁻² of atoms per pulse:

| E (µJ) | I_peak (W/cm²) | z_compress | UV frac | RDW | O₂ src (21 eV) | N₂ src (21 eV) | max f(O₃) |
|---|---|---|---|---|---|---|---|
| 2.0 | 1.40e14 | 7.84 cm | 0.174 | 268.2 nm | 1.89e14 | 7.09e14 | 0.403 |
| 2.2 | 1.43e14 | 7.39 cm | 0.185 | 268.5 nm | 2.20e14 | 8.28e14 | 0.449 |
| 2.5 | 1.48e14 | 7.04 cm | 0.193 | 268.5 nm | 4.46e14 | 1.68e15 | 0.520 |

At 2.5 µJ the old barrier gave an O₂ source of 1.78e18, so the recalibration is a factor
**3990** here. Naive steady state, source over per-pulse pump destruction:
4.46e14 / 0.52 = **8.6e14 cm⁻²**, versus ~3.4e18 before. Spread over the whole 22.5 cm
that is a fraction 2.6e-7; concentrated into ~1 cm around the compression point it is
~6e-6. Compare the switch below — this may now be *too little* ozone to modulate the
RDW, which would be the opposite failure to the one just fixed. `rdw_energy_match.jl`
is the test.

**Red flag on `:N2_diss`.** With both barriers at 21.0 eV, N₂ is 3.76× more abundant
than O₂, so the nitrogen source is 3.8× the oxygen source — the model now makes N atoms
faster than O atoms. That is almost certainly wrong, and it points the opposite way to
the user's suggestion of lowering `:N2_diss` to 15.5 eV, which would make it worse by
another ~1e4. Physically N₂ should sit *above* O₂: stronger bond (9.79 vs 5.12 eV) and
higher ionisation potential (15.58 vs 12.07 eV). Left unchanged pending an observable.

### Ozone back-action on the pump is weak; the switch is pure absorption

Scanning uniform ozone through the 6 bar / 2.5 µJ propagation:

| O₃ fraction | I_peak | z_compress | UV frac | O₂ src ratio |
|---|---|---|---|---|
| 0 | 1.482e14 | 7.04 cm | 0.194 | 1.000 |
| 1e-5 | 1.482e14 | 7.04 cm | 0.177 | 0.999 |
| 1e-4 | 1.481e14 | 7.04 cm | 0.097 | 0.995 |
| 3e-4 | 1.480e14 | 7.04 cm | 0.068 | 0.985 |
| 1e-3 | 1.475e14 | 7.04 cm | 0.059 | 0.952 |

So ozone attenuates the RDW strongly (5.2 dB by 1e-3) but barely touches the pump:
the compression point does not move at all and the O₂ source falls only 5 % at 1e-3.
**There is no useful optical feedback from ozone onto its own source.** Any oscillator
has to close its loop through the chemistry, not through the compression.

**`:N2_diss` remains at 21.0 eV and is now UNCONSTRAINED by any measurement.** It was
21 before this change and is unchanged; the user has suggested it might be 15.5. It was
deliberately left alone here because nothing in the 255 nm data constrains it, and
guessing at it would repeat the mistake this section exists to correct. Nitrogen needs
its own calibration observable — the NO side-scattering lines (onset ~35 s) are the
obvious candidate.

## THE 255 nm CURVE IS MATCHED — and it took three corrections (2026-09-08)

`examples/oscillation/probe255_match.jl` now reproduces the measured
I/I₀(t) at 255 nm (2 bar air, 33 cm, 1.3 µJ, 30 fs) to within ~5 % over the full
8 s, with **no free multiplier anywhere** — `diss_yield = 1`:

| t (s) | model column | exp column | model I/I₀ | exp I/I₀ |
|---|---|---|---|---|
| 0.60 | 2.51e16 | 2.16e16 | 0.749 | 0.780 |
| 1.67 | 6.94e16 | 7.14e16 | 0.450 | 0.440 |
| 3.03 | 1.26e17 | 1.44e17 | 0.236 | 0.190 |
| 4.99 | 2.06e17 | 2.14e17 | 0.093 | 0.085 |
| 6.91 | 2.85e17 | 2.92e17 | 0.038 | 0.035 |
| 7.85 | 3.23e17 | 3.40e17 | 0.024 | 0.020 |

Settings: `diss_yield=1.0, n2_yield=0.0, wall=0 (WallParams()), o3_from_adk=false`,
`:O2_diss = 21.0 eV`. One UPPE solve (frozen optics, justified above).

**Three independent things were wrong, and all three had to be fixed together.**
Each was masked by the others, which is why the indirect arguments only ever
suggested "about 10×".

### 1. The source: `:O2_diss` 15.5 → 21.0 eV (factor 9200)

See the section above. Validated directly: with the two sinks below removed, the
model gives 2.12e16 cm⁻² at t = 0.5 s against a measured 2.16e16 — a **2 % match on
a number set by the source alone**, since nothing has had time to remove anything.

### 2. The silica wall is essentially inert to O atoms: γ(O) ≲ 1e-5, not 1e-3

`Rates.WallParams` had γ(O) = 1e-3. That is not a small correction here, because the
pump re-dissociates its own ozone every pulse and the released O atoms then race:

- reform ozone, O + O₂ + M at **3.19e5 s⁻¹** at 2 bar
- die on the wall, γ·v̄/(2a) = **2.09e4 s⁻¹** at a = 15 µm

So ~6 % take the wall each time, giving ~0.6 % column loss per pulse and equilibrium
in ~160 pulses. Measured column instead rises for 8000 pulses, bounding per-pulse loss
at <1.2e-4. Scan (nitrogen off, 8 s plateau):

| wall scale | plateau (cm⁻²) |
|---|---|
| 1.0 | 2.9e15 |
| 0.2 | 1.45e16 |
| 0.02 | 9.5e16 (still rising) |
| 0 | 1.52e17 (still rising) |

Plateau goes exactly as 1/γ, which confirms the mechanism. **Physically this is a
passivated surface** — a fibre that has been sitting in ozone and atomic oxygen for
minutes is not clean silica, and recombination coefficients drop orders of magnitude
on oxidised/passivated SiO₂. Note the existing docstring already bounded γ(O₃) < 1e-9
from the He-O₂ 700 s persistence; the same argument evidently applies to O.

### 3. `:O3_diss = 12.53 eV` is far too low — the pump eats its own ozone

Even with zero wall the model saturated at 1.52e17 against 3.40e17. The residual sink
is **Chapman, amplified by the pump**: each pulse re-dissociates ~30 % of the ozone
(Dissfrac_O3 peaks at 0.31 here) into O atoms sitting *inside* the ozone cloud. Almost
all reform ozone, but a fraction

    k₃[O₃] / (k₁[O₂][M]) = 8e-15 × 3e16 / 3.19e5 = 7.5e-4

instead runs O + O₃ → 2 O₂, killing **two** ozone each. That is ~4.5e-4 of the column
per pulse, i.e. equilibrium near 9e16 — the right order. Setting `o3_from_adk=false`
removes it and the curve matches.

`o3_from_adk=false` is the *limiting* case, not the fix. `:O3_diss` is set equal to
ozone's real ionisation potential (12.53 eV) on the argument that ionised ozone
dissociates anyway; that argument may be fine but the *number* is not, and it should be
calibrated the same way `:O2_diss` just was. Rough size of the correction: the sink must
shrink ~10-20×. Taking `Dissfrac_O3` from 0.31 to ~0.03 needs the rate down ~12×; the O₂
barrier scan measured **~5.9× per eV** at 20-23 eV, and the Ip^0.5 scaling of the ADK
exponent softens that to ~3.9× per eV at 12.5 eV, giving **~1.8 eV**. So the barrier would
go **12.53 → ~14.3 eV**. Treat this as an order-of-magnitude guide to the SIZE of the
correction, not a value to adopt: it extrapolates a sensitivity measured at a different
barrier and a different saturation level, and the "10-20×" target is itself a judgement.
Not yet
fitted — it needs the same barrier scan against this curve, which is cheap now that
`PropAir`'s `diss_scan` exists.

### What is still open

- `:N2_diss` remains unconstrained (see the red flag above). These runs used
  `n2_yield = 0`. Once (2) and (3) are fixed, nitrogen's effect on THIS curve is
  small, so this measurement does not pin it — it needs the NO side-scattering data.
- The ~27th-order energy scaling of the paper's 1.2/1.3/1.4 µJ simulations is still
  unexplained, and by the tunnelling bound above it cannot come from the source.

## The 6 bar RDW trace: the sweep reproduces, the SWITCH does not (2026-09-08)

`examples/oscillation/rdw_energy_match.jl`, 6 bar air / 2.5 µJ / 22.5 cm / 30 fs, run
with the calibrated source and the limiting sinks from the 255 nm fit
(`diss_yield=1, n2_yield=0, wall=0, o3_from_adk=false`, 101 z-points, one UPPE re-solve
per second of experiment time). **Abandoned at t = 10 s** — see the cost note below.
Deliberately NO nitrogen chemistry, so this is the ozone-only control, not an
oscillation test.

| t (s) | RDW centroid | band energy | O₃ column (cm⁻²) |
|---|---|---|---|
| 0.0 | 282.1 nm | 0.00 dB | 0 |
| 1.0 | 288.3 nm | −1.12 dB | 4.08e17 |
| 2.0 | 290.4 nm | −1.50 dB | 7.26e17 |
| 4.0 | 291.8 nm | −1.82 dB | 1.18e18 |
| 7.0 | 292.8 nm | −2.04 dB | 1.60e18 |
| 10.0 | 293.5 nm | −2.17 dB | 1.88e18 |

**What works.** The RDW is born at **282.1 nm** — the experiment says it starts at
280 nm — and sweeps monotonically to longer wavelengths as ozone accumulates, which is
exactly the early behaviour described from the 6 bar spectrograms. That is the first
time the model has reproduced this, and it comes out of the calibrated source with no
tuning.

**What does not.** The band energy falls only **2.2 dB in 10 s and is plainly
saturating** (per-second decrements 1.12, 0.38, 0.20, 0.12, 0.09, 0.07, 0.06, 0.05,
0.04), against an observed switch depth of ~9.7 dB. It is not for want of ozone: the
column reaches 1.88e18 cm⁻², 5.5× the entire 2 bar / 8 s column, and is still rising.

Two reasons, and they compound:

1. **The sweep is self-protecting.** σ(O₃) is 5.5e-18 cm² at 282 nm but 2.6e-18 at
   293 nm. As ozone builds, phase matching moves the RDW to where ozone is more
   transparent. A band that runs away from its own absorber cannot switch itself off.
2. **The ozone sits at the generation point, not downstream of it.** The compression
   point is at 7.04 cm of 22.5 cm and the ozone peaks there too, so the RDW is created
   inside the absorber rather than having to traverse it. This is consistent with the
   uniform-ozone back-action scan above (5.7e-4 uniform would give ~5 dB; localised
   gives 2.2), and it is the same "production and sinks are co-located" problem flagged
   earlier, seen now from the optical side.

**So absorption of a freely-sweeping RDW does not produce the observed switch.** Either
something pins the band at 280 nm while ozone builds (the experiment does show it
returning to 280 nm to oscillate, *after* the initial sweep and disappearance), or the
switch is not absorption at all.

### Cost note — this is the practical blocker for any long run

At 201 z-points the 6 bar case runs ~0.3 s of experiment time per minute of wall clock,
i.e. **~7 hours per 200 s**, and the observed phenomenon needs ~1400 s. The UPPE solves
are *not* the bottleneck: at one re-solve per second of experiment time they are ~13 s
each against ~500 s of chemistry per 4 s of experiment time, so the cost is the
reaction-diffusion integration itself — ~0.25 s per kick for 2211 unknowns (201 nodes ×
11 species). The `FastRHS` work sped up the residual; the remaining cost is the implicit
solver's Jacobian and linear algebra. Halving the grid to 101 points helps but not
enough. **Any plan that needs hundreds of seconds of experiment time needs this
addressed first** — sparse/banded Jacobian, or a reduced model with the switch
tabulated.

## THE PROBE IS AN NO2 MONITOR (2026-09-08) — read this before touching the chemistry

This is the result that reorganises everything above. It comes from the user's
own pump-probe figure (two panels, time 0–1400 s on the vertical axis), plus the
line plot of the two integrated traces.

**The experiment.** 6 bar air, 2.5 µJ pump, 22.5 cm fibre, 1 kHz. A PCF-derived
broadband probe co-propagates at the same rep rate, **~100 ns behind the pump**.
Panel (a) is the transmitted probe, panel (b) the pump-generated RDW. The two
line traces are those panels summed over wavelength. A SEPARATE probe run at
**2 bar / 1.3 µJ generates no RDW at all** — that one is the clean
ozone-formation case (`examples/oscillation/probe_2bar.jl`).

**The probe band is 345–440 nm, and in that band ozone is invisible while NO₂ is
near its peak.** From the tables in `PhotoChem.jl`:

| λ (nm) | 350 | 370 | 390 | 410 | 430 |
|---|---|---|---|---|---|
| σ(O₃) cm² | 1.95e-22 | 1.04e-23 | 8.18e-24 | 2.79e-23 | 7.71e-23 |
| σ(NO₂) cm² | 4.81e-19 | 5.78e-19 | 6.15e-19 | 6.05e-19 | 5.52e-19 |
| σ(NO₃) cm² | 0 | 0 | 0 | 1.18e-20 | 1.75e-19 |

**σ(NO₂)/σ(O₃) = 5.57e4 at 370 nm.** So panel (a) measures NO₂, not ozone. The
white region below ~340 nm in that panel is simply where the probe source has no
output; it is not an absorption feature. (The void near 400 nm in the earlier
version of the figure is the anti-resonant fibre's own non-guided resonance
band.)

**Therefore the observed anti-phase is a DIRECT measurement of NOx titration.**
Probe transmission is high exactly when the RDW is absent, i.e. NO₂ is high
exactly when the RDW transmits, and the RDW transmits when ozone is low. NO₂ and
O₃ anti-correlated is the titration signature itself. This is measured, not
inferred, and it settles the mechanism question the whole "oscillation hunt"
section above was circling.

**It also calibrates the NITROGEN source, which nothing else in this project
does.** Converting the observed probe swing to a column at σ(370 nm) = 5.78e-19:

| probe swing | NO₂ column (cm⁻²) | mean over 22.5 cm (cm⁻³) | as % of 1.5e20 |
|---|---|---|---|
| 1 dB | 3.99e17 | 1.77e16 | 0.0118 |
| 2 dB | 7.97e17 | 3.54e16 | 0.0236 |
| 3 dB | 1.20e18 | 5.32e16 | 0.0354 |
| 5 dB | 1.99e18 | 8.86e16 | 0.0591 |

The observed swing is a factor ~2 (≈3 dB), so **target NO₂ column ≈ 4e17–1.2e18
cm⁻²**. The model brackets it:

| run | NO₂ column (cm⁻²) | dB at 370 nm | O₃ column (% of density) |
|---|---|---|---|
| f2D = 0.5, η = 0.1 (150 s) | 6.32e16 | 0.16 | 1.3 % |
| f2D = 0.9, η = 0.3 (at t = 12 s) | 6.25e18 | 15.7 | 0.50 % |

So f2D = 0.5 / η = 0.1 is ~20x short and f2D = 0.9 / η = 0.3 is ~5x over.
**Interpolate to the measured swing rather than scanning further**; roughly
f2D ≈ 0.7, η ≈ 0.2. Note NO₂ responds superlinearly (99x for an 18x change in
the product of net-NOx-yield and η) because titration feeds back on itself.

**The full time sequence, from the user's 6 bar figure with a zoomed first-10 s
panel (this CORRECTS an earlier reading in this file that said the band never
moves):**

1. **t ≈ 0–1.5 s** — no RDW.
2. **t ≈ 1.5–5 s** — RDW appears at 280 nm and **sweeps to ~330 nm**, driven by
   ozone accumulating at the compression point and moving the phase-matching.
3. **t ≈ 5–130 s** — dark. Attributed to O₃ absorption.
4. **t ≈ 130 s onward** — oscillation, bursts **at a FIXED ~280 nm**, period
   ~80–100 s at 6 bar. (In the 1400 s record at slightly different conditions:
   bursts at 130, 250, 400, 570, 740, 900, 1090, 1270 s, period growing 120 → 190 s
   — a threshold oscillator riding a slow drift, square tops, dwell ~60–80 s high
   and ~80–110 s low.)

**The model reproduces phase 2 correctly and is ~10x too fast**, which is the
same factor the Chappuis null and the constant 800 nm spectrum give
independently. Model centroid at φ = 0.05, 6 bar, f2D = 0.5, η = 0.1:

| t (s) | 0.00 | 0.01 | 0.21 | 0.37 | 1.0 | 3.5 | 5.1 |
|---|---|---|---|---|---|---|---|
| centroid (nm) | 272.8 | 314.1 | 324.7 | 332.2 | 326.6 | 330.2 | 328.1 |
| UV out (dB) | 0 | −12.8 | −19.8 | −21.3 | −23.7 | −19.3 | −17.9 |

i.e. 273 → 332 nm in 0.37 s against the experiment's ~3.5 s. **So the relocation
is NOT a model defect — it is the experiment's own phase 2, played 10x fast.**
Three independent measurements now agree on the same factor of ~10, which makes
it a calibration correction, not a mystery: either φ ≈ 0.005 rather than 0.05, or
a missing sink removing 90 % of the ozone.

**What the model does NOT reproduce is phases 3 and 4** — the return of
phase-matching to 280 nm and the sustained oscillation there. Returning to 280 nm
requires the ozone AT THE COMPRESSION POINT to be cleared back to near zero. The
125 s dark stretch (phase 3) is plausibly the NOx accumulation time before
titration can do that, which ties directly to the NO₂ calibration above and to
the f2D analysis (net NOx yield only 0.135 per N atom at f2D = 0.5).

**The probe's ~100 ns delay is itself a measurement of the per-pulse photolysis
fraction.** O + O₂ + M → O₃ reformation vs pressure, and how much has reformed
when the probe arrives:

| P (bar) | 2 | 3 | 5 | 6 | 8 | 12 |
|---|---|---|---|---|---|---|
| k₁[O₂][M] (s⁻¹) | 3.19e5 | 7.19e5 | 2.00e6 | 2.87e6 | 5.11e6 | 1.15e7 |
| 1/e (ns) | 3131 | 1392 | 501 | 348 | 196 | 87 |
| reformed at 100 ns | 3.1 % | 6.9 % | 18 % | 25 % | 40 % | 68 % |

At 6 bar a quarter has reformed, so the probe sees the gas with ~75 % of the
photolysed odd oxygen still atomic.

## Ozone sinks: the pump beats the RDW, and both are localized (2026-09-08)

- **The pump destroys 44.6 % of the ozone PER PULSE at the compression point**,
  by strong-field dissociation at 12.53 eV (`:O3_diss`). Measured from
  `Dissfrac_O3` at 6 bar / 2.2 µJ: peak 0.446 at z = 7.2 cm, against
  `Dissfrac_O2` = 1.10e-2 at the same z. An e-folding in 2 pulses. This is
  **3.7x larger than the RDW's 12 % photolysis** and is already wired into the
  chemistry through `apply_lunaion(...; o3_from_adk=true)` (the default).
- **`diss_yield` (φ) does NOT scale the O₃ sink**, only the O₂/N₂ source. So the
  local balance at the compression point is
  `2·φ·f_O2·[O₂] = f_O3·[O₃]`, giving
  `[O₃] = 2 × 0.05 × 0.011 × 0.21 / 0.446 = 5.2e-4` of total density = **0.052 %**
  — squarely inside the range the experiments allow. **The model gets the
  compression point about right and is wrong everywhere else.**
- **The real defect is that both sinks are violently localized while the source
  is not.** Strong-field dissociation follows a steep intensity law so it exists
  only at the compression point; the RDW is born there and travels forward so it
  only sweeps downstream. Production happens along the whole fibre. **Upstream of
  the compression point the model has NO ozone sink at all.** Even an entrance
  `Dissfrac_O2` of 1e-5 accumulates to percent levels over 1e5 shots. That is
  exactly the entrance plug (0–3 cm) that never clears in the 150 s runs and
  drags the band to 330 nm. The user's measurements are COLUMN measurements, so
  they say ozone is low everywhere including upstream — something sweeps that
  region and the model does not have it. **This is the open problem.**
- **Withdrawn: the per-shot photon-budget bound.** I argued the RDW could not
  destroy ozone because 65 nJ at 270 nm is 8.83e10 photons against 3.1e14 ozone
  molecules in the beam. That is per-shot reasoning and the process takes tens of
  seconds. Burn-through time = molecules / (photons × 1 kHz):

  | O₃ fraction | 1.3e-2 | 1e-3 | 3e-4 | 1e-4 |
  |---|---|---|---|---|
  | burn-through | 3.5 s | 0.3 s | 0.1 s | 0.03 s |

  Nothing is photon-limited on the oscillation timescale. **The Chappuis null and
  the constant 800 nm spectrum remain valid and independent**; they still say the
  ozone column is ~10x below the model.
- **Caveat on the φ = 0.05 calibration.** It was fitted to "9.3e24 m⁻³ after ~5 s"
  quoted from the paper in `examples/he_o2.jl`. That may be the PAPER'S MODEL
  output rather than a measurement — the same file says the paper's model
  reproduces the He-O₂ case. If so, φ was calibrated model-against-model, while
  the user's Chappuis/IR data are real and disagree by ~100x. **Unresolved; ask
  before relying on φ = 0.05.** The 2 bar / 1.3 µJ no-RDW probe run is the clean
  replacement calibration.
- Bug found on the way: `setup_propair`'s `densityfun` called CoolProp live with
  whatever the partial-pressure spline returned, and CoolProp's root-find throws
  on denormal-scale pressures (`P = 4.4e-71`) although `density(gas, 0.0)` is
  fine and 1e-30 bar is fine. `PropAir.safe_partial_pressure` now snaps anything
  below 1e-20 bar to zero. Only appeared at 2 bar; the 6 bar runs never drove a
  species that low.

**Current best-guess parameter set, and the run testing it (2026-09-08).**
Combining every constraint above: ozone accumulates ~10x too fast (three
measurements) and NOx accumulates too slowly (the probe bracketing), so
**φ = 0.005, f2D = 0.7, η = 0.2**. Running
`air_1d_live.jl 400 0.7 0.2 1.0 0.005 2 0.1 201`. The test is the full sequence:
a 280 → 330 nm sweep taking ~3.5 s rather than 0.37, a dark stretch of order
100 s, then a return to 280 nm. Anything that reproduces phases 1–3 is already
past every previous run; phase 4 is the real target. Note a 400 s run at
~150 s/s is ~17 h, so read the Monitor PNG/JLD2 as it goes rather than waiting.

---
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

## Dispersion: how the 11-species chemistry maps onto tracked gases

`State.effective_densities` collapses the 11-species chemistry onto the
gases Luna actually has dispersion data for, **conserving atoms**:

```
N2_eff = N2 + 0.5*(N + NO + NO2 + NO3) + N2O + N2O5
O2_eff = O2 + 0.5*(O + O1D) + 0.5*NO + 1.0*NO2 + 1.5*NO3 + 0.5*N2O + 2.5*N2O5
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
