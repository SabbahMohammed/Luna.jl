# Luna fork (SabbahMohammed/Luna.jl, branch `agent-work`) — working context

This fork exists for **LupoAirOsc** (`~/.julia/dev/LupoAirOsc`), pinned there by path in its `Manifest.toml`. **The control document is `LupoAirOsc/STATE.md`** (one page: configuration, the one question being worked, runs in flight, decision log); parameter values are in `LupoAirOsc/params/baseline.toml`. Read those first. The 1,900-line history this file used to hold is in `.github/archive/CLAUDE_2026-09-17.md` — the record of every bug found and every result up to 2026-09-17; consult it before re-deriving anything.

## Standing rules (user; details in STATE.md and baseline.toml)
- The pump pulse is the **measured FROG field** (`LupoAirOsc/data/frog_13-03-2019/…/Ek.dat`, 28.7 fs, negatively chirped); every energy quoted is the **COUPLED** energy, of which 20 % is ASE (the pulse carries 0.80 × coupled). Labels 3/4/5 µJ are incident: coupled = 0.62 × (2019 setup, 26 µm core) or 0.84 × (2020 setup, 30 µm core).
- Ionisation is Luna's PPT rate with the Talebpour Z_eff (O₂ 0.53, N₂ 0.9), rate scale 1; Raman (N₂ + O₂, added in this fork) on. The O₂/N₂ pseudo-ADK dissociation channels are retired; the O/N source is a yield × the mode-averaged ionisation fraction.
- The O₃ strong-field channel is on in propagation **and** chemistry (`dissociate=(:O3,)` ↔ `o3_from_adk = true`; `LupoAirOsc.State.check_o3_pairing` enforces it). Its barrier (12.53 eV) is provisional and is being calibrated.

## What this fork changes relative to upstream Luna (all in `src/`)
`Capillary.jl`: cached `neff`/`neff_wg` for `loss=:gas`, per-z `gas_mixture` cache, pressure-spline clamping, `losslabel` for `:gas`; `PhysData.jl`: ozone complex index (Hartley + Chappuis) and O₂/N₂ Raman parameters, the `:O2_diss`/`:N2_diss`/`:O3_diss` channel constants; tabulated `GasMixtureIndex` so a z-dependent linear operator costs ~35 % not ~200 %. See the archive for the why of each.

## Lessons that must survive (from the archive)
- When a comparison isolates one case and that case is the odd one out, suspect the code path only it exercises before inventing physics (the frozen-linop bug).
- A simulation figure in the paper is not a measurement; record for every number whether it was measured or simulated.
- Pulse energy in a file name is a label; the convention differs per dataset and an error is silent (the source runs at ~8th order in intensity).

Keep this file short. New facts go to `LupoAirOsc/STATE.md`; history goes to the archive.
