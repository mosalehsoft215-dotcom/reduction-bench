# Ring-down campaign folders

One folder per measuring day. Every folder holds the same four files.

## campaign.json

- `cavity.mirror_separation_cm`: distance between the two ring-down mirrors.
- `cavity.heated_length_cm`: length of the tube furnace around the reactor section, centred between the mirrors. The furnace holds its whole length at one temperature. The arms that join the furnace to the two mirror mounts are water-cooled, and the gas in them is at the run's `arm_temperature_K`.
- `cavity.mirror_purge`: the gas that keeps the mirrors clean, or `"off"`.
  - With a purge, a small flow enters next to each mirror and fills the end of that arm, so only part of the mirror separation holds the run's gas. Where the purge gas meets the run's gas has not been measured. The purge flows were the same in every run of a day.
  - `"off"`: the whole mirror separation holds the run's gas.
- `cavity.purge_o2_mole_fraction`: O2 mole fraction of the purge gas (null when the purge is off).
- `laser`: the dye laser's time-averaged spectral profile, measured each morning with a scanning Fabry-Perot.
  - `spectral_profile`: its shape, `"gaussian"` or `"sech2"` (intensity proportional to sech²).
  - `fwhm_cm`: its full width at half maximum, in cm-1; null on days when the Fabry-Perot was not available.
- `detector.time_constant_ns`: the amplified photodiode behind the output mirror has a single-pole (first-order low-pass) response with this time constant.
- `digitiser`: sample interval, samples per trace, and the number of laser shots averaged into each trace. Sample `a0` is the laser trigger. The pulse lasts about 5 ns.
- `runs`: the runs of the day. Each has its `pressure_torr`, the reading of the absolute capacitance manometer on the cell, and its `arm_temperature_K`, read by a thermocouple on the cooled arms.
  - The `vacuum` run is scanned with the cavity evacuated: gas and purge valves closed, turbo pump at its base pressure (below 1e-5 Torr). Some days have no vacuum run.
  - A `reference` run is a certified gas cylinder, scanned with the furnace off, so the whole cell is at the given gas temperature. Its O2 mole fraction is given. A day can have more than one reference run.
  - `sample` runs are the reactor feed. A sample run's temperature is the gas temperature inside the furnace; it is not known (`temperature_K` null). All sample runs of one day use the same feed gas.

## ringdowns.csv

One row per averaged trace.

| column | meaning |
|---|---|
| `run` | run name from campaign.json |
| `step` | step number within the run |
| `wavenumber_cm` | laser wavenumber read by the wavemeter, vacuum, cm-1 (±0.001) |
| `a0` … | raw digitiser counts, dark level included |

The laser is stepped across a window around each scanned line. The cavity length is dithered by a piezo during averaging. The cavity is mode-matched to its TEM00 mode.

## o2_aband_lines.csv

The 16O2 lines of the b-X (0,0) band, in HITRAN conventions.

| column | meaning |
|---|---|
| `intensity_296K_cm_per_molecule` | intensity at 296 K, per molecule of O2 |
| `lower_energy_cm` | E'' |
| `gamma_air_cm_per_atm` | air-broadened half width at 296 K |
| `n_air` | its temperature exponent |

- The line shape is a Voigt profile.
- N2 broadening equals air broadening.
- Pressure shifts are negligible at these pressures.
- Other isotopologues have no lines inside the scan windows.

## o2_partition_sums.csv

Total internal partition sum Q(T) of 16O2, on the temperature scale used by the line intensities.
