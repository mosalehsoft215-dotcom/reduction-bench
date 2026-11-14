# Trial-blow records

Each folder under `blows/` is one HRSG high-pressure circuit after chemical cleaning, with the records of its trial steam blows. During a blow, steam leaves the HP drum, passes SH1, the attemperator, SH2 and the main steam line, and is discharged to atmosphere through a temporary pipe, an exit nozzle and a silencer. Pressures are absolute.

Every unit folder holds the four files below and uses only the fields, values and event types listed here. A unit may use any of them, including ones the example units do not show.

## circuit.json

- `unit`, `title`.
- `drum`: inner diameter and length of the horizontal HP drum, and the level reference.
- `gas_side`: exhaust-gas flow and specific heat at the blowing gas-turbine load, the gas temperature entering SH2, and the order in which the gas meets the banks.
- `steam_path`: the sections in flow order:
  - `section`, `heated`, `description` (for heated banks, also how the bank exchanges heat with the gas);
  - for SH1, SH2 and MS, the maximum continuous rating (MCR) duty from the HRSG heat balance: `mcr_flow_kg_s` and the steam state (`P_bar`, `T_C`) at the section inlet and outlet. A `T_C` of null means saturated steam.
- `mcr_note`: how the MCR state is taken between the section ends.
- `instruments`: one entry per column of `trial_blows.csv`, with what it measures and where, in words.
- `instrument_setup`: the same instrument facts in fixed form:
  - `drum_pressure` and `sh1_outlet_pressure`: each `{"mounting": "close_coupled"}`, or `{"mounting": "below_tapping", "depth_m": h, "impulse_line": "condensate", "leg_C": t, "elevation_correction": c}`: the transmitter sits h metres below its tapping on a line full of condensate at t. With `c` false the DCS applies no correction, so the reading includes the head of that condensate column; with `c` true the DCS already subtracts it.
  - `feedwater_nozzle`: `transmitters` (`FT-A`, `FT-B`, `FT-C`, all on one flow nozzle upstream of the spray take-off), their `columns`, and `density`, which is either `{"basis": "live"}` (the DCS uses the density at the measured feedwater temperature and pressure) or `{"basis": "fixed", "P_bar": p, "T_C": t}` (the DCS uses the one density of feedwater at p and t). This is the configuration on the instrument data sheet; a `density_compensation` event in `shift_log.csv` records a later change for the transmitters it lists.
  - `feedwater_clamp_on`: the `column` of the temporary clamp-on meter, its `density` (always live), `uncertainty_pct_of_reading`, and `location`: `upstream_of_spray_takeoff` (between the flow nozzle and the spray take-off, so it carries the same flow as the nozzle) or `downstream_of_spray_takeoff` (it carries only the feedwater going on to the economiser and drum; the spray branch, including whatever passes the shut spray valve, is not in it).
  - `drum_level`: `selected_at_start` (the transmitter the DCS level signal comes from at the start of every trial) and `transmitters`, one entry each for `LT-A` and `LT-B` (one selected at a time in the DCS). Each gives the tap positions `tap_above_centreline_mm` and `tap_below_centreline_mm`, `reference_leg_C` (temperature of the water in its reference leg) and `pressure_compensation`. With `pressure_compensation` false the DCS converts with the saturated water and steam densities at `ranged_for_saturation_at_bar` and applies no other conversion; with true it converts with the saturated water and steam densities at the measured drum pressure (and the reference-leg density at `reference_leg_C`). The two transmitters need not agree in zero.
  - `spray_flow`: the spray meter's `element`, `range_kg_s`, `square_root_extraction`, `low_flow_cutoff_kg_s` (the DCS shows 0 for any flow below it) and `uncertainty_kg_s` (above the cut-off the reading is the flow through the orifice within this).
  - `blowdown_flow`: `source`, here `calculated_from_valve_position` (the DCS figure is calculated from the CBD control-valve position and the drum pressure; it is not a flow measurement).
  - `exit_nozzle`: `discharge` (to atmosphere through the silencer) and `isentropic_exponent` (treat the nozzle as an ideal-gas nozzle with this exponent).

## trial_blows.csv

One row every 10 s of steady blowing:
- `trial`: trial number;
- `t_s`: seconds from the start of steady blowing in that trial;
- one column per instrument listed in `circuit.json`. `L_drum_mm` is the level of whichever transmitter the DCS has selected at that moment, as that transmitter shows it; every trial starts with the transmitter `selected_at_start`.

## shift_log.csv

The shift log's plant and instrument entries for the trials, one per row. Columns: `event`, `trial` (blank when it applies to all trials), `t_s` (seconds into that trial; blank when not timed), `items` (instrument tags or columns, separated by `;`), `value`. The event types:

| event | trial, t_s | items | value | meaning |
|---|---|---|---|---|
| `barometric_pressure` | blank | blank | bar | barometric pressure at site during the trials |
| `clamp_on_fitted` | blank | column | blank | the temporary clamp-on feedwater meter was in service for the trials |
| `impulse_lines_moved` | blank | transmitter tags | text | the impulse lines of these transmitters were re-terminated before the trials (the value says where) |
| `density_compensation` | blank | transmitter tags | `live` | the DCS was changed before the trials to compute the listed transmitters' mass flow with the live water density, whatever the data sheet in `instrument_setup` says; transmitters not listed keep the data-sheet configuration |
| `blowdown_isolated` | first trial affected, 0 | `CBD` | blank | the continuous-blowdown isolation valve was shut before this trial; from this trial on no blowdown leaves the drum, while the DCS keeps calculating `F_cbd_kg_s` from the control-valve position |
| `level_selection` | trial, t_s | `LT-A>LT-B` or `LT-B>LT-A` | blank | the DCS level selection changed from the first to the second transmitter at this time; a trial may have several |

## commissioning_notes.txt

The commissioning engineer's notes on the trials in prose: dates, gas-turbine load, barometric pressure, plant line-up, the shift log's entries in words, and the drum pressure held in each trial.
