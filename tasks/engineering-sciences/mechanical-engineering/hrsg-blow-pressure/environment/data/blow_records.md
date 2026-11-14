# Trial-blow records

Each folder under `blows/` is one HRSG high-pressure circuit after chemical cleaning, with the records of its trial steam blows. During a blow, steam leaves the HP drum, passes SH1, the attemperator, SH2 and the main steam line, and is discharged to atmosphere through a temporary pipe and a silencer. Pressures are absolute.

## circuit.json

- `unit`, `title`.
- `drum`: inner diameter and length of the horizontal HP drum, and the level reference.
- `gas_side`: exhaust-gas flow and specific heat at the blowing gas-turbine load, the gas temperature entering SH2, and the order in which the gas meets the banks.
- `steam_path`: the sections in flow order:
  - `section`, `heated`, `description` (for heated banks, also how the bank exchanges heat with the gas);
  - for SH1, SH2 and MS, the maximum continuous rating (MCR) duty from the HRSG heat balance: `mcr_flow_kg_s` and the steam state (`P_bar`, `T_C`) at the section inlet and outlet. A `T_C` of null means saturated steam.
- `instruments`: one entry per column of `trial_blows.csv`, with what it measures and where.

## trial_blows.csv

One row every 10 s of steady blowing:
- `trial`: trial number;
- `t_s`: seconds from the start of steady blowing in that trial;
- one column per instrument listed in `circuit.json`.

## commissioning_notes.txt

The commissioning engineer's notes on the trials: dates, gas-turbine load, plant line-up and the drum pressure held in each trial.
