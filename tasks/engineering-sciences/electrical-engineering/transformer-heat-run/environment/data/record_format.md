# Temperature-rise test records

Each folder under `heat_runs/` is the record of one temperature-rise test of one oil-immersed power transformer, made by the short-circuit method: total losses first, then the test current for one hour, then switch-off and winding-resistance readings. All clocks in a record (the logger, the resistance meter and the test engineer) are the same clock. A folder holds four files.

## transformer.json

- `works`, `serial`: identification.
- `rated_power_MVA`, `rated_voltage_kV` (`HV`, `LV`), `rated_current_HV_A`.
- `cooling`: `ONAN`, `ONAF` or `ODAF`.
- `cooler_fans`: `null` for units without fans, otherwise a sentence giving the direction in which the fans blow air through the radiators or coolers.
- `windings`: for `HV` and `LV`, `conductor` is `copper` or `aluminium`.
- `guaranteed_losses_kW`: `no_load`, `load_at_75C` and `total` (rated total losses).
- `sensors`: one entry per temperature sensor in `heat_run.csv`:
  - `id`, `measures` (`top oil`, `bottom oil`, `cooler inlet oil` or `ambient air`), `where` (mounting and position in words);
  - `grade`: `test` (calibrated test sensor, with a `certificate_correction_K` to be added to every reading), `check` (a contact thermocouple on a pipe surface, under insulation) or `station` (an instrument of the plant or test bay, with an `accuracy_K`).
  - Ambient `test` sensors sit in oil cups whose time constant matches the transformer's.

## heat_run.csv

One row per minute from before the supply is switched on until 25 minutes after switch-off; each row is the average over the minute that starts at the stated `time` (`HH:MM`). Columns: `P_kW` (power fed into the transformer), `I_A` (HV line current), then one column per sensor id. `OTI` is read to 0.5 K.

## winding_resistance.csv

The readings of the winding-resistance meter after switch-off, both windings on the same terminals as the cold measurement (HV H1-H2, LV X1-X2): `time` (`HH:MM:SS`), `winding` (`HV` or `LV`), `current_A` (the meter's measuring current), `voltage_mV` (voltage across the winding terminals) and `resistance_mOhm` (the meter's quotient of the two).

## test_record.txt

The test engineer's record, with lines in these forms:

- `Cold resistance at HH:MM, ...` followed by `  top oil (TOT-1) reads <x> C, bottom oil (BOT) reads <y> C` (instrument readings, before certificate correction) and `  HV (line-to-line, H1-H2) <r> mOhm, LV (line-to-line, X1-X2) <r> mOhm`.
- `HH:MM  supply on; ...`, `HH:MM  end of total-loss period; ...`, `HH:MM:SS  supply off; ...`.
- `Test current in the last hour: <f> of rated HV current.`
- An optional `Remarks:` section in free text.

## Corrections to rated conditions

Oil temperature rises scale with the losses as (rated total losses / test losses)^x and winding gradients with the current as (rated current / test current)^y, with:

| cooling | x | y |
|---|---|---|
| ONAN | 0.8 | 1.6 |
| ONAF | 0.9 | 1.6 |
| ODAF | 1.0 | 2.0 |
