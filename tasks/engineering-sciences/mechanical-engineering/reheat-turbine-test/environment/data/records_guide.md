# Test campaign records

Each folder under `campaigns/` is one HP/IP section test campaign on one reheat unit: several test points run on one day, at one load, with the turbine inlet valves held in one position. Every campaign folder holds three files.

## unit.json

- `station`, `unit`, `rating_mw`: identification.
- `site_elevation_m`: elevation of the turbine hall above sea level.
- `turbine`: plain-text description of the steam path, with these keys:
  - `casings`
  - `mid_span_packing`: what the shaft packing between the HP and IP inlets separates. Either `shaft packing between the HP inlet chamber and the IP inlet chamber` or `shaft packing between the HP control-stage shell and the IP inlet chamber`.
  - `hp_exhaust_end_packing`
  - `hp_admission`: `full arc: all inlet valves open together, no control stage` or `partial arc: nozzle-group control valves and a control stage`.
  - `reheat_spray`
  - `cold_reheat_thermowells`: `upstream of the reheat desuperheater` or `downstream of the reheat desuperheater`. Downstream means far enough downstream for the spray water to be fully mixed and evaporated.
  - `cold_reheat_line`
  - `superheater_spray`
- `instruments`: one entry per tag in `trend.csv`.
  - `tag`, `measures`, `location`, `unit`.
  - Steam pressure transmitters have a `sensing_line` with:
    - `filled`: `condensate` (the line is full of water at the turbine-hall air temperature, `TE-AMB`) or `steam` (the weight of the line contents is negligible);
    - `transmitter_below_tap_m`: height of the tapping point above the transmitter in metres (negative when the transmitter is above the tap).
  - Test transmitters and test thermocouples have a `certificate_correction`, in the unit of the reading, to be added to every reading.
  - `TI-IPX` is the station thermocouple. It has `grade` `station` and an `accuracy_K` instead of a certificate.
- Pressure units: `bar(a)` absolute; `bar(g)` relative to the atmosphere at the operating floor, which `PB` records in `kPa(a)`.

The two main steam legs and the two hot reheat legs carry equal flows. `FT-FW` is the final feedwater flow from the test flow nozzle, already in kg/s. Steam leaves the boiler only through the main steam line, and only the reheat spray adds water to the reheat system. Casings, valves, packings and pipes exchange no heat with their surroundings.

## trend.csv

One-minute averages for the whole test day. The first column is `time` (`HH:MM`, local time on the test date), followed by one column per tag. An empty field means no reading.

## test_log.txt

The test engineer's log. Each test point is on one line:

`Test <n> <HH:MM>-<HH:MM>: <description>.`

The window is the steady period of that test point; samples whose time lies inside it, both ends included, belong to it.

A `Remarks:` section may follow. Lines in these forms carry instrument information:

- `<tag> post-test check: reads <x> K high` (or `low`, or `bar high` / `bar low`): the error of that instrument during the whole campaign, left after its certificate correction.
- `<tag> out of service: <reason>` or `<tag> out of service from <HH:MM>: <reason>`: readings of that tag are not to be used (from that time on, where a time is given).

Other remark lines are free text.

## Steam properties

IAPWS-IF97.
