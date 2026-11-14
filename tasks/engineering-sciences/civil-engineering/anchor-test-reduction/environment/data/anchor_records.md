# Records of a ground-anchor acceptance test

Each anchor folder holds:
- `anchor.json`: the anchor, its tendon and its test loads;
- `instruments.json`: the instruments;
- `logger.csv`: the data logger;
- `dial_sheet.csv`: the dial gauge readings written down on site;
- `levels.csv`: the levelling of the reference beam;
- `events.csv`: the test log.

The example anchors `A-07` and `A-12` also hold `accepted.json`, with their two results.

## The test

The anchor is loaded in cycles by a hollow-ram jack bearing on the anchor head. The load goes up from the datum load
to each fraction of the proof load listed in `cycle_levels`, is held, and goes back to the datum load, where it is
held again.
- **Holds below the proof load:** 1 min at the top.
- **The proof load:** held for 15 min. Where `hold_extended` is logged, the hold lasts the number of minutes given.
- **After the proof load:** the datum load is held for 5 min before the test ends.

The anchor axis is inclined at `inclination_deg` below the horizontal. Displacements are head movements along that
axis, outward positive.

## anchor.json

- **`tendon`:** the number of strands, the cross-section of one strand, the strands' modulus, and the free and bond
  lengths of the design.
- **Loads:** `proof_load_kn` and `datum_load_kn`.
- **`cycle_levels`:** the cycle levels as fractions of the proof load.

## instruments.json

**`jack`.** The ram area, as the jack maker gives it. PT-1 measures the oil pressure in bar. The ram's seal friction
makes the force on the anchor smaller than pressure times ram area while the ram moves out, and larger by the same
share after it has moved back. Holds keep the share of the last movement.

**`jack_table`.** The factor, in kN per bar, that turns the PT-1 pressure into force. It was found by calibrating the
jack against the named proving ring on the named branch.

**`load_cell`.** LC-1 sits between the jack and the anchor head and was calibrated against the named proving ring.

**`displacement`.** LVDT-1 is an electrical transducer whose reading in mm rests on its own calibration; it is logged.
DG-1 is a mechanical dial gauge whose readings are written on the dial sheet; `hold_clock` says from which moment the
sheet's hold minutes count. Both gauges are clamped to the extension bar and bear on the anchor head.

**`extension_bar`.** A bar that runs along the anchor axis from the reference beam to the anchor head.
- Its material is steel, which lengthens by 12e-6 of its length per K.
- TB-1 is its thermometer.

**How the gauges read.** They measure the gap between the bar's head end and the anchor head along the axis; a
smaller gap reads as an outward head movement. Anything that carries the bar's head end towards the head therefore
reads as an outward head movement as well: the bar lengthening, and the beam end of the bar moving inward along the
axis.

**`reference_beam`.** It holds the bar's outer end. Its supports stand at the given distance from the head and are
levelled.

## logger.csv

Every 10 s: `PT-1_bar`, `LC-1_kN`, `LVDT-1_mm` and `TB-1_C`.

## dial_sheet.csv

One row per reading, with:
- the time it was read;
- the stage (`datum`, `peak` or `proof`);
- `hold_min`, the hold minutes by the clock that `instruments.json` names;
- `DG-1_mm`, to the gauge's resolution.

## levels.csv

Staff readings on the two beam supports, in mm, at the times listed. A staff reading grows as the point it stands on
sinks.

## events.csv

`time,event,value`. The events:
- `test_start`;
- `increment_start` (value: the load being approached, kN);
- `unload_start` (value: the datum load, kN);
- `hold_extended` (value: the proof hold's new length, min);
- `test_end`.

## Other anchors

Anchors beyond the two examples may contain events and instrument conditions that the examples do not show.
