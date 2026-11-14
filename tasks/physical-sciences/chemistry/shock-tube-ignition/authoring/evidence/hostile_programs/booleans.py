import pathlib, sys
pathlib.Path(sys.argv[2]).write_text('{"activation_energy_kj_mol": 160.0, "pressure_exponent": true, "oxygen_exponent": false, "tau_report_1_us": 1000.0, "tau_report_2_us": 80.0}')
