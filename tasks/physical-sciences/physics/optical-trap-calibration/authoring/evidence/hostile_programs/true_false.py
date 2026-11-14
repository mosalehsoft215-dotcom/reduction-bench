import pathlib, sys
pathlib.Path(sys.argv[-1]).write_text('{"stiffness_per_watt_pn_nm_w": 0.25, "reference_stiffness_pn_nm": true, "reference_sensitivity_nm_v": 120.0, "bead_radius_nm": false}')
