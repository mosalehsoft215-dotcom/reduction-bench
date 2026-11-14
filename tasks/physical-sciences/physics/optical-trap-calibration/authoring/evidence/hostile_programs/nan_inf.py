import pathlib, sys
pathlib.Path(sys.argv[-1]).write_text('{"stiffness_per_watt_pn_nm_w": NaN, "reference_stiffness_pn_nm": 0.05, "reference_sensitivity_nm_v": Infinity, "bead_radius_nm": 500.0}')
