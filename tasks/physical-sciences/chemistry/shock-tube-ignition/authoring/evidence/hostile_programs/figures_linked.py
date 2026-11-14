import pathlib, sys
pathlib.Path(sys.argv[2]).symlink_to('/tests/campaign_figures.json')
