import source.Obstrategy as obs

tmp = obs.Obstrategy("./data/test_list.dat")
tmp.make_script_grid("script/grid_plan.csv")
tmp.make_script_gb("script/gb_plan.csv",band=["H","J","Y"])
tmp.print_time()