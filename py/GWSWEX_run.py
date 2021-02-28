import GWSWEX
import numpy as np
from tqdm import tqdm

#%% Initialization of the class/objects
times = GWSWEX.Timing("2020-05-01 00:00:00", "2020-06-01 00:00:00", 3600, 15) # Hint: times are stored as attributes. try viewing times.str
Fort = GWSWEX.Fort("../fortran/", times)
PET = GWSWEX.PET("../data/", times, Fort)
RES = GWSWEX.ResNC(times, Fort)

#%% Initialization of the GWSWEX model
PET.prep()
PET.get(throttle=True)

Fort.Ini.elems = 2*2 # number of elements; here a simple example of a 2*2 grid i.e. 4 elements

Fort.Ini.k = np.full(Fort.Ini.elems, 1e-3) # hydraulic conductivity
Fort.Ini.n = 0.4 # porosity
Fort.Ini.m = 0.1 # pwp
Fort.Ini.beta = 0.95 # model parameter
Fort.Ini.alpha = 0.90 # model parameter
Fort.Ini.sw_th = 0.1 # model parameter

Fort.Ini.gok = np.array([10, 11, 10.5, 11.5]) # the ground-surface elevation
Fort.Ini.chd_cells = np.array([0, 1, 0, 0]) # array indicating wheather the element is imposed by a constant-head boundary condition
Fort.Ini.gws = Fort.Ini.gok - 3 # initial groundwater elevation (here 3 under the ground surface)
Fort.Ini.epv = (Fort.Ini.gok - Fort.Ini.gws)*Fort.Ini.n  # initial effective pore volume
Fort.Ini.sm = Fort.Ini.epv*0.75 # initial soil moisture
Fort.Ini.sws = np.array([0.2, 0.3, 0, 0]) # initial surface water content

RES.build()

# Hint: Use the keys() function to see all the attributes

#%%
Fort.build()
Fort.Run() # Hint: the success state is stored in Fort.success
Fort.load() # Hint: Results are stored in Fort.Res
# Fort.plot(1)

#%%
RES.update(save=True)

#%%
for ts in tqdm(range(times.nTS["max"])):
	times.update()
	PET.get(throttle=True)
	Fort.update(run=True) # updates, runs, and loads
	RES.update(save=True)
RES.close()

#%% Visualization hints
# load the results netcdf file after the model run
nc = RES.loadNC()

# plot the desired data
RES.plot1D("gws", fromNC="mean") # mean of groundwater levels along the timesteps

# view the stored variables in the file to plot more. use:
nc.variables.keys()