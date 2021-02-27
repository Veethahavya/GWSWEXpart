import GWSWEX
import numpy as np
from tqdm import tqdm

#%%
times = GWSWEX.timing("2020-05-01 00:00:00", "2020-06-01 00:00:00", 3600, 15)

Fort = GWSWEX.Fort("../fortran/", times)

PET = GWSWEX.PET("../data/", times, Fort)

RES = GWSWEX.resNC(times, Fort)

#%%
PET.prep()
PET.get(throttle=True)

Fort.Ini.elems = 2*2

Fort.Ini.k = np.full(Fort.Ini.elems, 1e-3)
Fort.Ini.n = 0.4
Fort.Ini.m = 0.1
Fort.Ini.beta = 0.95
Fort.Ini.alpha = 0.90
Fort.Ini.sw_th = 0.1

Fort.Ini.gok = np.array([10, 11, 10.5, 11.5])
Fort.Ini.chd_cells = np.array([0, 1, 0, 0])
Fort.Ini.gws = Fort.Ini.gok - 3
Fort.Ini.epv = (Fort.Ini.gok - Fort.Ini.gws)*Fort.Ini.n
Fort.Ini.sm = Fort.Ini.epv*0.75
Fort.Ini.sws = np.array([0.2, 0.3, 0, 0])

RES.build()

#%%
Fort.build(restart=False, res=RES, run=True)
# Fort.plot(0)

#%%
RES.update(save=True)

#%%
for ts in tqdm(range(times.nTS["max"])):
	times.update()
	PET.get(throttle=True)
	Fort.update(run=True)
	RES.update(save=True)
RES.close()