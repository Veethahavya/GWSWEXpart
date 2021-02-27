import time
import numpy as np
import pandas as pd
import sys, os
import subprocess as sp
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from datetime import datetime, timedelta
import netCDF4 as nc
import logging
import shutil

logging.basicConfig(filename="../output/GWSWEX.log", format="%(asctime)s; %(levelname)s: %(message)s", level=logging.DEBUG)

#%%
class bundleIT(object):
	def __init__(self, **argd):
		self.__dict__.update(argd)
	def keys(self):
		for property, value in vars(self).items():
			print(property)
	def __main__(self):
		self.keys()


#%%
class Timing():
	def __init__(self, start, end, dt_exchange, ts_exchange):
		self.str = {}
		self.unix = {}
		self.rel = {}
		self.stmp = {}
		self.dt = {}
		self.nTS = {}
		self.ts = {}
		self.run = {}

		self.run["PET.prep"] = []
		self.run["PET.get"] = []
		self.run["Fort.build"] = []
		self.run["Fort.Run"] = []
		self.run["Fort.update"] = []
		self.run["Fort.load"] = []
		self.run["start"] = time.time()
		self.runtimes =[]

		self.str["global_start"] = start
		self.unix["global_start"] = time.mktime(datetime.strptime(self.str["global_start"], "%Y-%m-%d %H:%M:%S").timetuple())
		self.rel["global_start"] = 0.0
		self.stmp["global_start"] = datetime.strptime(self.str["global_start"], "%Y-%m-%d %H:%M:%S")
		self.str["global_end"] = end
		self.unix["global_end"] = time.mktime(datetime.strptime(self.str["global_end"], "%Y-%m-%d %H:%M:%S").timetuple())
		self.rel["global_end"] = self.unix["global_end"] - self.unix["global_start"]
		self.stmp["global_end"] = datetime.strptime(self.str["global_end"], "%Y-%m-%d %H:%M:%S")

		self.nTS["exchange"] = ts_exchange
		self.dt["exchange"] = dt_exchange
		self.dt["FortRun"] = self.dt["exchange"]/self.nTS["exchange"]
		self.nTS["ran"] = 0
		self.nTS["run_num"] = 0
		self.nTS["max"] = int(self.rel["global_end"]/self.dt["exchange"] - 1)

		self.str["local_start"] = start
		self.unix["local_start"] = time.mktime(datetime.strptime(self.str["local_start"], "%Y-%m-%d %H:%M:%S").timetuple())
		self.rel["local_start"] = self.unix["local_start"] - self.unix["global_start"]
		self.stmp["local_start"] = datetime.strptime(self.str["local_start"], "%Y-%m-%d %H:%M:%S")
		self.unix["local_end"] = self.unix["local_start"] + self.dt["exchange"]
		self.str["local_end"] = datetime.fromtimestamp(self.unix["local_end"]).strftime("%Y-%m-%d %H:%M:%S")
		self.rel["local_end"] = self.unix["local_end"] - self.unix["global_start"]
		self.stmp["local_end"] = datetime.strptime(self.str["local_end"], "%Y-%m-%d %H:%M:%S")

		self.ts["local"] = 0
		self.ts["global"] = self.ts["local"]

	def update(self, update_by=None, internal=False):
		self.run["start"] = time.time()
		if update_by is None:
			dt_update = self.dt["exchange"]
		else:
			dt_update = update_by
		self.unix["local_start"] = self.unix["local_start"] + dt_update
		self.str["local_start"] = datetime.fromtimestamp(self.unix["local_start"]).strftime("%Y-%m-%d %H:%M:%S")
		self.rel["local_start"] = self.unix["local_start"] - self.unix["global_start"]
		self.stmp["local_start"] = datetime.strptime(self.str["local_start"], "%Y-%m-%d %H:%M:%S")
		self.unix["local_end"] = self.unix["local_end"] + dt_update
		self.str["local_end"] = datetime.fromtimestamp(self.unix["local_end"]).strftime("%Y-%m-%d %H:%M:%S")
		self.rel["local_end"] = self.unix["local_end"] - self.unix["global_start"]
		self.stmp["local_end"] = datetime.strptime(self.str["local_end"], "%Y-%m-%d %H:%M:%S")
		self.ts["local"] = np.arange(self.rel["local_start"], self.rel["local_end"]+1, self.dt["FortRun"])
		self.ts["global"] = np.append(self.ts["global"], self.ts["local"][1:])
		if internal:
			logging.info("TIMING: Updated times to suit nTS change")
		else:
			logging.info("*TIMING: Updated. Start: {0}   End: {1}".format(self.str["local_start"], self.str["local_end"]))


#%%
class PET:
	def __init__(self, data_path, times, Fort, p_name="p.dat", et_name="et.dat"):
		self.data_path = data_path
		self.p_path = os.path.join(data_path, p_name)
		self.et_path = os.path.join(data_path, et_name)
		self.times = times
		self.Fort = Fort
		self.pDF = None
		self.etDF = None
		self.p = None
		self.et = None
		self.petDF_sel = None

	def prep(self):
		strt_time = time.time()
		pDF = pd.read_table(self.p_path)
		pDF[pDF.columns[0]] = pd.to_datetime(pDF[pDF.columns[0]])
		pDF.index = pDF[pDF.columns[0]]
		pDF = pDF.drop(pDF.columns[0], 1)
		etDF = pd.read_table(self.et_path)
		etDF[etDF.columns[0]] = pd.to_datetime(etDF[etDF.columns[0]])
		etDF.index = etDF[etDF.columns[0]]
		etDF = etDF.drop(etDF.columns[0], 1)
		start = datetime(self.times.stmp["global_start"].year, self.times.stmp["global_start"].month, self.times.stmp["global_start"].day,\
		self.times.stmp["global_start"].hour)
		end = datetime(self.times.stmp["global_end"].year, self.times.stmp["global_end"].month, self.times.stmp["global_end"].day,\
		self.times.stmp["global_end"].hour)
		p_mask = (pDF.index >= start) & (pDF.index <= end)
		et_mask = (etDF.index >= start) & (etDF.index <= end)
		pDF = pDF.loc[p_mask]
		etDF = etDF.loc[et_mask]
		self.pDF = pDF
		self.etDF = etDF
		self.times.run["PET.prep"].append(time.time() - strt_time)
		logging.info("PET: Prepped")

	def get(self, throttle=False, new_dt_exchange=None):
		strt_time = time.time()
		start = datetime(self.times.stmp["local_start"].year, self.times.stmp["local_start"].month, self.times.stmp["local_start"].day,\
		self.times.stmp["local_start"].hour)
		end = datetime(self.times.stmp["local_end"].year, self.times.stmp["local_end"].month, self.times.stmp["local_end"].day,\
		self.times.stmp["local_end"].hour)
		p_mask = (self.pDF.index >= start) & (self.pDF.index <= end)
		pDF_sel = self.pDF.loc[p_mask]
		p = pDF_sel[pDF_sel.columns[0]].to_numpy()[0]
		if throttle:
			if not new_dt_exchange:
				new_dt_exchange = self.times.dt["exchange"]
			if p.max() == 0:
				self.times.nTS["exchange"] = 1
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 0.1:
				self.times.nTS["exchange"] = 2
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 0.2:
				self.times.nTS["exchange"] = 4
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 0.3:
				self.times.nTS["exchange"] = 6
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 0.5:
				self.times.nTS["exchange"] = 8
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 1:
				self.times.nTS["exchange"] = 10
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			else:
				self.times.nTS["exchange"] = 15
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			self.times.dt["FortRun"] = self.times.dt["exchange"]/self.times.nTS["exchange"]
			self.times.update(update_by=self.times.dt["exchange"]-old_dt_exch, internal=True)
			logging.warning("PET: !Throttling! nTS: {0}. ts_exch: {1}s".format(self.times.nTS["exchange"], self.times.dt["exchange"]))
		p = np.repeat(p, (self.times.nTS["exchange"])+1)
		p = p/3600
		logging.info("PET: max P for this ts is {} mm/s".format(p.max()))
		et_mask = (self.etDF.index >= start) & (self.etDF.index <= end)
		etDF_sel = self.etDF.loc[et_mask]
		et = etDF_sel[etDF_sel.columns[0]].to_numpy()[0]
		et = np.repeat(et, (self.times.nTS["exchange"]/et.size)+1)
		et = et/3600
		logging.info("PET: max ET for this ts is {} mm/s".format(et.max()))
		if pDF_sel.index.equals(etDF_sel.index):
			petDF_sel = pDF_sel
			petDF_sel[etDF_sel.columns[0]] = etDF_sel[etDF_sel.columns[0]]
		else:
			logging.errorr("WARNING: P and ET indexes do not match")
			sys.exit("WARNING: P and ET indexes do not match")
		self.Fort.Ini.p = p
		self.Fort.Ini.et = et
		self.Fort.Ini.ts = p.size
		self.Fort.Ini.ts_size = self.times.dt["FortRun"]
		self.p = p
		self.et = et
		self.petDF_sel = petDF_sel
		self.times.run["PET.get"].append(time.time() - strt_time)
		logging.info("PET: Fetched and Passed P and ET values to FORT")


#%%
class Fort:
	def __init__(self, path, times, exe_path="GWSWEX.exe"):
		self.path = path
		self.exe_path = os.path.join(path, exe_path)
		self.times = times
		self.Ini = bundleIT(elems=None, ts=None, ts_size=None, n=None, m=None, beta=None, alpha=None, gok=None, k=None, sw_th=None,\
		p=None, et=None, gws=None, sws=None, epv=None, sm=None, chd=None)
		self.Res = bundleIT(gws=None, sws=None, sm=None, epv=None, Qin=None, Qout=None, Qdiff=None, gw_dis=None, sw_dis=None, sm_dis=None,\
		max_err_abs=None, max_err_perc=None)

	def build(self, run=False, restart=False, res=None):
		strt_time = time.time()
		if restart:
			fort_file = os.path.join(res.fort_path, "wasenmoos_fort_"+(self.times.stmp["local_start"]-timedelta(seconds=\
			self.times.dt["exchange"])).strftime("%Y%m%d_%H%M%S")+".npz")
			fort_rst_file = np.load(fort_file)
			self.Ini.sm = fort_rst_file["sm"][:,-1]
			self.Ini.epv = fort_rst_file["epv"][:,-1]
		self.Ini.chd = np.zeros(self.Ini.elems, dtype=np.int8)
		for elem in range(self.Ini.elems):
			if elem in self.Ini.chd_cells:
				self.Ini.chd[elem] = 1
			else:
				self.Ini.chd[elem] = 0
		missing = []
		for attr in self.Ini.__dict__.keys():
			if getattr(self.Ini, attr) is None:
				missing.append(attr)
		if len(missing) != 0:
			logging.error("FORT: Incomplete set of FortIni Objects. Missing values: ", missing[:])
			raise ValueError("Incomplete set of FortIni Objects. Missing values: ", missing[:])
		p_file = FortranFile(os.path.join(self.path,"p.ip"), "w")
		p_file.write_record(self.Ini.p)
		p_file.close()
		et_file = FortranFile(os.path.join(self.path,"et.ip"), "w")
		et_file.write_record(self.Ini.et)
		et_file.close()
		args = [self.Ini.elems, self.Ini.ts, self.Ini.ts_size, self.Ini.n, self.Ini.m, self.Ini.beta, self.Ini.alpha, self.Ini.sw_th]
		np.savetxt(os.path.join(self.path,"args.ip"), args)
		gws_ini_file = FortranFile(os.path.join(self.path,"gws_ini.ip"), "w")
		gws_ini_file.write_record(self.Ini.gws)
		gws_ini_file.close()
		sws_ini_file = FortranFile(os.path.join(self.path,"sws_ini.ip"), "w")
		sws_ini_file.write_record(self.Ini.sws)
		sws_ini_file.close()
		epv_ini_file = FortranFile(os.path.join(self.path,"epv_ini.ip"), "w")
		epv_ini_file.write_record(self.Ini.epv)
		epv_ini_file.close()
		sm_ini_file = FortranFile(os.path.join(self.path,"sm_ini.ip"), "w")
		sm_ini_file.write_record(self.Ini.sm)
		sm_ini_file.close()
		gok_file = FortranFile(os.path.join(self.path,"gok.ip"), "w")
		gok_file.write_record(self.Ini.gok)
		gok_file.close()
		k_file = FortranFile(os.path.join(self.path,"k.ip"), "w")
		k_file.write_record(self.Ini.k)
		k_file.close()
		chd_file = FortranFile(os.path.join(self.path,"chd.ip"), "w")
		chd_file.write_record(self.Ini.chd)
		chd_file.close()
		self.times.run["Fort.build"].append(time.time() - strt_time)
		logging.info("FORT: Built")
		if run:
			self.Run(load=True)

	def update(self, run=False):
		strt_time = time.time()
		self.Ini.sws = self.Res.sws[:,-1]
		self.Ini.gws = self.Res.gws[:,-1]
		self.Ini.sm = self.Res.sm[:,-1]
		self.Ini.epv = self.Res.epv[:,-1]
		p_file = FortranFile(os.path.join(self.path,"p.ip"), "w")
		p_file.write_record(self.Ini.p)
		p_file.close()
		et_file = FortranFile(os.path.join(self.path,"et.ip"), "w")
		et_file.write_record(self.Ini.et)
		et_file.close()
		args = [self.Ini.elems, self.Ini.ts, self.Ini.ts_size, self.Ini.n, self.Ini.m, self.Ini.beta, self.Ini.alpha, self.Ini.sw_th]
		np.savetxt(os.path.join(self.path,"args.ip"), args)
		gws_ini_file = FortranFile(os.path.join(self.path,"gws_ini.ip"), "w")
		gws_ini_file.write_record(self.Ini.gws)
		gws_ini_file.close()
		sws_ini_file = FortranFile(os.path.join(self.path,"sws_ini.ip"), "w")
		sws_ini_file.write_record(self.Ini.sws)
		sws_ini_file.close()
		epv_ini_file = FortranFile(os.path.join(self.path,"epv_ini.ip"), "w")
		epv_ini_file.write_record(self.Ini.epv)
		epv_ini_file.close()
		sm_ini_file = FortranFile(os.path.join(self.path,"sm_ini.ip"), "w")
		sm_ini_file.write_record(self.Ini.sm)
		sm_ini_file.close()
		self.Ini.chd = np.zeros(self.Ini.elems, dtype=np.int8)
		for elem in range(self.Ini.elems):
			if elem in self.Ini.chd_cells:
				self.Ini.chd[elem] = 1
			else:
				self.Ini.chd[elem] = 0
		chd_file = FortranFile(os.path.join(self.path,"chd.ip"), "w")
		chd_file.write_record(self.Ini.chd)
		chd_file.close()
		self.times.run["Fort.update"].append(time.time() - strt_time)
		logging.info("FORT: Updated")
		if run:
			self.Run(load=True)

	def Run(self, load=False):
		strt_time = time.time()
		logging.info("FORT: Run initialized")
		owd = os.getcwd()
		if os.path.exists(self.exe_path):
			dir_path = os.path.split(self.exe_path)[0]
			exe = os.path.split(self.exe_path)[1]
			os.chdir(dir_path)
			fort_run = sp.Popen(exe, shell=True, stdout = sp.PIPE)
			stdout_fort, stderr_fort = fort_run.communicate()
			if not fort_run.returncode == 0:
				logging.error("FORT: Failed!")
				raise RuntimeError("FortRun failed")
			else:
				Fort.success = True
		else:
			logging.error("FORT: executable not found")
			raise NameError("FortRun executable missing")
		os.chdir(owd)
		self.times.run["Fort.Run"].append(time.time() - strt_time)
		logging.info("FORT: Run complete")
		if load:
			self.load()

	def load(self):
		strt_time = time.time()
		shape = (self.Ini.elems, self.Ini.ts)
		Qin_file = FortranFile(os.path.join(self.path,"Qin.op"), "r")
		self.Res.Qin = (Qin_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		Qin_file.close()
		Qout_file = FortranFile(os.path.join(self.path,"Qout.op"), "r")
		self.Res.Qout = (Qout_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		Qout_file.close()
		Qdiff_file = FortranFile(os.path.join(self.path,"Qdiff.op"), "r")
		self.Res.Qdiff = (Qdiff_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		Qdiff_file.close()
		gws_file = FortranFile(os.path.join(self.path,"gws.op"), "r")
		self.Res.gws = gws_file.read_reals().reshape(shape, order="F")
		gws_file.close()
		sws_file = FortranFile(os.path.join(self.path,"sws.op"), "r")
		self.Res.sws = sws_file.read_reals().reshape(shape, order="F")
		sws_file.close()
		sm_file = FortranFile(os.path.join(self.path,"sm.op"), "r")
		self.Res.sm = sm_file.read_reals().reshape(shape, order="F")
		sm_file.close()
		epv_file = FortranFile(os.path.join(self.path,"epv.op"), "r")
		self.Res.epv = epv_file.read_reals().reshape(shape, order="F")
		epv_file.close()
		gw_dis_file = FortranFile(os.path.join(self.path,"gw_dis.op"), "r")
		self.Res.gw_dis = (gw_dis_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		gw_dis_file.close()
		sw_dis_file = FortranFile(os.path.join(self.path,"sw_dis.op"), "r")
		self.Res.sw_dis = (sw_dis_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		sw_dis_file.close()
		sm_dis_file = FortranFile(os.path.join(self.path,"sm_dis.op"), "r")
		self.Res.sm_dis = (sm_dis_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		sm_dis_file.close()
		self.Res.max_err_abs = np.amax(self.Res.Qdiff.sum(1))
		self.times.run["Fort.load"].append(time.time() - strt_time)
		logging.info("FORT: Loaded results")

	def plot(self, elem, plotPrec=False, plotDis=False, savefig=False, dDPI=90, pDPI=1600, alpha_scatter=0.7, scatter_size=3, fromMF6=None):
		pal = ["#E74C3C", "#2ECC71", "#5EFCA1", "#E7AC3C", "#2980B9", "#1A3BFF", "#FF6600"] #[gw, sm, epv, sv, sw, p, et]

		if plotDis:
			plt.figure(dpi=dDPI)
			plt.xlabel("Time Steps")
			plt.ylabel("Discharges in Storage")
			plt.scatter(range(0,self.Ini.ts), self.Res.gw_dis[elem,:], label="GW_dis", color=pal[0],\
			alpha=alpha_scatter, s=scatter_size)
			plt.scatter(range(0,self.Ini.ts), self.Res.sm_dis[elem,:], label="SM_dis", color=pal[1],\
			alpha=alpha_scatter, s=scatter_size)
			plt.scatter(range(0,self.Ini.ts), self.Res.sw_dis[elem,:], label="SW_dis", color=pal[4],\
			alpha=alpha_scatter, s=scatter_size)
			plt.legend(loc="best", fontsize="small")
			if savefig:
				plt.savefig("prec_dist.jpg", format="jpeg", dpi=pDPI)

		plt.figure(dpi=dDPI)
		plt.xlabel("Time Steps")
		plt.ylabel("Water Levels in mm")
		plt.ylim([self.Res.gws[elem,:].min()-50, self.Res.sws[elem,:].max()+50+self.Ini.gok[elem]])
		plt.stackplot(range(0,self.Ini.ts), self.Res.gws[elem,:], self.Res.sm[elem,:],\
		self. Res.epv[elem,:]-self.Res.sm[elem,:], (np.full(self.Ini.ts,self.Ini.gok[elem])-self.Res.gws[elem,:])*(1-self.Ini.n),\
		self.Res.sws[elem,:], labels=["Groundwater","Soil Moisture", "Field Capacity", "Soil Volume", "Surface Water"], colors=pal)
		plt.plot(range(0,self.Ini.ts), np.full(self.Ini.ts,self.Ini.gok[elem]), "k", linewidth=0.5, label="Ground Level")
		plt.legend(loc="best", fontsize="small")
		if savefig:
			plt.savefig("water_levels.jpg", format="jpeg", dpi=pDPI)

		if plotPrec:
			plt.figure(dpi=dDPI)
			plt.xlabel("Time Steps")
			plt.ylabel("Precipitation in mm")
			plt.scatter(range(0,self.Ini.ts), self.Ini.p, s=0.1, color=pal[6])
			if savefig:
				plt.savefig("prec.jpg", format="jpeg", dpi=pDPI)

	def delIP(self):
		dat_files = [f for f in os.listdir(self.path) if f.endswith(".ip") ]
		for f in dat_files:
			os.remove(os.path.join(self.path, f))
	def delOP(self):
		dat_files = [f for f in os.listdir(self.path) if f.endswith(".op") ]
		for f in dat_files:
			os.remove(os.path.join(self.path, f))


#%%
class ResNC:
	def __init__(self, times, Fort, op_path="../output/", op_name="GWSWEX"):
		self.path = op_path
		if not os.path.exists(self.path):
			os.mkdir(op_path)
		self.delft_path = os.path.join(self.path, "delft/")
		self.mf_path = os.path.join(self.path, "modflow/")
		self.fort_path = os.path.join(self.path, "fortran/")
		self.name = op_name
		self.times = times
		self.Fort = Fort
		self.save = None

	def build(self):
		self.OP = nc.Dataset(os.path.join(self.path, self.name+".nc"), "w")
		self.OP.source = "GWSWEX: Groundwater-Surfacewater Exchange Model"
		self.OP.createDimension("element", self.Fort.Ini.elems)
		self.OP.createDimension("time", None)

		timeV = self.OP.createVariable("time", "i", ("time", ))
		timeV.units = "seconds since " + self.times.str["global_start"]
		time_unixV = self.OP.createVariable("time_unix", "i", ("time", ))
		time_unixV.units = "unix time, seconds"

		pV = self.OP.createVariable("p", "d", ("time", ), fill_value=np.nan)
		pV.long_name = "Precipitation rate"
		pV.units = "mm/s"
		pV[0] = 0
		eV = self.OP.createVariable("et", "d", ("time", ), fill_value=np.nan)
		eV.long_name = "Evapotranspiration rate"
		eV.units = "mm/s"
		eV[0] = 0

		gwsV = self.OP.createVariable("gws", "d", ("time", "element"), fill_value=np.nan)
		gwsV.long_name = "Groundwater Levels (Layer 1: Unsaturated) from reference sea level"
		gwsV.units = "m"
		swsV = self.OP.createVariable("sws", "d", ("time", "element"), fill_value=np.nan)
		swsV.long_name = "Surfacewater Levels above ground level"
		swsV.units = "m"
		smV = self.OP.createVariable("sm", "d", ("time", "element"), fill_value=np.nan)
		smV.long_name = "Soil Moisture"
		smV.units = "mm"
		epvV = self.OP.createVariable("epv", "d", ("time", "element"), fill_value=np.nan)
		epvV.long_name = "Field Capacity"
		epvV.units = "mm"

		fort_residualV = self.OP.createVariable("fort_residual", "d", ("time", ), fill_value=np.nan)
		fort_residualV.long_name = "Residual from FortRun"
		fort_residualV.units = "mm"

		sw_flux_fortV = self.OP.createVariable("sw_dis_fort", "d", ("time", "element"), fill_value=np.nan)
		sw_flux_fortV.long_name = "SW Storage Discharge as reported by Fort"
		sw_flux_fortV.standard_name = "SW Discharge: Fort"
		sw_flux_fortV.units = "m3/dt"
		sw_flux_fortV[0] = 0
		gw_flux_fortV = self.OP.createVariable("gw_dis_fort", "d", ("time", "element"), fill_value=np.nan)
		gw_flux_fortV.long_name = "GW Storage Discharge as reported by Fort"
		gw_flux_fortV.standard_name = "GW Discharge: Fort"
		gw_flux_fortV.units = "m3/dt"
		gw_flux_fortV[0] = 0
		sm_fluxV = self.OP.createVariable("sm_dis", "d", ("time", "element"), fill_value=np.nan)
		sm_fluxV.long_name = "Soil Moisture Discharge"
		sm_fluxV.units = "m3/dt"
		sm_fluxV[0] = 0

		timeV[0] = 0
		time_unixV[0] = self.times.unix["global_start"]
		gwsV[0] = self.Fort.Ini.gws
		swsV[0] = self.Fort.Ini.sws
		smV[0] = self.Fort.Ini.sm
		epvV[0] = self.Fort.Ini.epv

		logging.info("RES: NC file initialized")

	def update(self, save=False):
		self.save = save
		idx = self.times.nTS["run_num"] + 1
		self.OP["time"][idx] = self.times.rel["local_end"]
		self.OP["time_unix"][idx] = self.times.unix["local_end"]
		self.OP["gws"][idx] = self.Fort.Res.gws[:,-1]
		self.OP["sws"][idx] = self.Fort.Res.sws[:,-1]
		self.OP["sm"][idx] = self.Fort.Res.sm[:,-1].T
		self.OP["epv"][idx] = self.Fort.Res.epv[:,-1].T
		self.OP["fort_residual"][idx] = np.abs(self.Fort.Res.Qdiff).sum()
		self.OP["p"][idx] = self.Fort.Ini.p.mean()
		self.OP["et"][idx] = self.Fort.Ini.et.mean()
		logging.info("*RES: Wrote OP to NC file\n")

		if save:
			self.saveOPs()
		else:
			self.times.runtimes.append(time.time() - self.times.run["start"])
			logging.info("*RUNTIME{}: {}\n\n".format(self.times.nTS["run_num"], self.times.runtimes[int(self.times.nTS["run_num"]-1)]))
		self.times.nTS["ran"] = self.times.nTS["ran"] + self.times.nTS["exchange"]

		self.times.nTS["run_num"] = self.times.nTS["run_num"] + 1

	def saveOPs(self):
		if not os.path.exists(self.fort_path):
			os.mkdir(self.fort_path)
		fort_file = os.path.join(self.fort_path, "wasenmoos_fort_"+self.times.stmp["local_start"].strftime("%Y%m%d_%H%M%S"))
		np.savez(fort_file, gws=self.Fort.Res.gws, sws=self.Fort.Res.sws, sm=self.Fort.Res.sm, epv=self.Fort.Res.epv,\
		Qin=self.Fort.Res.Qin, Qout=self.Fort.Res.Qout, Qdiff=self.Fort.Res.Qdiff,\
		gw_dis=self.Fort.Res.gw_dis, sw_dis=self.Fort.Res.sw_dis, sm_dis=self.Fort.Res.sm_dis)
		logging.info("*RES: Dumped FORT outputs to output folder")
		self.times.runtimes.append(time.time() - self.times.run["start"])
		logging.info("*RUNTIME{}: {}\n\n".format(self.times.nTS["run_num"], self.times.runtimes[int(self.times.nTS["run_num"]-1)]))

	def close(self):
		self.OP.close()
		logging.info("*RES: Closed NC file\nMaking a copy of the log")
		shutil.copy(os.path.join(self.path,"GWSWEX.log"), os.path.join(self.path,"GWSWEX_"+self.times.stmp["global_start"].\
		strftime("%Y%m%d_%H%M%S")+"-"+self.times.stmp["local_end"].strftime("%Y%m%d_%H%M%S")+".log"))
		logging.info("*LOG END*\n\n\n")
		logging.shutdown()

	def loadNC(self, nc_path=None):
		if not nc_path:
			nc_path = os.path.join(self.path, self.name+".nc")
		self.OP = nc.Dataset(nc_path)
		return self.OP


	def loadFort(self, load_time=None):
		if not load_time:
			load_time = (self.times.stmp["local_start"] - timedelta(seconds=self.times.dt["exchange"])).strftime("%Y%m%d_%H%M%S")
		fort_path = os.path.join(self.fort_path, "wasenmoos_fort_"+load_time+".npz")
		fort_file = np.load(fort_path)
		self.Fort.Res.gws = fort_file["gws"]
		self.Fort.Res.gw_dis = fort_file["gw_dis"]
		self.Fort.Res.sws = fort_file["sws"]
		self.Fort.Res.sw_dis = fort_file["sw_dis"]
		self.Fort.Res.sm = fort_file["sm"]
		self.Fort.Res.sm_dis = fort_file["sm_dis"]
		self.Fort.Res.epv = fort_file["epv"]
		self.Fort.Res.Qin = fort_file["Qin"]
		self.Fort.Res.Qout = fort_file["Qout"]
		self.Fort.Res.Qdiff = fort_file["Qdiff"]

	def plot1D(self, data, fromNC=None, items=None, labels=None, title=None, to_idx=None):
		pal = {"gw":"#E74C3C", "sm":"#2ECC71", "fc":"#5EFCA1", "sv":"#E7AC3C", "sw":"#2980B9", "p":"#1A3BFF", "et":"#FF6600"\
		, "red":"red", "green":"green", "grey":"grey", "cyan":"c", "pink":"m", "default":"black"}
		if not to_idx and type(data) != str:
			to_idx = data[0].size
		else:
			to_idx = self.OP["time_unix"][:].size
		fig, ax = plt.subplots()
		ax.ticklabel_format(useOffset=False)
		if title:
			fig.suptitle(title)
		c = 0
		posix_times = self.OP["time_unix"][:]
		x = [datetime.fromtimestamp(x) for x in posix_times]
		if items and not labels:
			labels = items
		if not items:
			items = ["default", ]
			labels = [data, ]
		if fromNC == "key":
			data = [self.OP[data][:], ]
		elif fromNC == "sum":
			data = [self.OP[data][:].sum(axis=1), ]
		elif fromNC == "mean":
			data = [self.OP[data][:].mean(axis=1), ]
		for lin in data:
			plt.plot(x[:to_idx], lin[:to_idx], color=pal[items[c]], label=labels[c])
			plt.legend()
			c += 1