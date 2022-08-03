# GWSWEX [obselete readme; to be updated from private repo]
A 0-dimensional model to facilitate exchange between groundwater and surface-water models. This is a part of the module forked and trimmed from [here](https://github.com/Veethahavya/GWSWEX).
An visual illustration of the numerical scheme of this model is as follows:
[obselete illustration; to be updated from private repo]
![0D Cell](https://user-images.githubusercontent.com/14050804/109405064-c85e8a80-796c-11eb-953f-9ba27bfe245e.png)
The model itself is written in Fortran and is compiled as an executable (windows) found [here](https://github.com/Veethahavya/GWSWEXpart/blob/main/fortran/).
[This](https://github.com/Veethahavya/GWSWEXpart/blob/main/py/GWSWEX.py) is a python wrapper module to handle the model and ensure ease of use.

The pivotal class Fort helps initialize, run, and load the results of the model. 
Since the model uses an explicitly coupling numerical scheme, a robust timing control is required to synchronize the three models between global and descrete local time periods. This is implemented in the timing class of the python module.
The class ResNC manages the accumulation of results for the global simulation period and stores it into a netCDF folder. This class also offers loading and visualization options.

The main python executable to define, initialize, and run the model is found [here](https://github.com/Veethahavya/GWSWEXpart/blob/main/py/GWSWEX_run.py).

## Dependencies
* numpy
* pandas
* netCDF4
* matplotlib
* scipy

## Getting started
* Pull the repositry
* Install any missing dependencies
* To modify the model input variables and parametres open the file [GWSWEX_run.py](https://github.com/Veethahavya/GWSWEXpart/blob/main/py/GWSWEX_run.py)
* Make changes if desired
* Run the file
* Use hints in the main executable file to view results
