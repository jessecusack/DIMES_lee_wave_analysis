EM-APEX Tools
=============
This reposity contains modules of oceanographic data analysis tools and scripts
that use them for EM-APEX observations.

Modules
------
Certain modules may require numpy, scipy, the gibbs seawater toolbox (gsw),
matplotlib, basemap and various other bits and bobs.

* utils.py: contains functions to convert MATLAB datenumbers to python datetime 
objects and to calculate the distance between lat lon points as well as the
Bunch class.
* sandwell.py: contains functions to read areas or tracks of bathymetry from the 
Smith and Sandwell binary file.
* window.py: contains functions that apply moving windows and binning to data.
* finescale.py: contains functions for applying finescale parameterisations of 
turbulent dissipation and analysing spectra of shear, strain. (Work in progress)
* GM79.py: the 1979 verion of the Garrett-Munk internal wave spectrum in pythonic 
form.
* emapex.py: contains a class for EM-APEX floats with a host of methods. 
* gravity_waves.py: all things related to internal gravity waves i.e. dispersion 
relations, polarisation relations etc. 

