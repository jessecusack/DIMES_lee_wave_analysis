EM-APEX Tools
=============
I put code here that I use to analyse data from EM-APEX floats.

Python
------
Certain modules may require numpy, scipy, the gibbs seawater toolbox (gsw),
matplotlib, basemap and various other bits and bobs.

Less specific bits of code:

* utils.py: contains functions to convert MATLAB datenumbers to python datetime 
objects and to calculate the distance between lat lon points as well as the
Bunch class.
* sandwell.py: contains functions to read areas or tracks of bathymetry from the 
Smith and Sandwell binary file. 

More specific bits of oceanography code:

* finescale.py: contains functions for applying finescale parameterisations of 
turbulent dissipation and analysing spectra of shear, strain. (Work in progress)
* GM79.py: the 1979 verion of the Garrett-Munk internal wave spectrum in pythonic 
form.(Work in progress)
* emapex.py: contains a class for EM-APEX floats with a host of methods. 
* gravity_waves.py: all things related to internal gravity waves i.e. dispersion 
relations, polarisation relations etc. (Work in progress)

Everything else is related to my own analysis of the data and is messy.


MATLAB
------

I've stopped using MATLAB completely but I leave the mess of functions here for 
posterity. Scripts require the gibbs seawater toolbox (gsw) and m_map toolbox.
