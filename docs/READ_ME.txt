Variables in new EM-APEX Matlab velocity files
John Dunlap, dunlap@apl.washington.edu
August 5, 2009

A typical file name is "ema-3763a-0065-vel.mat" where,
	3763a -- float serial number with a suffix to indicate run number.
	0065  -- "hpid", half-profile ID, odd numbers descend, even numbers ascend.

On CTD grid:
	ctd_mlt -- UTC of CTD samples (Matlab datenum)
	Pctd -- pressure (dbar)
	T -- temperature (degrees C)
	S -- salinity (PSU)

Note: The CTD data is better quality ascending than descending because
the CTD inlet is on the top.

On Processed electric field (EFP) grid:
	efp_mlt -- UTC of EFP estimates (Matlab datenum)
	Pef   -- pressure (dbar)
	U1,U2 -- east component of velocity from electrode sets 1,2 (m/s)
	V1,V2 -- north component of velocity from electrode sets 1,2 (m/s)
	Verr1 -- standard error of U1,V1 (m/s)
	Verr2 -- standard error of U2,V2 (m/s)
	Wef   -- vertical velocity derived from Pctd (m/s)
	RotP  -- rotation period (s)
	e1mean, e2mean -- mean value of EF (microvolts)
	e1sdev, e2sdev -- standard deviation of residuals of EF fit (microvolts)

Note: When RotP becomes more than about EmaProcessNvals / 2 the
velocity processing becomes noisy.  Also, the variables Verr1 and
Verr2 can be used to screen out noisy U1,V1 and U2,V2.

GPS data as sampled just before phoning home with data.
	MLT_GPS -- UTC of GPS sample (Matlab datenum)
	NSAT    -- number of satellites used
	LAT     -- latitude of this profile (degrees)
	LON     -- longitude of this profile (degrees)
	STAT    -- status from GPGGA NMEA message, 1: OK, 0: no fix
	HDOP    -- horizontal dilution of precision

Typically GPS data is sampled once every several profiles so the following
is interpolated position:
	lat     -- latitude interpolated from other profiles (degrees)
	lon     -- longitude interpolated from other profiles (degrees)

More:
	MLT_ref -- UTC of beginning of half-profile. (Matlab datenum)
	           For ascents the reference is at the deepest point.
	magvar  -- magnetic variation at lat/lon (degrees)
	fh,fz   -- earth's horizontal and vertical magnetic field (nT)
	
	hpid    -- half-profile number:  odd numbers descend, even numbers ascend.

