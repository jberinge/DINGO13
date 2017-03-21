def FFP_climatology(zm=None, z0=None, umean=None, h=None, ol=None, sigmav=None, ustar=None, 
                     domain=None, wind_dir=None, rs=None, smooth_data=None):
	"""
	Derive a flux footprint estimate based on the simple parameterisation FFP
	See Kljun, N., P. Calanca, M.W. Rotach, H.P. Schmid, 2015: 
	The simple two-dimensional parameterisation for Flux Footprint Predictions FFP.
	Geosci. Model Dev. 8, 3695-3713, doi:10.5194/gmd-8-3695-2015, for details.
	contact: n.kljun@swansea.ac.uk

	This function calculates footprints within a fixed physical domain for a series of
	time steps, rotates footprints into the corresponding wind direction and aggregates 
	all footprints to a footprint climatology. The percentage of source area is
	calculated for the footprint climatology.
	For determining the optimal extent of the domain (large enough to include footprints)
        use calc_footpritn_FFP.py.

	FFP Input
	zm       = Measurement height above displacement height (i.e. z-d) [m]
               usually a scalar, but can also be a vector 
	z0       = Roughness length [m]; enter None if not known 
               usually a scalar, but can also be a vector 
	umean    = Mean wind speed at zm [m/s]; enter None if not known 
	           Either z0 or umean is required. If both are given,
	           umean is selected to calculate the footprint
	h        = Vector of boundary layer height [m]
	ol       = Vector of Obukhov length [m]
	sigmav 	 = Vector of standard deviation of lateral velocity fluctuations [ms-1]
	ustar    = Vector of friction velocity [ms-1]
        domain   = Domain size as an array of xmin, xmax, ymin, ymax [m]
               Footprint will be calculated for a measurement at [0 0 zm] m
               Enter None for default [-1000 1000 -1000 1000]m
	wind_dir = Vector of wind direction in degrees (of 360) for rotation of the footprint     

	optional inputs:
	rs       	= Percentage of source area, i.e. a value between 10 and 90. 
				  Can be either a single value (e.g., "80") or an array of 
				  increasing percentage values (e.g., "[10:10:80]") 
        smooth_data = Apply convolution filter to smooth footprint climatology if smooth_data=1

	FFP output
	FFP      = Structure array with footprint climatology data for measurement at [0 0 zm] m
	x_2d	 = x-grid of 2-dimensional footprint [m]
	y_2d	 = y-grid of 2-dimensional footprint [m]
	fclim_2d = Normalised footprint function values of footprint climatology [m-2]
        rs       = Percentage of footprint as in input, if provided
        fr       = Footprint value at r, if r is provided
        xr_yr    = x and y-array for contour line of r, if r is provided
	n        = Number of footprints calculated and included in footprint climatology
	flag_err = 1 in case of error, 2 if not all contour plots (r%) within specified domain,
			   0 otherwise

    Created: 19 May 2016 natascha kljun
    Translated to python, 26 May, together with gerardo fratini, LI-COR Biosciences Inc.
    version: 1.0
    last change: 26/05/2016 natascha kljun
    Copyright (C) 2015, Natascha Kljun
	"""
	
	import numpy as np
	import sys
	#pip install scipy
	from scipy import signal as sg

	## Check existence of required pars
	#if None in [zm, h, ol, sigmav, ustar] or (z0 is None and umean is None):
#		raise_ffp_exception(1)

#	# Check passed values
#	if zm <= 0.: raise_ffp_exception(2)
#	if z0 is not None and umean is None and z0 <= 0.: raise_ffp_exception(3)
#	if h <= 10.: raise_ffp_exception(4)
#	if zm > h : raise_ffp_exception(5)
#	if z0 is not None and umean is None and zm < 12.5*z0: raise_ffp_exception(6)
#	if float(zm)/ol <= -15.5: raise_ffp_exception(7)
#	if sigmav <= 0: raise_ffp_exception(8)

	#Resolve ambiguity if both z0 and umean are passed (defaults to using z0)
	if None not in [z0, umean]:
		print '''Both z0 and umean were provided. Proceeding using z0.\n
			 Provided umean value will be ignored'''

        #!!Need to create the following values by checking input - to be done
        flag_error =0
	len_series = len(ustar)
        valid = np.ones(len_series)
         	
	#Model parameters
	a = 1.4524
	b = -1.9914
	c = 1.4622
	d = 0.1359
	ac = 2.17 
	bc = 1.66
	cc = 20.0   
	oln = 5000 #limit to L for neutral scaling
	k = 0.4 #von Karman
 
	#===========================================================================
	#Define domain
	#default
	nx = 1000
	ny = nx
	if domain is None:
		xmin = -nx
		xmax = nx
		ymin = -ny
		ymax = ny
	else:
		xmin, xmax, ymin, ymax = domain 
		dx = (xmax-xmin) / nx
		dy = (ymax-ymin) / ny
		
	#Define physical domain in cartesian and polar coordinates
	# Cartesian coordinates
	x = np.linspace(xmin, xmax, nx + 1)
	y = np.linspace(ymin, ymax, ny + 1)
	x_2d, y_2d = np.meshgrid(x, y)
	# Polar coordinates
	# Set theta such that North is pointing upwards and angles increase clockwise
	rho = np.sqrt(x_2d**2 + y_2d**2)
	theta = np.arctan2(x_2d, y_2d)

	rho = np.sqrt(x_2d**2 + y_2d**2)
	theta = np.arctan2(y_2d, x_2d)
	#JB added
	theta = theta - wind_dir * np.pi / 180.


	# initialize raster for footprint climatology
	fclim_2d = np.zeros(x_2d.shape)
	
	#===========================================================================
	# Start loop on time series
	for foot_loop in range (0, len_series):
	
		if foot_loop % 100 == 1:
			print 'Calculating footprint ', repr(foot_loop), ' of ', repr(len_series)
		
		if valid[foot_loop]:
		
	#===========================================================================
	# Rotate coordinates into wind direction
			if wind_dir[foot_loop] is not None: 
				wind_dir_rad = wind_dir[foot_loop] * np.pi / 180.
				theta_loop = theta - wind_dir_rad
		
	#===========================================================================
	# Create real scale crosswind integrated footprint and dummy for
	# rotated scaled footprint
			fstar_ci_dummy = np.zeros(x_2d.shape)
			f_ci_dummy = np.zeros(x_2d.shape)
			if z0 is not None:
				#Use z0
				if ol[foot_loop] <= 0 or ol[foot_loop] >= oln:
					xx = (1 - 19.0 * zm[foot_loop]/ol[foot_loop])**0.25
					psi_f = (np.log((1 + xx**2) / 2.) + 2. * np.log((1 + xx) / 2.) - 2. * np.arctan(xx) + 
					        np.pi/2)
				elif ol[foot_loop] > 0 and ol[foot_loop] < oln:
					psi_f = -5.3 * zm[foot_loop] / ol[foot_loop]
				if (np.log(zm[foot_loop] / z0[foot_loop])-psi_f)>0:
					xstar_ci_dummy = (rho * np.cos(theta_loop) / zm[foot_loop] * (1. - (zm[foot_loop] / 
					                 h[foot_loop])) / (np.log(zm[foot_loop] / z0[foot_loop]) - psi_f))
					px = np.where(xstar_ci_dummy > d)
					fstar_ci_dummy[px] = a * (xstar_ci_dummy[px] - d)**b * np.exp(-c / (xstar_ci_dummy[px] - d))
					f_ci_dummy[px] = (fstar_ci_dummy[px] / zm[foot_loop] * (1. - (zm[foot_loop] / 
					                 h[foot_loop])) / (np.log(zm[foot_loop] / z0[foot_loop]) - psi_f))
				else:
                                        flag_err = 1

			else:
				#Use umean if z0 is not available
				xstar_ci_dummy = (rho * np.cos(theta_loop) / zm[foot_loop] * (1. - (zm[foot_loop] / 
				                  h[foot_loop])) / (umean[foot_loop] / ustar[foot_loop] * k))
				px = np.where(xstar_ci_dummy > d)
				fstar_ci_dummy[px] = a * (xstar_ci_dummy[px] - d)**b * np.exp(-c / (xstar_ci_dummy[px] - d))
				f_ci_dummy[px] = (fstar_ci_dummy[px] / zm[foot_loop] * (1. - (zm[foot_loop] / 
				                 h[foot_loop])) / (umean[foot_loop] / ustar[foot_loop] * k))

	#===========================================================================
	# Calculate dummy for scaled sig_y* and real scale sig_y
			sigystar_dummy = np.zeros(x_2d.shape)
			sigystar_dummy[px] = (ac * np.sqrt(bc * np.abs(xstar_ci_dummy[px])**2 / (1 + 
			                      cc * np.abs(xstar_ci_dummy[px]))))

			if abs(ol[foot_loop]) > oln:
				ol[foot_loop] = -1E6
			if ol[foot_loop] <= 0:   #convective
				scale_const = 1E-5 * abs(zm[foot_loop] / ol[foot_loop])**(-1) + 0.80
			elif ol[foot_loop] > 0:  #stable
				scale_const = 1E-5 * abs(zm[foot_loop] / ol[foot_loop])**(-1) + 0.55
			if scale_const > 1:
				scale_const = 1.0

			sigy_dummy = np.zeros(x_2d.shape)
			sigy_dummy[px] = (sigystar_dummy[px] / scale_const * zm[foot_loop] * sigmav[foot_loop] 
			                 / ustar[foot_loop])
			sigy_dummy[sigy_dummy < 0] = np.nan

	#===========================================================================
	# Calculate real scale f(x,y)
			f_2d = np.zeros(x_2d.shape)
			f_2d[px] = (f_ci_dummy[px] / (np.sqrt(2 * np.pi) * sigy_dummy[px]) *
           			   np.exp(-(rho[px] * np.sin(theta_loop[px]))**2 / ( 2. * sigy_dummy[px]**2)))

	#===========================================================================
	# Add to footprint climatology raster
			fclim_2d = fclim_2d + f_2d;
			
	#===========================================================================
	# Normalize and smooth footprint climatology	
	fclim_2d = fclim_2d / np.sum(valid);

	if smooth_data is not None:
		skernel  = np.matrix('0.05 0.1 0.05; 0.1 0.4 0.1; 0.05 0.1 0.05')
		fclim_2d = sg.convolve2d(fclim_2d,skernel,mode='same');
		fclim_2d = sg.convolve2d(fclim_2d,skernel,mode='same');

	#===========================================================================
	# Derive footprint ellipsoid incorporating R% of the flux, if requested,
	# starting at peak value
	v = None
	clevs = None
	if rs is not None:
		clevs = get_contour_levels(fclim_2d, dx, dy, rs)
		v = []
		for clev in clevs:
			v.append(get_contour_vertices(x_2d, y_2d, fclim_2d, clev[2]))

	#===========================================================================
	# Fill output structure
	n = np.sum(valid)
	return {'x_2d': x_2d, 'y_2d': y_2d, 'fclim_2d': fclim_2d,
                'rs': rs, 'fr': clevs, 'xr_yr': v,'n':n,'flag_error':flag_error}

#===============================================================================
#===============================================================================
def get_contour_levels(f, dx, dy, rs=None):
	'''Contour levels of f at percentages of f-integral given by rs'''

	import numpy as np
	from numpy import ma

	#Check input and resolve to default levels in needed
	if not isinstance(rs, (int, float, list)):
		rs = list(np.linspace(0.10, 0.90, 9))
	if isinstance(rs, (int, float)): rs = [rs]
	
	if np.min(rs) > 1: rs = [x/100. for x in rs]

	#Sort levels in ascending order
	rs = np.sort(rs)

	#Levels
	pclevs = np.empty(len(rs))
	pclevs[:] = np.nan
	ars = np.empty(len(rs))
	ars[:] = np.nan

	sf = np.sort(f, axis=None)[::-1]
	msf = ma.masked_array(sf, mask=(np.isnan(sf) | np.isinf(sf))) #Masked array for handling potential nan
	csf = msf.cumsum().filled(np.nan)*dx*dy
	for ix, r in enumerate(rs):
		dcsf = np.abs(csf - r)
		pclevs[ix] = sf[np.nanargmin(dcsf)]
		ars[ix] = csf[np.nanargmin(dcsf)]

	return [(round(r, 3), ar, pclev) for r, ar, pclev in zip(rs, ars, pclevs)]

#===============================================================================
def get_contour_vertices(x, y, f, lev):
	import matplotlib._cntr as cntr
	c = cntr.Cntr(x, y, f)
	nlist = c.trace(lev, lev, 0)
	segs = nlist[:len(nlist)//2]
	N = len(segs[0][:, 0])
	return [(segs[0][ix, 0], segs[0][ix, 1]) for ix in range(N)] # x,y coords of contour points.	

#===============================================================================
def plot_footprint(x_2d, y_2d, fs, clevs, show_footprint, normalize, colormap, line_width, group_x,night_day_label,year_label,month_label,myBaseforResults,Site_ID):
	'''Plot footprint function and contours if request'''

	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.colors import LogNorm

	#If input is a list of footprints, don't show footprint but only contours, 
	#with different colors
	if isinstance(fs, list): 
		show_footprint = False
	#else:
		#fs = [fs]

	if colormap is None: colormap = cm.jet
	#Define colors for each contour set
	cs = [colormap(ix) for ix in np.linspace(0, 1, len(fs))]

	# Initialize figure
	fig, ax = plt.subplots(figsize=(12, 10))
	# fig.patch.set_facecolor('none')
	# ax.patch.set_facecolor('none')

	if clevs is not None:
		#Plot isopleth
		levs = [clev[2] for clev in clevs]
                for f, c in zip(fs, cs):
			cc = [c]*len(clevs) 
			cp = ax.contour(y_2d, x_2d, f, levs, colors = c, linewidths=line_width) 
			#Isopleth Labels
			pers = [str(int(clev[0]*100))+'%' for clev in clevs]
			fmt = {}
			for l,s in zip(cp.levels, pers):
				fmt[l] = s
			plt.clabel(cp, cp.levels[:], inline=1, fmt=fmt, fontsize=7)
		#Isopleth Labels
		pers = [str(int(clev[0]*100))+'%' for clev in clevs]
		fmt = {}
		for l,s in zip(cp.levels, pers):
			fmt[l] = s
		plt.clabel(cp, cp.levels[:], inline=1, fmt=fmt, fontsize=7)

	#plot footprint function
	if show_footprint:
		if normalize == 'log': 
			norm = LogNorm()
		else:
			norm = None 
		
		xmin = np.nanmin(x_2d)
		xmax = np.nanmax(x_2d)
		ymin = np.nanmin(y_2d)
		ymax = np.nanmax(y_2d)
		im = ax.imshow(fs[:, :].T, cmap=colormap, extent=(xmin, xmax, ymin, ymax), 
			norm=norm, origin='lower', aspect=1) 


		#Colorbar
		cbar = fig.colorbar(im, shrink=1.0, format='%.8f')
		cbar.set_label('Flux contribution', color = 'k')
	
	# JB added text box with mean quantities
	# JB then added save figure
	zm_mean = group_x['zm'].mean()
	z0_mean = group_x['z0'].mean()
	umean_mean = group_x['umean'].mean()
	h_mean= group_x['h'].mean()
	ol_mean = group_x['ol'].mean()
	sigmav_mean = group_x['vv'].mean()
	ustar_mean = group_x['ustar'].mean()
	n_mean=len(group_x)     	
	
	graphtext1=str('Zm          ' + str("{0:.2f}".format(zm_mean)) +'\n' +
	               'Zo          ' + str("{0:.2f}".format(z0_mean)) +'\n' +
	               'u mean      ' + str("{0:.2f}".format(umean_mean)) +'\n' +
	               'h           ' + str("{0:.2f}".format(h_mean)) +'\n' +
	               'MOL mean    ' + str("{0:.2f}".format(ol_mean)) +'\n' +
	               'sigma v mean' + str("{0:.2f}".format(sigmav_mean)) +'\n' +
	               'n samples   ' + str("{0:.2f}".format(n_mean)) +'\n' +
	               'ustar mean  ' + str("{0:.2f}".format(ustar_mean))       )  
	
	plt.figtext(0.2,0.2,graphtext1, color='w', bbox=dict())
	plt.suptitle('Footprint climatology'+Site_ID + ' '+year_label+' '+month_label + ' '+night_day_label)
	plt.savefig(myBaseforResults+'/Footprint_climatology_Kljun_'+Site_ID + '_'+year_label+'_'+month_label + '_'+night_day_label)
	print 'Saving Footprint_climatology_Kljun_'+Site_ID + '_'+year_label+'_'+month_label + '_'+night_day_label
	#return fig, ax




#===============================================================================
#===============================================================================
exTypes = {'message': 'Message',
			 'alert': 'Alert',
			 'error': 'Error',
			 'fatal': 'Fatal error'}

exceptions = [
	{'code': 1,
	 'type': exTypes['fatal'],
	 'msg': 'At least one required parameter is missing. Please enter all '
			'required inputs. Check documentation for details.'},				 
	{'code': 2,
	 'type': exTypes['fatal'],
	 'msg': 'zm (measurement height) must be larger than zero.'},				 
	{'code': 3,
	 'type': exTypes['fatal'],
	 'msg': 'z0 (roughness length) must be larger than zero.'},				 
	{'code': 4,
	 'type': exTypes['fatal'],
	 'msg': 'h (BPL height) must be larger than 10 m.'},				 
	{'code': 5,
	 'type': exTypes['fatal'],
	 'msg': 'zm (measurement height) must be smaller than h (PBL height).'},				 
	{'code': 6,
	 'type': exTypes['fatal'],
	 'msg': 'zm (measurement height) must be larger than 12.5*z0 (roughness sub-layer) '},
	{'code': 7,
	 'type': exTypes['fatal'],
	 'msg': 'zm/ol (measurement height to Obukhov length ratio) must be equal or larger than -15.5'},
	{'code': 8,
	 'type': exTypes['fatal'],
	 'msg': 'sigmav (standard deviation of crosswind) must be larger than zero'},
	]

def raise_ffp_exception(code):
	'''Raise exception or prints message according to specified code'''
	
	ex = [it for it in exceptions if it['code'] == code][0]
	string = ex['type'] + '(' + str(ex['code']).zfill(4) + '):\n '+ ex['msg'] 

	print('')
	if ex['type'] == exTypes['fatal']:
		string = string + '\n FFP_fixed_domain execution aborted.'
		raise Exception(string)
	else:
		print(string)
