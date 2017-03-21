def FFP(zm=None, z0=None, umean=None, h=None, ol=None, sigmav=None, ustar=None, 
		wind_dir=None, rs=None):
	"""
	Derive a flux footprint estimate based on the simple parameterisation FFP
	See Kljun, N., P. Calanca, M.W. Rotach, H.P. Schmid, 2015: 
	The simple two-dimensional parameterisation for Flux Footprint Predictions FFP.
	Geosci. Model Dev. 8, 3695-3713, doi:10.5194/gmd-8-3695-2015, for details.
	contact: n.kljun@swansea.ac.uk

	FFP Input
	zm		= Measurement height above displacement height (i.e. z-d) [m]
	z0		= Roughness length [m]; enter None if not known 
	umean	= Mean wind speed at zm [m/s]; enter None if not known 
			  Either z0 or umean is required. If both are given,
			  umean is selected to calculate the footprint
	h		= Boundary layer height [m]
	ol		= Obukhov length [m]
	sigmav	= standard deviation of lateral velocity fluctuations [ms-1]
	ustar   = friction velocity [ms-1]

	optional inputs:
	wind_dir = wind direction in degrees (of 360) for rotation of the footprint     
	rs       = percentage of source area, i.e. a value between 10 and 90. 

	FFP output
	x_ci_max  = x location of footprint peak (distance from measurement) [m]
	x_ci	  = x array of crosswind integrated footprint [m]
	f_ci	  = array with footprint function values of crosswind integrated footprint [m-1] 
	x_2d	  = x-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
	y_2d	  = y-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
	f_2d	  = footprint function values of 2-dimensional footprint [m-2]
        rs        = percentage of footprint as in input, if provided
        fr        = footprint value at r, if r is provided
        xr_yr     = x and y-array for contour line of r, if r is provided
	flag_err  = 1 in case of error, 0 otherwise

    created: 15 April 2015 natascha kljun
    translated to python, December 2015 gerardo fratini, LI-COR Biosciences Inc.
    version: 1.1
    last change: 14/12/2015 natascha kljun
    Copyright (C) 2015, Natascha Kljun
	"""
	
	import numpy as np
	import sys


	# Check existence of required pars
	if None in [zm, h, ol, sigmav, ustar] or (z0 is None and umean is None):
		raise_ffp_exception(1)

	# Check passed values
	if zm <= 0.: raise_ffp_exception(2)
	if z0 is not None and umean is None and z0 <= 0.: raise_ffp_exception(3)
	if h <= 10.: raise_ffp_exception(4)
	if zm > h and ol > 0.: raise_ffp_exception(5)
	if zm > 0.8*h and ol < 0.: raise_ffp_exception(6)
	if z0 is not None and umean is None and zm < 12.5*z0: raise_ffp_exception(7)
	if float(zm)/ol <= -15.5: raise_ffp_exception(8)
	if sigmav <= 0: raise_ffp_exception(9)

	#Resolve ambiguity if both z0 and umean are passed (defaults to using z0)
	if None not in [z0, umean]:
		print '''Both z0 and umean were provided. Proceeding using umean.\n 
			Provided z0 value will be ignored'''

	#Grid resolution		
    #Selection of nx has impact on accuracy and on output file size, 
	#decrease for speed, increase for accuracy (nx=3000 ideal but slow, 
	#nx=600 fast but may not resolve details)
	nx = 800
	xstar_end=30
	
	#Model parameters
	a = 1.4524
	b = -1.9914
	c = 1.4622
	d = 0.1359
	ac = 2.17 
	bc = 1.66
	cc = 20.0

	#limit to L for neutral scaling
	oln = 5000

	#von Karman
	k = 0.4

	#Create scaled X* for crosswind integrated footprint
	xstar_ci_param = np.linspace(d, xstar_end, nx+2)
	xstar_ci_param = xstar_ci_param[1:]

	#Calculate crosswind integrated scaled F* 
	fstar_ci_param = a * (xstar_ci_param-d)**b * np.exp(-c/ (xstar_ci_param-d))
	ind_notnan	 = ~np.isnan(fstar_ci_param)
	fstar_ci_param = fstar_ci_param[ind_notnan]
	xstar_ci_param = xstar_ci_param[ind_notnan]

	#Calculate scaled sig_y*
	sigystar_param = ac * np.sqrt(bc * xstar_ci_param**2 / (1 + cc * xstar_ci_param))

	#Calculate real scale x and f_ci
	if umean is not None:
		#Use umean if available
		x = xstar_ci_param * zm / (1. - zm / h) * (umean / ustar * k)
		if umean / ustar > 0:
			x_ci = x
			f_ci = fstar_ci_param / zm * (1. - zm / h) / (umean / ustar * k)
		else:
			x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d = None
			flag_err = 1
	else:
		#Use z0 if umean is not available
		if ol <= 0 or ol >= oln:
			xx  = (1 - 19.0 * zm/ol)**0.25
			psi_f = np.log((1 + xx**2) / 2.) + 2. * np.log((1 + xx) / 2.) - 2. * np.arctan(xx) + np.pi/2
		elif ol > 0 and ol < oln:
			psi_f = -5.3 * zm / ol

		x = xstar_ci_param * zm / (1. - (zm / h)) * (np.log(zm / z0) - psi_f)
		if np.log(zm / z0) - psi_f > 0:
			x_ci = x
			f_ci = fstar_ci_param / zm * (1. - (zm / h)) / (np.log(zm / z0) - psi_f)
		else:
			x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d = None
			flag_err = 1

	#Calculate maximum location of influence (peak location)
	xstarmax = -c / b + d
	if umean is not None:
		x_ci_max = xstarmax * zm / (1. - (zm / h)) * (umean / ustar * k)
	else:
		x_ci_max = xstarmax * zm / (1. - (zm / h)) * (np.log(zm / z0) - psi_f)

	#Calculate real scale sig_y
	if abs(ol) > oln:
		ol = -1E6
	if ol <= 0:   #convective
		scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.80
	elif ol > 0:  #stable
		scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.55

	if scale_const > 1:
		scale_const = 1.0

	sigy = sigystar_param / scale_const * zm * sigmav / ustar
	sigy[sigy < 0] = np.nan

	#Calculate real scale f(x,y)
	dx = x_ci[2] - x_ci[1]
	y_pos = np.arange(0, (len(x_ci) / 2.) * dx * 1.5, dx)
	#f_pos = np.full((len(f_ci), len(y_pos)), np.nan)
	f_pos = np.empty((len(f_ci), len(y_pos)))
	f_pos[:] = np.nan
	for ix in range(len(f_ci)):
		f_pos[ix,:] = f_ci[ix] * 1 / (np.sqrt(2 * np.pi) * sigy[ix]) * np.exp(-y_pos**2 / ( 2 * sigy[ix]**2))

	#Complete footprint for negative y (symmetrical)
	y_neg = - np.fliplr(y_pos[None, :])[0]
	f_neg = np.fliplr(f_pos)
	y = np.concatenate((y_neg[0:-1], y_pos))
	f = np.concatenate((f_neg[:, :-1].T, f_pos.T)).T

	#Matrices for output
	x_2d = np.tile(x[:,None], (1,len(y)))
	y_2d = np.tile(y.T,(len(x),1))

	# Rotate 3d footprint if requested
	if 0. < wind_dir < 360.:			
		wind_dir = wind_dir * np.pi / 180.
		dist = np.sqrt(x_2d**2 + y_2d**2)
		angle = np.arctan2(y_2d, x_2d)
		x_2d = dist * np.sin(wind_dir - angle)
		y_2d = dist * np.cos(wind_dir - angle)

	#Calculate contours for given percentages
	#Check input and resolve to default levels in needed
	v = None
	clevs = None
	if rs is not None:
		clevs = get_contour_levels(f, dx, dx, rs)
		v = []
		for clev in clevs:
			v.append(get_contour_vertices(x_2d, y_2d, f, clev[2]))

	return {'x_ci_max': x_ci_max, 'x_ci': x_ci, 'f_ci': f_ci,
                'x_2d': x_2d, 'y_2d': y_2d, 'f_2d': f,
		'rs': rs, 'fr': clevs, 'xr_yr': v}

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

	#Calculate levels
	pclevs = np.empty(len(rs))
	pclevs[:] = np.nan
	ars = np.empty(len(rs))
	ars[:] = np.nan
	#Flatten and sort in descending order the 2D footprint function
	sf = np.sort(f, axis=None)[::-1]
	#Masked array for handling potential nan
	msf = ma.masked_array(sf, mask=(np.isnan(sf) | np.isinf(sf)))
	#Cumulative sum
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
	return [(segs[0][ix, 0], segs[0][ix, 1]) for ix in range(N)]   # x,y coords of contour points.	

#===============================================================================
def plot_footprint(X, Y, fs, clevs=None, 
				   show_footprint=True,
				   normalize=None, 
				   colormap=None, line_width=0.3):
	'''Plot footprint function and contours if request'''

	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.colors import LogNorm

	#If input is a list of footprints, don't show footprint but only contours, 
	#with different colors
	if isinstance(fs, list): 
		show_footprint = False
	else:
		fs = [fs]

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
			cc = [c]*len(levs) 
			cp = ax.contour(Y, X, f, levs, colors = cc, linewidths=line_width) 
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
		
		xmin = np.nanmin(X)
		xmax = np.nanmax(X)
		ymin = np.nanmin(Y)
		ymax = np.nanmax(Y)
		
		data = fs[0][:,:].T
		im = ax.imshow(data, cmap=colormap, extent=(xmin, xmax, ymin, ymax), 
			norm=norm, origin='lower', aspect=1) 

		#Colorbar
		cbar = fig.colorbar(im, shrink=1.0, format='%.8f')
		cbar.set_label('Flux contribution', color = 'k')
	plt.show()

	return fig, ax




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
	 'msg': 'h (BPL height) must be larger than 10m.'},				 
	{'code': 5,
	 'type': exTypes['fatal'],
	 'msg': 'zm (measurement height) must be smaller than h (PBL height).'},				 
	{'code': 6,
	 'type': exTypes['fatal'],
	 'msg': 'For convective stratifications, zm (measurement height) must '
			' be smaller than 0.8*h (the entrainment height).'}, 
	{'code': 7,
	 'type': exTypes['fatal'],
	 'msg': 'zm (measurement height) must be larger than 12.5*z0 (roughness sub-layer) '},
	{'code': 8,
	 'type': exTypes['fatal'],
	 'msg': 'zm/ol (measurement height to Obukhov length ratio) must be equal or larger than -15.5'},
	{'code': 9,
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
