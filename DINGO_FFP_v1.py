import calc_footprint_FFP_climatology as Footprint_climatology
import calc_footprint_FFP as Footprint
import numpy as np
import meteorologicalfunctions as metfunc
from configobj import ConfigObj
import pandas as pd
import os
from windrose import plot_windrose
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import math


def dataframe_check(Dataframe, FluxFreq):
    #Check dataframe for duplicates, pad as necessary and sort
    Dataframe.sort(inplace=True)
    Dataframe["index"] = Dataframe.index
    Dataframe.drop_duplicates('index', take_last=True, inplace=True)
    del Dataframe["index"]	
    Dataframe=Dataframe.asfreq(FluxFreq, method=None)      
    return Dataframe

def CoordRotation2D(ds):
    #Taken directly from Isaac OzFluxQC v2.9.5
    """
        2D coordinate rotation to force v = w = 0.  Based on Lee et al, Chapter
        3 of Handbook of Micrometeorology.  This routine does not do the third
        rotation to force v'w' = 0.
        
        Usage qcts.CoordRotation2D(ds)
        ds: data structure
        """
    # get the raw wind velocity components
    Ux = ds.Ux          # longitudinal component in CSAT coordinate system
    Uy = ds.Uy          # lateral component in CSAT coordinate system
    Uz = ds.Uz          # vertical component in CSAT coordinate system
    # get the raw covariances
    UxUz = ds.UxUz # covariance(Ux,Uz)
    UyUz = ds.UyUz # covariance(Uy,Uz)
    UxUy = ds.UxUy      # covariance(Ux,Uy)
    UyUy = ds.UyUy      # variance(Uy)
    UxUx = ds.UxUx      # variance(Ux)
    UzUz = ds.UzUz      # variance(Ux)
    UzC = ds.UzC    # covariance(Uz,C)
    UzA = ds.UzA    # covariance(Uz,A)
    UzT = ds.UzT    # covariance(Uz,T)
    UxC = ds.UxC        # covariance(Ux,C)
    UyC = ds.UyC       # covariance(Uy,C)
    UxA = ds.UxA        # covariance(Ux,A)
    UyA = ds.UyA        # covariance(Ux,A)
    UxT = ds.UxT        # covariance(Ux,T)
    UyT = ds.UyT        # covariance(Uy,T)
    nRecs = len(ds)     # number of records

    print ' Applying 2D coordinate rotation (components and covariances)'
    # get the 2D and 3D wind speeds
    ws2d = np.sqrt(Ux**2 + Uy**2)
    ws3d = np.sqrt(Ux**2 + Uy**2 + Uz**2)
    # get the sine and cosine of the angles through which to rotate
    #  - first we rotate about the Uz axis by eta to get v = 0
    #  - then we rotate about the v axis by theta to get w = 0
    ce = Ux/ws2d          # cos(eta)
    se = Uy/ws2d          # sin(eta)
    ct = ws2d/ws3d        # cos(theta)
    st = Uz/ws3d          # sin(theta)
    # get the rotation angles
    theta = np.rad2deg(np.arctan2(st,ct))
    eta = np.rad2deg(np.arctan2(se,ce))
    # do the wind velocity components first
    u = Ux*ct*ce + Uy*ct*se + Uz*st           # longitudinal component in natural wind coordinates
    v = Uy*ce - Ux*se                         # lateral component in natural wind coordinates
    w = Uz*ct - Ux*st*ce - Uy*st*se           # vertical component in natural wind coordinates
    # do the variances
    uu = UxUx*ct**2*ce**2 + UyUy*ct**2*se**2 + UzUz*st**2 + 2*UxUy*ct**2*ce*se + 2*UxUz*ct*st*ce + 2*UyUz*ct*st*se
    vv = UyUy*ce**2 + UxUx*se**2 - 2*UxUy*ce*se
    ww = UzUz*ct**2 + UxUx*st**2*ce**2 + UyUy*st**2*se**2 - 2*UxUz*ct*st*ce - 2*UyUz*ct*st*se + 2*UxUy*st**2*ce*se
    # now do the scalar covariances
    wT = UzT*ct - UxT*st*ce - UyT*st*se       # covariance(w,T) in natural wind coordinate system
    wA = UzA*ct - UxA*st*ce - UyA*st*se       # covariance(w,A) in natural wind coordinate system
    wC = UzC*ct - UxC*st*ce - UyC*st*se       # covariance(w,C) in natural wind coordinate system
    # now do the momentum covariances
    # full equations, Wesely PhD thesis via James Cleverly and EddyPro
    uw = UxUz*ce*(ct*ct-st*st) - 2*UxUy*ct*st*ce*se + UyUz*se*(ct*ct-st*st) - \
         UxUx*ct*st*ce*ce - UyUy*ct*st*se*se + UzUz*ct*st # covariance(w,x) in natural wind coordinate system
    uv = UxUy*ct*(ce*ce-se*se) + UyUz*st*ce - UxUz*st*se - \
         UxUx*ct*ce*se + UyUy*ct*ce*se                    # covariance(x,y) in natural wind coordinate system
    vw = UyUz*ct*ce - UxUz*ct*se - UxUy*st*(ce*ce-se*se) + \
         UxUx*st*ce*se - UyUy*st*ce*se                    # covariance(w,y) in natural wind coordinate system

    ds['vv']=vv
    return ds


def distance_for_r(x,r):
    print "Doing FFP for r% " + str(r) + " timedate " + str(x.name) 
    print "zm" + str(x['zm']) + " and  OL " + str(x['ol']) + " and ratio " + str(x['ol']/x['zm'])
    print x
    a = np.array(Footprint.FFP(x['zm'], x['z0'], x['umean'], x['h'], x['ol'], x['vv'], x['ustar'], x['Wd'], r)['xr_yr'])[0]
    b= np.sqrt((a[:,0]**2) + (a[:,1]**2))
    c=max(b)
    
    return c            

def calculate_h(x):
    #For stable and neutral conditions there are simple diagnostic
    #relations with which the boundary layer height can be
    #estimated. Nieuwstadt (1981) proposed an interpolation formula
    #for neutral to stable conditions    
    omega =  7.2921*(10.0**-5.0) #being the angular velocity of the Earth’s rotation 7.2921 × 10-5 rad/s
    theta =  latitude * np.pi /180  #being  latitude in radians so convert from degree
    f = abs(2*omega*np.sin(theta))
    h = abs(((x['ol'])/3.8) * ( -1.0 + (1.0 + 2.28 * ( x['ustar']/ (f * x['ol']) ) )**0.5    ))
    if h >1500: h=1500
    if h <50: h=50
    if math.isnan(h): 
        h=500
    if (x['ol']) < 0: h=1500
    return h

def main():
    
    #At least for testing and stand alone get the required information from the config file
    #Later just run from DINGO and pass info as required
    #cf = ConfigObj('E:/My Dropbox/Dropbox/Data_flux_data/controlfiles/Whroo/Advanced/AdvancedProcessing_config_Whroo_v13.txt') 
    cf = ConfigObj('E:/My Dropbox/Dropbox/Data_flux_data/controlfiles/HowardSprings/Advanced/AdvancedProcessing_config_HowardSprings_v13.txt') 
    
    
    #Input site details
    versionID="v12a"    
    Site_ID=cf['Site']['Site_ID']    
    Canopy_height=float(cf['Site']['Canopy_height'])
    Instrument_height =float(cf['Site']['Instrument_height'])
    mypathforFluxdata=cf['Files']['mypathforFluxdata']+Site_ID+ "/Advanced_"+versionID
    myBaseforResults_FFP=cf['Files']['mypathforFluxdata']+Site_ID+ "/Advanced_"+versionID+"/Diagnostics/Footprint"
    myBaseforResults_windrose=cf['Files']['mypathforFluxdata']+Site_ID+ "/Advanced_"+versionID+"/Diagnostics/Windrose"
    Ws_variable_name = cf['Options']['Ws_variable_name']+'_Con'
    global latitude
    latitude= float(cf['Site']['Tower_Lat'])
    
    #Toggle what you want to do
    do_ffp_distances = True
    do_ffp_climatology_plots = True
    do_windrose_climatology_plots = True
    
    #Read dataframe from site
    FLUXDataframe_orig= pd.read_pickle(mypathforFluxdata+'/Advanced_processed_data_'+Site_ID+'_'+versionID+'.df')
    
    print "Starting Footprint analysis"
    #Check for place to put results - does it exist? If not create    
    #Then subdirectories
    if not os.path.isdir(myBaseforResults_FFP):
        os.mkdir(myBaseforResults_FFP) 
    if not os.path.isdir(myBaseforResults_windrose):
        os.mkdir(myBaseforResults_windrose) 
    if not os.path.isdir(myBaseforResults_windrose+'/Annual'):
        os.mkdir(myBaseforResults_windrose+'/Annual')        
    if not os.path.isdir(myBaseforResults_windrose+'/Monthly'):
        os.mkdir(myBaseforResults_windrose+'/Monthly') 
                        
    #Prepare data for footprint analysis

    # Apply 2D coordinate rotation from Isaac OzFluxQC v2.9.5
    # We use this to get the sigma v required for input to FFP
    FLUXDataframe_orig = CoordRotation2D(FLUXDataframe_orig)
    
    # To calculate a single FFP flux footprint with calc_footprint_FFP.
    # Call calc_footprint_FFP(zm,z0,umean,h,ol,sigmav,ustar) using inputs
    # zm = Measurement height above displacement height z minus d [m]    
    FLUXDataframe_orig['zm']=Instrument_height-(Canopy_height*0.666)
    # z0 = Roughness length [m]  enter [NaN] if not known
    FLUXDataframe_orig['z0'] = Canopy_height * 0.1
    # umean = Mean wind speed at zm [ms-1] enter [NaN] if not known
    FLUXDataframe_orig['umean'] = FLUXDataframe_orig[Ws_variable_name]
    
    #Subset data if required for testing
    FFP_Dataframe=FLUXDataframe_orig["2015-01-30 00:00":"2015-02-01 23:30"]
    #FFP_Dataframe=FLUXDataframe_orig
    
    #Check for missing values and drop them
    FFP_Dataframe= FFP_Dataframe[['Ta_Con','Ah_Con','ps_Con','ustar','Fh_Con', 'wT','vv', Ws_variable_name, 'Wd','zm','z0','umean','day_night']].dropna(how='any')    
    #Calculate Monin Obukov Length calling Isaac OxFluxQC Metfunction
    print "Applying ML function"
    FFP_Dataframe['ol'] = FFP_Dataframe.apply(lambda x: metfunc.molen(x['Ta_Con'],x['Ah_Con'],x['ps_Con'],x['ustar'],x['Fh_Con'],fluxtype='sensible'),axis=1)
    # h = Boundary layer height [m] set differently for time of day/night
    
    FFP_Dataframe['h'] = FFP_Dataframe.apply(lambda x: calculate_h(x),  axis=1)
    
    print 'test'
    #FFP_Dataframe['h'][FFP_Dataframe['day_night']==2] = 1500.0
    #FFP_Dataframe['h'][FFP_Dataframe['day_night']==3] = 750.0


    #In the FFP code there is are exception for numerous items so catch them here
    #	if sigmav <= 0: raise_ffp_exception(8)
    #    zm/ol (measurement height to Obukhov length ratio) must be equal or larger than -15.5'}
    #    So lets catch that here
    
    print "Number rows total                " + str(len(FFP_Dataframe))
    print "Number rows with zmol exceptions " + str(len(FFP_Dataframe[FFP_Dataframe['ol']/FLUXDataframe_orig['zm'] < -15.5]))
    print "Number rows with sigma v exceptions " + str(len(FFP_Dataframe[FFP_Dataframe['vv']==0]))
    
    FFP_Dataframe=FFP_Dataframe[FFP_Dataframe['ol']/FLUXDataframe_orig['zm'] >= -15.5]
    FFP_Dataframe=FFP_Dataframe[FFP_Dataframe['vv']!=0]
    
    #First run FFP for each time step then later do footprint climatology    
    if do_ffp_distances == True:
        #Apply FFP to dataframe as function
        #The following Function returns the max distance away from the tower for any given r%
	 
        FFP_Dataframe['FFP_x_10%'] = FFP_Dataframe.apply(lambda x: distance_for_r(x,10),  axis=1)
        FFP_Dataframe['FFP_x_30%'] = FFP_Dataframe.apply(lambda x: distance_for_r(x,30),  axis=1)
        FFP_Dataframe['FFP_x_50%'] = FFP_Dataframe.apply(lambda x: distance_for_r(x,50),  axis=1)
        FFP_Dataframe['FFP_x_70%'] = FFP_Dataframe.apply(lambda x: distance_for_r(x,70),  axis=1)
        FFP_Dataframe['FFP_x_90%'] = FFP_Dataframe.apply(lambda x: distance_for_r(x,90),  axis=1)
        FFP_Dataframe['FFP_x_peak'] = FFP_Dataframe.apply(lambda x: Footprint.FFP(x['zm'], x['z0'], x['umean'], x['h'], x['ol'], x['vv'], x['ustar'], x['Wd'], None)['x_ci_max'],axis=1)
        #FFP_Dataframe['FFP_flag'] = FFP_Dataframe.apply(lambda x: Footprint.FFP(x['zm'], x['z0'], x['umean'], x['h'], x['ol'], x['vv'], x['ustar'], x['Wd'], None)['flag_error'],axis=1)
        
        #Do a check on the dataframe including resampling to 30 minutes before merging
        FFP_Dataframe=dataframe_check(FFP_Dataframe, '30T')
        
        #Output the FFP input and output
        FFP_Dataframe[['zm', 'z0', 'umean', 'h', 'ol', 'vv', 'ustar', 'Wd','FFP_x_10%','FFP_x_30%','FFP_x_50%','FFP_x_70%','FFP_x_90%','FFP_x_peak']].to_csv(myBaseforResults_FFP+'/FFP_output_'+Site_ID+'.csv')
    
        #Next take the results and merge with the original dataframe
        FLUXDataframe_orig_plus_FFP = pd.merge(FLUXDataframe_orig, FFP_Dataframe[['ol', 'vv', 'FFP_x_10%','FFP_x_30%','FFP_x_50%','FFP_x_70%','FFP_x_90%','FFP_x_peak']],left_index=True,right_index=True,how="left")   #Now join
    


    if do_ffp_climatology_plots == True:
        # Next Do the climatology footprints.  First groupby year and month
        by = lambda x: lambda y: getattr(y, x)
        FFP_Dataframe_grouped=FFP_Dataframe.groupby([by('year'),by('month')])    
        
        for group_index, group_x in FFP_Dataframe_grouped:
            year_label=str (group_index[0])
            month_label=str(group_index[1])
            print "Doing footrint for year " + year_label + "and " + month_label + " month"
            #Do for Night and then day seperately
            nightdata=group_x[group_x['day_night']==3]
            night_day_label = 'Night'
            zm_list = nightdata['zm'].values
            z0_list = nightdata['z0'].values
            umean_list = nightdata['umean'].values
            h_list = nightdata['h'].values
            ol_list = nightdata['ol'].values
            sigmav_list = nightdata['vv'].values
            ustar_list = nightdata['ustar'].values
            wind_dir_list = nightdata['Wd'].values        
            # Call Kljun footprint module passing lists of required parameters and then call the plotting routine
            FFP_output1=Footprint_climatology.FFP_climatology(zm_list,z0_list,umean_list,h_list,ol_list,sigmav_list,ustar_list,[-250 ,250, -250 ,250],wind_dir_list,None,None) 
            Footprint_climatology.plot_footprint(FFP_output1['x_2d'], FFP_output1['y_2d'], FFP_output1['fclim_2d'], None , True, None, None, 0.3, nightdata,night_day_label,year_label,month_label,myBaseforResults_FFP,Site_ID)        
            
                  
    
            #Do for Night and then day seperately
            daydata=group_x[group_x['day_night']==1]
            night_day_label = 'Day'
            zm_list = daydata['zm'].values
            z0_list = daydata['z0'].values
            umean_list = daydata['umean'].values
            h_list = daydata['h'].values
            ol_list = daydata['ol'].values
            sigmav_list = daydata['vv'].values
            ustar_list = daydata['ustar'].values
            wind_dir_list = daydata['Wd'].values        
            # Call Kljun footprint module passing lists of required parameters and then call the plotting routine
            FFP_output1=Footprint_climatology.FFP_climatology(zm_list,z0_list,umean_list,h_list,ol_list,sigmav_list,ustar_list,[-250 ,250, -250 ,250],wind_dir_list)     
            Footprint_climatology.plot_footprint(FFP_output1['x_2d'], FFP_output1['y_2d'], FFP_output1['fclim_2d'], None , True, None, None, 0.3, daydata,night_day_label,year_label,month_label,myBaseforResults_FFP,Site_ID)    
        
    #Now do Wind rose plots
    #Create a small DF for Wind rose plots
    DF_windrose=FFP_Dataframe[['umean','Wd','day_night']]  
    DF_windrose.columns = ['speed', 'direction','day_night']    
    
    by = lambda x: lambda y: getattr(y, x)
    

    if do_windrose_climatology_plots == True:
        #Define which type of plot(s) to do. Create a list and plot different types as you want
        #kind : kind of plot (might be either, 'contour', 'contourf', 'bar', 'box', 'pdf')
        plotkinds=['box','contourf','pdf']
        #Colourmap type (i.e. hot, summer, winter, plasma, etc.)     
        cmtype = cm.plasma
        
        for kindrose in plotkinds:
                
            #Do ALL data plot first
            DF_windrose_temp=DF_windrose
            #Do Windrose for all
            startyear=str(DF_windrose_temp.index[0].year)
            endyear=str(DF_windrose_temp.index[-1].year)
            
            print "Doing Windrose for year " + startyear + '_to_' + endyear
            fig = plt.figure(1, figsize=(6, 6), dpi=300, facecolor='w', edgecolor='k')
            plot_windrose(DF_windrose_temp, kind=kindrose, bins=np.arange(0.01,8,1), cmap=cmtype, lw=2)
            plt.suptitle('Windrose climatology'+Site_ID + ' '+  startyear + '_to_' + endyear,size=16)
            #plt.show() 
            plt.savefig(myBaseforResults_windrose+'/Annual'+'/Windrose climatology '+Site_ID + '_for_' + startyear + '_to_' + endyear +'_' + kindrose+'_' + versionID)
            print 'Saving Windrose_climatology_'+Site_ID + '_' + startyear + '_to_' + endyear    +'_' + kindrose       
            
            #Do annual plots first
            DF_windrose_temp=DF_windrose.groupby([by('year')])
            for group_index, group_x in DF_windrose_temp:
                #Do Windrose for all
                year_label=str (group_index)
    
                print "Doing Windrose for year " + year_label
                fig = plt.figure(1, figsize=(6, 6), dpi=300, facecolor='w', edgecolor='k')
                plot_windrose(group_x, kind=kindrose, bins=np.arange(0.01,8,1), cmap=cmtype, lw=2)
                plt.suptitle('Windrose climatology'+Site_ID + ' '+year_label,size=16)
                #plt.show() 
                plt.savefig(myBaseforResults_windrose+'/Annual'+'/Windrose climatology '+Site_ID + ' '+year_label+'_' + kindrose+'_'+versionID)
                print 'Saving Windrose_climatology_'+Site_ID + '_'+year_label +'_' + kindrose   +'_'+versionID 
            
            #Then do monthly plots
            DF_windrose_temp=DF_windrose.groupby([by('year'),by('month')])    
            for group_index, group_x in DF_windrose_temp:
                #Do Windrose for night
                year_label=str (group_index[0])
                month_label=str(group_index[1])
                print "Doing Windrose for year " + year_label + "and " + month_label + " month"
                #Do for Night and then day seperately
                nightdata=group_x[group_x['day_night']==3]
                night_day_label = 'Night'
                fig = plt.figure(1, figsize=(6, 6), dpi=300, facecolor='w', edgecolor='k')
                plot_windrose(nightdata, kind=kindrose, bins=np.arange(0.01,8,1), cmap=cmtype, lw=2)
                plt.suptitle('Windrose climatology'+Site_ID + ' '+year_label+' Month '+month_label + ' and '+night_day_label+'_'+versionID,size=16)
                #plt.show()
                plt.savefig(myBaseforResults_windrose+'/Monthly'+'/Windrose climatology '+Site_ID + ' '+year_label+' Month '+month_label + ' and '+night_day_label+'_'+kindrose+'_'+versionID)
                print 'Saving Windrose_climatology_'+Site_ID + '_'+year_label+'_month_'+month_label + '_'+night_day_label+'_'+kindrose+'_'+versionID
                
                
                #Do WIndrose for Day
                year_label=str (group_index[0])
                month_label=str(group_index[1])
                print "Doing Windrose for year " + year_label + "and " + month_label + " month"
                #Do for Night and then day seperately
                daydata=group_x[group_x['day_night']==1]
                night_day_label = 'Day'
                fig = plt.figure(1, figsize=(6, 6), dpi=300, facecolor='w', edgecolor='k')
                plot_windrose(daydata, kind=kindrose, bins=np.arange(0.01,8,1), cmap=cmtype, lw=2)
                plt.suptitle('Windrose climatology'+Site_ID + ' '+year_label+' Month '+month_label + ' and '+night_day_label+'_'+versionID,size=16)
                #plt.show()
                plt.savefig(myBaseforResults_windrose+'/Monthly'+'/Windrose climatology '+Site_ID + ' '+year_label+' Month '+month_label + ' and '+night_day_label+'_'+kindrose+'_'+versionID)
                print 'Saving Windrose_climatology_'+Site_ID + '_'+year_label+'_month_'+month_label + '_'+night_day_label+'_'+kindrose+'_'+versionID

    
    print "Finished"   
        
    #return FLUXDataframe_orig_plus_FFP

main()
