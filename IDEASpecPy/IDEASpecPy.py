import numpy as np
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit #see the following link for info on curve_fit


wavelengths = [475,488,505,520,535,545]


"""
old values
ecs=np.array([-0.86,-0.62,0.23,1.,0.51, 0.2])
qE=np.array([-.01,-.5,-.2,2.9,5.2,3.2])
qE=qE/np.max(qE)
#Zx=np.array([-0.45,-0.24,1.3,0.45,0.31,0.16])
#Zx=Zx-.3
Zx = np.array([-0.01708832,  0.02287851,  0.56242522,  0.51708838,  0.2996292 , 0.15525777])
scatter=np.array([1.1,.6,.25,.1,.07,.05])
scatter = scatter + .3
drift=np.array([.2,.4,.6,.75,.65,.55])
"""


# New values??
# Define the spectral signatures (the set of effective extinction coefficients) for the species 
# to deconvolute. This case, the electrochromic shift ('ecs'), the 535nm 'qE' signal, the 
# signal associated with violaxantin-->zeaxanthin conversion ('Zx') and signals associated with 
# changing leaf orientation ('drift') and light scattering ('scatter').

#Define the ecsc oeffficients

def defineFitCoefficients():
    ecs=2*np.array([-1.18558135, -0.72516973,  0.13613229,  1. ,  0.21642832,-0.46697259])
    ecs=ecs+.5

    #Define the qE coeffficients

    qE=np.array([-.01,-.5,-.2,2.9,5.2,3.2])
    qE=10*qE/np.max(qE)

    #Define the qE coeffficients

    Zx=1000*np.array([-0.00179332, -0.00144022,  0.00184943, -0.00044366, -0.00158732,-0.00190934])
    Zx=Zx+2.0

    # Zx=np.array([-0.15,-0.24,0.8,0.45,0.11,0.03])
    # Zx=Zx-2

    #Define the 'scatter' coeffficients

    scatter=10*np.array([1.1,.6,.25,.1,.07,.05])
    #scatter = np.array([0,0,0,0,0,0])
    #drift=-10*np.array([.2,.4,.6,.75,.65,.55])

    #Define the 'scatter' coeffficients
    drift=np.array([-0.65471441,  0.35155469,  2.03802372,  3.25785718,  3.482087 , 1.36204459])
    return (ecs,qE,Zx,drift,scatter)

ecs,qE,Zx,drift,scatter=defineFitCoefficients()

import glob

#baseFileName='/Users/davidkramer/Dropbox/Data/atsuko/jsontest 050919/col-0_3g-simple_1'

def generateCombinedJSON (baseFileName):
    oldJSON = glob.glob(baseFileName + '.datcombined' '*.json')

    print(oldJSON)
    for oldJSONfile in oldJSON:
        os.rename(oldJSONfile, oldJSONfile + 'x') 


    filesToCombine = glob.glob(baseFileName + '.dat' '*.json')
    f= open(baseFileName+'_combined.json',"w+")
    f.write('[')
    f.close()
    f= open(baseFileName+'_combined.json',"a+")

    for index, file in enumerate(filesToCombine):
        c = open(file,"r+")
        f.write(c.read())
        if (index<len(filesToCombine)-1):
            f.write(', ')
    f.write(']')        
    f.close()

   
# """
#         Recalculates the Abs data using a particular set of I0 points
#         input the dataframe, and optionally a list of 
#         measuring_light_names and the suffix to use, 
#         and a list of points to average for the I0 value to use in the calculation.

#         for example, 

#         recalcAbsUsingI0 (sample_df, ['475', '488', '505', '520', '535', '545'], I_suffix='_I', [0,1]):
#         the '_I' will be attached to the end of the wavelength when searching the dataframe
#         in this 

# """
    
# def recalcAbsUsingI0 (df, list_of_measuring_lights, I_suffix='_I', I0_points=[0,1], trim_points=[0,-1], newSufficx='_rDA'):
#     if (I0_points[0] == I0_points[1]):
#         I0_points[1]=I0_points[0]+1
#     for wl in df['measuring_light_names'][0]:
#         df[wl + newSufficx] = 0
#         df[wl + newSufficx] = df[wl + newSufficx].astype('object')

#     for t in range(len(df)): #['475']:
#         for wl in df['measuring_light_names'][t]:
            
#             I0=np.mean(np.array(df[wl+ '_I'][t][I0_points[0]:I0_points[1]]))
#             daTrace = -1*np.log10(np.array(df[wl+ '_I'][t])/I0) #np.array(df[wl+ '_I'][t][0]))
#             df[wl + newSufficx][t] = daTrace #[trim_points[0]:trim_points[1]]
            

def recalcAbsUsingBurnInProfile (df, wavelengths,  I_suffix, burn_in_profile, burn_in_suffix, burn_in_index, I0_trace=0, I0_points=[-1,-1], newSuffix='_rDA'):
    """
    #                 Recalculates the Abs data using a particular set of I0 points
    #                 and a separately measured "burn in" profile. 
    #                 for example, 
    #                 recalcAbsUsingBurnInProfile (sample_df, ['475', '488', '505', '520', '535', '545'], '_I', 
    #                     burn_in_profile, '_I', [0,1], '_rDA'):
    #                 burn_in_profile points to a location in a dataframe containing traces taken with the same protocol, but 
    #                 without a photosynthetically active sample, e.g. a sheet of green paper. 
    #                 New columns with named wavelength + newSuffix (in this case '_rDA') will be generated and attached to the 
    #                 dataframe                
    """
    for wli in wavelengths:  # generate new columns to hold objects
        wl=str(wli)
        df[wl + newSuffix] = 0
        df[wl + newSuffix] = df[wl + newSuffix].astype('object')
        
    for t in range(len(df)): 
        for wli in wavelengths:
            wl=str(wli)
            I0=np.array(burn_in_profile.loc[burn_in_index][wl+ burn_in_suffix])
            I0=I0/np.mean(I0[I0_points[0]:I0_points[1]])
            I = np.array(df[wl+ I_suffix][t])
            Iz = np.array(df[wl+ I_suffix][I0_trace])
            if (I0_points[0]>-1):    # if negative number is input, do not do 
                I = I/np.mean(Iz[I0_points[0]:I0_points[1]])
    #             elif (I[0]<0):
    #                 I=-1*I
            daTrace = -1*np.log10(I/I0) #np.array(df[wl+ '_I'][t][0]))
            
            df[wl + newSuffix][t] = daTrace #[trim_points[0]:trim_points[1]]
    for wli in wavelengths:
            wl=str(wli)
            df[wl + newSuffix][t] = np.array(df[wl + newSuffix][t])-np.mean(np.array(df[wl + newSuffix][I0_trace])[I0_points[0]:I0_points[1]])

 

# The subtraceBaseline funciton adds a set of traces where the baseline is subtracted. The baseline MUST BE
# the same trace EXCEPT that no catinic was changed AND it is saved with a specific label, e.g. "baseline". 
# The wavelengths paramter is a list of the wavelengths upon which traces to subtract.
# The suffix paramter is the suffix for the traces of interest (e.g. '_calc') 


# def subtractBaseline(sample_df, wavelengths, traceSuffix, baselineSuffix, newSuffix, baselineLabel = 'baseline'):
#     """
#         The subtraceBaseline funciton adds a set of traces where the baseline is subtracted. The baseline MUST BE
#         the same trace EXCEPT that no catinic was changed AND it is saved with a specific trace_label, e.g. "baseline". 
#         The wavelengths paramter is a list of the wavelengths upon which traces to subtract.
#         The suffix paramter is the suffix for the traces of interest (e.g. '_calc') 
         
#          fpor example: 
#          wavelenths = [475,488,505,520,535,545]
#          Ipy.subtractBaseline(sample_df, wavelengths, "_calc", "_calc_m_b", 'baseline)

#          where sample_df is the dataframe holding the data sets, wavelenths are the wavelengths of interest
#          the columns you wish to subtract the baseline from are '475_calc', '488_calc'...
#          the trace with the baseline data has trace_label of 'baseline'
#          and you wish to name the new columns '475_calc_m_b', '488_calc_m_b'...

#     """

#     baseline=sample_df[sample_df['trace_label']==baselineLabel]
#     bindex=baseline.index[0]
#     traces=sample_df[sample_df['trace_label']=='fluct']
    
#     for  wli in range(len(wavelengths)):
        
#         wl=str(wavelengths[wli])
#         #print(wl)
#         sample_df[wl + newSuffix] = 0
#         sample_df[wl + newSuffix] = sample_df[wl + newSuffix].astype('object')
#         for i in range(0, len(sample_df[wl+traceSuffix])): # cyclie through the traces for each wl
#             #print(i)
#             sample_df[wl + newSuffix][i] = np.array(sample_df[wl + traceSuffix][i]) - np.array(baseline[wl + baselineSuffix][bindex])
    
# #len(sample_df['475_time'][0]), len(sample_df['475_I'][0])

def subtractBaseline(sample_df,wavelengths, suffix, AvPoints):

    """
    The subtraceBaseline funciton adds a set of traces where the baseline is subtracted. The baseline MUST BE
    the same trace EXCEPT that no catinic was changed AND it is saved with a specific label, e.g. "baseline". 
    The wavelengths paramter is a list of the wavelengths upon which traces to subtract.
    The suffix paramter is the suffix for the traces of interest (e.g. '_calc') 
    """


    #plt.figure()
    for wli in wavelengths:  # generate new columns to hold objects
        wl=str(wli)
        for t in range(len(sample_df)):
            avbl=np.mean(sample_df[str(wl) + suffix][t][AvPoints[0]:AvPoints[1]])
            sample_df[str(wl)+suffix][t] = np.array(sample_df[str(wl)+suffix][t])-avbl



from scipy.optimize import curve_fit #see the following link for info on curve_fit


def fit_burn(x,a,b,c):
    
    """
        Function for fitting absorbance traces to a funciton to account for 
        I will fit the baseline with a hyperbolic or exponential
        
        a is the amplitude of the burn-in
        b is the time constant
        c is the offset

    """

    fit_burn = c + a/(1+x/b)
    return (fit_burn)
def burnCorrection(sample_df,wavelengths,suffix,baselineArray,newSuffix = '_bcor'):
    fit_b=baselineArray[0]
    fit_e=baselineArray[1]
    
    for i in range(len(wavelengths)):
        wl = str(wavelengths[i])
        sample_df[str(wl) + newSuffix] = 0
        sample_df[str(wl) + newSuffix] = sample_df[str(wl) + newSuffix].astype('object')
        tName = wl + suffix
        xName = wl + '_time'
        #corTraces=[]
        #print(len(sample_df[tName]))
        for traceNum in range (0, len(sample_df[tName])):
            #print(tName)
            y_data=sample_df[tName][traceNum]
            x_data=np.array(sample_df[xName][traceNum])
            x_data=x_data-x_data[0]
            popt, pcov = curve_fit(fit_burn, x_data[fit_b:fit_e], y_data[fit_b:fit_e], p0=[.001,1,.001], bounds=([-np.inf, 0, -np.inf], np.inf))
            corTrace=[]
            #print(len(y_data))
            for i, x in enumerate(x_data):
                vc = y_data[i] - fit_burn(x, popt[0],popt[1],popt[2])
                corTrace.append(vc)
            #print(str(wl) + newSuffix)
            sample_df[str(wl) + newSuffix][traceNum]=corTrace
                            
from scipy.signal import savgol_filter

def smoothTraces(sample_df,wavelengths,suffix,smooth_window=50,newSuffix = '_smooth'):

    for i in range(len(wavelengths)):
        
        wl = str(wavelengths[i])
        tName = wl + suffix
        newColName= tName + newSuffix
        sample_df[newColName] = 0
        sample_df[newColName] = sample_df[newColName].astype('object')
        for traceNum in range (0, len(sample_df[tName])):  # cycle through each trace in experiment 
            subtrace=sample_df[tName][traceNum]
            smoothedTrace=savgol_filter(subtrace , smooth_window, 3)
            sample_df[newColName][traceNum]=smoothedTrace


def generateDAS(sample_df, wavelengths, suffix='_rDA', newColName='das'):

    #sample_df[colName]= pd.Series([0]*len(sample_df))
    #sample_df[colName] = [[0]*len(sample_df)]
    sample_df[newColName] = 0
    sample_df[newColName] = sample_df[newColName].astype('object')

    #wavelengths=['475', '488', '505', '520', '535', '545']
    for t in range(0,len(sample_df)):
        dass=[]
        wl = str(wavelengths[0]) #start with the first wl
        for pt in range(len(sample_df[wl + suffix][t])): #cycle through all points in the trace
            das=[]
            for wl in range(len(wavelengths)):
                das.append(sample_df[str(wavelengths[wl]) + suffix][t][pt])
            dass.append(das)
        sample_df[newColName][t] = dass


# def fit_spec_5(x,a,b,c,d,e): #function to fit the spectroscopic data
#     global ecs
#     global qE
#     global Zx
#     global scatter
#     global drift
    
#     # Expecting absorbance data from 5 wavelengths per time point
#     # The component spectra also have 7 wavelength
#     #print(len(components[0]),len(components[1]),len(components[2]),len(components[3]))
#     fit_spectrum = a*ecs + b*qE + c*Zx + d*scatter + e*drift
#     return (fit_spectrum)

def fit_spec_5(x,a,b,c,d,e): #function to fit the spectroscopic data
    global ecs
    global qE
    global Zx
    global scatter
    global drift
    
    # Expecting absorbance data from 5 wavelengths per time point
    # The component spectra also have 7 wavelength
    #print(len(components[0]),len(components[1]),len(components[2]),len(components[3]))
    fit_spectrum = a*ecs + b*qE + c*Zx + d*scatter + e*drift
    return (fit_spectrum)



# def fitDAS_5(sample_df, wavelengths, dasColName, newSuffix=''):
#     components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
#     for c in components:
#         newColName= c + newSuffix
#         sample_df[newColName] = 0
#         sample_df[newColName] = sample_df[newColName].astype('object')

#     for traceNum in range (0, len(sample_df[dasColName])):  # cycle through each trace in experiment 
#         comp={}
#         for c in components:
#             comp[c]=[]
       
#         for time_index in range(len(sample_df[dasColName][0])):
#             y_data=sample_df[dasColName][traceNum][time_index]
#             popt, pcov = curve_fit(fit_spec_5, wavelengths, y_data, p0=[0,0,0,0,0])
#             for i, c in enumerate(components):
#                 comp[c].append(popt[i])
            
#         for c in components:
#             newColName= c + newSuffix
#             sample_df[newColName][traceNum] = comp[c]

# this cell only defines a function, see the next cell for an example of how you would use this function

# add/replace one row to the given dataframe, by computing averages for certain rows/columns
def compute_average_row( sample_df, rows_to_average, columns_to_average, label_for_new_averaged_row):
    # build a list of cell values for the new row, including blanks for non-relevenat columns
    values_for_new_row = []
    for column in sample_df.columns:
        if column in columns_to_average:

            # collect all values to be averaged
            array_to_average_2d = []
            for row in rows_to_average:
                array_to_average_2d.append(np.array(sample_df.loc[row,column]))

            # compute average
            values_for_new_row.append(np.mean(np.array(array_to_average_2d), axis=0))
        else:
            values_for_new_row.append("")

    # add the new row
    if label_for_new_averaged_row in sample_df.index:
        sample_df = sample_df.drop(label_for_new_averaged_row,axis=0)
    sample_df.loc[label_for_new_averaged_row] = values_for_new_row



def fitDAS_5(sample_df, wavelengths, dasColName, newSuffix=''):
    components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
    for c in components:
        newColName= c + newSuffix
        sample_df[newColName] = 0
        sample_df[newColName] = sample_df[newColName].astype('object')

    for traceNum in range (0, len(sample_df[dasColName])):  # cycle through each trace in experiment 
        comp={}
        for c in components:
            comp[c]=[]
       
        for time_index in range(len(sample_df[dasColName][0])):
            y_data=sample_df[dasColName][traceNum][time_index]
            popt, pcov = curve_fit(fit_spec_5, wavelengths, y_data, p0=[0,0,0,0,0])
            for i, c in enumerate(components):
                comp[c].append(popt[i])
            
        for c in components:
            newColName= c + newSuffix
            sample_df[newColName][traceNum] = comp[c]


def calcActinicLightProfile(sample_df):

    """
        Makes new columns in the sample dataframe with the calculated LED settings for 
        actinic, blue and far red lights:
        actinicIntensityTime contains the time axis information
        actinicIntensity contains the main actinic_intenisty information
        blueActinic contains the main blue_actinic information (just 0 or 1)
        farRed contains the main far_red information (just 0 or 1)

    """
    sample_df['actinicIntensity'] = 0
    sample_df['actinicIntensity'] = sample_df['actinicIntensity'].astype('object')
    sample_df['actinicIntensityTime'] = 0
    sample_df['actinicIntensityTime'] = sample_df['actinicIntensityTime'].astype('object')
    sample_df['blueActinic'] = 0
    sample_df['blueActinic'] = sample_df['blueActinic'].astype('object')
    sample_df['farRed'] = 0
    sample_df['farRed'] = sample_df['farRed'].astype('object')
    
    for traceNum in range(len(sample_df)): 
        current_time = 0
        times=[]
        #times.append(0)
        actinic_intensity=[]
        blue_actinic=[]
        far_red=[]
        #actinic_intensity.append(sample_df['actinic_intensity'][traceNum][0])
        for loop in range(len(sample_df['number_pulses'][traceNum])):
            times.append(current_time)
            actinic_intensity.append(sample_df['actinic_intensity'][traceNum][loop])
            blue_actinic.append(sample_df['blue_actinic'][traceNum][loop])
            far_red.append(sample_df['far_red'][traceNum][loop])
        #     print(sample_df['actinic_intensity'][0][loop])
            for c in range(len(sample_df['number_pulses'][traceNum])):
                cycle = sample_df['number_pulses'][traceNum][c]
                for mpc in range(0,cycle):
                    for mp in range(len(sample_df['measuring_light_names'][traceNum])):
                        current_time = current_time + sample_df['measuring_interval'][traceNum][mp] + sample_df['measuring_pulse_duration'][traceNum][mp]
                        #print(sample_df['measuring_interval'][0][mp])
                times.append(current_time)
                actinic_intensity.append(sample_df['actinic_intensity'][traceNum][loop])
                blue_actinic.append(sample_df['blue_actinic'][traceNum][loop])
                far_red.append(sample_df['far_red'][traceNum][loop])
                
        sample_df['actinicIntensityTime'][traceNum] = times
        sample_df['actinicIntensity'][traceNum] = actinic_intensity
        sample_df['blueActinic'][traceNum] = blue_actinic
        sample_df['farRed'][traceNum] = far_red

        