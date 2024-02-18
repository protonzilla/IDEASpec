# from winreg import DisableReflectionKey
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.signal import savgol_filter

from scipy.optimize import curve_fit #see the following link for info on curve_fit
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime
import os.path as path
#import matplotlib.colors 
import matplotlib.cm as cm
import glob 
import re
import os, fnmatch
import json
from pathlib import Path
import hashlib
from IPython.display import HTML
import pprint 
from sklearn.decomposition import FastICA, PCA
import scipy.optimize
from scipy.stats import linregress

#test for update synchronization

global s_offset

wavelengths = [475,488,505,520,535,545]
wls = wavelengths

driftTraces=np.array([-0.97930457, -0.84510825, -0.76824325, -0.68557158, -0.65721521,
                     -0.71411414, -0.80057651, -0.9111386 , -0.99341127, -1.        ,
                     -0.90233654, -0.7229873 , -0.50413059, -0.30143799, -0.150493  ,
                     -0.05590374,  0.00378046,  0.03427544,  0.04805041,  0.05288039,
                     0.05259516,  0.033782  ,  0.        ])
    
driftTracesWLs=np.array([460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])

    
print("I am renewed!")

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

global fudge 
global scatter_var
global drift_var 

# New values??
# Define the spectral signatures (the set of effective extinction coefficients) for the species 
# to deconvolute. This case, the electrochromic shift ('ecs'), the 535nm 'qE' signal, the 
# signal associated with violaxantin-->zeaxanthin conversion ('Zx') and signals associated with 
# changing leaf orientation ('drift') and light scattering ('scatter').

#Define the ecs coeffficients
def closestToDf(df, param, v):# finds the index of the value for param closes to the value v
    x=(df[param]-v)**2
    xx=x.sort_values(axis=0)
    return df.iloc[xx.index[0]]
    
def closestToList(alist, v):# finds the index of the value for param closes to the value v
    xx=(np.array(alist)-v)**2
    return xx.argmin()

def hx():
    print("hello")
    
def defineFitCoefficients():
    ecs=2*np.array([-1.18558135, -0.72516973,  0.13613229,  1. ,  0.21642832,-0.46697259])
    ecs=ecs+.5

    #Define the qE coeffficients

    qE=np.array([-.01,-.5,-.2,2.9,5.2,3.2])
    qE=10*qE/np.max(qE)

    #Define the qE coeffficients

    Zx=np.array([-0.00179332, -0.00144022,  0.00184943, -0.00044366, -0.00158732,-0.00190934])
    Zx=Zx/np.max(Zx)
    # Zx=Zx+2.0

    # Zx=np.array([-0.15,-0.24,0.8,0.45,0.11,0.03])
    # Zx=Zx-2

    #Define the 'scatter' coeffficients

    scatter=10*np.array([1.1,.6,.25,.1,.07,.05])
    #scatter = np.array([0,0,0,0,0,0])
    #drift=-10*np.array([.2,.4,.6,.75,.65,.55])

    #Define the 'drift' coeffficients
    drift=np.array([-0.65471441,  0.35155469,  2.03802372,  3.25785718,  3.482087 , 3.36204459])
    # drift=np.array([1.2,  1.35155469,  2.03802372,  3.25785718,  3.482087 , 3.36204459])
    print("Component spectra reset.")

    return (ecs,qE,Zx,drift,scatter)


ecs,qE,Zx,drift,scatter=defineFitCoefficients()

def cols_with (df, part, ignore_case = True):
    cols = []
    for col in list(df.columns):
        if ignore_case:
            col_lc = str.lower(col)
            part = str.lower(part)
        else:
            col_lc = col
        if part in col_lc:
            cols.append(col)
    return cols

#baseFileName='/Users/davidkramer/Dropbox/Data/atsuko/jsontest 050919/col-0_3g-simple_1'


# def hi():
#     print("hey there!")
# print("Again, I am renewed!")


# def fixIt (baseFileName):
# #     if (destinationBaseFileName==""):
# #         destinationBaseFileName=baseFileName
#     print(baseFileName)
#     filesToCombine = glob.glob(baseFileName + '.dat' '*.json')
#     print(filesToCombine)
# #     f= open(destinationBaseFileName+'_combined.json',"w+")
# #     f.write('[')
# #     f.close()
# #     f= open(destinationBaseFileName+'_combined.json',"a+")

#     #print("combining files:")
#     for index, file in enumerate(filesToCombine):
#         print(file)
#         #print(file)    
# #         c = open(file,"r+")
# #         f.write(c.read())
        
# #         if (index<len(filesToCombine)-1):
# #             f.write(', ')
# #     f.write(']')        
# #     f.close()   
# #     print('Combined ' + str(index + 1) + " json files.")
#         with open(file) as myfile:
#             content = myfile.read()

#         text = re.search(r'"I_File_Name".*.dat,', content, re.DOTALL).group()
        
#         with open("result.txt", "w") as myfile2:
#             myfile2.write(text)
    
    
def CharlieDAS(df, traceColName):
    """
    CharlieDAS(df = dataframe containing the experimental results, 
               traceColName = a string for the column that contains the kinetic trace.
    returns a new dataframe with the transposed data set with columns:
    
    time   wl1  wl2 wl3...
    
    CharlieDAS takes in a decay associated spectra experiment from Charlie-type IDEASpec 
    intruments. These contain a range of traces, each at different wavelengths. The fumnction
    combines the traces into a matrix and transposes them so that each row is a different time 
    point and each column is the absorbance change at a given wavelength.
    
    The traces in col traceColName must have the same length!
    
    """
    df=df.sort_values('wl')
    wls = list(df['wl'])
    times = df['1_time'][0]  # make a list of the times col values in the first dataFrame to use in contructing 
                            # the new transposed dataFrame
    allData = []
    for i, wl in enumerate(wls):
        t = df[df['wl']==wl][traceColName]

        trace=list(t.loc[t.index[0]])
        allData.append(trace)
    allData=np.asarray(allData)
    allDataT = allData.transpose()
    timePoint=0
    tdf = pd.DataFrame([np.append(times[timePoint], allDataT[timePoint])], columns = ["time"] + wls)
    for timePoint in range(1,allDataT.shape[0]):
        tdf = tdf.append(pd.DataFrame([np.append(times[timePoint], allDataT[timePoint])], columns = ["time"] + wls), ignore_index=True)
    return tdf

import json

def generateCombinedJSON (baseFileName, destinationBaseFileName="", append_bases = [''], verbose = False):
    if (destinationBaseFileName==""):
        destinationBaseFileName=baseFileName

    oldJSON = glob.glob(baseFileName + '.datcombined' '*.json')

    # print("removed old JSON files:" + str(oldJSON))
    # for oldJSONfile in oldJSON:
    #     os.rename(oldJSONfile, oldJSONfile + 'x') 
    # append_bases = ['']
    for append_base in append_bases:
        if verbose:
            print(append_base)
        addendum = baseFileName + '.dat' + append_base
        filesToCombine = glob.glob(addendum  + '*.json')
        if (len(filesToCombine)<1):
            print("No JSON files with base file name = " + baseFileName + "found ")
            return 
        f= open(destinationBaseFileName+'_combined.json',"w+")
        f.write('[')
        f.close()
        f= open(destinationBaseFileName+'_combined.json',"a+")
        if verbose:
            print("combined file to: " + destinationBaseFileName)
        for index, file in enumerate(filesToCombine):
            #print(file)
            succeeded = False    
            try:
                with open(file) as json_file:
                    jd = json.load(json_file)
                succeeded = True
            except:
                print("failed:" + file)
                succeeded = False
            if(succeeded == True):
                c = open(file,"r+")
                f.write(c.read())
                if (index<len(filesToCombine)-1):
                    f.write(', ')
        f.write(']')        
        f.close()   
        if verbose:
            print('Combined JSON = ' + destinationBaseFileName)

    
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

print('combineExperiment exists')
def combineExperiment(baseFileNames, destinationFolder):
    #dictionary holding all experiments
    experiment = {}  

    for baseFileName in baseFileNames:
        experiment_name = baseFileName.split('/')[-1]

        if (destinationFolder == ''):
            destinationBaseFileName=baseFileName #'/Volumes/D_KRAMER/Thekla/Data FL and LL grown plants/FL'
        else:
            destinationBaseFileName = destinationFolder + experiment_name
        #test if combined_json already exists
        destinationCombinedJSONfileName = destinationBaseFileName + '_combined.json'
        combined_json_exists = path.isfile(destinationCombinedJSONfileName)

        # combined_json_exists= False

        if (combined_json_exists == True):
            print("Found existing combined JSON, called:" + destinationCombinedJSONfileName)
        else:
            print("Generating new combined JSON, called")
            # each experiment should contain a gob of json files. We will combine them
            # go through all baseFileNames and generate the ombined json.



            # call the function that combines the individual JSON files into one big one and saves to 
            print(baseFileName)
            if (combined_json_exists == False):
                generateCombinedJSON(baseFileName, destinationBaseFileName)

#             folder_name='' 

            #Define the file name for the main program to open.

        file_name= destinationBaseFileName + '_combined.json'

        # open the combined JSON into a dictionary of data frames
        experiment[experiment_name] = {}
        df = pd.read_json(file_name) 

        #remove spaces from column names

        cols = df.columns
        newCols = []
        for k in cols:
            newCols.append(k.strip())
        df.columns = newCols
        df.columns
        experiment[experiment_name]['data']=df

        #reorder the data fram by the start_time, so the traces are in chronological order

        experiment[experiment_name]['data'] = experiment[experiment_name]['data'].sort_values('start_time')

        # if notes were taken, add them to the combined json
        try:
            notes_file = open(baseFileName + '.dat_notes.txt', 'r')
            experiment[experiment_name]['notes'] = notes_file.read()
        except:
            experiment[experiment_name]['notes'] = ""
    #make a list that contains the names of all the experiments
    allExperiments = list(experiment.keys())
    return experiment, allExperiments

def recalcAbsUsingPointsFromI_no_offset (df, list_of_measuring_lights, I_suffix='_I', I0_points=[10,20], trim_points=[0,-1], newSufficx='_rDA'):
    if (I0_points[0] == I0_points[1]):
        I0_points[1]=I0_points[0]+1
    for wl in list_of_measuring_lights:
        df[str(wl) + newSufficx] = 0
        df[str(wl) + newSufficx] = df[str(wl) + newSufficx].astype('object')

    
    for wl in list_of_measuring_lights:
        I0 = np.mean(df.iloc[0][str(wl)+ '_I'][I0_points[0]:I0_points[1]])
        # print(I0)
        # I0=np.mean(np.array(df.loc[t, str(wl)+ '_I'][I0_points[0]:I0_points[1]]))
        for t in df.index: #range(len(df)): #['475']:
            
            daTrace = -1*np.log10(np.array(df[str(wl)+ '_I'][t])/I0) #np.array(df[wl+ '_I'][t][0]))
            df[str(wl) + newSufficx][t] = daTrace #[trim_points[0]:trim_points[1]]
            
def recalcAbsUsingPointsFromI (df, list_of_measuring_lights, I_suffix='_I', I0_points=[0,1], trim_points=[0,-1], newSufficx='_rDA'):
    if (I0_points[0] == I0_points[1]):
        I0_points[1]=I0_points[0]+1
    for wl in list_of_measuring_lights:
        df[str(wl) + newSufficx] = 0
        df[str(wl) + newSufficx] = df[str(wl) + newSufficx].astype('object')

    for t in df.index: #range(len(df)): #['475']:
        for wl in list_of_measuring_lights:
            # print(wl)
            I0=np.mean(np.array(df[str(wl)+ '_I'][t][I0_points[0]:I0_points[1]]))
            daTrace = -1*np.log10(np.array(df[str(wl)+ '_I'][t])/I0) #np.array(df[wl+ '_I'][t][0]))
            df[str(wl) + newSufficx][t] = daTrace[trim_points[0]:trim_points[1]]

def recalcAbsUsingPointsFromI_selected_indexes (df, list_of_measuring_lights, I_suffix='_I', 
                    I0_points=[0,1], trim_points=[0,-1], newSufficx='_rDA',
                    indexes = []):
    if (len(indexes) == 0):
        indexes = df.index
    if (I0_points[0] == I0_points[1]):
        I0_points[1]=I0_points[0]+1
    for wl in list_of_measuring_lights:
        df[str(wl) + newSufficx] = 0
        df[str(wl) + newSufficx] = df[str(wl) + newSufficx].astype('object')

    for t in indexes: #range(len(df)): #['475']:
        for wl in list_of_measuring_lights:
            # print(wl)
            I0=np.mean(np.array(df[str(wl)+ '_I'][t][I0_points[0]:I0_points[1]]))
            daTrace = -1*np.log10(np.array(df[str(wl)+ '_I'][t])/I0) #np.array(df[wl+ '_I'][t][0]))
            df[str(wl) + newSufficx][t] = daTrace[trim_points[0]:trim_points[1]]

def recalcAbsUsingI0_no_offset (df, list_of_measuring_lights, I_suffix='_I', I0_suffix='_I0', newSuffix='_rDA'):

    """
        Recalculate the absorbance traces using the reference channel as the _I0 values. Do not subtrace off any
        offsets. This is useful for cases where series of traces are to be combined, among which the relative changes
        in the absorbance are desired. 

    """
    for wl in list_of_measuring_lights: #df['measuring_light_names'][0]:
        newColName=str(wl) + newSuffix
        df[newColName] = 0
        df[newColName] = df[newColName].astype('object')

    for t in df.index: #['475']:
        for wl in list_of_measuring_lights: #df['measuring_light_names'][t]:
            newColName=str(wl) + newSuffix
            I0=np.array(df.loc[t,str(wl)+ I0_suffix])
            I=np.array(df.loc[t, str(wl)+ I_suffix])
            daTrace = -1*np.log10(I/I0) #np.array(df[wl+ '_I'][t][0]))
            df[newColName][t] = daTrace #[trim_points[0]:trim_points[1]]

def recalcAbsUsingI0_pt0_offset (df, list_of_measuring_lights, I_suffix='_I', I0_suffix='_I0', newSuffix='_rDA'):
    for wl in list_of_measuring_lights: #df['measuring_light_names'][0]:
        newColName=str(wl) + newSuffix
        df[newColName] = 0
        df[newColName] = df[newColName].astype('object')
    
    for t in df.index: #['475']:
        for wl in list_of_measuring_lights: #df['measuring_light_names'][t]:
            newColName=str(wl) + newSuffix
            I0=np.array(df.loc[t,str(wl)+ I0_suffix])
            I=np.array(df.loc[t, str(wl)+ I_suffix])
            daTrace = np.array(-1*np.log10(I/I0)) #np.array(df[wl+ '_I'][t][0]))
            daTrace = daTrace - daTrace[0]
            df[newColName][t] = daTrace #[trim_points[0]:trim_points[1]]

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
#         for i in range(0, len(sample_df[wl+traceSuffix])): # cycle through the traces for each wl
#             #print(i)
#             sample_df[wl + newSuffix][i] = np.array(sample_df[wl + traceSuffix][i]) - np.array(baseline[wl + baselineSuffix][bindex])
    
# #len(sample_df['475_time'][0]), len(sample_df['475_I'][0])

def subtractBaseline(sample_df, wavelengths, suffix, AvPoints, new_suffix = '_sub'):

    """
    The subtraceBaseline function adds a set of traces where the baseline is subtracted. The baseline MUST BE
    the same trace EXCEPT that no atinic was changed AND it is saved with a specific label, e.g. "baseline". 
    The wavelengths paramter is a list of the wavelengths upon which traces to subtract.
    The suffix paramter is the suffix for the traces of interest (e.g. '_calc') 
    """


    #plt.figure()
    for wavelength in wavelengths:  # generate new columns to hold objects
        # wl=str(wli)
        data_col = str(wavelength) + suffix

        if (new_suffix == ''):
            newColName = data_col
        else:
            newColName = str(wavelength) + new_suffix #+ '_sub'
            sample_df[newColName] = 0
            sample_df[newColName] = sample_df[newColName].astype('object')
        
        for index in sample_df.index:
                # print(AvPoints)
                avbl = np.mean(sample_df[data_col][index][AvPoints[0]:AvPoints[1]])
                # print(avbl)
                dif = np.array(sample_df[data_col][index])-avbl
                # print(dif)
                sample_df[newColName][index] = dif
    return sample_df


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


def burnCorrection(sample_df,wavelengths,suffix,baselineArray,newSuffix = '_bcor', retainOffset = False):
    fit_b=baselineArray[0]
    fit_e=baselineArray[1]
    print("de-burning")
    for i in range(len(wavelengths)):
        wl = str(wavelengths[i])
        sample_df[str(wl) + newSuffix] = 0
        sample_df[str(wl) + newSuffix] = sample_df[str(wl) + newSuffix].astype('object')
        # print(str(wl) + newSuffix)
        tName = wl + suffix
        # print("response (x) colummn = " + tName)
        xName = wl + '_time'
        # print("time column = " + xName)
        #corTraces=[]
        #print(len(sample_df[tName]))
        # sample_df_s = sample_df[tName].select_dtypes(include=['object'])

        mask = sample_df[tName].apply(lambda x: isinstance(x, list))

        for traceNum in sample_df.index: #range (0, len(sample_df[tName])):
            # print('trace index = ' + str(traceNum))
            if (1): #mask[traceNum] == True):
                # print('trace Name = ' + str(tName))
                y_data = sample_df.loc[traceNum, tName]
                # print(len(y_data))
                # x_data=np.array(sample_df.loc[traceNum, xName])
                # x_data=x_data-x_data[0]
                x_data = np.linspace(0,len(y_data)-1, len(y_data))

                popt, pcov = curve_fit(fit_burn, x_data[fit_b:fit_e], y_data[fit_b:fit_e], p0=[.1,.1,.001], maxfev=100000 ) #bounds=([-np.inf, 0, -np.inf], np.inf), )
                # print(popt)
                corTrace=[]
                offset = 0.0
                if (retainOffset == True):
                    offset = np.mean(y_data[fit_b:fit_e])

                # print(len(x_data), len(y_data))
                for i, x in enumerate(y_data):
                    vc = y_data[i] - fit_burn(x_data[i], popt[0],popt[1],popt[2]) + offset
                    corTrace.append(vc)
                # print(str(wl) + newSuffix)
                # print(corTrace)
                
                sample_df[str(wl) + newSuffix].loc[traceNum]=np.array(corTrace)
                # print(sample_df[str(wl) + newSuffix].loc[traceNum])
                # sample_df[str(wl) + newSuffix][traceNum]=np.array(corTrace)
                
                # sample_df.loc[traceNum, str(wl) + newSuffix]=corTrace


def burnCorrectionTwoRegions(sample_df,wavelengths,suffix,baselineArray,newSuffix = '_bcor', retainOffset = False):

    fit_b=baselineArray[0]
    fit_e=baselineArray[1]
    print("de-burning")
    for i in range(len(wavelengths)):
        wl = str(wavelengths[i])
        sample_df[str(wl) + newSuffix] = 0
        sample_df[str(wl) + newSuffix] = sample_df[str(wl) + newSuffix].astype('object')
        # print(str(wl) + newSuffix)
        tName = wl + suffix
        # print("response (x) colummn = " + tName)
        xName = wl + '_time'
        # print("time column = " + xName)
        #corTraces=[]
        #print(len(sample_df[tName]))
        

        for traceNum in sample_df.index: #range (0, len(sample_df[tName])):
            #print('trace index = ' + str(traceNum))

            #print('trace Name = ' + str(tName))
            y_data=np.array(sample_df.loc[traceNum, tName])
            # x_data=np.array(sample_df.loc[traceNum, xName])
            # x_data=x_data-x_data[0]
            x_data = np.linspace(0,len(y_data)-1, len(y_data))
            popt, pcov = curve_fit(fit_burn, x_data[fit_b:fit_e], y_data[fit_b:fit_e], p0=[.1,.1,.001], maxfev=100000 ) #bounds=([-np.inf, 0, -np.inf], np.inf), )
            corTrace=[]
            offset = 0.0
            if (retainOffset == True):
                offset = np.mean(y_data[fit_b:fit_e])

            # print(len(x_data), len(y_data))
            for i, x in enumerate(y_data):
                vc = y_data[i] - fit_burn(x_data[i], popt[0],popt[1],popt[2]) + offset
                corTrace.append(vc)
            # print(str(wl) + newSuffix)
            # print(corTrace)
            sample_df[str(wl) + newSuffix].loc[traceNum]=np.array(corTrace)

def burnCorrectionMultipleRegions(sample_df, wavelengths, suffix, baselineArrays, newSuffix = '_bcor', retainOffset = False):


    """
    baselineArrays is a list of lists containing segments to use in the
    burn-in fit. 
    Each list with a begin and end point

    """

    print("de-burning")
    for i in range(len(wavelengths)):
        wl = str(wavelengths[i])
        sample_df[str(wl) + newSuffix] = 0
        sample_df[str(wl) + newSuffix] = sample_df[str(wl) + newSuffix].astype('object')
        tName = wl + suffix
        xName = wl + '_time'
        

        for traceNum in sample_df.index: #range (0, len(sample_df[tName])):
            #print('trace index = ' + str(traceNum))

            #print('trace Name = ' + str(tName))
            y_data=np.array(sample_df.loc[traceNum, tName])
            # x_data=np.array(sample_df.loc[traceNum, xName])
            # x_data=x_data-x_data[0]
            x_data = np.linspace(0,len(y_data)-1, len(y_data))

            fit_x = np.array([])
            fit_y = np.array([])

            for baselineArray in baselineArrays:
                fit_b=baselineArray[0]
                fit_e=baselineArray[1]
                fit_x = np.concatenate((fit_x[fit_b:fit_e], x_data[fit_b:fit_e]))                
                fit_y = np.concatenate((fit_y[fit_b:fit_e], y_data[fit_b:fit_e]))

            popt, pcov = curve_fit(fit_burn, fit_x, fit_y, p0=[.1,.1,.001], maxfev=100000 ) #bounds=([-np.inf, 0, -np.inf], np.inf), )
            corTrace=[]
            offset = 0.0
            if (retainOffset == True):
                offset = np.mean(fit_y)

            # print(len(x_data), len(y_data))
            for i, x in enumerate(y_data):
                vc = y_data[i] - fit_burn(x_data[i], popt[0],popt[1],popt[2]) + offset
                corTrace.append(vc)
            # print(str(wl) + newSuffix)
            # print(corTrace)
            sample_df[str(wl) + newSuffix].loc[traceNum]=np.array(corTrace)
         


def smoothTracesWL(sample_df, wavelengths, suffix, smooth_window=50, newSuffix = '_smooth'):

    for i in range(len(wavelengths)):
        
        wl = str(wavelengths[i])
        tName = wl + suffix
        newColName= tName + newSuffix
        sample_df[newColName] = 0
        sample_df[newColName] = sample_df[newColName].astype('object')
        for traceNum in sample_df.index:  # cycle through each trace in experiment 
            subtrace=sample_df.loc[traceNum, tName]
            smoothedTrace=savgol_filter(subtrace , smooth_window, 3)
            sample_df[newColName][traceNum]=smoothedTrace

def smoothTraces(df,tName,smooth_window=25,newSuffix = '_smooth'):
    """
     smoothTraces(dataFrame, name of the column holding the traces to be smoothed, smooth wondow (must be odd), suffix for the new column name)
    """
    if (smooth_window % 2) == 0:
        smooth_window = smooth_window+1
    newColName= tName + newSuffix
    df[newColName] = 0
    df[newColName] = df[newColName].astype('object')
    for traceNum in range(len(df[tName])):  # cycle through each trace in experiment
        subtrace=df[tName].iloc[traceNum]
        smoothedTrace=savgol_filter(subtrace , smooth_window, 3)
        df[newColName].iloc[traceNum]=smoothedTrace
    
def subtractStraightLine(df,y_column_name, beg_points_array,end_points_array, newSuffix = '_sub',x_column_name=None):
    
    """
     subtractStraightLine(dataFrame, name of the column holding the x-values of the traces, 
                name of the column holding the y-values of the traces,
                [beg, end time point indexes for the left-hand point for subtraced line],
                [beg, end time point indexes for the RIGHT hand point of line. Values within range will be averaged], 
                suffix for the new column name)
    """
    newColName= y_column_name + newSuffix
    df[newColName] = 0
    df[newColName] = df[newColName].astype('object')
    for traceNum in range(len(df[y_column_name])):  # cycle through each trace in experiment
        
        subtrace_y=np.array(df[y_column_name].iloc[traceNum])
        subtrace_y = subtrace_y - np.mean(subtrace_y[beg_points_array[0]:beg_points_array[1]])    #subtract the first point so that the y-intercept will be zero
        if (x_column_name==None):
            subtrace_x=np.linspace(0,len(subtrace_y), len(subtrace_y))
        else:
            subtrace_x=df[x_column_name].iloc[traceNum]

        x1 = np.mean(subtrace_x[beg_points_array[0]:beg_points_array[1]])
        x2 = np.mean(subtrace_x[end_points_array[0]:end_points_array[1]])

        y1 = np.mean(subtrace_y[beg_points_array[0]:beg_points_array[1]])
        y2 = np.mean(subtrace_y[end_points_array[0]:end_points_array[1]])
        slope = (y2-y1)/(x2-x1)
        subtracted_trace = []
        for i in range(len(subtrace_y)):
            y_offset = slope*subtrace_x[i]
            subtracted_trace.append(subtrace_y[i]-y_offset)
        df[newColName].iloc[traceNum]=subtracted_trace

def generateDAS(sample_df, wavelengths, suffix='_rDA', newColName='das', selected_indexes = []):
    """
    Generate a series of decay associated spectra from a series of kinetic traces, each taken at a differnet wavelength.
    
    This is for traces with multiple measuring LEDs. 
    IMPORTANT: It is NOT for data with traces taken with different m easuring LEDs at different wavelengths, e.g. as in Charlie results.

    generateDAS(sample_df (the dataFrame), wavelengths (a list of wavelengths), suffix='_rDA' 
                    (the name of the columns that contain the wavelength with the given suffix e.g. 505_rDA, newColName='das'),
                    a new suffix for the DAS series.

    selected_indexes is a list of indexes to use

    """

    sample_df[newColName] = 0
    sample_df[newColName] = sample_df[newColName].astype('object')
    # global test_wl 
    # test_wl = []
    selected_indexes = list(selected_indexes)
    if (selected_indexes == []):
        use_index = sample_df.index
        print("using all indexes: " + str(use_index))

    else:
        use_index = selected_indexes
        print("using selected_indexes :" + str(use_index))
    for t_index in use_index: 
        # print(t)
        dass=[]
        wl = str(wavelengths[0]) #start with the first wl
        # print(len(sample_df[wl + suffix][t_index]))
        # try:
        if (type(sample_df[str(wl) + suffix][t_index]) == np.nan):
            print("Nan")
        else:
            for pt in range(len(sample_df[wl + suffix][t_index])): #cycle through all points in the trace
                das=[]
                #test_wl = []
                for wl in range(len(wavelengths)):
                    das.append(sample_df[str(wavelengths[wl]) + suffix][t_index][pt])
                    #test_wl.append(wavelengths[wl])
                dass.append(das)
        # except:
        #     pass
        # print(test_wl)
        sample_df[newColName][t_index] = dass
    return sample_df


def generateDAS_baseline_subtracted(sample_df, wavelengths, suffix='_rDA', newColName='das'):
    """
    Generate a series of decay associated spectra from a series of kinetic traces, each taken at a differnet wavelength.
    
    In this version, a baseline is subtraced between the points with the lowest and highest wavelenghth

    This is for traces with multiple measuring LEDs. 
    IMPORTANT: It is NOT for data with traces taken with different m easuring LEDs at different wavelengths, e.g. as in Charlie results.

    generateDAS(sample_df (the dataFrame), wavelengths (a list of wavelengths), suffix='_rDA' 
                    (the name of the columns that contain the wavelength with the given suffix e.g. 505_rDA, newColName='das'),
                    a new suffix for the DAS series.


    """

    sample_df[newColName] = 0
    sample_df[newColName] = sample_df[newColName].astype('object')

    for t in sample_df.index: 
        dass=[]
        wl = str(wavelengths[0]) #start with the first wl
        for pt in range(len(sample_df[wl + suffix][t])): #cycle through all points in the trace
            das=[]
            for wl in range(len(wavelengths)):
                das.append(sample_df[str(wavelengths[wl]) + suffix][t][pt])
            das = np.array(das)  # convert to an array
            das = das - das[0]   # subtract off the value of the first point
            das_at_max_wl = das[-1]  #assuming that the wavelengths are in ascending order!!
            # print(das_at_max_wl)
            max_wl_change = wavelengths[-1] - wavelengths[0]
            slope = das_at_max_wl/max_wl_change
            # print(slope)
            for i in range(1,len(das)): 
                das[i] = das[i] - slope * (wavelengths[i]- wavelengths[0])
            
            dass.append(das)

        sample_df[newColName][t] = dass

def fit_spec_5_old(x,a,b,c,d,e): #function to fit the spectroscopic data
    global ecs
    global qE
    global Zx
    global scatter
    global drift
    
    # Expecting absorbance data from 5 wavelengths per time point
    # The component spectra also have 7 wavelength
    #print(len(components[0]),len(components[1]),len(components[2]),len(components[3]))
    fit_spectrum = a*ecs + b*qE + c*Zx + d*scatter #+ e*drift

    return (fit_spectrum)

def fit_spec_5_4_old(x,a,b,c,d,e): #function to fit the spectroscopic data
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

def fit_spec_5_4_genotype(xxspecs,a,b,c,d,e): #function to fit the spectroscopic data
    
    # xxecs, xxqE, xxZx, xxscatter, xxdrift = xxspecs

    # Expecting absorbance data from 5 wavelengths per time point
    # The component spectra also have 5 wavelength
    fit_spectrum = a*xxspecs[0] + b*xxspecs[1] + c*xxspecs[2] + d*xxspecs[3] + e*xxspecs[4]

    return (fit_spectrum)

def fit_spec_5_4_var(x,a,b,c,d,e): #function to fit the spectroscopic data
    global ecs
    global qE
    global Zx
    global scatter_var
    global drift
    global drift_var
    # Expecting absorbance data from 5 wavelengths per time point
    # The component spectra also have 7 wavelength
    #print(len(components[0]),len(components[1]),len(components[2]),len(components[3]))
    fit_spectrum = a*ecs + b*qE + c*Zx + d * scatter_var + e*drift_var
    return (fit_spectrum)


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
    ecs, qE, Zx, scatter, drift = x
    # print(x,a,b,c,d,e)
    # print(type(ecs))
    # print(type(qE))
    # print(type(Zx))
    # print(type(scatter))
    # print(type(drift))
    # Expecting absorbance data from 5 wavelengths per time point
    # The component spectra also have 7 wavelength
    #print(len(components[0]),len(components[1]),len(components[2]),len(components[3]))
    fit_spectrum = a*ecs + b*qE + c*Zx + d*scatter + e*drift
    # print(fit_spectrum)  
    return (fit_spectrum)

def fit_spec_5_Thekla(x,a,b,c,d,e): #function to fit the spectroscopic data
    global ecs
    global qE
    global Zx
    global scatter
    global drift
    
    # Expecting absorbance data from 5 wavelengths per time point
    # The component spectra also have 7 wavelength
    #print(len(components[0]),len(components[1]),len(components[2]),len(components[3]))
    fit_spectrum = a*ecs + b*qE + c*Zx + d*scatter #+ e*drift
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



def fitDAS_5(sample_df, wavelengths, dasColName, specs, newSuffix='', zeroEnds = False, ):
    components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
    # global ecs
    # global qE
    # global Zx
    # global drift
    # global scatter

    if (zeroEnds == True):
        ecs=ecs-np.linspace(ecs[0],ecs[-1],len(ecs))
        qE=qE-np.linspace(qE[0],qE[-1],len(qE))
        Zx=Zx-np.linspace(Zx[0],ecs[-1],len(Zx))
        drift=drift-np.linspace(drift[0],drift[-1],len(drift))
        scatter=scatter-np.linspace(scatter[0],scatter[-1],len(scatter))

    #drift = np.linspace(-1,15,len(wavelengths)) #
    #scatter = np.linspace(10,10,len(wavelengths))
    #ecs = ecs_spec(wavelengths)
    # qE = qE_spec(wavelengths)
    #Zx = Zx_spec(wavelengths)
    # # print(drift)
    for c in components:
        newColName= c + newSuffix
        sample_df[newColName] = 0
        sample_df[newColName] = sample_df[newColName].astype('object')

    for traceNum in sample_df.index: #range (0, len(sample_df[dasColName])):  # cycle through each trace in experiment 
        comp={}
        for c in components:
            comp[c]=[]
        # plt.figure()
        p0=[0,0,0,0,0]
        for time_index in range(len(sample_df[dasColName][sample_df.index[0]])):
            y_data=sample_df[dasColName][traceNum][time_index]
            if (zeroEnds == True):
                y_data=y_data-np.linspace(y_data[0],y_data[-1],len(y_data))
            # plt.plot(wavelengths, y_data)
             
            popt, pcov = curve_fit(fit_spec_5, specs, y_data, p0=[0,0,0,0,0])
            p0=[popt[0],popt[1],popt[2],popt[3] ]
            for i, c in enumerate(components):
                comp[c].append(popt[i])
        # plt.show()
        for c in components:
            newColName= c + newSuffix
            sample_df[newColName][traceNum] = comp[c]

# def fitDAS_6(sample_df, wavelengths, dasColName, newSuffix=''):
#     components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
#     global ecs
#     global qE
#     global Zx
#     global drift
#     global scatter
#     # drift = drift_spec(wavelengths)
#     # scatter = scatter
#     #np.linspace(1,5,len(wavelengths))
#     # ecs = ecs_spec(wavelengths)
#     # qE = qE_spec(wavelengths)
#     # Zx = Zx_spec(wavelengths)
#     # # print(drift)
#     for c in components:
#         newColName= c + newSuffix
#         sample_df[newColName] = 0
#         sample_df[newColName] = sample_df[newColName].astype('object')

#     for traceNum in sample_df.index: #range (0, len(sample_df[dasColName])):  # cycle through each trace in experiment 
#         comp={}
#         for c in components:
#             comp[c]=[]
       
#         for time_index in range(len(sample_df[dasColName][sample_df.index[0]])):
#             y_data=sample_df[dasColName][traceNum][time_index]
#             popt, pcov = curve_fit(fit_spec_, wavelengths, y_data, p0=[0,0,0,0,0])
#             for i, c in enumerate(components):
#                 comp[c].append(popt[i])
            
#         for c in components:
#             newColName= c + newSuffix
#             sample_df[newColName][traceNum] = comp[c]

def fitDAS_5_old(sample_df, wavelengths, dasColName, newSuffix=''):
    components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
    global ecs
    global qE
    global Zx
    global drift

    scatter = np.array([11. ,  6. ,  2.5,  1. ,  0.7,  0.5]) 
    drift = np.linspace(0,1,len(wavelengths)) +1

    ecs = np.array([ 0.16831426, -0.06476744,  0.51011528,  1.0,  0.52315683,0.30481469]) #ecs_spec(wavelengths)
    qE = np.array([2.86328688e-04, 5.47754353e-03, 1.14799114e-01, 5.33707544e-01, 9.38366272e-01, 7.96077155e-01])
        #qE_spec(wavelengths)
    qE = qE + (drift * .3)
    Zx = np.array([0.2785128, 0.40347622, 1.0, 0.46877969, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)
    # print(drift)
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
            popt, pcov = curve_fit(fit_spec_5_old, wavelengths, y_data, p0=[0,0,0,0,0])
            
            for i, c in enumerate(components):
                comp[c].append(popt[i])
            
        for c in components:
            newColName= c + newSuffix
            sample_df[newColName][traceNum] = comp[c]

def fitDAS_5_4_old(sample_df, wavelengths, dasColName, newSuffix=''):
    components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
    global ecs
    global qE
    global Zx
    global drift
    global scatter
    global fudge
    #fudge = 1
    # scatter = np.array([11. ,  6. ,  2.5,  1. ,  0.7,  0.5]) 
    scatter = np.array([.8 ,  .7 ,  0.6,  0.5 ,  .4,  .3]) 
    drift = np.linspace(0, 1, len(wavelengths)) + fudge
    # drift = 10*np.array([0.01009748, 0.0576072 , 0.05214349, 0.04538453, 0.03669655,0.0339339 ])
    # drift = drift + np.array([0.01265002, 0.0553461 , 0.04891906, 0.04261369, 0.03844194,0.02890673])

    ecs = np.array([ 0.16831426, -0.06476744,  0.51011528,  1.0,  0.52315683,0.30481469]) #ecs_spec(wavelengths)
    qE = np.array([2.86328688e-04, 5.47754353e-03, 1.14799114e-01, 4.33707544e-01, 9.38366272e-01, 7.96077155e-01])
        #qE_spec(wavelengths)
    # qE = qE + (drift * .3)


    Zx = np.array([0.2785128, 0.40347622, 1.0, 0.56877969, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)
    # if (fudge == 1):
    #     qE = np.array([2.86328688e-04, 5.47754353e-03, 1.14799114e-01, 4.33707544e-01, 9.38366272e-01, 7.96077155e-01])
        # Zx = np.array([0.2785128, 0.70347622, 1.0, 0.36877969, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)
        
    # print(drift)
    # qE=qE+0.3
    
    # Zx[3] =  Zx[3] + fudge
    arabid = 1

    # Zx[4] = Zx[4] + 0.1
    if (arabid == 1):
        qE[2] = 0 #fudge
        qE[1] = 0.05
        qE[3] = 0.4

    # ecs= ecs + fudge   #-.02
    # drift=drift-fudge
    # scatter=scatter+fudge

    zero_475 = 1

    if (zero_475 == 1):
        qE[0] = 0
        ecs[0]=0
        Zx[0]=0
        scatter[0]=0
        drift[0]=0


    for c in components:
        newColName= c + newSuffix
        sample_df[newColName] = 0
        sample_df[newColName] = sample_df[newColName].astype('object')

    for traceNum in sample_df.index:  # cycle through each trace in experiment 
        comp={}
        for c in components:
            comp[c]=[]
        p0=[0,0,0,0,0]
        for time_index in range(0, len(sample_df.loc[traceNum, dasColName])):
            y_data=sample_df.loc[traceNum, dasColName][time_index]
            if (zero_475 == 1):
                y_data[0]=0
            popt, pcov = curve_fit(fit_spec_5_4_old, wavelengths, y_data, p0=p0)
            p0 = popt
            for i, c in enumerate(components):
                comp[c].append(popt[i])
            
        for c in components:
            newColName= c + newSuffix
            sample_df[newColName][traceNum] = comp[c]

def fitDAS_5_Thekla(sample_df, wavelengths, dasColName, newSuffix=''):
    components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
    global ecs
    global qE
    global Zx
    global drift

    drift = np.linspace(1,1,len(wavelengths))
    ecs = ecs_spec_Thekla(wavelengths)
    qE = qE_spec_Thekla(wavelengths)
    Zx = Zx_spec_Thekla(wavelengths)
    # print(drift)
    for c in components:
        newColName= c + newSuffix
        sample_df[newColName] = 0
        sample_df[newColName] = sample_df[newColName].astype('object')

    for traceNum in sample_df[dasColName].index:  # cycle through each trace in experiment 
        comp={}
        for c in components:
            comp[c]=[]
        
        # sample_df.sort_values('start_time')

        for time_index in range(len(sample_df[dasColName][traceNum])):
            # print(time_index)
            try:
                y_data=sample_df[dasColName][traceNum][time_index]
                popt, pcov = curve_fit(fit_spec_5_Thekla, wavelengths, y_data, p0=[0,0,0,0,0])
                for i, c in enumerate(components):
                    comp[c].append(popt[i])
            except:
                pass
        for c in components:
            newColName= c + newSuffix
            sample_df[newColName][traceNum] = comp[c]

def fitDAS_5_4_arabidopsis(sample_df, wavelengths, dasColName, newSuffix='', thekla = 1, arabid = 1, s_offset = 0.0, 
                            zero_475 = 1, das_baseline = -1, d_offset = 0.0, genotype = ''):
    components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
    global ecs
    global qE
    global Zx
    global drift
    global scatter
    global fudge
    # global s_offset

    #fudge = 0.0
    # scatter = np.array([11. ,  6. ,  2.5,  1. ,  0.7,  0.5]) 
    scatter = np.array([.8 ,  .7 ,  0.6,  0.5 ,  .4,  .3]) 
    drift = np.linspace(0,1,len(wavelengths)) + s_offset
    # drift = 10*np.array([0.01009748, 0.0576072 , 0.05214349, 0.04538453, 0.03669655,0.0339339 ])
    # drift = drift + np.array([0.01265002, 0.0553461 , 0.04891906, 0.04261369, 0.03844194,0.02890673])

    ecs = np.array([ 0.16831426, -0.06476744,  0.51011528,  1.0,  0.52315683,0.30481469]) #ecs_spec(wavelengths)
    qE = np.array([2.86328688e-04, 5.47754353e-03, 1.14799114e-01, 4.33707544e-01, 9.38366272e-01, 7.96077155e-01])
        #qE_spec(wavelengths)
    # qE = qE + (drift * .3)


    Zx = np.array([0.2785128, 0.40347622, 1.0, 0.56877969, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)
    # if (fudge == 1):
    #     qE = np.array([2.86328688e-04, 5.47754353e-03, 1.14799114e-01, 4.33707544e-01, 9.38366272e-01, 7.96077155e-01])
        # Zx = np.array([0.2785128, 0.70347622, 1.0, 0.36877969, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)
        
    # print(drift)
    # qE=qE+0.3
    
    # Zx[3] =  Zx[3] + fudge
    
# [ 0.         -0.20850231  0.38283129  1.          0.55103526  0.20525645]
# [0.         0.00547754 0.11479911 0.53370754 0.93836627 0.79607715]
# [ 0.          0.12928294  1.          0.44407102  0.05743463 -0.03206746]
# [0.   1.   0.85 0.7  0.6  0.45]
# [0. 1. 1. 1. 1. 1.]

    # Zx[4] = Zx[4] + 0.1
    if (arabid == 1):
        # qE = np.array([2.86328688e-04, 5.47754353e-03, 1.14799114e-01, 4.33707544e-01, 9.38366272e-01, 7.96077155e-01])
        qE = np.array([0.002, 0.05, 0.1, 0.45, 1.0, 0.85])
        drift = np.array([0.6, 0.60, 0.7, 0.8, 0.9, 1.0])
        scatter = np.array([2.0, 1.5, 0.85, 0.7, 0.6, 0.45])
        ecs = np.array([ 0.16831426, -0.40,  0.4,  1.0,  0.35, 0.1]) #ecs_spec(wavelengths)
        Zx = np.array([0.2785128, 0.40347622, 1.0, 0.4, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)

        # qE[2] = 0.1 #fudge
        # qE[1] = -0.15
        # qE[3] = 0.4
        # qE[5] = 0.8
        # Zx[3] = 0.4
        # ecs[1] = -0.4
        # ecs[2] = 0.4
        # ecs[4] = 0.5
        # ecs[5] = 0.2
        
        # scatter = scatter + 0.1

    # ecs= ecs + fudge   #-.02
        # drif[5] = 
    # scatter=scatter+fudge
    if (thekla == 1):
        print('t:', end='')
        # drift = np.linspace(1,1,len(wavelengths))
        # ecs = ecs_spec_Thekla(wavelengths)
        # qE = qE_spec_Thekla(wavelengths)
        # Zx = Zx_spec_Thekla(wavelengths)

        # S: 0.0, 1.0, 0.7,0.60,0.45
        # D: 0.0, 0.6, 0.7, 0.85, 1.0
        # drift = np.array([0.6, 0.60, 0.7, 0.8, 0.9, 1.0]) + fudge
        scatter = np.array([2.0, 1.5, 1.04, 0.7, 0.5, 0.45]) #+ fudge

        # JUST CHANGED:
        scatter = np.array([0.5, 2.0, 0.7, 0.6, 0.5, 0.45]) + s_offset

        # drift = np.array([0.6, 0.60, 0.7, 0.8, 0.9, 1.0]) + fudge
        # scatter = np.array([1.0, 1.0, 0.7, 0.6, 0.5, 0.45]) + fudge

        # qE = np.array([0.002, 0.05, 0.1, 0.45, 1.0, 0.85])
        # ecs = np.array([ 0.16831426, -0.40,  0.4,  1.0,  0.35, 0.1]) #ecs_spec(wavelengths)
        # Zx = np.array([0.2785128, 0.40347622, 0.8, 0.3, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)

        drift = np.linspace(1.0,3.0,len(wavelengths)) + d_offset
        # drift = drift 
        # drift = np.array([0.01265002, 0.0553461 , 0.04891906, 0.04261369, 0.03844194,0.02890673])
        
        print(s_offset, drift)
        ecs = ecs_spec_Thekla(wavelengths)
        # ecs = np.array([ 0.16831426, -0.40,  0.4,  1.0,  0.35, 0.1]) #ecs_spec(wavelengths)
        
        # for npq2
        # genoptype = 'npq2'

        
        # for Col

        qE = qE_spec_Thekla(wavelengths)
        #qE = np.array([0.002, 0.05, 0.1, 0.45, 1.0, 0.85])
        Zx = Zx_spec_Thekla(wavelengths)
        
        #from Col:
        if (genotype == 'Col'):
            ecs = 1000*np.array([-0.00112671, -0.00202359,  0.00027301,  0.00227476,  0.00129693,
            0.00047521])
            Zx= 1000*np.array([-0.00012045,  0.00108412,  0.00418004,  0.00189427,  0.00124919,
            0.00125713])

        elif (genotype == 'npq2'):
            ecs = 1000*np.array([-0.000181465, -0.00023744 ,  0.000010198,  0.00184786,  0.00112505,
            0.0003081 ])
            Zx= 1000*np.array([-0.000012045,  0.000108412,  0.00418004,  0.00189427,  0.00124919,
            0.00125713])

        elif (genotype == 'npq1'):
            ecs= np.array([-0.00157078, -0.00172374,  0.00085025,  0.00240309,  0.00156114,
                    0.00074263])

        elif (genotype == 'npqx'):
            Zx= np.array([-3.80720439e-05, 1.48578924e-03,  8.10177959e-03,  4.86061694e-04,
        3.89543452e-05, 2.12303333e-05])
            ecs = 1000*np.array([-0.00112671, -0.00202359,  0.00027301,  0.00227476,  0.00129693,0.00047521])
            
    new_ext = 0
    if (new_ext == 1):
        best_spec = {'ecs': np.array([0, -0.20661602,  0.37254349,  1.04767565,  0.42276887,  0.29680173]),
        'qE': np.array([ 0, 1.96338 , -3.888,  2.262,  6.645,
                -1.151]),
        'Zx': np.array([ 0, 0.1312866 ,  0.99978185,  0.44760946,  0.03740546, -0.01714042]),
        'scat': np.array([0, -4.62484149, -0.16461381,  4.88319402, -1.74890509, -3.3586597 ]),
        'drift': np. array([0, -1.45153232, -0.77336078, -1.09298723,  0.11316727,  0.48537003])}

        qE = best_spec['qE']
        ecs = best_spec['ecs']
        Zx = best_spec['Zx']
        scatter = best_spec['scat']
        drift = best_spec['drift']

    if (das_baseline > -1):
        
        drift = drift - drift[das_baseline]
        ecs = ecs - ecs[das_baseline]
        qE = qE - qE[das_baseline]
        Zx = Zx - Zx[das_baseline]
        scatter = scatter - scatter[das_baseline]

    zero_475 = 0

    if (zero_475 == 1):
        qE[0] = 0
        ecs[0]=0
        Zx[0]=0
        scatter[0]=0
        drift[0]=0

    print("____________________")
    print(repr(ecs), repr(qE), repr(Zx), repr(drift), repr(scatter))
    print("____________________")

    for c in components:
        newColName= c + newSuffix
        sample_df[newColName] = 0
        sample_df[newColName] = sample_df[newColName].astype('object')

    for traceNum in sample_df.index:  # cycle through each trace in experiment 
        comp={}
        for c in components:
            comp[c]=[]
       
        for time_index in range(0, len(sample_df.loc[traceNum, dasColName])):
            y_data=sample_df.loc[traceNum, dasColName][time_index]
            if (das_baseline > -1):
                y_data = y_data  - y_data[das_baseline]
            if (zero_475 == 1):
                y_data[0]=0
            popt, pcov = curve_fit(fit_spec_5_4_old, wavelengths, y_data, p0=[0,0,0,0,0])
            
            for i, c in enumerate(components):
                comp[c].append(popt[i])
            
        for c in components:
            newColName= c + newSuffix
            sample_df[newColName][traceNum] = comp[c]




def fitDAS_5_4_genotype(sample_df, wavelengths, dasColName, genotype_specs, genotype, newSuffix='',  zero_475 = 1):

    components = ['ecsfit', 'qEfit', 'Zxfit', 'scatterfit', 'drift']


    # das_baseline = 0
    # if (das_baseline > -1):
        
    #     drift = drift - drift[das_baseline]
    #     ecs = ecs - ecs[das_baseline]
    #     qE = qE - qE[das_baseline]
    #     Zx = Zx - Zx[das_baseline]
    #     scatter = scatter - scatter[das_baseline]

    # print(genotype)

    xspecs = [ genotype_specs[genotype]['ecs'],
            genotype_specs[genotype]['qE'],
            genotype_specs[genotype]['Zx'],
            genotype_specs[genotype]['scatter'],
            genotype_specs[genotype]['drift']]

    if (zero_475 == 1):
        for ii in range(len(xspecs)):
            xspecs[ii][0] = 0.0
        # qE[0] = 0
        # ecs[0]=0
        # Zx[0]=0
        # scatter[0]=0
        # drift[0]=0
    print(xspecs)
    for c in components:
        newColName= c + newSuffix
        sample_df[newColName] = 0
        sample_df[newColName] = sample_df[newColName].astype('object')

    for traceNum in sample_df.index:  # cycle through each trace in experiment 
        print (type(sample_df[dasColName][traceNum]))
        if isinstance(sample_df[dasColName][traceNum],(list,pd.core.series.Series,np.ndarray)):
            comp={}
            for c in components:
                comp[c]=[]  
        
            
            for time_index in range(0, len(sample_df.loc[traceNum, dasColName])):
                y_data=sample_df.loc[traceNum, dasColName][time_index]

                if (zero_475 == 1):
                    y_data[0]=0
                popt, pcov = curve_fit(fit_spec_5_4_genotype, xspecs, y_data, p0=[0,0,0,0,0])
                
                for i, c in enumerate(components):
                    comp[c].append(popt[i])
                
            for c in components:
                newColName= c + newSuffix
                sample_df[newColName][traceNum] = comp[c]
        else:
            print("not a list or array")




def set_scatter (scat):
    global scatter_var
    scatter_var = scat

def fitDAS_5_4_var_ext(sample_df, wavelengths, dasColName, newSuffix='', thekla = 1, arabid = 1, s_offset = 0.0, 
                            zero_475 = 1, das_baseline = -1):
    components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
    global ecs
    global qE
    global Zx
    global drift
    global scatter_var
    global drift_var 

    global fudge


    drift = np.linspace(1.0,3.0,len(wavelengths)) 
    ecs = ecs_spec_Thekla(wavelengths)
    qE = qE_spec_Thekla(wavelengths)
    Zx = Zx_spec_Thekla(wavelengths)
    scatter_var = np.array([0.32723259, 2.18691524, 0.31623998, 0.70914631, 0.59932656, 0.42981241])
    #drift_var = np.array([ 0.11880038,  1.32202375, -0.6938962 ,  0.8793105 ,  1.74722079, 2.29522145])

    if (das_baseline > -1):
        
        drift = drift - drift[das_baseline]
        ecs = ecs - ecs[das_baseline]
        qE = qE - qE[das_baseline]
        Zx = Zx - Zx[das_baseline]
        scatter_var = scatter_var - scatter_var[das_baseline]

    # zero_475 = 0

    if (zero_475 == 1):
        qE[0] = 0
        ecs[0]=0
        Zx[0]=0
        scatter[0]=0
        drift[0]=0


    for c in components:
        newColName= c + newSuffix
        sample_df[newColName] = 0
        sample_df[newColName] = sample_df[newColName].astype('object')
    
    
    # print(str(drift_var), end=', ')

    for traceNum in sample_df.index:  # cycle through each trace in experiment 
        comp={}
        for c in components:
            comp[c]=[]
       
        for time_index in range(0, len(sample_df.loc[traceNum, dasColName])):
            y_data=sample_df.loc[traceNum, dasColName][time_index]
            if (das_baseline > -1):
                y_data = y_data  - y_data[das_baseline]
            if (zero_475 == 1):
                y_data[0]=0
            popt, pcov = curve_fit(fit_spec_5_4_var, wavelengths, y_data, p0=[0,0,0,0,0])
            
            for i, c in enumerate(components):
                comp[c].append(popt[i])
            
        for c in components:
            newColName= c + newSuffix
            sample_df[newColName][traceNum] = comp[c]






# def fitDAS_5_4_old(sample_df, wavelengths, dasColName, newSuffix=''):
#     components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
#     global ecs
#     global qE
#     global Zx
#     global drift
#     global scatter
#     global fudge
#     fudge = 0
#     # scatter = np.array([11. ,  6. ,  2.5,  1. ,  0.7,  0.5]) 
#     scatter = np.array([.8 ,  .7 ,  0.6,  0.5 ,  .4,  .3]) 
#     drift = np.linspace(0,1,len(wavelengths)) +1
#     # drift = 10*np.array([0.01009748, 0.0576072 , 0.05214349, 0.04538453, 0.03669655,0.0339339 ])
#     # drift = drift + np.array([0.01265002, 0.0553461 , 0.04891906, 0.04261369, 0.03844194,0.02890673])

#     ecs = np.array([ 0.16831426, -0.06476744,  0.51011528,  1.0,  0.52315683,0.30481469]) #ecs_spec(wavelengths)
#     qE = np.array([2.86328688e-04, 5.47754353e-03, 1.14799114e-01, 4.33707544e-01, 9.38366272e-01, 7.96077155e-01])
#         #qE_spec(wavelengths)
#     # qE = qE + (drift * .3)


#     Zx = np.array([0.2785128, 0.40347622, 1.0, 0.56877969, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)
#     # if (fudge == 1):
#     #     qE = np.array([2.86328688e-04, 5.47754353e-03, 1.14799114e-01, 4.33707544e-01, 9.38366272e-01, 7.96077155e-01])
#         # Zx = np.array([0.2785128, 0.70347622, 1.0, 0.36877969, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)
        
#     # print(drift)
#     # qE=qE+0.3
    
#     Zx[3] =  Zx[3] + fudge
#     arabid = 1

#     # Zx[4] = Zx[4] + 0.1
#     if (arabid == 1):
#         qE[2] = 0 #fudge
#         qE[1] = 0.05
#         qE[3] = 0.4

#     # ecs= ecs + fudge   #-.02
#     # drift=drift-fudge
#     # scatter=scatter+fudge

#     zero_475 = 1

#     if (zero_475 == 1):
#         qE[0] = 0
#         ecs[0]=0
#         Zx[0]=0
#         scatter[0]=0
#         drift[0]=0


#     for c in components:
#         newColName= c + newSuffix
#         sample_df[newColName] = 0
#         sample_df[newColName] = sample_df[newColName].astype('object')

#     for traceNum in range (0, len(sample_df[dasColName])):  # cycle through each trace in experiment 
#         comp={}
#         for c in components:
#             comp[c]=[]
       
#         for time_index in range(len(sample_df[dasColName][traceNum])):
#             y_data=sample_df[dasColName][traceNum][time_index]
#             if (zero_475 == 1):
#                 y_data[0]=0
#             popt, pcov = curve_fit(fit_spec_5_4_old, wavelengths, y_data, p0=[0,0,0,0,0])
            
#             for i, c in enumerate(components):
#                 comp[c].append(popt[i])
            
#         for c in components:
#             newColName= c + newSuffix
#             sample_df[newColName][traceNum] = comp[c]

def fitDAS_5_4_old_tobacco(sample_df, wavelengths, dasColName, newSuffix=''):
    components = ['ecs', 'qE', 'Zx', 'scatter', 'drift']
    global ecs
    global qE
    global Zx
    global drift
    global scatter
    global fudge

    # scatter = np.array([11. ,  6. ,  2.5,  1. ,  0.7,  0.5]) 
    scatter = np.array([0 ,  .7 ,  0.6,  0.5 ,  .4,  .3]) 
    # drift = np.linspace(0,1,len(wavelengths)) +1
    drift = np.array([0.00056166, 0.0520648 , 0.04416302, 0.03697218, 0.03030906, 0.02142226])

    drift = drift + fudge #np.array([1,1,1,1,1,1])


    ecs = np.array([ 0.16831426, -0.06476744,  0.51011528,  1.0,  0.52315683,0.30481469]) #ecs_spec(wavelengths)
    qE = np.array([2.86328688e-04, 5.47754353e-03, 1.14799114e-01, 4.33707544e-01, 9.38366272e-01, 7.96077155e-01])
        #qE_spec(wavelengths)
    # qE = qE + (drift * .3)
    Zx = np.array([0.2785128, 0.40347622, 1.0, 0.56877969, 0.12686525, 0.0360355 ] )#Zx_spec(wavelengths)
    # print(drift)
    
    qE[1] =qE[1] + 0.2
    
    zero_to_545 = 0
    if (zero_to_545 == 1):
        ecs=ecs-ecs[-1]
        qE = qE - qE[-1]
        Zx=Zx-Zx[-1]
        scatter=scatter-scatter[-1]
        drift=drift-drift[-1]
    
    # Zx[4] = 0.6 #
    # Zx[3] = Zx[3] + fudge
    zero_475 = 0
    if (zero_475 == 1):
        qE[0] = 0
        Zx[0]=0
        # ecs[2] = ecs[2] + fudge
        ecs[0]=0
        # drift = drift + fudge
        drift[0]=0
        # scatter = scatter + fudge
        scatter[0]=0
        # scatter[1] = scatter[1] + fudge

    
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
            # y_data=y_data-y_data[-1] # subtract higest wavelength value from all values
            # y_data[0]=0
            popt, pcov = curve_fit(fit_spec_5_4_old, wavelengths, y_data, p0=[0,0,0,0,0])
            
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

    
#cyt b spectrum

def cytB_spec(wl):
    cytbWL=[500,505,510,515,520,525,530,535,540,545,550,552.5, 555,560,565,570,575,580]
    cytbWL=np.array(cytbWL)
    cytbWL=cytbWL+6
    cytbSpec = np.array([19,18.5,17,15.3,13.7,13,15.3,17,18,15.5,7.5,2.5, 1.5,13,20,21.5,21,20.5])
    cytbSpec=cytbSpec-cytbSpec[closestToList(cytbWL, 545)]
    cytbSpec=cytbSpec/cytbSpec[closestToList(cytbWL, 563)]
    spec = []
    for i in range(len(wl)):
        spec.append(np.interp(wl[i],cytbWL,cytbSpec))
    spec=np.array(spec)
    return spec
    


#cyt f spectrum

def cytF_spec(wl):
    # cytfWL=[450, 500,505,510,515,520,525,530,535,540,545,550,552.5, 555,560,565,570,575,580]
    # cytfWL=np.array(cytfWL)
    # #cytfWL=cytfWL
    # cytfSpec = np.array([22, 19,18.5,17,15.3,13.7,13,15.3,17,18,15.5,7.5,2.5, 1.5,13,20,21.5,21,20.5])
    # cytfSpec=cytfSpec-cytfSpec[closestToList(cytfWL, 545)]
    # cytfSpec=cytfSpec/cytfSpec[closestToList(cytfWL, 555)]
    cytfWL=np.array([450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510,515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
    cytfSpec=np.array([0.01361136, 0.01988276, 0.02824803, 0.03903336, 0.05245899,0.06857097, 0.08717653, 0.10780354, 
    0.12977418, 0.15289982, 0.18051922, 0.22618189, 0.31385579, 0.44725321, 0.56574545, 0.57799012, 0.47110094, 0.32856016, 0.23591982, 0.35417867,1.04405923, 1.40069565, 0.62199098, 0.13665604, 0.05894475])

    spec = []
    for i in range(len(wl)):
        spec.append(np.interp(wl[i],cytfWL,cytfSpec))
    spec=np.array(spec)
    return spec

#     driftTraceWL=np.array([590, 580, 578, 576, 574, 572, 570, 568, 566, 564, 562, 560, 558, 556, 554, 552, 550, 
#                            548, 546, 544, 542, 540, 538, 536, 534, 532, 530, 528, 526, 524, 522, 520, 518, 516, 
#                            514, 512, 510, 508, 506, 504, 502, 500, 498, 496, 494, 492, 490, 488, 486, 484, 482, 
#                            480, 478, 476, 474, 472, 470, 468, 466, 464, 462] )
#     driftTrace=np.array([-3, -2.23259625, -2.08900144, -1.90676129, -1.675778,   -1.43736771, -1.28648502,
#      -1.15140822, -0.9901862,  -0.8268231,  -0.67734669, -0.55846313, -0.4926361,
#      -0.46835363, -0.45765929, -0.46030131, -0.45289486, -0.43868051, -0.40763071,
#      -0.38775439, -0.35447138, -0.34874509, -0.34048107, -0.34603687, -0.37985102,
#      -0.41839697, -0.48259808, -0.59052032, -0.72087349, -0.92584169, -1.19200849,
#      -1.5119132,  -1.90420981, -2.36021374, -2.88309839, -3.50197487, -4.13192967,
#      -4.82742742, -5.60321121, -6.39996584, -7.09626046, -7.74269712, -8.28250766,
#      -8.68809366, -8.86180847, -8.94101401, -8.87006093, -8.68051578, -8.4522491,
#      -8.26082746, -7.90264827, -7.67191769, -7.37697394, -7.23126009, -7.10643427,
#      -7.24440586, -7.45629293, -7.80926795, -8.14807294, -8.46294016, -8.7906025 ])

#     driftTrace=np.array([ 0.        , -0.10397884,  0.21093104, -0.73842213, -0.16607826,
#                          0.26878166,  0.14159894,  1.        ,  0.52912004, -0.04063416,
#                          -1.02106684, -2.50854863, -2.69910322, -3.26596692, -4.08074846,
#                          -4.29461607, -4.43580393, -4.49831963, -4.35877979, -4.37127129,
#                          -4.3599782 , -4.16608803, -3.76691498, -3.04949697, -2.86542408,
#                          -2.37663864])
#     driftTraceWL=np.array([460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570, 575, 580, 585])

    #     driftTraceWL=np.flip(driftTraceWL)
    #     driftTrace=np.flip(driftTrace)

def drift_spec_from_trace(wl, driftTraceWL, driftTrace):
#     driftTrace=np.array([ 0.        , -0.09785897, -0.086966  , -0.03990309,  0.00588594,0.02904874, -0.07482107, -0.21361617, -0.32485364, -0.45968053,
#                          -0.46807275, -0.29403091, -0.02755454,  0.33847267,  0.67510016, 0.84530737,  0.97552712,  0.99503461,  1.        ,  0.84823227,
#                          0.85294959,  0.92271972,  0.91801522,  0.91141255,  0.91054229, 0.86078758,  0.83688201,  0.73295667,  0.6493151 ,  0.61022771,
#                          0.57411674,  0.56226763,  0.56638073,  0.56789946,  0.53810999,0.42679676,  0.28241103,  0.15982337,  0.06346248, -0.00105235,
#                          -0.06746963, -0.04847341, -0.14903963])
#     driftTraceWL=np.array([464, 467, 470, 473, 476, 479, 482, 485, 488, 491, 494, 497, 500, 503, 506, 509, 512, 515, 518, 521, 524, 527, 530, 533, 536, 539,
#                        542, 545, 548, 551, 554, 557, 560, 563, 566, 569, 572, 575, 578, 581, 584, 587, 590])
    
#     driftTrace=np.array([-1., -0.96627725, -0.90695315, -0.88016554, -0.83421386, 
#                          -0.83003877, -0.85046984, -0.863734  , -0.88017335, -0.86681812,
#                          -0.83408534, -0.7532907 , -0.64931702, -0.54754943, -0.43186477,
#                          -0.33026628, -0.27214651, -0.25390547, -0.24432818, -0.24877121,
#                          -0.255649  , -0.27150303, -0.29383148, -0.31180011, -0.32601945,
#                          -0.33942055, -0.35029493, -0.36257122, -0.37053509, -0.37777462,
#                          -0.37886552, -0.37664543, -0.37134035, -0.36172252, -0.35611495,
#                          -0.3555346 , -0.35997958, -0.36298778, -0.36666854, -0.37249348,
#                          -0.37897204, -0.38831918, -0.40072305])
#     driftTraceWL=np.array([464, 467, 470, 473, 476, 479, 482, 485, 488, 491, 494, 497, 500,
#                             503, 506, 509, 512, 515, 518, 521, 524, 527, 530, 533, 536, 539,
#                             542, 545, 548, 551, 554, 557, 560, 563, 566, 569, 572, 575, 578,
#                             581, 584, 587, 590])
    #driftTrace=driftTrace-np.linspace(driftTrace[0],driftTrace[-1],len(driftTrace))
#     driftTrace=np.array([-0.97930457, -0.84510825, -0.76824325, -0.68557158, -0.65721521,
#                          -0.71411414, -0.80057651, -0.9111386 , -0.99341127, -1.        ,
#                          -0.90233654, -0.7229873 , -0.50413059, -0.30143799, -0.150493  ,
#                          -0.05590374,  0.00378046,  0.03427544,  0.04805041,  0.05288039,
#                          0.05259516,  0.033782  ,  0.        ])
    
#     driftTraceWL=np.array([460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
#     global driftTraces
#     global driftTraceWLs
#     driftTrace=driftTraces
#     driftTracesWL=driftTracesWLs
    driftTrace=driftTrace/np.max(np.abs(driftTrace))
    print(driftTrace)
    spec = []
    for  i in range(len(wl)):
        print(i)
        spec.append(np.interp(wl[i],driftTraceWL,driftTrace))
    spec=np.array(spec)
    return spec

def drift_spec(wl):
#     driftTrace=np.array([ 0.        , -0.09785897, -0.086966  , -0.03990309,  0.00588594,0.02904874, -0.07482107, -0.21361617, -0.32485364, -0.45968053,
#                          -0.46807275, -0.29403091, -0.02755454,  0.33847267,  0.67510016, 0.84530737,  0.97552712,  0.99503461,  1.        ,  0.84823227,
#                          0.85294959,  0.92271972,  0.91801522,  0.91141255,  0.91054229, 0.86078758,  0.83688201,  0.73295667,  0.6493151 ,  0.61022771,
#                          0.57411674,  0.56226763,  0.56638073,  0.56789946,  0.53810999,0.42679676,  0.28241103,  0.15982337,  0.06346248, -0.00105235,
#                          -0.06746963, -0.04847341, -0.14903963])
#     driftTraceWL=np.array([464, 467, 470, 473, 476, 479, 482, 485, 488, 491, 494, 497, 500, 503, 506, 509, 512, 515, 518, 521, 524, 527, 530, 533, 536, 539,
#                        542, 545, 548, 551, 554, 557, 560, 563, 566, 569, 572, 575, 578, 581, 584, 587, 590])
    
#     driftTrace=np.array([-1., -0.96627725, -0.90695315, -0.88016554, -0.83421386, 
#                          -0.83003877, -0.85046984, -0.863734  , -0.88017335, -0.86681812,
#                          -0.83408534, -0.7532907 , -0.64931702, -0.54754943, -0.43186477,
#                          -0.33026628, -0.27214651, -0.25390547, -0.24432818, -0.24877121,
#                          -0.255649  , -0.27150303, -0.29383148, -0.31180011, -0.32601945,
#                          -0.33942055, -0.35029493, -0.36257122, -0.37053509, -0.37777462,
#                          -0.37886552, -0.37664543, -0.37134035, -0.36172252, -0.35611495,
#                          -0.3555346 , -0.35997958, -0.36298778, -0.36666854, -0.37249348,
#                          -0.37897204, -0.38831918, -0.40072305])
#     driftTraceWL=np.array([464, 467, 470, 473, 476, 479, 482, 485, 488, 491, 494, 497, 500,
#                             503, 506, 509, 512, 515, 518, 521, 524, 527, 530, 533, 536, 539,
#                             542, 545, 548, 551, 554, 557, 560, 563, 566, 569, 572, 575, 578,
#                             581, 584, 587, 590])
    #driftTrace=driftTrace-np.linspace(driftTrace[0],driftTrace[-1],len(driftTrace))
    driftTrace=np.array([-0.97930457, -0.84510825, -0.76824325, -0.68557158, -0.65721521,
                         -0.71411414, -0.80057651, -0.9111386 , -0.99341127, -1.        ,
                         -0.90233654, -0.7229873 , -0.50413059, -0.30143799, -0.150493  ,
                         -0.05590374,  0.00378046,  0.03427544,  0.04805041,  0.05288039,
                         0.05259516,  0.033782  ,  0.        ])
    
    driftTraceWL=np.array([460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
    # global driftTraces
#     global driftTraceWLs
#     driftTrace=driftTraces
#     driftTracesWL=driftTracesWLs
    # driftTrace=driftTrace/np.max(np.abs(driftTrace))
    # print(driftTrace)
    spec = []
    for  i in range(len(wl)):
        print(i)
        spec.append(np.interp(wl[i],driftTraceWL,driftTrace))
    spec=np.array(spec)
    return spec

def ecs_spec(wl):
    
#     pulseECSWL=[590, 585,580,575,570,565,560,555,550,545,540,535,530,525,520,515,510,505,500,495,490,485,480,475,470,465]
#     pulseECS = np.array([ 0.        ,  0.03344159,  0.10252101,  0.17498088,  0.28436124,
#              0.37270553,  0.39466636,  0.41549049,  0.52507503,  0.64408459,
#              0.73683989,  0.82697959,  0.93609262,  1.        ,  0.97564624,
#              0.83643035,  0.59682901,  0.42985326,  0.25337081,  0.10998487,
#             -0.02828996, -0.07548181, -0.14295679, -0.13838168,  0.00708065,
#              0.        ])
    
#     ecsWL=np.flip(np.array([590, 585,580,575,570,565,560,555,550,545,540,535,530,525,520,515,510,505,500,495,490,485,480,475,470,465]))
#     ecsWL=np.flip(np.array([580,575,570,565,560,555,550,545,540,535,530,525,520,515,510,505,500,495,490,485,480,475,470,465,460,455]))
#     ecsSpec=np.flip(np.array([ 0.        ,  0.00026393,  0.0005899 ,  0.0011091 ,  0.00146355,
#          0.00225297,  0.00202352,  0.00205593,  0.00257003,  0.0032522 ,
#          0.00366017,  0.00408537,  0.00447852,  0.0049811 ,  0.00479048,
#          0.00389282,  0.00247211,  0.00121837,  0.00052815,  0.00049321,
#         -0.00072343, -0.00070766, -0.00167152, -0.00073296, -0.00063791,
#          0.        ]))
#     ecsSpec=ecsSpec/np.max(ecsSpec)
# #     pulseECSWL=np.flip(pulseECSWL)
# #     pulseECS=np.flip(pulseECS)
#     ecsWL=ecsWL+3
#     ecsSpec=np.array([-0.71277281, -0.43899874, -0.10387849,  0.09269294,  0.13036653,
#                       0.22194588,  0.46461142,  0.80149811,  0.98423798,  1.        ,
#                       0.85163492,  0.65056971,  0.45509381,  0.29812798,  0.1620286 ,
#                       0.05349879,  0.01913389,  0.0399582 ,  0.03149595,  0.        ])
#     ecsSpec=np.array([ 0.00000000e+00,  8.78528548e-02,  2.26741113e-01,  2.40000754e-01,
#                       -3.31505254e-01, -1.16754523e-01, -3.62251730e-01, -3.51848102e-01,
#                       -1.60837185e-01,  2.68709078e-02,  1.76662957e-01,  3.82831285e-01,
#                       6.98142034e-01,  9.27353062e-01,  1.00000000e+00,  8.96634243e-01,
#                       7.39476823e-01,  5.51035261e-01,  3.54006517e-01,  2.05256450e-01,
#                       8.48566119e-02,  3.77208458e-05,  3.90792755e-02,  3.99415059e-02, 5.89143128e-03])
    # ecsSpec=np.array([ 0.00000000e+00,  8.78528548e-02,  2.26741113e-01,  1.40000754e-01,
    #                       -3.01505254e-01, -3.6754523e-01, -3.8e-01, -2.8e-01,
    #                       -1.60837185e-01,  2.68709078e-02,  1.76662957e-01,  3.82831285e-01,
    #                       6.98142034e-01,  9.27353062e-01,  1.00000000e+00,  8.96634243e-01,
    #                       7.39476823e-01,  5.51035261e-01,  3.54006517e-01,  2.05256450e-01,
    #                       9.8e-02,  7.7e-02,  4.3e-02,  3.99415059e-02,5.89143128e-03])

    # ecsWL=np.array([450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
    #print('hello')
    ecsWL=np.array([450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
    # ecsWL = ecsWL - 2
    ecsSpec=np.array([ 0.        ,  0.08785285,  0.22674111,  0.14000075, -0.30150525,
    -0.36754523, -0.38      , -0.28      , -0.16083718,  0.02687091,
    0.17666296,  0.38283129,  0.69814203,  0.92735306,  1.        ,
    0.89663424,  0.73947682,  0.55103526,  0.35400652,  0.20525645,
    0.098     ,  0.077     ,  0.043     ,  0.03994151,  0.00589143])
    # print(len(wl), len(ecsWL), len(ecsSpec))
    spec = []
    for i in range(len(wl)):
        spec.append(np.interp(wl[i],ecsWL,ecsSpec))
    spec=np.array(spec)
    return spec

def Zx_spec(wl):
    #     ZxTraceStd5 = np.array([ 0.00000000e+00, -1.79886312e-04, -3.25289036e-05, -1.15619522e-04,
    #             0.00000000e+00, -1.06225262e-04,  7.33846784e-05,  1.38237125e-04,
    #             4.88201525e-04,  4.73997165e-04,  5.98396301e-04,  1.03501985e-03,
    #             1.21480343e-03,  1.34991721e-03,  1.22544035e-03,  1.25055240e-03,
    #             1.38857424e-03,  1.40426691e-03,  1.76330198e-03,  2.08156458e-03,
    #             2.54884787e-03,  3.17621349e-03,  3.54474194e-03,  3.87742627e-03,
    #             4.03303620e-03,  4.09680338e-03,  4.26200253e-03,  4.38726155e-03,
    #             4.79697190e-03,  5.14672030e-03,  5.38124169e-03,  5.46379428e-03,
    #             5.73546494e-03,  6.24380746e-03,  6.97425469e-03,  7.44846748e-03,
    #             8.11431609e-03,  8.66275501e-03,  8.89903164e-03,  9.04655284e-03,
    #             9.69320636e-03,  1.00677371e-02,  9.93035609e-03,  9.36338487e-03,
    #             8.81188930e-03,  7.52567521e-03,  5.94754468e-03,  4.79249617e-03,
    #             3.39647970e-03,  1.77938303e-03,  9.12544997e-04, -8.01840151e-05,
    #            -5.84768300e-04, -1.46712032e-03, -1.47155542e-03, -1.13722256e-03,
    #            -1.16044516e-03, -5.59467011e-04,  8.28921536e-04,  3.41894540e-04,
    #             1.04732335e-03,  3.48071832e-04,  8.47159471e-04,  6.40058052e-04,
    #             0.00000000e+00])
    #     ZxTraceStd5=ZxTraceStd5-np.linspace(ZxTraceStd5[0],ZxTraceStd5[-1],len(ZxTraceStd5))
    #     ZxTraceStd5=ZxTraceStd5/np.max(ZxTraceStd5)
    #     ZxTraceStd5=np.flip(ZxTraceStd5)

    #     wlZxStd5=[590, 588,586,584,582,580,578,576,574,572,570,568,566,564,562,560,558,556,554,552,550,
    #      548,546,544,542,540,538,536,534,532,530,528,526,524,522,520,518,516,514,512,510,508,
    #      506,504,502,500,498,496,494,492,490,488,486,484,482,480,478,476,474,472,470,468,466,464,462]

    #     wlZxStd5=np.flip(wlZxStd5)
    ZxSpec=np.array([ 0.        , -0.03599246,  0.02244886,  0.00645265, -0.03741914,
         0.00254552,  0.21377456,  0.44296526,  0.80471238,  1.        ,
         0.89261038,  0.7168024 ,  0.44407102,  0.23854102,  0.14500318,
         0.05743463,  0.00804195, -0.03206746, -0.05262834, -0.06410261,
        -0.06272308, -0.04705852, -0.01548376, -0.02279776,  0.00120181,
         0.00533344])
    ZxWL=np.array([460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520,
        525, 530, 535, 540, 545, 550, 555, 560, 565, 570, 575, 580, 585])
    # ZxWL=ZxWL+2
    spec = []
    for i in range(len(wl)):
        spec.append(np.interp(wl[i],ZxWL,ZxSpec))
    spec=np.array(spec)
    return spec

def qE_spec(wl):
    #qEWL=[460, 475,488,505,520,530, 535, 545, 560, 590]
    #qESpectrum=np.array([0, -.01,-.05,.5,2.9,4.9, 5.2, 4.2, 2.5, 0])
    #qESpectrum=qESpectrum/np.max(qESpectrum)
    # b=453
    # e=580
    # qEWL=np.linspace(b,e,25)
    # qESpectrum = [0,-.1,-.5,-4,0,0,.5,4,3,-2,-4,2,5,12,18,22,22,19.5,15,11,8,5,5,4,0]
    # #qESpectrum = [0,-.1,-.5,-.4,0,0,.5,.4,.3,-.2,-.4,2,5,12,18,22,22,19.5,15,11,8,5,5,4,0]
    # qESpectrum=np.array(qESpectrum)
    # qESpectrum=qESpectrum/np.max(qESpectrum)
    qEWL=np.array([450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
    qESpectrum=np.array([-5.78085357e-04, -6.78553127e-05, -5.04221075e-06,  
    6.01235341e-06, 3.10455355e-05,  1.33747482e-04,  5.15200497e-04,  1.77585925e-03,
    5.47754353e-03,  1.51184821e-02,  3.73401108e-02,  8.25254823e-02, 
    1.63209561e-01,  2.88833911e-01,  4.57399644e-01,  6.48169395e-01,
    8.21913809e-01,  9.32629376e-01,  9.46971617e-01,  8.60418801e-01,
    6.99564686e-01,  5.08968521e-01,  3.31359193e-01,  1.93042167e-01,
    1.00635313e-01])
    qEWL=qEWL-2
    spec = []
    for i in range(len(wl)):
        spec.append(np.interp(wl[i],qEWL,qESpectrum))
    spec=np.array(spec)
    return spec


# def get_extintion_coefficients(component_name, wl):
    

def calculateAbsorbanceChangesFromScans(sample_df, colName, baselineTraceIndex, newColName, offset=0):

    #calculate the deltaA values and store in a new column
    #use baselineTraceIndex to indicate which trace is the baseline
    
    # set up the new column to hold results.
    
    sample_df[newColName] = 0
    sample_df[newColName] = sample_df[newColName].astype('object')

    for i in range(0,len(sample_df)):
        #print(sample_df[colName].iloc[i][0], sample_df[colName].iloc[baselineTraceIndex][0])
        sample_df[newColName].iloc[i] = -1*np.log10((sample_df[colName].iloc[i]-offset)/(sample_df[colName].iloc[baselineTraceIndex]-offset))

# def calculate_phi2(sample_df, fluorescence_col, Fs_begin, Fs_end,Fmp_beg, Fmp_end):
#     for index in sample_df.index:
#         trace = sample_df.loc[index,fluorescence_col]
#         Fs = np.mean(trace[Fs_begin:Fs_end])
#         Fmp = np.mean(trace[Fs_begin:Fs_end])


def actinicProfile(sample_df, colName = 'actinic_intensity'):
    
    """
    Return the segments with the same actinic intensity, store in a list of lists designating the start, end indexes and the intensities
    The actinic intensity values are usually in the column with colName = 'actinic_intensity', but this can be overridden
    """

    
    actinic_segments = []
    intensities = []
    actTime=[]
    
    current_actinic_intensity=np.sum(sample_df[colName].iloc[0])
    #print(current_actinic_intensity)
    actinic_segments.append(0)            
    intensities.append(current_actinic_intensity)
    actTime.append(sample_df["start_time"].iloc[0])
    for i in range(1,len(sample_df)):
        intensities.append(np.sum(sample_df[colName].iloc[i]))
        actTime.append(sample_df["start_time"].iloc[i])
        if (np.sum(sample_df[colName].iloc[i]) != current_actinic_intensity):
            current_actinic_intensity=np.sum(sample_df[colName].iloc[i])
            actinic_segments.append(i)
            
            
    
    actinic_segments.append(len(sample_df))  #append the end of the segment
    #intensities.append(sample_df[colName].iloc[len(sample_df)-1])  #append the end of the segment
    
    actinicProfile=[np.array(intensities),actinic_segments,np.array(actTime) ]
#     for i in range(0,len(actinic_segments)-1):
#         actinicProfile.append([actinic_segments[i]+1, intensities[i]])
#         actinicProfile.append([actinic_segments[i+1]+1, intensities[i]])
    return actinicProfile

def avPulsesPerWL(sample_df, colName, newColName, moveEvery, settleSteps):
    #     newColName = "meas0_I_av"
    sample_df[newColName] = 0
    sample_df[newColName] = sample_df[newColName].astype('object')


    for i in range(0,len(sample_df[colName])):
        #print(i)
        avTrace=[]  # holds new trace consisting of average of I values at each wavelength
        for j in range(0, int(len(sample_df[colName].iloc[0])/moveEvery)):  #cycle through the WLs to get average at each
            avTrace.append(np.mean(sample_df[colName].iloc[i][j*moveEvery+settleSteps:(j+1)*moveEvery]))

        sample_df[newColName][i] = np.array(avTrace)  #store avTrace in the dataFrame

def fit_spec_4(wl,a,b,c,d,e,f): #function to fit the spectroscopic data


    global ecs
    global qE
    global Zx
    global drift
    global cytF
    global cytB
    #print(cytF)
    global bwli
    global ewli
    bwli=0
    ewli=len(wl)

    xecs=ecs[bwli:ewli]
    xqE=qE[bwli:ewli]
    xZx=Zx[bwli:ewli]
    xdrift=drift[bwli:ewli]
    xcytF=np.array(cytF[bwli:ewli])
    xcytB=np.array(cytB[bwli:ewli])
    
    xecs=xecs-np.linspace(xecs[0],xecs[-1],len(xecs))
    xqE=xqE-np.linspace(xqE[0],xqE[-1],len(xqE))
    xZx=xZx-np.linspace(xZx[0],xZx[-1],len(xZx))
    #xdrift=xdrift-np.linspace(xdrift[0],xdrift[-1],len(xdrift)) + .01
    xcytF=xcytF-np.linspace(xcytF[0],xcytF[-1],len(xcytF))
    xcytB=xcytB-np.linspace(xcytB[0],xcytB[-1],len(xcytB))
    #a=0
    b=0
    #c=0
    d=0
    #e=0
    #f=0
    # Expecting absorbance data from 5 wavelengths per time point
    # The component spectra also have 7 wavelength
    #print(len(components[0]),len(components[1]),len(components[2]),len(components[3]))
    fit_spectrum = a*xecs + b*xqE + c*xZx + d*xdrift + e*xcytF + f*xcytB
    return (fit_spectrum)



def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def gaussianSpec(wl,center,width):
    # wl is a an input, list or array of wavelength values for which data is desired

    gSpec=[]
    for w in wl: #w will be set, sequentially, to each wavelength in the list of wavelegngths, wl
        #print(wl)
        gSpec.append(gaussian(w,center,width))
    return np.array(gSpec)

# def model3(wl, coeffs):
#     global ecs
#     global Zx
#     global unkSpec
    
#     fitSpec = coeffs[0] *  unkSpec #gaussianSpec(wl, coeffs[1], coeffs[2] )
#     fitSpec = fitSpec + coeffs[3] * ecs 
#     fitSpec = fitSpec + coeffs[4] * Zx 
#     return fitSpec # * np.exp( - ((t-coeffs[2])/coeffs[3])**2)



def model4(wl, coeffs):
    global ecs
    global Zx
    global qE
    global unkSpec
    
    fitSpec = coeffs[0] *  unkSpec #gaussianSpec(wl, coeffs[1], coeffs[2] )
    fitSpec = fitSpec + coeffs[3] * ecs 
    fitSpec = fitSpec + coeffs[4] * Zx 
    fitSpec = fitSpec + coeffs[5] * qE 
    return fitSpec # * np.exp( - ((t-coeffs[2])/coeffs[3])**2)

# def residualsLocal4(coeffs, y, t):
#     global unkSpec
#     global globSpecSet
# #     unkSpec = coeffs[0] * gaussianSpec(wl, coeffs[1], coeffs[2] )
#     r = y - model4(t, coeffs)
#     return r

# def residualsGlobal4(coeffs, y, wl):
#     global unkSpec
#     global globSpecSet
#     unkSpec = gaussianSpec(wl, coeffs[1], coeffs[2] )
#     unkSpec = unkSpec - unkSpec[-1]
#     residualVals = np.array([0] * len(wl))
#     for index, y in enumerate(globSpecSet):
#         lx0=[coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4], coeffs[5]]
#         x = scipy.optimize.least_squares(residualsLocal, lx0, args=(y, wl))
#         r = np.array(y - model4(wl, x.x))
#         residualVals = residualVals + r #np.abs(r)/len(globSpecSet)
#     return residualVals



# def residualsGlobal3_ecs_f_unk(coeffs, y, wl, ecsY, cytFY, unkSpecY):

#     unkSpecY = gaussianSpec(wl, coeffs[1], coeffs[2] )

#     residualVals = np.array([0] * len(wl))
#     for index, y in enumerate(globSpecSetY):
#         lx0=[coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4]]
#         x = scipy.optimize.least_squares(residualsLocal3_ecs_f_unk, lx0, args=(y, wl, ecsY, cytFY, unkSpecY))
#         r = np.array(y - model3_ecs_f_unk(wl, x.x, ecsY, cytFY, unkSpecY))
#         residualVals = residualVals + r #np.abs(r)/len(globSpecSet)
#     return residualVals

# def residualsLocal3(coeffs, y, t):
#     global unkSpec
#     global globSpecSet
# #     unkSpec = coeffs[0] * gaussianSpec(wl, coeffs[1], coeffs[2] )
#     r = y - model3(t, coeffs)
#     return r


def residualsLocal3_ecs_f_unk(coeffs, y, wl, ecsY, cytFY, unkSpecY):
    #global unkSpec
    #global globSpecSet
#     unkSpec = coeffs[0] * gaussianSpec(wl, coeffs[1], coeffs[2] )
    r = y - model3_ecs_f_unk(wl, coeffs, ecsY, cytFY, unkSpecY)
    return r


def model3_ecs_f_unk(wl, coeffs, ecsY, cytFY, unkSpec):
    fitSpec = coeffs[0] *  unkSpec #gaussianSpec(wl, coeffs[1], coeffs[2] )
    fitSpec = fitSpec + coeffs[3] * ecsY 
    fitSpec = fitSpec + coeffs[4] * cytFY 
    return fitSpec # * np.exp( - ((t-coeffs[2])/coeffs[3])**2)



# def residualsGlobal_ecs_f_qE_unk(coeffs, yY, wl, globSpecSetY, ecsY, cytFY, qEY):
def residualsGlobal_ecs_f_qE_unk(coeffs, wl, globSpecSetY, ecsY, cytFY, qEY):
    #print(wl)
    unkSpecY = gaussianSpec(wl, coeffs[1], coeffs[2] )
    #print(len(unkSpecY[0]))
    residualVals = np.array([0] * len(wl))
    for y in globSpecSetY:
        lx0=[coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5]]
        #print(lx0)
        #print(y)
        x = scipy.optimize.least_squares(residualsLocal_ecs_f_qE_unk, lx0, args=(y, ecsY, cytFY, qEY, unkSpecY))
        
        r = y - np.array(model_ecs_f_qE_unk(x.x, ecsY, cytFY, qEY, unkSpecY))
        # r = y - model_ecs_f_qE_unk(coeffs, ecsY, cytFY, qEY, unkSpecY)
        residualVals = residualVals + r #np.abs(r)/len(globSpecSet)
    return residualVals

def residualsLocal_ecs_f_qE_unk(coeffs, y, ecsY, cytFY, qEY, unkSpecY):
    
    r = y - model_ecs_f_qE_unk(coeffs, ecsY, cytFY, qEY, unkSpecY)
    return r


def model_ecs_f_qE_unk(coeffs, ecsY, cytFY, qEY, unkSpecY):
    
    fitSpec = coeffs[0] *  unkSpecY #gaussianSpec(wl, coeffs[1], coeffs[2] )
    fitSpec = fitSpec + coeffs[3] * ecsY 
    fitSpec = fitSpec + coeffs[4] * cytFY
    fitSpec = fitSpec + coeffs[5] * qEY
    return fitSpec #

def generateShiftSpec(wlVals, amp1,cwl1,hwidth1,amp2,cwl2,width2): #function to generate a ECS type spectrum
    """
    generateSOGSpec generates a spectrum from the sum of two Gaussian curves. The curves can have positive
    or negative amplitudes so that a "shift" spectrum can be generated.
    wlVals, an array of wavelength values for which the spectrum will be generated
    amp1, amplitude of the first component
    cwl1, center wavelength of the first component
    hwidth1, the half width at half height for component 1
    amp2,cwl2,hwidth2, same as above for component 2

    """
    hwidth1 = 10
    hwidth2 = 20
    componentSpec1 = amp1 * gaussianSpec(wlVals, cwl1, hwidth1 )
    componentSpec2 = amp2 * gaussianSpec(wlVals, cwl2, hwidth2 )
    return componentSpec2-componentSpec1

def shiftSpecResiduals(coefs, fitToThisSpec, wlVals):
    """
    shiftSpecResiduals returns the residuals between a given spectrum called fitToThisSpec, 
    wlVals, array of wavelengths for the given spectrum
    coefs = inputs for generateShiftSpec, i.e. [amp1,cwl1,hwidth1,amp2,cwl2,hwidth2]

    """
    amp1,cwl1,hwidth1,amp2,cwl2,width2 = coefs
    shiftSpec=generateShiftSpec(wlVals, amp1,cwl1,hwidth1,amp2,cwl2,width2)
    return fitToThisSpec - shiftSpec

#

# test if a simple sum-of-Gaussians is a reasonable mimic of the red-ox cyt f dif spectrum 

#
def generateSOGSpec(wlVals, amp1,cwl1,hwidth1,amp2,cwl2,hwidth2,amp3,cwl3,hwidth3): #function to generate a ECS type spectrum
    """
    generateSOGSpec generates a spectrum from the sum of three Gaussian curves. The curves can have positive
    or negative amplitudes so that a "shift" spectrum can be generated.
    wlVals, an array of wavelength values for which the spectrum will be generated
    amp1, amplitude of the first component
    cwl1, center wavelength of the first component
    hwidth1, the half width at half height for component 1
    amp2,cwl2,hwidth2, same as above for component 2
    amp3,cwl3,hwidth3, same as above for component 3

    """
    componentSpec1 = amp1 * gaussianSpec(wlVals, cwl1, hwidth1 )
    componentSpec2 = amp2 * gaussianSpec(wlVals, cwl2, hwidth2 )
    componentSpec3 = amp3 * gaussianSpec(wlVals, cwl3, hwidth3 )
    return componentSpec1 + componentSpec2 + componentSpec3

def shiftSOGResiduals(coefs, fitToThisSpec, wlVals):

    """
    shiftSOGResiduals returns the residuals between a givesn spectrum called fitToThisSpec, 
    wlVals, array of wavelengths for the given spectrum
    coefs = inputs for generateSOGSpec, i.e. [amp1,cwl1,hwidth1,amp2,cwl2,hwidth2,amp3,cwl3,hwidth3]

    """

    amp1,cwl1,hwidth1,amp2,cwl2,hwidth2,amp3,cwl3,hwidth3 = coefs
    shiftSpec=generateSOGSpec(wlVals, amp1,cwl1,hwidth1,amp2,cwl2,hwidth2,amp3,cwl3,hwidth3)
    return fitToThisSpec - shiftSpec


# --------------------
def residualsGlobal_ecs_f_qE_Zx_unk(coeffs, wl, globSpecSetY, ecsY, cytFY, qEY, ZxY):
    #print(wl)
    unkSpecY = gaussianSpec(wl, coeffs[1], coeffs[2] )
    #print(len(unkSpecY[0]))
    residualVals = np.array([0] * len(wl))
    for y in globSpecSetY:
        lx0=[coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5], coeffs[6]]
        #print(lx0)
        #print(y)
        x = scipy.optimize.least_squares(residualsLocal_ecs_f_qE_Zx_unk, lx0, args=(y, ecsY, cytFY, qEY, unkSpecY, ZxY))
        
        r = y - np.array(model_ecs_f_qE_Zx_unk(x.x, ecsY, cytFY, qEY, unkSpecY, ZxY))
        # r = y - model_ecs_f_qE_unk(coeffs, ecsY, cytFY, qEY, unkSpecY)
        residualVals = residualVals + r #np.abs(r)/len(globSpecSet)
    return residualVals

def residualsLocal_ecs_f_qE_Zx_unk(coeffs, y, ecsY, cytFY, qEY, unkSpecY, ZxY):
    
    r = y - model_ecs_f_qE_Zx_unk(coeffs, ecsY, cytFY, qEY, unkSpecY, ZxY)
    return r


def model_ecs_f_qE_Zx_unk(coeffs, ecsY, cytFY, qEY, unkSpecY, ZxY):
    
    fitSpec = coeffs[0] *  unkSpecY #gaussianSpec(wl, coeffs[1], coeffs[2] )
    fitSpec = fitSpec + coeffs[3] * ecsY 
    fitSpec = fitSpec + coeffs[4] * cytFY
    fitSpec = fitSpec + coeffs[5] * qEY
    fitSpec = fitSpec + coeffs[6] * ZxY
    return fitSpec #


# -------------------- using a 3-Gaussian unknown spectrum -----------------------

def residualsGlobal_ecs_f_qE_Zx_unk3(coeffs, wl, globSpecSetY, ecsY, cytFY, qEY, ZxY):
    #print(wl)
    #unkSpecY = gaussianSpec(wl, coeffs[1], coeffs[2] )
    #coeffs = [unkAmpl,cwl1,hwidth1,amp2,cwl2,hwidth2,amp3,cwl3,hwidth3,ecsAmp, cytFAmpl, qEAmpl, ZxAmpl]
    unkSpecY = generateSOGSpec(wl, 1.0,coeffs[1],coeffs[2],
                                    coeffs[3],coeffs[4],coeffs[5],
                                    coeffs[6], coeffs[7],coeffs[8])
    #Note, I pass 1.0 for the first coeff in the above because I will also adjust the 
    # amplitude of the unkSpecY in the global fit, below
    unkSpecY=unkSpecY-unkSpecY[-1]

    #print(len(unkSpecY[0]))
    residualVals = np.array([0] * len(wl))
    lx0=[coeffs[0], # unk component amplitude
    coeffs[9], #ecs amplitude
    coeffs[10], #cyt f amplitude
    coeffs[11], #qE amplitude
    coeffs[12]] #Zx amplitude

    for y in globSpecSetY:
        # lx0=[coeffs[0],coeffs[1],coeffs[2], #unk component 1
        #     coeffs[3],coeffs[4],coeffs[5], #unk component 2
        #     coeffs[6],coeffs[7],coeffs[8], #unk component 3
        #     coeffs[9], #ecs amplitude
        #     coeffs[10], #cyt f amplitude
        #     coeffs[11], #qE amplitude
        #     coeffs[12]] #Zx amplitude
        #print(lx0)
        #print(y)
        x = scipy.optimize.least_squares(residualsLocal_ecs_f_qE_Zx_unk3, lx0, args=(y, ecsY, cytFY, qEY, ZxY, unkSpecY))
        lx0=x.x  

        # boundsLow=[]
        # boundsHigh=[]
        # for i in range(len(lx0)):
        #     boundsLow.append(lxo[i]-.01)
        #     boundsHigh.append(lxo[i]+.01)

        r = y - np.array(model_ecs_f_qE_Zx_unk3(x.x, ecsY, cytFY, qEY, ZxY, unkSpecY))
        # r = y - model_ecs_f_qE_unk(coeffs, ecsY, cytFY, qEY, unkSpecY)
        residualVals = residualVals + r #np.abs(r)/len(globSpecSet)
    return residualVals

def residualsLocal_ecs_f_qE_Zx_unk3(coeffs, y, ecsY, cytFY, qEY, ZxY, unkSpecY):
    unkSpecY=unkSpecY-unkSpecY[-1]
    r = y - model_ecs_f_qE_Zx_unk3(coeffs, ecsY, cytFY, qEY, ZxY, unkSpecY)
    return r

def model_ecs_f_qE_Zx_unk3(coeffs, ecsY, cytFY, qEY, ZxY, unkSpecY):
    
    fitSpec = coeffs[0] *  unkSpecY #gaussianSpec(wl, coeffs[1], coeffs[2] )
    fitSpec = fitSpec + coeffs[1] * ecsY 
    fitSpec = fitSpec + coeffs[2] * cytFY
    fitSpec = fitSpec + coeffs[3] * qEY
    fitSpec = fitSpec + coeffs[4] * ZxY
    return fitSpec #


# --------------------
def residualsGlobal_ecs_f_qE_Zx_unk(coeffs, wl, globSpecSetY, ecsY, cytFY, qEY, ZxY):
    #print(wl)
    unkSpecY = gaussianSpec(wl, coeffs[1], coeffs[2] )
    #print(len(unkSpecY[0]))
    residualVals = np.array([0] * len(wl))
    for y in globSpecSetY:
        lx0=[coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5], coeffs[6]]
        #print(lx0)
        #print(y)
        x = scipy.optimize.least_squares(residualsLocal_ecs_f_qE_Zx_unk, lx0, args=(y, ecsY, cytFY, qEY, unkSpecY, ZxY))
        
        r = y - np.array(model_ecs_f_qE_Zx_unk(x.x, ecsY, cytFY, qEY, unkSpecY, ZxY))
        # r = y - model_ecs_f_qE_unk(coeffs, ecsY, cytFY, qEY, unkSpecY)
        residualVals = residualVals + r #np.abs(r)/len(globSpecSet)
    return residualVals

def residualsLocal_ecs_f_qE_Zx_unk(coeffs, y, ecsY, cytFY, qEY, unkSpecY, ZxY):
    
    r = y - model_ecs_f_qE_Zx_unk(coeffs, ecsY, cytFY, qEY, unkSpecY, ZxY)
    return r


def model_ecs_f_qE_Zx_unk(coeffs, ecsY, cytFY, qEY, unkSpecY, ZxY):
    
    fitSpec = coeffs[0] *  unkSpecY #gaussianSpec(wl, coeffs[1], coeffs[2] )
    fitSpec = fitSpec + coeffs[3] * ecsY 
    fitSpec = fitSpec + coeffs[4] * cytFY
    fitSpec = fitSpec + coeffs[5] * qEY
    fitSpec = fitSpec + coeffs[6] * ZxY
    return fitSpec #


# -------------------- using cyt b and 1 3-Gaussian unknown spectrum -----------------------

def residualsGlobal_ecs_f_qE_Zx_unk_b(coeffs, wl, globSpecSetY, ecsY, cytFY, qEY, ZxY, cytBY):
    #print(wl)
    #unkSpecY = gaussianSpec(wl, coeffs[1], coeffs[2] )
    #coeffs = [unkAmpl,cwl1,hwidth1,amp2,cwl2,hwidth2,amp3,cwl3,hwidth3,ecsAmp, cytFAmpl, qEAmpl, ZxAmpl]
    unkSpecY = generateSOGSpec(wl, 1.0,coeffs[1],coeffs[2],
                                    coeffs[3],coeffs[4],coeffs[5],
                                    coeffs[6], coeffs[7],coeffs[8])
    #Note, I pass 1.0 for the first coeff in the above because I will also adjust the 
    # amplitude of the unkSpecY in the global fit, below
    unkSpecY=unkSpecY-unkSpecY[-1]

    #print(len(unkSpecY[0]))
    residualVals = np.array([0] * len(wl))
    lx0=[coeffs[0], # unk component amplitude
    coeffs[9], #ecs amplitude
    coeffs[10], #cyt f amplitude
    coeffs[11], #qE amplitude
    coeffs[12], #Zx amplitude
    coeffs[13]] #cyt b  amplitude

    for y in globSpecSetY:
        # lx0=[coeffs[0],coeffs[1],coeffs[2], #unk component 1
        #     coeffs[3],coeffs[4],coeffs[5], #unk component 2
        #     coeffs[6],coeffs[7],coeffs[8], #unk component 3
        #     coeffs[9], #ecs amplitude
        #     coeffs[10], #cyt f amplitude
        #     coeffs[11], #qE amplitude
        #     coeffs[12]] #Zx amplitude
        #print(lx0)
        #print(y)
        x = scipy.optimize.least_squares(residualsLocal_ecs_f_qE_Zx_unk_b, lx0, args=(y, ecsY, cytFY, qEY, ZxY, unkSpecY, cytBY))
        lx0=x.x  

        # boundsLow=[]
        # boundsHigh=[]
        # for i in range(len(lx0)):
        #     boundsLow.append(lxo[i]-.01)
        #     boundsHigh.append(lxo[i]+.01)

        r = y - np.array(model_ecs_f_qE_Zx_unk_b(x.x, ecsY, cytFY, qEY, ZxY, unkSpecY, cytBY))
        # r = y - model_ecs_f_qE_unk(coeffs, ecsY, cytFY, qEY, unkSpecY)
        residualVals = residualVals + r #np.abs(r)/len(globSpecSet)
    return residualVals

def residualsLocal_ecs_f_qE_Zx_unk_b(coeffs, y, ecsY, cytFY, qEY, ZxY, unkSpecY, cytBY):
    unkSpecY=unkSpecY-unkSpecY[-1]
    r = y - model_ecs_f_qE_Zx_unk_b(coeffs, ecsY, cytFY, qEY, ZxY, unkSpecY, cytBY)
    return r

def model_ecs_f_qE_Zx_unk_b(coeffs, ecsY, cytFY, qEY, ZxY, unkSpecY, cytBY):
    
    fitSpec = coeffs[0] *  unkSpecY #gaussianSpec(wl, coeffs[1], coeffs[2] )
    fitSpec = fitSpec + coeffs[1] * ecsY 
    fitSpec = fitSpec + coeffs[2] * cytFY
    fitSpec = fitSpec + coeffs[3] * qEY
    fitSpec = fitSpec + coeffs[4] * ZxY
    fitSpec = fitSpec + coeffs[5] * cytBY
    return fitSpec #



def spectral_color(wl): #/ RGB <0,1> <- lambda l <400,700> [nm]
    """
    Returns the apparent RGB color as a tuple of a monochromatic light source at wavelenth = wl.

    """
    
    l=wl
    t=0.0
    r=0.0
    g=0.0
    b=0.0
    if ((l>=400.0) and (l<410.0)):
        t=(l-400.0)/(410.0-400.0)
        r=    +(0.33*t)-(0.20*t*t)
    elif ((l>=410.0) and (l<475.0)):
        t=(l-410.0)/(475.0-410.0)
        r=0.14-(0.13*t*t)
    elif ((l>=545.0) and (l<595.0)):
        t=(l-545.0)/(595.0-545.0)
        r=+(1.98*t)-(t*t)
        
    elif ((l>=595.0)  and (l<650.0)):
        t=(l-595.0)/(650.0-595.0)
        r=0.98+(0.06*t)-(0.40*t*t)
        
    elif ((l>=650.0) and(l<700.0)): 
        t=(l-650.0)/(700.0-650.0)
        r=0.65-(0.84*t)+(0.20*t*t)
    if ((l>=415.0) and (l<475.0)): 
        t=(l-415.0)/(475.0-415.0)
        g=+(0.80*t*t)
    elif ((l>=475.0)  and (l<590.0)):
        t=(l-475.0)/(590.0-475.0)
        g=0.8 +(0.76*t)-(0.80*t*t)
    elif ((l>=585.0) and (l<639.0)): 
        t=(l-585.0)/(639.0-585.0)
        g=0.84-(0.84*t)
    if ((l>=400.0) and (l<475.0)): 
        t=(l-400.0)/(475.0-400.0)
        b=+(2.20*t)-(1.50*t*t)
    elif ((l>=475.0) and (l<560.0)):
        t=(l-475.0)/(560.0-475.0)
        b=0.7 -(t)+(0.30*t*t)
    return (r,g,b)




# Function to combine multiple project JSON files into one big one.

# {
#     "tags": [
#         "hide-input",
#     ]
# }

def generateCombinedJSON (baseFileName, destinationBaseFileName="", append_bases = [''], verbose = False, overwrite = False):
    if (destinationBaseFileName==""):
        destinationBaseFileName=baseFileName

    # oldJSON = glob.glob(baseFileName + '.datcombined' '*.json')

    
    my_file = Path(destinationBaseFileName+'_combined.json')
    if ((my_file.is_file()) and (overwrite == False)):
        # print ('combined exists')
        if verbose:
            print ('combined json file already exists:', destinationBaseFileName+'_combined.json')
    else:
        if verbose:
            print( 'generating new' )
        for append_base in append_bases:
            if verbose:
                print(append_base)
            addendum = baseFileName + '.dat' + append_base
            filesToCombine = glob.glob(addendum  + '*.json')
            if (len(filesToCombine)<1):
                if verbose:
                    print("No JSON files with base file name = " + baseFileName + "found ")
                # return 
            else:
                f= open(destinationBaseFileName+'_combined.json',"w+")
                f.write('[')
                f.close()
                f= open(destinationBaseFileName+'_combined.json',"a+")
                if verbose:
                    print("combining files into: " + destinationBaseFileName)
                for index, file in enumerate(filesToCombine):
                    #print(file)
                    succeeded = False    
                    try:
                        with open(file) as json_file:
                            jd = json.load(json_file)
                        succeeded = True
                    except:
                        print("failed:" + file)
                        succeeded = False
                    if(succeeded == True):
                        c = open(file,"r+")
                        f.write(c.read())
                        if (index<len(filesToCombine)-1):
                            f.write(', ')
                f.write(']')        
                f.close()   
            if verbose:
                print('Combined JSON = ' + destinationBaseFileName+'_combined.json')

# Test if there is already a notes CSV, if so, open it and see what's there. Use the columns and data
# in any existing file. 

# Note: This is important if the user is adding additional files to the project and does not 
# want to over-write the previously recorded stuff. 

def generate_project_csv (project_name, user_data_folder = '', master_folder = '', notes_file = '.csv', over_write = False):

    notes_file_name = user_data_folder + master_folder + project_name + notes_file
    print('notes_file_name = ' + notes_file_name)
    new_dataframe = False

    if ((os.path.exists(notes_file_name)) and (over_write == False)):
        print("Previous notes file found and will be modified. \n Re-run with over_write = True to overwrite all previous edits.")

        project_dataframe = pd.read_csv(notes_file_name)
        new_dataframe = True
    else:
        print("Generating new dataframe")
        project_dataframe = pd.DataFrame(columns =  ['name', 'subfolder', 'notes', 'checksum', 'added notes'])
        project_dataframe['protocol_labels'] = 0
        project_dataframe['protocol_labels'] = project_dataframe['protocol_labels'].astype('object')
    return project_dataframe
    # Find all the files in the indicated master folder and its sub-folders.

    all_project_files = glob.glob(user_data_folder + master_folder + '/**', recursive = True)

# Funciton to get a list of all subdirectories in and below the master-folder


 
def listdirs(rootdir):
    sub_folders = []
    for path in Path(rootdir).iterdir():
        if path.is_dir():
            listdirs(path)
            sub_folders.append(str(path))
    return sub_folders

# The following is used to determine if a file has been changed. It calculated the 
# check sum for all the data in the file. 

def generate_file_hash( path) :
    with open(path, "rb") as f:
        file_hash = hashlib.md5()
        chunk = f.read(8192)
        while chunk:
            file_hash.update(chunk)
            chunk = f.read(8192)

    return file_hash.hexdigest()

def find_all_project_subfolders():
    project_sub_folders_with_local = listdirs( '.') #+ user_data_folder + master_folder)

    project_sub_folders = []

    print("The following sub folders are included in the project :")
    for  project_sub_folder in  project_sub_folders_with_local:
        project_sub_folders.append(project_sub_folder) #.split(master_folder)[-1])

    return(project_sub_folders)

def combined_jsons_and_add_to_csv(project_dataframe, project_sub_folders, user_data_folder = '', 
master_folder = '', verbose = False, over_write = False):
    

    project_checksums = list(project_dataframe['checksum'])

    unique_project_jsons = []  # any combined_json with a unique sub_folder and file name

    index = 0

    all_notes = {}
    for  project_sub_folder_short in  project_sub_folders:
        project_sub_folder = user_data_folder + master_folder + project_sub_folder_short
        print("-------------------------------------------")
        print("scanning folder for JSON files: " + str(project_sub_folder))
        print("-------------------------------------------")
        sub_folder_project_files = glob.glob(project_sub_folder + '/**', recursive = True)
        all_jsons = {}
        all_combinated_jsons = []
        for fn in sub_folder_project_files:
            full_file_name = os.path.basename(fn)
            subfolder = os.path.dirname(fn)

            if full_file_name.endswith('_combined.json'):
                if (verbose == True):
                    print("found _combined.json file for: " + str(full_file_name))
                all_combinated_jsons.append(full_file_name.split('_combined.json')[0])
            elif full_file_name.endswith('.json'):
                base_file_name = full_file_name.split('.dat')[0]
                all_jsons[base_file_name] = subfolder + '/' + base_file_name
        unique_jsons = set(all_jsons.keys())
        if (verbose == True):
            print(unique_jsons, all_combinated_jsons)
        for unique_json in unique_jsons:

            if (unique_json in all_combinated_jsons):
                if (verbose == True):
                    print("already combined : " + unique_json)
            else:
                if (verbose == True):
                    print("**** generating combined JSON for: " + unique_json)
                generateCombinedJSON(all_jsons[unique_json], all_jsons[unique_json])

            if (verbose == True):
                print(all_jsons[unique_json])
            write_time = os.path.getctime(all_jsons[unique_json] + '_combined.json')

            checksum = generate_file_hash(all_jsons[unique_json] + '_combined.json')
            if ((checksum in project_checksums) and (over_write == False)):
                print("Checksum matched for: " + unique_json)
            else:
                print("new file")
                unique_project_jsons.append(all_jsons[unique_json]) #.split(master_folder)[-1])

                project_dataframe.loc[index, 'name'] = all_jsons[unique_json].split('/')[-1]
    #             print("all_jsons[unique_json] = " + all_jsons[unique_json])

                test = all_jsons[unique_json] # the full path and file
#                 test = test[0:test.rfind('/')] # find the last path indicator
#                 test = test.split(master_folder)[-1] # strip off the file name by plitting at the master_folder name

                project_dataframe.loc[index, 'subfolder'] = test #subfolder.split(master_folder)[-1]

                project_dataframe.loc[index, 'checksum'] = checksum

                notes_file_name = all_jsons[unique_json] + '.dat_notes.txt'
                if (verbose == True):
                    print (notes_file_name)
                if (os.path.exists(notes_file_name)):

                    with open(notes_file_name) as f:
                        notes_text = f.read()

                    all_notes[index] = notes_text

                    protocol = notes_text.split("script = ")[-1]
                    protocol = protocol.split('\n****')[0]
                    protocol_name = protocol.split('\n')[0]
                    print ("protocol_name = " + protocol_name)
                    project_dataframe.loc[index, 'notes_beginning'] = protocol_name

                    notes_text_before_protocol =  notes_text.split('*********')[0]
                    notes_text_before_protocol = notes_text_before_protocol.replace('\n', ' ')
                    notes_text_before_protocol = notes_text_before_protocol.replace('"', '')

                    notes_text_after_protocol =  notes_text.split('*********')[0]

                    project_dataframe.loc[index, 'notes'] = notes_text_before_protocol



                else:
                    project_dataframe.loc[index, 'notes'] = "none"
                    notes_text_after_protocol =  "none"
                    all_notes[index] = "none"
                    project_dataframe.loc[index, 'protocol_name'] = "unknown"

            # Load the data file to extract certain data 
                file_name = all_jsons[unique_json] + '_combined.json'
                df = pd.read_json(file_name)
                df = df.sort_values('start_time') # Get them in the right order!
                #project_dataframe.loc[index, 'set_temperature'] = np.mean(df['set_temperature']) 
                #project_dataframe.loc[index, 'ambient_temperatures'] = str(list(df['ambient_temperature'])) 
                #project_dataframe.loc[index, 'av_ambient_temperature'] = np.round(np.mean(df['ambient_temperature']),2)
                project_dataframe.loc[index, 'actinic_intensity'] = str(list(df['actinic_intensity'])) 
                project_dataframe.loc[index, 'start_times'] = str(list(df['start_time'])) 
                all_unique_protocol_labels = list(set(df['protocol_label']))
                print(all_unique_protocol_labels)
                for j, protocol_label in enumerate(all_unique_protocol_labels):
                    print(j, protocol_label)
                project_dataframe['protocol_labels'][index] = all_unique_protocol_labels
            index += 1
            
def export_project_dataframe (project_name, project_dataframe, user_data_folder = '', master_folder = '', notes_file = '.csv'):
    project_dataframe.to_csv(user_data_folder + master_folder + project_name + notes_file)

def open_and_display_generated_project_dataframe(project_name, user_data_folder = '', master_folder = '', notes_file = '.csv'):
    # Open the modified version to see what's in it. 
    return (pd.read_csv(user_data_folder + master_folder + project_name +  notes_file))
    

def set_up_project_dataframe (project_name, over_write = False):
    """
    Perform this operation in a folder with all the 
    data in sub-folders.

    """

    project_dataframe = generate_project_csv(project_name, over_write = over_write)
    project_sub_folders = find_all_project_subfolders()
    combined_jsons_and_add_to_csv(project_dataframe, project_sub_folders, over_write = over_write)
    export_project_dataframe(project_name, project_dataframe)
    # Ipy.open_and_display_generated_project_dataframe(project_name)
    return project_dataframe



# The saved time array can be truncated. This function recalculates the 
# values based on those in the protocol.
# Note that the starting time for each trace is set to zero

def recalculate_time_scale_from_protocol (project_DF, allExperiments):
    for experiment_index in allExperiments: 
        dataSet = project_DF['data'][experiment_index]
        dataSet['trace_times'] = 0
        dataSet['trace_times'] = dataSet['trace_times'].astype('object')
        for index in dataSet.index:
            total_trace_time = 0
            trace_times = []
            for i, pulses in enumerate(dataSet.loc[index, 'number_pulses']):
                pulse_time = dataSet.loc[index, 'measuring_interval'][i]
                for pulse in range(0, pulses):
                    total_trace_time += pulse_time
                    trace_times.append(total_trace_time)
            dataSet['trace_times'][index] = trace_times
    project_DF['data'][experiment_index] = dataSet
    return project_DF

def recalculate_time_information (combined_df):
    """
    Recalculates the times for each trace based on the protocol information. 
    input the dataframe as combined_df.

    The saved time array can be truncated. This function recalculates the 
    values based on those in the protocol.
    Note that the starting time for each trace is set to zero


    """
    for index in combined_df.index:
        current_time = 0.0

        timing_information = {} # set up a dict to contain the timing information for each measuring_light_name in measuring_light_names
        for measuring_light_name in combined_df['measuring_light_names'][index]:
            timing_information[measuring_light_name] = []

        # the 'measuring_interval' columns has timing information for each LED for each cycle.
        # The following reshales measuring_interval to make an array of arrays, one for each cycle

        timings_for_each_segment = np.reshape(np.array(combined_df['measuring_interval'][index]), [len(combined_df['number_pulses'][index]),-1])

        for pulse_segment_index in range(len(combined_df['number_pulses'][index])):  # go through each pulse segment (loops)
            for pulses in range(0, combined_df['number_pulses'][index][pulse_segment_index]): # go through each set of pulses
                for measuring_light_index, measuring_light in enumerate(combined_df['measuring_light_names'][index]):  # this is done for each meauring light
                    current_time = current_time + timings_for_each_segment[pulse_segment_index][measuring_light_index]
                    timing_information[measuring_light].append(current_time)

        for key in timing_information.keys():
            combined_df[key + '_time'][index] = np.array(timing_information[key])

    
# project_DF = recalculate_time_scale_from_protocol (project_DF)

def exp_decay(x, amplitude, tau):
    return amplitude * np.exp(-x/tau)





def drift_spec_Thekla(wl):
#     driftTrace=np.array([ 0.        , -0.09785897, -0.086966  , -0.03990309,  0.00588594,0.02904874, -0.07482107, -0.21361617, -0.32485364, -0.45968053,
#                          -0.46807275, -0.29403091, -0.02755454,  0.33847267,  0.67510016, 0.84530737,  0.97552712,  0.99503461,  1.        ,  0.84823227,
#                          0.85294959,  0.92271972,  0.91801522,  0.91141255,  0.91054229, 0.86078758,  0.83688201,  0.73295667,  0.6493151 ,  0.61022771,
#                          0.57411674,  0.56226763,  0.56638073,  0.56789946,  0.53810999,0.42679676,  0.28241103,  0.15982337,  0.06346248, -0.00105235,
#                          -0.06746963, -0.04847341, -0.14903963])
#     driftTraceWL=np.array([464, 467, 470, 473, 476, 479, 482, 485, 488, 491, 494, 497, 500, 503, 506, 509, 512, 515, 518, 521, 524, 527, 530, 533, 536, 539,
#                        542, 545, 548, 551, 554, 557, 560, 563, 566, 569, 572, 575, 578, 581, 584, 587, 590])
    
#     driftTrace=np.array([-1., -0.96627725, -0.90695315, -0.88016554, -0.83421386, 
#                          -0.83003877, -0.85046984, -0.863734  , -0.88017335, -0.86681812,
#                          -0.83408534, -0.7532907 , -0.64931702, -0.54754943, -0.43186477,
#                          -0.33026628, -0.27214651, -0.25390547, -0.24432818, -0.24877121,
#                          -0.255649  , -0.27150303, -0.29383148, -0.31180011, -0.32601945,
#                          -0.33942055, -0.35029493, -0.36257122, -0.37053509, -0.37777462,
#                          -0.37886552, -0.37664543, -0.37134035, -0.36172252, -0.35611495,
#                          -0.3555346 , -0.35997958, -0.36298778, -0.36666854, -0.37249348,
#                          -0.37897204, -0.38831918, -0.40072305])
#     driftTraceWL=np.array([464, 467, 470, 473, 476, 479, 482, 485, 488, 491, 494, 497, 500,
#                             503, 506, 509, 512, 515, 518, 521, 524, 527, 530, 533, 536, 539,
#                             542, 545, 548, 551, 554, 557, 560, 563, 566, 569, 572, 575, 578,
#                             581, 584, 587, 590])
    #driftTrace=driftTrace-np.linspace(driftTrace[0],driftTrace[-1],len(driftTrace))
    driftTrace=np.array([-0.97930457, -0.84510825, -0.76824325, -0.68557158, -0.65721521,
                         -0.71411414, -0.80057651, -0.9111386 , -0.99341127, -1.        ,
                         -0.90233654, -0.7229873 , -0.50413059, -0.30143799, -0.150493  ,
                         -0.05590374,  0.00378046,  0.03427544,  0.04805041,  0.05288039,
                         0.05259516,  0.033782  ,  0.        ])
    
    driftTraceWL=np.array([460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
    # global driftTraces
#     global driftTraceWLs
#     driftTrace=driftTraces
#     driftTracesWL=driftTracesWLs
    # driftTrace=driftTrace/np.max(np.abs(driftTrace))
    # print(driftTrace)
    spec = []
    for  i in range(len(wl)):
        print(i)
        spec.append(np.interp(wl[i],driftTraceWL,driftTrace))
    spec=np.array(spec)
    return spec

def ecs_spec_Thekla(wl):
    
#     pulseECSWL=[590, 585,580,575,570,565,560,555,550,545,540,535,530,525,520,515,510,505,500,495,490,485,480,475,470,465]
#     pulseECS = np.array([ 0.        ,  0.03344159,  0.10252101,  0.17498088,  0.28436124,
#              0.37270553,  0.39466636,  0.41549049,  0.52507503,  0.64408459,
#              0.73683989,  0.82697959,  0.93609262,  1.        ,  0.97564624,
#              0.83643035,  0.59682901,  0.42985326,  0.25337081,  0.10998487,
#             -0.02828996, -0.07548181, -0.14295679, -0.13838168,  0.00708065,
#              0.        ])
    
#     ecsWL=np.flip(np.array([590, 585,580,575,570,565,560,555,550,545,540,535,530,525,520,515,510,505,500,495,490,485,480,475,470,465]))
#     ecsWL=np.flip(np.array([580,575,570,565,560,555,550,545,540,535,530,525,520,515,510,505,500,495,490,485,480,475,470,465,460,455]))
#     ecsSpec=np.flip(np.array([ 0.        ,  0.00026393,  0.0005899 ,  0.0011091 ,  0.00146355,
#          0.00225297,  0.00202352,  0.00205593,  0.00257003,  0.0032522 ,
#          0.00366017,  0.00408537,  0.00447852,  0.0049811 ,  0.00479048,
#          0.00389282,  0.00247211,  0.00121837,  0.00052815,  0.00049321,
#         -0.00072343, -0.00070766, -0.00167152, -0.00073296, -0.00063791,
#          0.        ]))
#     ecsSpec=ecsSpec/np.max(ecsSpec)
# #     pulseECSWL=np.flip(pulseECSWL)
# #     pulseECS=np.flip(pulseECS)
#     ecsWL=ecsWL+3
#     ecsSpec=np.array([-0.71277281, -0.43899874, -0.10387849,  0.09269294,  0.13036653,
#                       0.22194588,  0.46461142,  0.80149811,  0.98423798,  1.        ,
#                       0.85163492,  0.65056971,  0.45509381,  0.29812798,  0.1620286 ,
#                       0.05349879,  0.01913389,  0.0399582 ,  0.03149595,  0.        ])
#     ecsSpec=np.array([ 0.00000000e+00,  8.78528548e-02,  2.26741113e-01,  2.40000754e-01,
#                       -3.31505254e-01, -1.16754523e-01, -3.62251730e-01, -3.51848102e-01,
#                       -1.60837185e-01,  2.68709078e-02,  1.76662957e-01,  3.82831285e-01,
#                       6.98142034e-01,  9.27353062e-01,  1.00000000e+00,  8.96634243e-01,
#                       7.39476823e-01,  5.51035261e-01,  3.54006517e-01,  2.05256450e-01,
#                       8.48566119e-02,  3.77208458e-05,  3.90792755e-02,  3.99415059e-02, 5.89143128e-03])
    # ecsSpec=np.array([ 0.00000000e+00,  8.78528548e-02,  2.26741113e-01,  1.40000754e-01,
    #                       -3.01505254e-01, -3.6754523e-01, -3.8e-01, -2.8e-01,
    #                       -1.60837185e-01,  2.68709078e-02,  1.76662957e-01,  3.82831285e-01,
    #                       6.98142034e-01,  9.27353062e-01,  1.00000000e+00,  8.96634243e-01,
    #                       7.39476823e-01,  5.51035261e-01,  3.54006517e-01,  2.05256450e-01,
    #                       9.8e-02,  7.7e-02,  4.3e-02,  3.99415059e-02,5.89143128e-03])

    # ecsWL=np.array([450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
    #print('hello')
    ecsWL=np.array([450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
    # ecsWL = ecsWL - 2
    ecsSpec=np.array([ 0.        ,  0.08785285,  0.22674111,  0.14000075, -0.30150525,
    -0.36754523, -0.38      , -0.28      , -0.16083718,  0.02687091,
    0.17666296,  0.38283129,  0.69814203,  0.92735306,  1.        ,
    0.89663424,  0.73947682,  0.55103526,  0.35400652,  0.20525645,
    0.098     ,  0.077     ,  0.043     ,  0.03994151,  0.00589143])
    # print(len(wl), len(ecsWL), len(ecsSpec))
    spec = []
    for i in range(len(wl)):
        spec.append(np.interp(wl[i],ecsWL,ecsSpec))
    spec=np.array(spec)
    return spec

def Zx_spec_Thekla(wl):
    #     ZxTraceStd5 = np.array([ 0.00000000e+00, -1.79886312e-04, -3.25289036e-05, -1.15619522e-04,
    #             0.00000000e+00, -1.06225262e-04,  7.33846784e-05,  1.38237125e-04,
    #             4.88201525e-04,  4.73997165e-04,  5.98396301e-04,  1.03501985e-03,
    #             1.21480343e-03,  1.34991721e-03,  1.22544035e-03,  1.25055240e-03,
    #             1.38857424e-03,  1.40426691e-03,  1.76330198e-03,  2.08156458e-03,
    #             2.54884787e-03,  3.17621349e-03,  3.54474194e-03,  3.87742627e-03,
    #             4.03303620e-03,  4.09680338e-03,  4.26200253e-03,  4.38726155e-03,
    #             4.79697190e-03,  5.14672030e-03,  5.38124169e-03,  5.46379428e-03,
    #             5.73546494e-03,  6.24380746e-03,  6.97425469e-03,  7.44846748e-03,
    #             8.11431609e-03,  8.66275501e-03,  8.89903164e-03,  9.04655284e-03,
    #             9.69320636e-03,  1.00677371e-02,  9.93035609e-03,  9.36338487e-03,
    #             8.81188930e-03,  7.52567521e-03,  5.94754468e-03,  4.79249617e-03,
    #             3.39647970e-03,  1.77938303e-03,  9.12544997e-04, -8.01840151e-05,
    #            -5.84768300e-04, -1.46712032e-03, -1.47155542e-03, -1.13722256e-03,
    #            -1.16044516e-03, -5.59467011e-04,  8.28921536e-04,  3.41894540e-04,
    #             1.04732335e-03,  3.48071832e-04,  8.47159471e-04,  6.40058052e-04,
    #             0.00000000e+00])
    #     ZxTraceStd5=ZxTraceStd5-np.linspace(ZxTraceStd5[0],ZxTraceStd5[-1],len(ZxTraceStd5))
    #     ZxTraceStd5=ZxTraceStd5/np.max(ZxTraceStd5)
    #     ZxTraceStd5=np.flip(ZxTraceStd5)

    #     wlZxStd5=[590, 588,586,584,582,580,578,576,574,572,570,568,566,564,562,560,558,556,554,552,550,
    #      548,546,544,542,540,538,536,534,532,530,528,526,524,522,520,518,516,514,512,510,508,
    #      506,504,502,500,498,496,494,492,490,488,486,484,482,480,478,476,474,472,470,468,466,464,462]

    #     wlZxStd5=np.flip(wlZxStd5)
    ZxSpec=np.array([ 0.        , -0.03599246,  0.02244886,  0.00645265, -0.03741914,
         0.00254552,  0.21377456,  0.44296526,  0.80471238,  1.        ,
         0.89261038,  0.7168024 ,  0.44407102,  0.23854102,  0.14500318,
         0.05743463,  0.00804195, -0.03206746, -0.05262834, -0.06410261,
        -0.06272308, -0.04705852, -0.01548376, -0.02279776,  0.00120181,
         0.00533344])
    ZxWL=np.array([460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520,
        525, 530, 535, 540, 545, 550, 555, 560, 565, 570, 575, 580, 585])
    # ZxWL=ZxWL+2
    spec = []
    for i in range(len(wl)):
        spec.append(np.interp(wl[i],ZxWL,ZxSpec))
    spec=np.array(spec)
    return spec

def qE_spec_Thekla(wl):
    #qEWL=[460, 475,488,505,520,530, 535, 545, 560, 590]
    #qESpectrum=np.array([0, -.01,-.05,.5,2.9,4.9, 5.2, 4.2, 2.5, 0])
    #qESpectrum=qESpectrum/np.max(qESpectrum)
    # b=453
    # e=580
    # qEWL=np.linspace(b,e,25)
    # qESpectrum = [0,-.1,-.5,-4,0,0,.5,4,3,-2,-4,2,5,12,18,22,22,19.5,15,11,8,5,5,4,0]
    # #qESpectrum = [0,-.1,-.5,-.4,0,0,.5,.4,.3,-.2,-.4,2,5,12,18,22,22,19.5,15,11,8,5,5,4,0]
    # qESpectrum=np.array(qESpectrum)
    # qESpectrum=qESpectrum/np.max(qESpectrum)
    qEWL=np.array([450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570])
    qESpectrum=np.array([-5.78085357e-04, -6.78553127e-05, -5.04221075e-06,  
    6.01235341e-06, 3.10455355e-05,  1.33747482e-04,  5.15200497e-04,  1.77585925e-03,
    5.47754353e-03,  1.51184821e-02,  3.73401108e-02,  8.25254823e-02, 
    1.63209561e-01,  2.88833911e-01,  4.57399644e-01,  6.48169395e-01,
    8.21913809e-01,  9.32629376e-01,  9.46971617e-01,  8.60418801e-01,
    6.99564686e-01,  5.08968521e-01,  3.31359193e-01,  1.93042167e-01,
    1.00635313e-01])
    qEWL=qEWL-2
    spec = []
    for i in range(len(wl)):
        spec.append(np.interp(wl[i],qEWL,qESpectrum))
    spec=np.array(spec)
    return spec

# from scipy.stats import linregress

def subtract_multi_line_from_single_trace ( x, y, selected_ranges):
# all_genotypes
    """
        select out regions of x and y to fit a straight line
        selected_ranges is s list of lists, e.g., [[0,10], [100, 110]],
        which will select the points between 0 and 9 and 10 and 109.
        Will return the trce - baseline
    """

    # make lists to hold the selected data
    selected_x = np.array([])
    selected_y = np.array([])
    
    for selected_range in selected_ranges:
        selected_x = np.concatenate((selected_x, x[selected_range[0]:selected_range[1]]))
        selected_y = np.concatenate((selected_y, y[selected_range[0]:selected_range[1]]))
    
    results = linregress(selected_x, selected_y) 
    fit_line = x * results.slope + results.intercept
    corrected_trace = y - fit_line
    return corrected_trace


# CODE FOR SETTING UP NEW DATAFRAMES

import matplotlib.colors as mcolors

all_color_names = list(mcolors.XKCD_COLORS.keys())

def getDirectoryList(path):
    directoryList = []

    #return nothing if path is a file
    if os.path.isfile(path):
        return []

    #add dir to directorylist if it contains .txt files
    if len([f for f in os.listdir(path) if f.endswith('.json')])>0:
        directoryList.append(path)

    for d in os.listdir(path):
        new_path = os.path.join(path, d)
        if os.path.isdir(new_path):
            directoryList += getDirectoryList(new_path)

    return directoryList

def getAllExperimentBasefiles(path, verbose = False):
    all_non_combined_jsons = []
    for directory in getDirectoryList(path):
        for file in os.listdir(directory):
            if verbose:
                print('*', directory, file)
            if file.endswith(".json"):
                if verbose:
                    print(file)
                jsonfile = os.path.join(directory, file)
                if ('_combined' not in jsonfile):
                    all_non_combined_jsons.append(jsonfile.split('.dat')[0])
    return list(set(all_non_combined_jsons))

def generateBaseFileNamesDF (paths, verbose = True):
    """

    if path is a string (one path) it will do just that. 
    If input is a list of paths, it will sum up all of them. 

    The function tries to find a "genotype name" from the base file name. This will be the first part of the 
    name separated by an underscore, for example:
    
    npq4_test_this_

    will find and use "npq4" as the "genotype". If no underscores are used, genotype will be set to the entire 
    basefilename.


    """
    if isinstance(paths, str):
        paths = [paths]
    print(paths)
    baseFileNames = pd.DataFrame(columns = ['color', 'group', 'intensity factor', 'genotype', 'files_folder'])
    for path in paths:
        allBF = getAllExperimentBasefiles(path, verbose)
        for i, baseFile in enumerate(allBF):
            name = baseFile.split('/')[-1]
            if (verbose == True):
                print(name)
            genotype = name.split('_')[0]
            group = name.split('_')[-1]
            color = all_color_names[i]
            baseFileNames.loc[name] = {'files_folder': baseFile,
                                                        'color': color,
                                                        'group' : group,
                                                        'genotype' : genotype,
                                                        'intensity factor': 1.0}
    return baseFileNames

# baseFileNames = generateBaseFileNamesDF (path)

def recalculate_time_information (combined_df):
    """
    Recalculates the times for each trace based on the protocol information. 
    input the dataframe as combined_df 
    """
    for index in combined_df.index:
        current_time = 0.0

        timing_information = {} # set up a dict to contain the timing information for each measuring_light_name in measuring_light_names
        for measuring_light_name in combined_df['measuring_light_names'][index]:
            timing_information[measuring_light_name] = []

        # the 'measuring_interval' columns has timing information for each LED for each cycle.
        # The following reshales measuring_interval to make an array of arrays, one for each cycle

        timings_for_each_segment = np.reshape(np.array(combined_df['measuring_interval'][index]), [len(combined_df['number_pulses'][index]),-1])

        for pulse_segment_index in range(len(combined_df['number_pulses'][index])):  # go through each pulse segment (loops)
            for pulses in range(0, combined_df['number_pulses'][index][pulse_segment_index]): # go through each set of pulses
                for measuring_light_index, measuring_light in enumerate(combined_df['measuring_light_names'][index]):  # this is done for each meauring light
                    current_time = current_time + timings_for_each_segment[pulse_segment_index][measuring_light_index]
                    timing_information[measuring_light].append(current_time)

        for key in timing_information.keys():
            combined_df[key + '_time'][index] = np.array(timing_information[key])
            

def generate_combo_dataframe(baseFileNames, verbose = True, recalculate_times = True):
    for i, experiment_name in enumerate(baseFileNames.index):

        if (verbose == True):
            print(experiment_name)


        # Enter a path and "base" file name to save the combined JSON. The file name for the combined JSON
        # will be: baseFileName_combined.json
        baseFileName = baseFileNames.loc[experiment_name, 'files_folder']
        destinationBaseFileName=baseFileName #'/Volumes/D_KRAMER/Thekla/Data FL and LL grown plants/FL'

        # call the function that combines the individual JSON files into one big one and saves to 
        if (verbose == True):
            print(baseFileName)

        generateCombinedJSON(baseFileName, destinationBaseFileName)

        folder_name='' 

        #Define the file name for the main program to open.
    #     print(file_name)
        file_name= destinationBaseFileName + '_combined.json'

        # open the combined JSON into a dictionary of data frames
    #     experiment[experiment_name] = {}
        if verbose:
            print(file_name)

        one_experiment = pd.read_json(file_name)
        one_experiment['experiment_name'] = experiment_name
        one_experiment.index = [str(i) + " " + experiment_name for i in one_experiment.index]
        one_experiment['genotype'] = baseFileNames.loc[experiment_name, 'genotype'] 
        one_experiment['plot_color'] = baseFileNames.loc[experiment_name, 'color'] 

        try:
            notes_file = open(baseFileName + '.dat_notes.txt', 'r')
            one_experiment['notes'] = notes_file.read()
        except:
            one_experiment['notes'] = ""


        if (i == 0): #combining files
            combo = one_experiment.copy()
        else:
            combo = pd.concat([combo, one_experiment])
    if recalculate_times:
        recalculate_time_information (combo)
    return combo

    # There is sometimes an issue with the stored timing information. 
# This will over-write the values with corrected ones.

from itertools import accumulate

def recalculate_measuring_times (selected_experinent, index):

    # total_times = len(selected_experinent['measuring_light_names'][index]) *  [selected_experinent['start_time'][index]]

    cumulative_time = 0 #[selected_experinent['start_time'][index]] #selected_experinent['start_time'][index]  # holds the accumulated time for sets of pulses

    total_times = np.array(len(selected_experinent['measuring_light_names'][index]) *  [0]) 


    for loop in range(0, selected_experinent['number_loops'][index]):  # go through each loop
        
        nmls = len(selected_experinent['measuring_light_names'][index])
        # find the timings for each meauring pulse and LED in the loop
        loop_times = np.array(list(accumulate(selected_experinent['measuring_interval'][index][loop * nmls : (loop + 1) * nmls]))) 
        
        # go through each measuring pulse in the loop
        
        for i in range(0,selected_experinent['number_pulses'][index][loop]):

            total_times = np.vstack((total_times, loop_times + cumulative_time)) # This adds all the times for each LED
            cumulative_time =  total_times[-1][-1]   # increment the cumulative_time holder

    for i, measuring_light_name in enumerate(selected_experinent['measuring_light_names'][index]):
        col_name = str(measuring_light_name) + '_time'
        selected_experinent[col_name][index] = total_times[1:,i]



def add_object_column (data_frame, column_name):

    """
    Add a new column to a data frame and set it to type object.
    This is used, for example, when one wants to add traces (arrays) to cells

    data_frame: the data_frame 
    column_name: string for the name of the column

    Returns:
        Modified data_frame with the added column

    """

    data_frame[column_name] = 0
    data_frame[column_name] = data_frame[column_name].astype(object)
    return data_frame



# Given a selected experimenet, index and a trace name, find the Fs and Fmp values.
# It is assumed that there will be an Fmp measurement when the sat_intensity_setting (default = 256) is found.
# The function returns a list of dictionaries (one entry per sat pulse found) containing the
# time that the sat pulse ofccured as well as the fluorscence parameters.  

def mean_top (arr, av_number, skip = 0):
    return np.mean(np.flip(np.sort(arr))[skip:av_number])

def mean_bottom (arr, av_number, skip = 0):
    return np.mean(np.sort(arr)[skip:av_number])

def capture_actinic_trace_segments (df, 
                                    col_name,
                                    capture_at_these_actinic_settings,
                                    time_col_name = '',
                                    capture_segment = 'during',
                                    av_method = 'top',
                                    sub_select = [0,-1],
                                    av_number = 3,
                                    skip = 0
                                   ):
        sub_traces = pd.Series()
        sub_traces = sub_traces.astype('object')
        
        sample_index = pd.Series()
        sample_index = sample_index.astype('object')

        sample_times = pd.Series()
        sample_times = sample_times.astype('object')

        sample_times_rel = pd.Series()
        sample_times_rel = sample_times_rel.astype('object')

        offset_index = {'during' : 0,
                        'before' : -1,
                        'after': 1}
        
        if (time_col_name == ''):
            time_col_name = col_name.split('_')[0] + '_time'
            
            if (time_col_name in df.columns):
                skip_rel_time = False
            else: 
                skip_rel_time = True
                
        for index in df.index:
            actinic_intensities = np.array(df['actinic_intensity'][index])
            
            # subtrace_indexes_with_actinic_settings that match any in capture_at_these_actinic_settings 
            subtrace_indexes = [i for i in range(len(actinic_intensities)) if actinic_intensities[i] in capture_at_these_actinic_settings]

            # get \number of pulses for each subtrace
            pulses = df['number_pulses'][index] 
            
            # Find cumulative number of pulses, i.e. at what index does an actinic light occur
            cumulative_index = [0] + list(accumulate(pulses))  # cumulative indexes
#             print(cumulative_index)
            sub_traces[index] = []
#             df['cumulative_index'][index] = cumulative_index
            sample_index[index] = []
            sample_times[index] = []
            sample_times_rel[index] =[]
#             t = []
            for j, i in enumerate(subtrace_indexes):
                
                capture_index = i + offset_index[capture_segment]
                if (capture_index < 0):
                    capture_index = 0
                if (capture_index > len(cumulative_index)):
                    capture_index = -1
                beg_index = cumulative_index[capture_index]
                end_index = cumulative_index[capture_index+1]
                sample_index[index].append(beg_index)
                sample_times[index].append(df[time_col_name][index][beg_index])
                
#                 sample_times[index].append(df[time_col_name][index][beg_index])
#                 print(beg_index, end_index)
                sub_trace = df[col_name][index][beg_index:end_index][sub_select[0]:sub_select[1]]
                if (av_method == 'top'):
                    
                    av_val =  mean_top (sub_trace, av_number, skip = skip)
                elif (av_method == 'bottom'):
                    av_val =  mean_bottom (sub_trace, av_number, skip = skip)
                else:
                    av_val = np.mean(sub_trace)
                sub_traces[index].append(av_val)
#             print(t)
            
            tt = np.array(sample_times[index])
#             print(tt)
            sample_times_rel[index] = tt-tt[0]
            
        return sub_traces, sample_index, sample_times, sample_times_rel

def trace_actinic_intentites (df, capture_at_these_actinic_settings ):                
    sample_index = {}
    for index in df.index:
        actinic_intensities = np.array(df['actinic_intensity'][index])
        
        # subtrace_indexes_with_actinic_settings that match any in capture_at_these_actinic_settings 
        subtrace_indexes = [i for i in range(len(actinic_intensities)) if actinic_intensities[i] in capture_at_these_actinic_settings]

        # get \number of pulses for each subtrace
        pulses = df['number_pulses'][index] 
        
        # Find cumulative number of pulses, i.e. at what index does an actinic light occur
        cumulative_index = [0] + list(accumulate(pulses))  # cumulative indexes
#             print(cumulative_index)
        # sub_traces[index] = []
#             df['cumulative_index'][index] = cumulative_index
        sample_index[index] = []
        # sample_times[index] = []
        # sample_times_rel[index] =[]
#             t = []
        for j, i in enumerate(subtrace_indexes):
            
            capture_index = i #+ offset_index[capture_segment]
            if (capture_index < 0):
                capture_index = 0
            if (capture_index > len(cumulative_index)):
                capture_index = -1
            beg_index = cumulative_index[capture_index]
            end_index = cumulative_index[capture_index+1]
            for ii in range(int(beg_index), int(end_index)):
                sample_index[index].append(ii)
            # sample_times[index].append(df[time_col_name][index][beg_index])
            
#                 sample_times[index].append(df[time_col_name][index][beg_index])
#                 print(beg_index, end_index)
#                 sub_trace = df[col_name][index][beg_index:end_index][sub_select[0]:sub_select[1]]
#                 if (av_method == 'top'):
                
#                     av_val =  mean_top (sub_trace, av_number, skip = skip)
#                 elif (av_method == 'bottom'):
#                     av_val =  mean_bottom (sub_trace, av_number, skip = skip)
#                 else:
#                     av_val = np.mean(sub_trace)
#                 sub_traces[index].append(av_val)
# #             print(t)
        
#             tt = np.array(sample_times[index])
# #             print(tt)
#             sample_times_rel[index] = tt-tt[0]
        
    return sample_index

# Dave's version of Pandas "groupby" function that works on arrays, lists, and text
# An automated averager/combiner for multiple rows in a dataframe containing floats, ints, text
# or array-like objects (lists and arrays).

#
# Pass:  all_traces: a dictionary of dataframes, called,
#              each of these will be selected for the appropriate rows to average.
#       selected_groups: which of the dataframes to include in the 
#              averageing/pompression

# Returns a dictionary containing the averaged/compressed rows, with
# all columns as keys.
# Rules:
#    if cells are all floats, return average of floats
#    if cells are all ints, return average of ints converted into int
#    if array-like, return numpy array of averages of each element
#    if str, return the first instance 

def average_compress_rows (df, selected_groups, exclude_indexes = []):
    # define what data types we are lookign for. These will be treated differently
    
    # separate by data types, aggegate those in the touble
    data_type_sets = {'strings': (str), 
                      'numbers': (float, int, np.floating, np.integer),
                      'array_like': (list, np.ndarray)}

    new_series_dict = {} # save everything in a dictionary, will be converted to dataframe later

#     for selected_wl_ranges in selected_groups: # Do this in two stages so we only average the correct data sets

#         df = all_traces[selected_wl_ranges] # select the data subset

    for data_type_set in data_type_sets.keys():  # now, we go through all the expect data types

        selected_types = data_type_sets[data_type_set]
        for col in df.columns: #= 'start_time' #'trace_label'

            # Only cells rows that have the desired type
            # and ignore rows with other types. This will
            # eliminate cases where the data is Nan, for instance. 
            rows_with_selected_type = df[df[col].apply(lambda x: isinstance(x, selected_types))][col]
            
            if (len(rows_with_selected_type) > 0):  # there are some values of the selected type
                if (len(rows_with_selected_type) == len(df)):  # all rows have the same type
                    if (data_type_set == 'array_like'):  # if the data is in an array form

                        if (len(np.array(rows_with_selected_type))>0): # don't try to average zero length arrays
                            # check if all members of the arrays are numbers.
                            if all(all(isinstance(e, (int, float, np.floating, float, int, np.integer)) for e in np.array(r)) for r in rows_with_selected_type):
                                try: #catch all errors, e.g. the arrays are not the same size
                                    new_series_dict[col] = np.mean(np.array([np.array(r ) for r in rows_with_selected_type]), axis = 0)
                                    new_series_dict[col + '_std'] = np.std(np.array([np.array(r ) for r in rows_with_selected_type]), axis = 0)
                                except:
                                    print('error averaging col', col)
                                    new_series_dict[col] = []
                            else:  # In some cases, e.g. wl names, the array contains text.
                                   # For the moment, assume these are all the same
                                new_series_dict[col] = np.array(rows_with_selected_type)[0]  # use the first instance
                    else:
                        if (len(set(rows_with_selected_type)) == 1):  # check to see if all values are the same
                            new_series_dict[col] = list(set(rows_with_selected_type))[0]
                        elif  any(isinstance(i, (float, np.floating)) for i in rows_with_selected_type): # if any value is a float average them!
                            new_series_dict[col] = np.mean(rows_with_selected_type)
                        else:  # if all are ints, we can take average and convert back to int
                            print("ints", col)        
                            new_series_dict[col] = int(np.mean(rows_with_selected_type))
    return new_series_dict

def subtract_line_from_trace (x, y, ranges):
        
    x_values = np.array([])
    y_values = np.array([])
#     print(ranges)
    for rr in ranges:
        x_values = np.concatenate((x_values, x[rr[0]:rr[1]]))
        y_values = np.concatenate((y_values, y[rr[0]:rr[1]]))
    lin_fit = linregress(x_values, y_values)

    fit_line_y = x * lin_fit.slope + lin_fit.intercept
    corrected_y_values = y - fit_line_y
    return corrected_y_values

def subtract_line (df, x_column, y_column, ranges, index):
        
    x_values = np.array([])
    y_values = np.array([])
#     print(ranges)
    for rr in ranges:
        # print(df[x_column]) #[index][rr[0]:rr[1]])
        x_values = np.concatenate((x_values, df[x_column][index][rr[0]:rr[1]]))
        y_values = np.concatenate((y_values, np.array(df[y_column][index])[rr[0]:rr[1]]))
    lin_fit = linregress(x_values, y_values)

    fit_line_y = df[x_column][index] * lin_fit.slope + lin_fit.intercept
    corrected_y_values = df[y_column][index] - fit_line_y
    return corrected_y_values

def subtract_line_multiple_wavelengths (df, x_column_suffix, y_column_suffix, wavelength_names, ranges, new_suffix = '', indexes = []):
    # print(wavelength_names)
    """

    returns dataframe with subtracted traces using multiple indexes to describe the baseline

    x_column_suffix which column to use for trace x data
    x_column_suffix which column to use for trace y data
    
    Works over a list of wavelengths in wavelength_names used as prefixes to the _suffix names

    ranges are a list of lists containing at least one region of indexes to use for the baseline.

    new_suffix is the new column name suffix to store data. If set to '' will overwrite the current data. 

    indexes is a list of row names. If zero length will do all rows

    """
    # plt.figure()
    if len(indexes) == 0: 
        indexes = df.index

    for wavelength_name in wavelength_names:
        # print(wavelength_name)
        y_column = str(wavelength_name) + y_column_suffix
        x_column = str(wavelength_name) + x_column_suffix
        if (new_suffix != ''):
            df[y_column + new_suffix] = 0
            df[y_column + new_suffix] = df[y_column + new_suffix].astype(object)
        for index in indexes:

            corrected_y_trace = subtract_line (df, x_column, y_column, ranges, index)
            df[y_column + new_suffix][index] = corrected_y_trace
        # for index in indexes:
        #     plt.plot(df[y_column + new_suffix][index])
        # plt.show()
    return df

def fit_burn(x, a, b, c):
    
    """
        Function for fitting absorbance traces to a funciton to account for 
        I will fit the baseline with a hyperbolic or exponential
        
        a is the amplitude of the burn-in
        b is the time constant
        c is the offset

    """

    fit_burn = c + a/(1+x/b)
    return (fit_burn)

    
def subtract_burn_region (x, y, ranges):
        
    x_values = np.array([])
    y_values = np.array([])
    # print(ranges)
    for rr in ranges:
        # print(rr)
        x_values = np.concatenate((x_values, x[rr[0]:rr[1]]))
        y_values = np.concatenate((y_values, y[rr[0]:rr[1]]))
#     plt.plot(x_values, y_values)
    return x_values, y_values
#     lin_fit = linregress(x_values, y_values)

#     fit_line_y = df[x_column][index] * lin_fit.slope + lin_fit.intercept
#     corrected_y_values = df[y_column][index] - fit_line_y
#     return corrected_y_values
# from scipy.interpolate import CubicSpline


def subtract_burn_artifact (df, xcol, ycol, index):
    """
    Returns the corrected trace removing the "burn-in" artifacts caused by LED heating.
    Uses the baseline points to estimate the burn-in
    parameters:
        df is the dataframe name
        xcol is the time column
        ycol is the absorbance values
        index is the row index
    """
    number_baseline_pulses = (df['number_pulses'][index][0])  # find the number of pulses used in the baseline,
                                                              # i.e., the first set of pulses
    trace = df[ycol][index] # get tje trace 
    trace = trace - trace[0]  # subtract the first point (makes it easier to calculate things)
    trace_x = df[xcol][index]  # get the time (x) values

    s_x, s_y = subtract_burn_region (trace_x, trace, [[0, int(0.7 * number_baseline_pulses)]])
    popt, pcov = curve_fit(fit_burn, s_x, s_y, p0=[.1,.1,.001], maxfev = 30000)
#     ,bounds=([-10, -10, -10], [10, 10, 10]), )
    fit = fit_burn(trace_x, *popt)
    return trace - fit

def plot_test_figure (df, xcol, ycol, std_col = None):
    plt.figure()
    for index in df.index:
        x = df[xcol][index]
        trace = df[ycol][index]
        plt.plot(x, trace)
        if std_col:
            std = df[std_col][index]
            plt.fill_between(x, trace + std, trace - std, alpha = 0.3)
    plt.show()

import re

def replace_set_temperature_with_values_from_file_name (df, set_as = 'int', 
                                                        new_col_name = 'set_temperature', 
                                                        temperature_end_chars = 'c'):
    """
    Take in dataframe of IDEASpec data, as df.
    Make new column list of all unique set temperatures based on the file names

    It is assumed that the file name contains the following information separated by _:

    genotype_XXc_experimentalConditionsetc.

    where XXc is the set temperature and follows the first _ and ends with temperature_end_char,
    which has a default values of lower case c, but could be changed to, for example, ''

    set_as determines what data type is written, options are 'int', 'float' and 'text'

    new_col_name holds the column name, with defaults as 'set_temperature'



    """ 
    
    if 'set_temperature' not in df.columns:
        df[new_col_name] = 'None'
    for index in df.index:
        experiment_name = df['experiment_name'][index]
        # split_strings = experiment_name.split('_')
        split_strings = re.split('_|,| ',experiment_name)

        set_temperature_value = ''
        for split_string in split_strings:
            # print(split_string)
            if split_string.endswith(temperature_end_chars):
                set_temperature_value = split_string 

        # set_temperature_value = [1]
        # if set_temperature_value.endswith(temperature_end_chars):
        if set_temperature_value == '':
            print('Did not find set temperature')
            df[new_col_name][index] = None
        else:
            if set_as == 'text':
                df[new_col_name][index] = set_temperature_value
            elif set_as == 'int':
                df[new_col_name][index] = int(set_temperature_value[0:-1])
            elif set_as == 'float':
                df[new_col_name][index] = float(set_temperature_value[0:-1])
    return df

# averaged_compressed_dataframe = replace_set_temperature_with_values_from_file_name(averaged_compressed_dataframe)

# averaged_compressed_dataframe['set_temperature']

def print_all_col_names (df):
    for c in df.columns:
        print(c)

# kinetic filters 
from scipy import signal
from scipy.signal import butter, lfilter, freqz

# from scipy.signal import butter, lfilter

def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y


def generate_filtered_specs (df, wls, source_x_col = '_time', source_y_col = '_rDA', 
                    destination_col_name = 'all_filtered_specs', 
                    bandpasses = [[100, 1000], [30, 300], [10, 100], [1, 10]],
                    threshold_fraction = 0.1, 
                    normalize_each_das = True,
                    save_kinetics_traces = False,
                    plot_it = False, 
                    verbose = False):

    df[destination_col_name] = df[destination_col_name].astype(object)
    
    for i, bandpass in enumerate(bandpasses):
        for wl in wls:
            df[destination_col_name + 'kinetics_' + str(wl)] = 0
            df[destination_col_name + 'kinetics_' + str(wl)] = df[destination_col_name + 'kinetics_' + str(wl)] .astype(object)

            df[destination_col_name + 'kinetics_time_' + str(wl)] = 0
            df[destination_col_name + 'kinetics_time_' + str(wl)] = df[destination_col_name + 'kinetics_time_' + str(wl)] .astype(object)

    if save_kinetics_traces:
        print('not implemented')

    for j, index in enumerate(df.index):
        for i, bandpass in enumerate(bandpasses):
            if verbose:
                print('bandpass = ', bandpass)
            if plot_it:
                plt.figure()
                experiment_name = np.array(df['experiment_name'][index])
                plt.title(experiment_name)
            filtered_traces = []
            for wl in wls:
                raw_y_data = np.array(df[str(wl) + source_y_col][index])
                raw_x_data = np.array(df[str(wl) + source_x_col][index])

                x_data = np.linspace(0, raw_x_data[-1], 3 * bandpass[1])
                data = np.interp(x_data, raw_x_data, raw_y_data)  
                if plot_it:
                    plt.plot(data, alpha = 0.2)
                # fs = 3000
                fs = np.max(bandpass) * 2.2
                sos = signal.butter(2, bandpass, 'bandpass', fs=fs, output='sos') #butter_highpass_filter(data, 0.99, 10)
                filtered = signal.sosfilt(sos, data)

                df[destination_col_name + 'kinetics_' + str(wl)][index] = filtered
                filtered_traces.append(filtered)
                if plot_it:
                    plt.plot(filtered, label = wl)
            if plot_it:
                plt.legend()
                plt.show()

            filtered_specs = np.array(filtered_traces).T
            if (i ==0):
                all_filtered_specs = filtered_specs
            else:
                all_filtered_specs = np.concatenate((all_filtered_specs, filtered_specs), axis = 0)
            
            threshold = threshold_fraction * mean_top([np.max(np.abs(m)) for m in all_filtered_specs], 5, 1)

            thresholded_spectra = [m for m in all_filtered_specs if np.max(np.abs(m)) > threshold]

            if normalize_each_das:
                thresholded_spectra = [m/np.max(np.abs(m)) for m in thresholded_spectra]

        df[destination_col_name][index] = thresholded_spectra
    return df

def select_data_by_criteria (df, groups={}, verbose = False):
    # Written by Sebastian 
    
    """
    Generate a query to select data from a dataframe, df,
    
    return a dictionary with selected indexes, query, and pass-through variables.
    
    Send df, and c dictionary of groups, 
        each groups entry has:
            name (key)
            'cols' is a dict containing:
                "name of column": [values to include]
           
           The function iterates through these cols and sub-selects the rows with the selected values.
           Thus, the final selected indexes are OR within each cols entry and AND begtween the different
           cols entries. 
    
    Sending a query without 'cols' will return all indexes in the dataframe.

    Example:
    
    #define a set of threee groups:
    
    groups = {}
    
    groups['one'] = {'cols': {'genotype': ['npq2','Col0'],
                 'set_temperature' : [10, 20]}, 
                'xKWARGS': {'use_filtered':True}}

    groups['two'] = {'cols': {'genotype': all_genotypes,
                 'set_temperature' : set_temperatures}
                }

    groups['three'] = {}

    
    select_data_by_criteria(averaged_compressed_dataframe, groups)
    
    # returns the following:
    
        {
            'one': {'params': {'genotype': ['npq2', 'Col0'], 'set_temperature': [10, 20]},
          'query': 'genotype in ["npq2","Col0"] & set_temperature in [10,20]',
          'KWARGS': {'set_temperature < ': [1e-05, 1e-06]},
          'idx': Int64Index([3, 10, 1, 8], dtype='int64')},


         'two': {'params': {'genotype': ['npq4', 'npq2', 'Col0', 'npq1'],
           'set_temperature': [10, 35, 20, 30]},
          'query': 'genotype in ["npq4","npq2","Col0","npq1"] & set_temperature in [10,35,20,30]',
          'KWARGS': {},
          'idx': Int64Index([3, 4, 10, 14, 0, 7, 9, 15, 1, 5, 8, 12, 2, 6, 11, 13], dtype='int64')},
         'three': {'params': {},
          'query': '',
          'KWARGS': {},
          'idx': Int64Index([3, 4, 10, 14, 0, 7, 9, 15, 1, 5, 8, 12, 2, 6, 11, 13], dtype='int64')}}
  
    
    """
    
    collect_idx = {}

    for group, keys in groups.items():
        query = ''
        if 'cols' in keys:
            
            for key, values in keys['cols'].items():
                # print(key)
                if query != '':
                    query = query + ' & '
                params = []
                for value in values:
                    # print(value)
                    if isinstance(value, str):
                        params.append('"%s"' % value)
                    else:
                        params.append(str(value))
                query = query + '{0} in [{1}]'.format(key, ','.join(params))
            if verbose:
                print(query)

        collect_idx[group] = {
            'params': keys['cols'] if 'cols' in keys else {},
            'query': query,
            'KWARGS': keys['KWARGS'] if 'KWARGS' in keys else {}, 
            'idx': df.query(query).index if query != '' else df.index
        }
    return collect_idx

def plot_ICA_components (ECS_A_, ICA_and_das_params, title = 'ICA components', verbose = True):
    
    # selected_indexes, n_components, threshold, sub_selected_wls
    
    plt.figure()
    
    plt.title(title)
    for i, a_fit in enumerate(ECS_A_.T):
        plt.plot(ICA_and_das_params['selected_wls'], a_fit/np.max(np.abs(a_fit)), label = '# ' + str(i))
    plt.xlabel('wavelength (nm)')
    plt.ylabel('relative extinction coefficients')
    plt.legend()
    plt.show()
    if verbose:
        pprint.pprint(ICA_and_das_params)

def compute_ICA (all_das, n_components, all_das_params = {}, whiten = "arbitrary-variance"):

    # set the desired number of components
    # whiten = 'unit-variance'
    algorithm = 'parallel' #, 'deflation'
    tol = 1e-5
    ECS_ica = FastICA(n_components = n_components, whiten = whiten, algorithm = algorithm,
            max_iter = 10000,
            fun = "exp", tol = tol) #True) #'unit-variance') #"arbitrary-variance")
    ECS_S_ = ECS_ica.fit_transform(all_das, all_das)  # Reconstruct signals
    ECS_A_ = ECS_ica.mixing_  # Get estimated mixing matrix
    params = ECS_ica.get_params()
#     print(all_das_params)
    z = params.copy()
    z.update(all_das_params)
    return ECS_S_, ECS_A_, z

import numpy as np
from scipy.spatial import procrustes

# a = np.array([[1, 6], [1, 2], [1, 1], [2, 1]], 'd')
# b = np.array([[4, -2], [4, -4], [4, -6], [2, -6]], 'd')
# mtx1, mtx2, disparity = procrustes(a, b)
# disparity

def procrustes_distance (ref_spec, test_spec):
    """
    Use procrustes to rotate and scale a test spectra to match as best as possible a reference spectrum.
    Return the "disparity", a measure of how similar the matched spectra are. 
    
    """
    wls = np.linspace(0, len(ref_spec)-1, len(ref_spec))
    results = spectral_disparity (wls, ref_spec, test_spec)
    
    return results['disparity']

def spectral_disparity (wls, ref_spec, test_spec):
    """
    Use procrustes to rotate and scale a test spectra to match as best as possible a reference spectrum.
    Return the rotated, scaled spectra as well as the "disparity", a measure of how similar the 
    matched spectra are. 
    
    """
    a = np.vstack([wls, ref_spec]).T
    b = np.vstack([wls, test_spec]).T

    mtx1, mtx2, disparity = procrustes(a, b)
    return {"mtx1": mtx1, "mtx2": mtx2, "disparity": disparity}

from scipy.cluster.hierarchy import fclusterdata


def procrusties_heirarchical_clusters (all_ECS_A_, desired_number_of_clusters, verbose = False, wls = []):
    #  fclust3 = fclusterdata(all_spectra, t, criterion = 'distance', metric = procrustes_distance)
    # Use the fclusterdata function to cluster by the procrustes_distance function.
    # Rotate the spectra in each cluster to be similar to the first spec in the cluster,
    # return a dict containing arrays of the spectra in each unique cluster
    # make an array of all spectra we want to cluster:

    all_spectra = np.array([])
    all_spectra_dict = {}
    for i, key in enumerate(all_ECS_A_.keys()):    
        for j, spectrum in enumerate(all_ECS_A_[key].T):
            if (len(all_spectra) == 0):
                    all_spectra = np.array(spectrum)
            else:
                all_spectra = np.vstack([all_spectra, np.array(spectrum)])
            all_spectra_dict[key + "_" + str(j)] = np.array(spectrum)
    
    fclust3 = fclusterdata(all_spectra, desired_number_of_clusters, 
                           criterion = 'maxclust', metric = procrustes_distance)
    distinct_clusters = list(set(fclust3))
    if verbose:
        print(distinct_clusters)
    spec_groups = {}
    for cluster in distinct_clusters:
        if verbose:
            print(np.where(fclust3 == cluster))
        spec_groups[cluster] = all_spectra[np.where(fclust3 == cluster)]
    spec_groups_rotated = {}
    for key in spec_groups:
        spec_groups_rotated[key] = []
        for i in range(len(spec_groups[key])):
            results = spectral_disparity (wls, spec_groups[key][0], spec_groups[key][i])
            spec = np.array(results['mtx2'][:,1])
            spec = spec/np.max(np.abs(spec))
            spec_groups_rotated[key].append(spec)

    return spec_groups_rotated, all_spectra_dict, all_spectra


def procrusties_heirarchical_clusters_df (df, desired_number_of_clusters, ECS_A_col = 'ECS_A_', verbose = False, wls = []):
    #  fclust3 = fclusterdata(all_spectra, t, criterion = 'distance', metric = procrustes_distance)
    # Use the fclusterdata function to cluster by the procrustes_distance function.
    # Rotate the spectra in each cluster to be similar to the first spec in the cluster,
    # return a dict containing arrays of the spectra in each unique cluster
    
    # make an array of all spectra we want to cluster:
    all_spectra = np.array([])

    all_spectra_dict = {}
    all_spectra_df = pd.DataFrame(columns = ['origin_index', 'spec_index', 'cluster', 'rotated_spectrum'])

    df['component_spectra'] = 0
    df['component_spectra'] = df['component_spectra'].astype(object)
    # experiment_components = []
    all_index = 0

    for i, index in enumerate(df.index):    
        spectra = df[ECS_A_col][index].T
        for j, spectrum in enumerate(spectra):
            if (len(all_spectra) == 0):
                    all_spectra = np.array(spectrum)
            else:
                all_spectra = np.vstack([all_spectra, np.array(spectrum)])
        # df['component_spectra'][index] = np.array(spectrum)
            all_spectra_dict[all_index]['origin_index'] = index

            all_spectra_dict[all_index]['spec_index'] = j  # record which spectrum this was
            #  index + "_" + str(j)] = np.array(spectrum)
            all_index = all_index + 1

        all_spectra_df[all_index]

    # get the desired number of distinct clusters for the entire data set

    fclust3 = fclusterdata(all_spectra, desired_number_of_clusters, 
                           criterion = 'maxclust', metric = procrustes_distance)

    # get a list of the distinct clusters for the above
    distinct_clusters = list(set(fclust3))

    if verbose:
        print(distinct_clusters)

    # get the sets of spectra and indexes for each distinct cluster
    # and place into the spec_groups dictionary

    spec_groups = {}
    for cluster in distinct_clusters:
        if verbose:
            print(np.where(fclust3 == cluster))
        spec_groups[cluster] = all_spectra[np.where(fclust3 == cluster)]
    spec_groups_rotated = {}

    for key in spec_groups:
        spec_groups_rotated[key] = []
        for i in range(len(spec_groups[key])):
            results = spectral_disparity (wls, spec_groups[key][0], spec_groups[key][i])
            spec = np.array(results['mtx2'][:,1])
            spec = spec/np.max(np.abs(spec))
            spec_groups_rotated[key].append(spec)

    return spec_groups_rotated, all_spectra_dict, all_spectra

    