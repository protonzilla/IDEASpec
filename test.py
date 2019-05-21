# Test the functions in IDEASpecPy
# These tests are designed to run as a part of the continuous integration pipeline
# These tests are designed to emulate the way the module is used in the example notebook
# These tests (and IDEASpecPy) do NOT follow the standard python package structure


# to start out, we have basic tests just assert that each function exists and can be run without raising an exception
# ideally, the output from each function would also be investigated during testing
# in a perfect world, we would follow Test-driven development (TDD)

import sys
import os
import pandas as pd

# "install" the module and import it
folder_that_contains_code="./IDEASpecPy/"
sys.path.append(os.path.abspath(folder_that_contains_code))
import IDEASpecPy as Ipy


# test the function generateCombinedJSON
# assert that outputs are compatible with pd.read_json
baseFileName = './IDEASpecPy/example_data/col-0_3g-simple_1'
Ipy.generateCombinedJSON(baseFileName)
sample_df=pd.read_json(baseFileName + '_combined.json') 

burnInFileName = './IDEASpecPy/example_data/card3_3g-simple_1'
Ipy.generateCombinedJSON(burnInFileName)
burn_in_profile=pd.read_json(burnInFileName + '_combined.json') 


# test the function smoothTraces
wavelengths=[475,488,505,520,535,545]
Ipy.smoothTraces(sample_df, wavelengths, '_I', 31, '_smooth')
Ipy.smoothTraces(burn_in_profile, wavelengths, '_I', 31, '_smooth')


# test the function compute_average_row
label_for_new_averaged_row = "avg"
rows_to_average = [0,1,2,3]
columns_to_average = ["475_I_smooth","488_I_smooth", "505_I_smooth", "520_I_smooth", "535_I_smooth", "545_I_smooth"]
Ipy.compute_average_row( burn_in_profile, rows_to_average, columns_to_average, label_for_new_averaged_row)

# test the function recalcAbsUsingBurnInProfile
Ipy.recalcAbsUsingBurnInProfile (sample_df, wavelengths,  '_I_smooth', burn_in_profile, 
                                 '_I_smooth', 'avg', I0_points=[4,50], newSuffix='_bDA')
								 
# test the function subtractBaseline
Ipy.subtractBaseline(sample_df, wavelengths, '_bDA', [100,500])

# test the function generateDAS
Ipy.generateDAS(sample_df, wavelengths, '_bDA', 'das')

# test the function defineFitCoefficients
ecs,qE,Zx,drift,scatter=Ipy.defineFitCoefficients()

# test the function fitDAS_5
Ipy.fitDAS_5(sample_df, wavelengths, 'das', newSuffix='')

