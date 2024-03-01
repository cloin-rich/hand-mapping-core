import os

script_dir = os.path.dirname(os.path.realpath('__file__'))
os.chdir(script_dir)

import hand_utils

AOI = os.path.join(script_dir, 'NC_Catch_MERIT_combined.shp') 

hand_utils.getHANDData(AOI, 'florence',rm_temp=False)




