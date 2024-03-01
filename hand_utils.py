import arcpy
import os
import sys
import time
import json
import glob
from numpy import nan
from tqdm import tqdm
from urllib.request import urlopen, urlretrieve
import zipfile
import ssl
import pandas
import subprocess
import xarray
import multiprocessing
from functools import partial
import task

def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r...%d%%, %d MB, %d KB/s, %d seconds passed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()


def whichHUC6(AOI, name, wdir=None, outdir=None):

    arcpy.env.overwriteOutput = True

    crs = arcpy.SpatialReference(4326)

    pdir = os.path.join(os.getcwd(), name)
    if not os.path.isdir(pdir):
        os.makedirs(pdir)

    if wdir is None:
        cwd = os.path.join(pdir,name+'-tempdata')
        if not os.path.isdir(cwd):
            os.makedirs(cwd)
    else:
        cwd = wdir

    if outdir is None:
        odir = os.path.join(pdir,name+'-output')
        if not os.path.isdir(odir):
            os.makedirs(odir)
    else:
        odir = outdir

    AOI_proj = arcpy.Project_management(AOI, odir + "\\" + name + '-AOI_proj.shp', crs)
    ext = str(arcpy.Describe(AOI_proj).extent).split()[0:4]
    ext_ord = ext[0] + ',' + ext[1] + ',' + ext[2] + ',' + ext[3]

    url = 'https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/3/query?where=&text=&objectIds=&time=&geometry=' + ext_ord + '&geometryType=esriGeometryEnvelope&inSR=4326&spatialRel=esriSpatialRelIntersects&distance=&units=esriSRUnit_Foot&relationParam=&outFields=huc6&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=4326&havingClause=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&featureEncoding=esriDefault&f=pjson'
    data_json = json.loads(urlopen(url).read())

    if 'error' in data_json:
        print('There may be an issue with the WBD server...please check the server status at https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer')
        return()

    huc6s = []

    for f in data_json['features']:
        huc6 = f['attributes']['huc6']
        huc6s.append(huc6)

    return huc6s


def getHANDData(AOI, name, wdir=None, outdir=None, rm_temp=True):

    arcpy.env.overwriteOutput = True
    arcpy.env.parallelProcessingFactor = str(multiprocessing.cpu_count()-1)

    ssl._create_default_https_context = ssl._create_unverified_context

    pdir = os.path.join(os.getcwd(), name)
    if not os.path.isdir(pdir):
        os.makedirs(pdir)

    if wdir is None:
        cwd = os.path.join(pdir,name+'-tempdata')
        if not os.path.isdir(cwd):
            os.makedirs(cwd)
    else:
        cwd = wdir

    if outdir is None:
        odir = os.path.join(pdir,name+'-output')
        if not os.path.isdir(odir):
            os.makedirs(odir)
    else:
        odir = outdir
    
    huc6s = whichHUC6(AOI, name, cwd, odir)


    # should use multiprocessing to speed this up
    for huc_i in range(len(huc6s)):
        huc = huc6s[huc_i]
        url = 'https://cfim.ornl.gov/data/HAND/20200601/'+huc+'.zip'
        filedir = os.path.join(cwd, huc + '.zip')
        if os.path.exists(filedir):
            print(filedir,' already exists!\n')
            continue
        else:
            print('############################################################\n\n'+'('+str(huc_i+1)+' of '+str(len(huc6s))+') '+'Retrieving data from...'+
                url+'\n\n############################################################\n')
            urlretrieve(url, filedir, reporthook)
        
    os.chdir(cwd)
    for file in glob.glob('*.zip'):
        print('############################################################\n\nExtracting...'+ file +
              '\n\n############################################################\n')
        with zipfile.ZipFile(file, 'r') as zf:
            for member in tqdm(zf.infolist(), desc='Extracting '):
                try:
                    zf.extract(member, cwd)
                except zipfile.error as e:
                    continue
        os.remove(file)

    
    folds = next(os.walk(os.path.join(cwd,'.')))[1]
    wbd_paths = []
    for fold in folds:
        wbd_path = os.path.join(cwd,fold,fold+'-wbdraw.shp')
        wbd_paths.append(wbd_path)
    
    arcpy.Merge_management(wbd_paths,cwd+r'\wbd_merged.shp')
    intersect_paths = [cwd+r'\wbd_merged.shp',odir + "\\" + name + '-AOI_proj.shp']
    arcpy.Intersect_analysis(intersect_paths,cwd+r'\wbd_AOI_intersect.shp')

    huc6_int = []

    fc = cwd+r'\wbd_AOI_intersect.shp'
    fields = ['HUC6']
    with arcpy.da.SearchCursor(fc, fields) as cursor:
        for row in cursor:
            huc6_int.append(row[0])

    extr_folds = list(set(folds) - set(huc6_int))

    if len(extr_folds) > 0:

        print('############################################################\n\nDeleting extraneous HAND folders...\n\n'+'############################################################\n')

        for fold in extr_folds:
            extr_path = os.path.join(cwd,fold)
            os.remove(extr_path)
    else: 
        pass
    
    folds = next(os.walk(os.path.join(cwd,'.')))[1]
    hand_paths = []
    src_paths = []
    mask_paths = []
    dem_paths = []
    fline_paths = []
    for fold in folds:
        hand_path = os.path.join(cwd,fold,fold+'hand.tif')
        src_path = os.path.join(cwd,fold,'hydrogeo-fulltable-' + fold + '.csv')
        mask_path = os.path.join(cwd,fold,fold+'catchmask.tif')
        raw_dem_path = os.path.join(cwd,fold,fold+'.tif')
        trim_dem = arcpy.sa.ExtractByMask(raw_dem_path, hand_path)
        dem_path = os.path.join(cwd,fold,fold+'_dem_trim.tif')
        trim_dem.save(dem_path)
        fline_path = os.path.join(cwd,fold,fold+'-flows.shp')
        hand_paths.append(hand_path)
        src_paths.append(src_path)
        mask_paths.append(mask_path)
        dem_paths.append(dem_path)
        fline_paths.append(fline_path)
        
    print('############################################################\n\nMerging and Cropping Rasters...\n\n'+'############################################################\n')
    
    arcpy.MosaicToNewRaster_management(hand_paths, cwd, name + '_hand_merge.tif',coordinate_system_for_the_raster=4326,number_of_bands=1,pixel_type= "32_BIT_FLOAT")
    arcpy.MosaicToNewRaster_management(mask_paths, cwd, name + '_catchmask_merge.tif',coordinate_system_for_the_raster=4326,number_of_bands=1, pixel_type="32_BIT_SIGNED")
    arcpy.MosaicToNewRaster_management(dem_paths, cwd, name + '_dem_merge.tif', coordinate_system_for_the_raster=4326,number_of_bands=1, pixel_type='32_BIT_FLOAT')
    arcpy.Merge_management(fline_paths, cwd + r'\flowlines_merged.shp')
    hand_clip = arcpy.sa.ExtractByMask(cwd+"\\"+name+r'_hand_merge.tif',odir + "\\" + name + '-AOI_proj.shp')
    hand_clip.save(odir+"\\"+name+r'_hand.tif')
    catchmask_clip = arcpy.sa.ExtractByMask(cwd+"\\"+name+r'_catchmask_merge.tif',odir + "\\" + name + '-AOI_proj.shp')
    catchmask_clip.save(cwd+"\\"+name+r'_catchmask.tif')
    dem_clip = arcpy.sa.ExtractByMask(cwd+"\\"+name+r'_dem_merge.tif', odir + "\\" + name + '-AOI_proj.shp')
    dem_clip.save(odir+"\\"+name+r'_dem.tif')
    arcpy.Clip_analysis(cwd + r'\flowlines_merged.shp',odir + "\\" + name + '-AOI_proj.shp',odir + "\\"+ name + '-flowlines.shp')

    print('############################################################\n\nMerging SRC Tables...\n\n'+'############################################################\n')

    frames = [0]*len(src_paths)    
    for i in range(len(src_paths)):
        frames[i] = pandas.read_csv(src_paths[i])

    comb_tbl = pandas.concat(frames, keys=folds)
    comb_tbl.to_csv(odir+"\\"+name+r'-hydrogeo-fulltable.csv', chunksize = 1000000)

    if rm_temp is True:
        print('############################################################\n\nDeleting temp data...\n\n'+'############################################################\n')
        os.remove(cwd)
    else:
        pass 


def getNWMData(dates, name, wdir=None, outdir=None, nco_pth=None, rm_temp=True):

    arcpy.env.overwriteOutput = True
    arcpy.env.parallelProcessingFactor = str(multiprocessing.cpu_count()-1)

    pdir = os.path.join(os.getcwd(), name)
    if not os.path.isdir(pdir):
        os.makedirs(pdir)

    if wdir is None:
        cwd = os.path.join(pdir,name+'-tempdata')
        if not os.path.isdir(cwd):
            os.makedirs(cwd)
    else:
        cwd = wdir

    if outdir is None:
        odir = os.path.join(pdir,name+'-output')
        if not os.path.isdir(odir):
            os.makedirs(odir)
    else:
        odir = outdir
    
    if len(dates) > 1:
        strt_date = dates[0]
        end_date = dates[1]
        date_range = pandas.date_range(strt_date, end_date)
        for date in date_range:
            year = str(date.year)
            month = str(date.strftime('%m'))
            day = str(date.strftime('%d'))
            date_text = year + month + day + '1200' 
            url = "https://noaa-nwm-retrospective-2-1-pds.s3.amazonaws.com/model_output/" + year + "/" + date_text + ".CHRTOUT_DOMAIN1.comp"
            filedir = os.path.join(cwd,'flows_' + date_text + '.nc')
            print('\n\n############################################################\n\nRetrieving data from...\n'+
              url+'\nto...\n'+filedir+'\n\n############################################################\n\n')
            urlretrieve(url, filedir)

        print("\n\n############################################################\n\nCalculating NWM ensemble maxima with NCO...\n\n############################################################\n\n")
        if nco_pth is None:
            os.chdir(r'C:\nco')
        else:
            os.chdir(nco_pth)
        input_str = cwd + r'\flows*.nc'
        ofile = odir + "\\" + name + "_flows.nc"
        list1 = ['nces', '-ymax', '--no_tmp_fl']
        list2 = glob.glob(input_str)
        command = list1 + list2 + [ofile]
        subprocess.run(command)

    else:
        year = str(date.year)
        month = str(date.strftime('%m'))
        day = str(date.strftime('%d'))
        date_text = year + month + day + '1200'
        url = "https://noaa-nwm-retrospective-2-1-pds.s3.amazonaws.com/model_output/" + year + "/" + date_text + ".CHRTOUT_DOMAIN1.comp"
        filedir = os.path.join(odir,name + '_flows.nc')
        print('\n\n############################################################\n\nRetrieving data from...\n'+
              url+'\nto...\n'+filedir+'\n\n############################################################\n\n')
        urlretrieve(url, filedir)
    
    if rm_temp is True:
        print('\n\n############################################################\n\nDeleting temp data...\n'+'########################################\n\n')
        os.remove(cwd)
    else:
        pass 


def mapFlood(name, wdir=None, outdir=None, rm_temp=None):
    
    arcpy.env.overwriteOutput = True
    arcpy.env.parallelProcessingFactor = str(multiprocessing.cpu_count()-1)
    
    pdir = os.path.join(os.getcwd(), name)
    if not os.path.isdir(pdir):
        os.makedirs(pdir)

    if wdir is None:
        cwd = os.path.join(pdir,name+'-tempdata')
        if not os.path.isdir(cwd):
            os.makedirs(cwd)
    else:
        cwd = wdir

    if outdir is None:
        odir = os.path.join(pdir,name+'-output')
        if not os.path.isdir(odir):
            os.makedirs(odir)
    else:
        odir = outdir

    flows = xarray.open_dataset(odir + "\\" + name + r'_flows.nc')
    flows_df = flows.to_dataframe()
    flows_df = flows_df.reset_index()
    
    src = pandas.read_csv(odir+"\\"+name+r'-hydrogeo-fulltable.csv')

    src['hand_flow'] = ''
    src['hand_stage'] = ''
    
    pool = multiprocessing.Pool()
    print('\n\n############################################################\n\nCalculating catchment stage values...\n\n'+'############################################################\n\n')
    results = list(tqdm(pool.imap(partial(task.task, src = src, flows_df = flows_df), range(len(src['CatchId'].unique())), chunksize=750),total = len(src['CatchId'].unique())))
    pool.close()
    pool.join()

    results_df = pandas.concat(results).reset_index()
    results_df.to_csv(odir+"\\"+name+r'-hydrogeo-fulltable-hand_vals.csv', chunksize=1000000)

    arcpy.RasterToPolygon_conversion(cwd + "\\" + name + r'_catchmask.tif',odir + "\\" + name +'-catchments.shp', "NO_SIMPLIFY")
    arcpy.JoinField_management(odir + "\\" + name +'-catchments.shp', 'gridcode', odir+"\\"+name+r'-hydrogeo-fulltable-hand_vals.csv', 'CatchId', ['hand_flow','hand_stage'])
    arcpy.JoinField_management(odir + "\\" + name + '-flowlines.shp', 'COMID', odir+"\\"+name+r'-hydrogeo-fulltable-hand_vals.csv', 'CatchId', ['hand_flow','hand_stage'])
    arcpy.env.snapRaster = odir+"\\"+name+r'_hand.tif'
    arcpy.env.extent = odir+"\\"+name+r'_hand.tif' 
    arcpy.PolygonToRaster_conversion(odir + "\\" + name +'-catchments.shp','hand_stage',cwd + "\\" + name + '-catchments-stage.tif', cellsize=odir+"\\"+name+r'_hand.tif')
    hand_tif = odir+"\\"+name+r'_hand.tif'
    stage_tif = cwd + "\\" + name + '-catchments-stage.tif'
    dem_tif = odir+"\\"+name+r'_dem.tif'
    ras_list = [hand_tif,stage_tif]
    ras_names = ['hand', 'stage']
    depth = arcpy.sa.RasterCalculator(ras_list,ras_names,'Con(hand - stage > 0, 0, stage-hand)',extent_type = 'FirstOf')
    depth = arcpy.sa.ExtractByMask(depth,odir + "\\" + name + '-AOI_proj.shp')
    depth.save(odir + "\\" + name + '-depth_map.tif')
    bound = arcpy.sa.RasterCalculator([depth], ['depth'], 'depth > 0')
    bound = arcpy.sa.Reclassify(bound, "Value", arcpy.sa.RemapValue([[1,1], [0,'NODATA']]))
    bound = arcpy.RasterToPolygon_conversion(bound, odir + "\\" + name + '-flood_bound.shp', "NO_SIMPLIFY")
    elevation = arcpy.sa.RasterCalculator([depth,dem_tif], ['depth','dem'], 'depth + dem')
    elevation = arcpy.sa.ExtractByMask(elevation, bound)
    elevation.save(odir + "\\" + name + '-elevation_map.tif')

    if rm_temp is True:
        print('\n\n############################################################\n\nDeleting temp data...\n\n'+'############################################################\n\n')
        os.remove(cwd)
    else:
        pass 





