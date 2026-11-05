## Load the pre-DR4 data, select one of the 12 targets, and convert its epoch astrometry to
## the format expected by Octofitter

import sys
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.time import Time

## Unit of scan_pos_angle in the output: Octofitter v9 takes degrees, the archive's own
## unit (v8 took radians). ../tests/test_scan_angle.jl checks the converted files against
## the installed Octofitter.
SCAN_ANGLE_UNIT = 'deg'


def gaia_data_converter(dr4_source_id, outfile='GaiaDR4_data.csv'):
    #dr4_source_id = int(sys.argv[1])
    keep_cols = ['transit_id', 
                 'obs_time_tcb', 'scan_pos_angle', 'centroid_pos_al', 
                 'centroid_pos_error_al', 'ipd_error_al', 'used_by_agis_al', 
                 'parallax_factor_al', 'ra0', 'dec0']
    array_cols = ['obs_time_tcb', 'scan_pos_angle', 'centroid_pos_al', 
                  'centroid_pos_error_al', 'ipd_error_al', 'used_by_agis_al']
    
    raw_data = Table.read("GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml").to_pandas()
    single_source_data = raw_data.query("source_id==@dr4_source_id")[keep_cols]
    
    ## Create an array of masks for each transit
    ssd = single_source_data.copy()

    ## Some columns have array values. Expand (aka 'explode') these 
    ## so that each value in an array gets its own row
    ssd["any_mask"] = ssd[array_cols].apply(
                       lambda row: np.any(
                          np.column_stack([
                             np.ma.getmaskarray(row[col])
                             for col in array_cols]),
                             axis=1),
                          axis=1,)
    ssd_expanded = ssd.explode(
                    array_cols + ["any_mask"],
                    ignore_index=True,)
    ssd_nomask = ssd_expanded.query("any_mask==False")
    
    ## Modify some cols to match Octofitter expectation
    ssd_nomask['epoch'] = ssd_nomask.obs_time_tcb /1e9/86400 + 2455197.5 - 2400000.5
    if SCAN_ANGLE_UNIT == 'rad':
        ssd_nomask['scan_pos_angle'] = ssd_nomask['scan_pos_angle'] * np.pi/180
    ssd_nomask['outlier_flag'] = ssd_nomask.used_by_agis_al.apply(lambda x: 0 if x else 1)
    
    
    final_columns = ['transit_id', 'epoch', 'scan_pos_angle', 'centroid_pos_al',
                     'centroid_pos_error_al', 'ipd_error_al', 'parallax_factor_al',
                     'used_by_agis_al', 'outlier_flag']
    
    clean_data = ssd_nomask[final_columns].sort_values(by='epoch')
    
    clean_data.to_csv(outfile)
    #import pdb; pdb.set_trace()
    
    return


def rv_data_converter(target_name):


    keep_cols = ['bjd', 'rv', 'σ_rv', 'tel']
    
    ## Target-specific choices
    if target_name=='HD114762':
        rv_df = pd.read_csv('RAW_RV_'+target_name+'.csv')
        rv_df.columns = keep_cols
        rv_df['mjd'] = Time(rv_df.bjd, format='jd').mjd
        #import pdb; pdb.set_trace()
        
        rv_df = rv_df.query("tel=='j'")
        #rv_df = rv_df.query("58299<mjd<58340") # 3 HIRES RVs in 2018
        rv_df = rv_df.query("57938<mjd<58364") # 10 RVs between Jul 2017-Jul2018
        
        
        #col_rename_dict = {'mjd':'epoch', 'mnvel':'rv', 'errvel':'σ_rv'}
        #rv_df = rv_df.rename(columns=col_rename_dict)
        
        rv_df = rv_df[['mjd', 'rv', 'σ_rv']]
        
    if target_name=='Gaia4':
        rv_df = pd.read_csv('RAW_RV_'+target_name+'.txt', 
                            skiprows=31, 
                            names=keep_cols,
                            sep='\s+')
        rv_df['mjd'] = Time(rv_df.bjd, format='jd').mjd
        
        #rv_df = rv_df.query("tel=='NEID'") # All 8 NEID RVs, between Dec 2023 - Sept 2024
        #rv_df = rv_df.query("tel=='NEID'").head(3) # First 3 NEID RVs in Dec. 2023
        #import pdb; pdb.set_trace()
        #rv_df = rv_df.query("tel=='HPF'").iloc[[0,4,9]] # 3 HPF RVs spaced over ~1.3 yr
        rv_df = rv_df.query("tel=='NEID'").iloc[[0,3,7]] # 3 NEID RVs spaced over ~6 months
        
        
        #rv_df = rv_df.query("tel=='HPF'") # All 10 HPF RVs, between Apr 2023 - Aug 2024
        #import pdb; pdb.set_trace()
        
        rv_df = rv_df[['mjd', 'rv', 'σ_rv', 'tel']]
        
    
    #import pdb; pdb.set_trace()
    
    #rv_df = rv_df[['mjd', 'rv', 'σ_rv']]
    rv_df = rv_df.rename(columns={'mjd':'epoch'})

    rv_df.to_csv(f'RV_data_{target_name}.csv', index=False)
    
    return




if __name__=="__main__":
    
    #gaia_data_converter(3937211745905473024)
    #rv_data_converter('HD114762')
    rv_data_converter('Gaia4')
















