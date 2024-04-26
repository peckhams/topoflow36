
#  May want to convert lat and lon in DMS to decimal degrees.
#  Are all stations in HCDN also in our collated TSV?

#  Copyright (c) 2023-2024, Scott D. Peckham
#
#  Jan 2024. Renamed all "ngen/utils" files to end in "_utils.py"
#            instead of "_tools.py"
#            Modified to use new data_utils.py. 
#  Nov 2023. Original version. 
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import hcdn_utils as hcdn
#  >>> hdcn.convert_dat_to_tsv()
#
#---------------------------------------------------------------------
#
#  convert_dat_to_tsv()
#
#---------------------------------------------------------------------

import numpy as np
from topoflow.utils.ngen import data_utils as dtu

#---------------------------------------------------------------------
def convert_dat_to_tsv( tsv_file='hcdn_stations.tsv'):

    #----------------------------------------------------------------
    # Note: Original file with ".dat" extension has data in columns
    #       of fixed widths, with spaces between.  The Name column
    #       can contain a comma.  stations.dat has no header line.
    #----------------------------------------------------------------
    data_dir = dtu.get_data_dir( 'USGS_HCDN' )
    new_dir  = dtu.get_new_data_dir( 'USGS_HCDN' )
    station_file = 'stations.dat'   
    station_unit = open( data_dir + station_file, 'r' )
    tsv_unit     = open( new_dir  + tsv_file, 'w')

    #----------------------
    # Write a header line
    #----------------------
    delim = '\t'
    header = ''
    header += 'usgs_stn_id' + delim
    header += 'usgs_stn_name' + delim
    header += 'huc' + delim
    header += 'area_sqmi' + delim
    header += 'time_scale_cd' + delim
    header += 'record_part_cd' + delim
    header += 'n_comment_lines' + delim
    header += 'district_code' + delim
    header += 'state_code' + delim
    header += 'county_code' + delim
    header += 'gage_lat_dms' + delim
    header += 'gage_lon_dms' + delim
    header += 'gage_elev_ft' + delim
    header += 'non_tca_pct'  + delim
    header += 'stn_type_code' + '\n'
    tsv_unit.write( header )
    
    #--------------------------------    
    # Process dat file line by line
    #--------------------------------
    print('Working...')
    while (True):
        line = station_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        usgs_stn_id     = line[0:9].strip()
        usgs_stn_name   = line[9:58].strip()
        huc             = line[58:67].strip()
        area_sqmi       = line[67:76].strip()
        time_scale_cd   = line[76:77].strip()
        record_part_cd  = line[77:79].strip()
        n_comment_lines = line[79:82].strip()
        district_code   = line[82:85].strip()
        state_code      = line[85:88].strip()
        county_code     = line[88:92].strip()
        gage_lat_dms    = line[92:100].strip()
        gage_lon_dms    = line[100:108].strip()
        gage_elev_ft    = line[108:117].strip()
        non_tca_pct     = line[117:121].strip()
        stn_type_code   = line[121:].strip()    # B=Benchmark, C=Current Conditions
                                                                 
        #-----------------------------         
        # Write new line to TSV file
        #-----------------------------
        line2 = ''
        line2 += usgs_stn_id   + delim
        line2 += usgs_stn_name + delim
        line2 += huc + delim
        line2 += area_sqmi       + delim
        line2 += time_scale_cd   + delim
        line2 += record_part_cd  + delim
        line2 += n_comment_lines + delim
        line2 += district_code   + delim
        line2 += state_code      + delim
        line2 += county_code     + delim
        line2 += gage_lat_dms    + delim
        line2 += gage_lon_dms    + delim
        line2 += gage_elev_ft    + delim
        line2 += non_tca_pct     + delim
        line2 += stn_type_code   + '\n'
        tsv_unit.write( line2 )
                                                                      
    #------------------       
    # Close the files
    #------------------             
    station_unit.close()
    tsv_unit.close()
    print('Finished converting DAT to TSV.')
    print()

#   convert_dat_to_tsv()
#---------------------------------------------------------------------

