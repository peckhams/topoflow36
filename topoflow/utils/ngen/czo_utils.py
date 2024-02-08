
# Copyright (c) 2024, Scott D. Peckham
#
# Jan 2024. Started.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import czo_utils as cz
#  >>> cz.get_usgs_site_ids()
#  >>> 
#---------------------------------------------------------------------
#
#  get_usgs_site_ids()
#  get_czo_site_id()
#  get_czo_long_name()
#  get_czo_site_url()
#
#---------------------------------------------------------------------

import numpy as np

from topoflow.utils.ngen import data_utils as dtu

#---------------------------------------------------------------------
def get_usgs_site_ids():

    #-------------------------------------------------------
    # Get the USGS site IDs for sites close to CZO basins.
    # Checked "Juniata" is not supposed to be "Juanita".
    #--------------------------------------------------------------
    # There are also 10 headwater catchments in the KREW (Kings
    # River Experimental Watersheds).  Some are within Providence
    # Creek and some are in Bull Creek.
    # KREW is also the location for the Southern Sierra CZO.
    #--------------------------------------------------------------
    # Reference:
    # https://www.hydroshare.org/resource/a60faac4c73649efb4736de7ec587f7c/
    #--------------------------------------------------------------    
    id_list = [
    '06730200',   # Boulder Creek at North 75TH St. Near Boulder, CO; Boulder
    '06730500',   # Boulder Creek at Mouth Near Longmont, CO; Boulder 
    '06727500',   # Fourmile Creek at Orodell, CO; Boulder
    '06727410',   # Fourmile Creek at Logan Mill Road Near Crisman, CO; Boulder
    '021564493',  # Broad River Below Neal Shoals Res NR Carlisle, SC; Calhoun
    '02156500',   # Broad River Near Carlisle, SC; Calhoun
    '02160700',   # Enoree River at Whitmire, SC; Calhoun
    '02160105',   # Tyger River Near Delta, SC; Calhoun
    '08324000',   # Jemez River Near Jemez, NM; Catalina-Jemez
    '09484000',   # Sabino Creek Near Tucson, AZ; Catalina-Jemez
    '11475560',   # ELDER C NR BRANSCOMB CA; Eel
    '11475800',   # SF EEL R A LEGGETT CA; Eel
    '05573540',   # Sangamon River at Route 48 at Decatur, IL; IML
    '05570910',   # Sangamon River at Fisher, IL; IML
    '05572000',   # Sangamon River at Monticello, IL; IML
    '50063440',   # Quebrada Sonadora NR EL Verde, PR; Luquillo
    '50076000',   # Rio Blanco NR Florida, PR; Luquillo
    '50065500',   # Rio Mameyes NR Sabana, PR; Luquillo
    '50067000',   # Rio Sabana at Sabana, PR; Luquillo 
    '01558000',   # Little Juniata River at Spruce Creek, PA; Shale Hills
    '01558500'    # Shaver Creek near Petersburg, PA; Shale Hills
    ]
    
    return id_list

#   get_usgs_site_ids()
#---------------------------------------------------------------------
def get_czo_site_id( values ):

    #----------------------------------------------
    # Note:  CZO sites don't have an official ID,
    #        so we make one from other info.
    #----------------------------------------------
    obs_name = values[0].strip()
    p1 = obs_name.find('(')
    obs_abbr = obs_name[p1:]
    obs_abbr = obs_abbr[1:-1]
    #--------------------------------
    catch_name = values[1].strip()
    p2 = catch_name.find('(')
    catch_abbr = catch_name[p2:]
    catch_abbr = catch_abbr[1:-1]  
    #----------------------------------------------
    czo_id = 'CZO-' + obs_abbr + '-' + catch_abbr
    return czo_id

#   get_czo_site_id()
#---------------------------------------------------------------------
def get_czo_long_name( values ):

    catch_name = values[1].strip()
    state_code = values[2].strip()
        
    if ('(P30' in catch_name):
        long_name = catch_name
    else:
        p = catch_name.find('(')
        long_name = catch_name[:p].strip()
    long_name += ', ' + state_code

    return long_name
    
#   get_czo_long_name() 
#---------------------------------------------------------------------
def get_czo_site_url( site_id, REPORT=False,
                      ARCHIVE_SITES=True, ARCHIVE_BASE=False,
                      HYDROSHARE=False ):

    #----------------------------------------------------------
    # Note: Santa Catalina Mountains (SCM) and Jemez River
    # Basin (JRB) are now treated as one CZO for some reason,
    # called "Catalina-Jemez CZO".
    #----------------------------------------------------------
    # See CZO Catalina-Jemez User Profile in HydroShare:
    # https://www.hydroshare.org/user/5407/
    #----------------------------------------------------------
    czo_archive_base_url = 'https://czo-archive.criticalzone.org/'
    hydroshare_base_url  = 'https://www.hydroshare.org/group/'
    # earthchem_base_url   = 'https://www.earthchem.org/' # has chem data?
 
    #-----------------------------------------------------------------    
    # Notes:  (1) in CH CZO, there are 8 research areas; which one?
    #         (2) In Eel CZO, could not find URL for Elder Creek.
    #         (3) In JRB CZO, could not find URL for Upper Jaramillo.
    #         (4) In JRB CZO, could not find URL for La Jara.
    #         (5) In SH CZO,  could not determine watershed.
    #------------------------------------------------------------------           
    site_id_map = {
    'CZO-BC-GG':     ['boulder',         'gordon-gulch', 140],
    'CZO-CH-WS4':    ['calhoun',         'calhoun-czo-research-area-1', 148], # there are 8?
    'CZO-CRB-CR':    ['christina',       'christina-river-basin', 146],
    'CZO-ER-EC':     ['eel',             'eel-river-watershed', 141],  # not for elder creek
    'CZO-IML-A-US':  ['iml',             'sangamon-river-basin', 150],
    'CZO-IML-A2-US': ['iml',             'sangamon-river-basin', 150],
    'CZO-IML-B-CC':  ['iml',             'clear-creek-watershed', 150], 
    'CZO-IML-C-MR':  ['iml',             'minnesota-river-basin', 150],
    'CZO-JRB-UJ':    ['catalina-jemez',  'santa-catalina-mountains', 142], ## better site?
    'CZO-JRB-LJ':    ['catalina-jemez',  'santa-catalina-mountains', 142], ## better site?
    'CZO-LUQ-MAM':   ['luquillo',        'rio-icacos',  144],
    'CZO-LUQ-IC':    ['luquillo',        'rio-mameyes', 144],
    'CZO-RC-JD':     ['reynolds',        'johnston-draw', 143],
    'CZO-SCM-OR':    ['catalina-jemez',  'oracle-ridge-mid-elevation', 142],
    'CZO-SCM-MG':    ['catalina-jemez',  'bigelow-tower-marshall-gulch-high-elevation', 142],
    'CZO-SS-P301':   ['sierra',          'providence-creek-subcatchment-p301', 145],
    'CZO-SS-P303':   ['sierra',          'providence-creek-subcatchment-p303', 145],
    'CZO-SS-P304':   ['sierra',          'providence-creek-subcatchment-p304', 145],
    'CZO-SH-SH':     ['shale-hills',     'cole-farm-agricultural-site', 147] }  # what river?
    
    # Note:  hs_group_num = HydroShare Group number
    czo_abbr, site_abbr, hs_group_num = site_id_map[ site_id ]
    czo_field_areas_str = 'infrastructure/field-areas-'
    czo_field_area_str  = 'infrastructure/field-area/'  # specific area
    #----------------------------------------------------------------------
    czo_archive_url1  = czo_archive_base_url + czo_abbr + '/'
    czo_archive_url2  = czo_archive_url1
    czo_archive_url2 += czo_field_area_str   + site_abbr + '/'
    ## czo_archive_url2 += czo_field_areas_str  + czo_abbr + '/'
    #----------------------------------------------------------------------
    hydroshare_url = hydroshare_base_url + str(hs_group_num)
    
    if (REPORT):
        print('czo_url1 =', czo_archive_url1)
        print('czo_url2 =', czo_archive_url12)
        print('hydroshare_url =', hydroshare_url)
        print()

    if (ARCHIVE_SITES):
        return czo_archive_url2   # includes url1
    elif (ARCHIVE_BASE):
        return czo_archive_url1
    if (HYDROSHARE):
        return hydroshare_url

    #--------------------------------------------------
    # Note that CZEN also has some CZO websites like:
    #--------------------------------------------------
    # https://www.czen.org/content/calhoun-czo-1
        
    #------------------------------------------------------------
    # Note that individual field areas have URLs that are like
    # the one below.  Could scrape these to check data values.
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/catalina-jemez/infrastructure/
    #         field-area/oracle-ridge-mid-elevation/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/national/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/boulder/infrastructure/field-area/betasso/
    # https://czo-archive.criticalzone.org/boulder/infrastructure/field-area/gordon-gulch/
    # https://czo-archive.criticalzone.org/boulder/infrastructure/field-area/green-lakes-valley/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-1/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-2/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-3/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-4/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-5/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-6/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-7/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-8/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/catalina-jemez/infrastructure/field-area/santa-catalina-mountains/
    # https://czo-archive.criticalzone.org/catalina-jemez/infrastructure/field-area/oracle-ridge-mid-elevation/
    # https://czo-archive.criticalzone.org/catalina-jemez/infrastructure/field-area/bigelow-tower-marshall-gulch-high-elevation/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/christina/about/
    # https://czo-archive.criticalzone.org/christina/infrastructure/field-area/christina-river-basin/
    #------------------------------------------------------------
    # THERE SEEMS TO BE NO PAGE FOR ELDER CREEK.
    # https://czo-archive.criticalzone.org/eel/infrastructure/field-area/eel-river-watershed/
    # https://czo-archive.criticalzone.org/eel/infrastructure/field-area/rivendell/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/iml/infrastructure/field-area/sangamon-river-basin/
    # https://czo-archive.criticalzone.org/iml/infrastructure/field-area/clear-creek-watershed/
    # https://czo-archive.criticalzone.org/iml/infrastructure/field-area/minnesota-river-basin/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/bisley/
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/puente-roto/
    # Puente Roto USGS station:  USGS 50065500 RIO MAMEYES NR SABANA, PR
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/rio-blanco-nr-florida/
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/rio-icacos/
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/rio-mameyes/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/reynolds-creek-experimental-watershed/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/murphy-creek/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/tollgate/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/upper-sheep-creek/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/johnston-draw/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/mountain-big-sage/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/lower-sage/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/wyoming-big-sage/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/flats/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/shale-hills/infrastructure/field-area/cole-farm-agricultural-site/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/sierra/infrastructure/field-area/providence-creek-subcatchment-p301/        
    # https://czo-archive.criticalzone.org/sierra/infrastructure/field-area/providence-creek-subcatchment-p303/
    # https://czo-archive.criticalzone.org/sierra/infrastructure/field-area/providence-creek-subcatchment-p304/

#   get_czo_site_url()
#---------------------------------------------------------------------

