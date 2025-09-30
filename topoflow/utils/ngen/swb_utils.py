
# Copyright (c) 2023-2024, Scott D. Peckham
#
# Jan 2024. Renamed all "ngen/utils" files to end in "_utils.py"
#           instead of "_tools.py".
#           Modified to use new data_utils.py.
# Oct 2023. Moved SWB routines here from gages2_utils.py.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import swb_utils as swb
#  >>> names = swb.get_swb_class_names()
#  >>> delta_p = 0.4
#  >>> f_s = 0.3
#  >>> phi = 2.5
#  >>> c = swb.get_swb_class(delta_p, f_s, phi)
#
#---------------------------------------------------------------------
#
#  get_swb_class_names()
#  get_swb_class()
#  get_seasonal_water_balance_class()   # alias
#
#---------------------------------------------------------------------
def get_swb_class_names( USE_B3=False ):

    if (USE_B3):
        names = ['A1','A2','A3','B1','B2','B3',
                   'C1','C2','D1','D2','D3','None']
    else:
        names = ['A1','A2','A3','B1','B2',
                   'C1','C2','D1','D2','D3','None']
    return names
    
#   get_swb_class_names()
#---------------------------------------------------------------------
def get_swb_class( delta_p, f_s, phi, ORIGINAL=False,
                   USE_B3=True, REPORT=False):

    #-------------------------------------------------------------
    # Note: Set ORIGINAL=True and USE_B3=False to use the
    #       original 10 class boundaries of the Seasonal Water
    #       Balance (SWB) classification system as defined by
    #       Berghuijs et al. (2014).
    #       This leaves 1008 of the 1947 basins unclassified.
    #-------------------------------------------------------------
    #       Set ORIGINAL=False and USE_B3=False to use expanded
    #       class boundaries that close the gaps between cuboids
    #       in the 3-parameter space, without any overlap.
    #       This results in all 1947 basins being classified.
    #       The 3D parameter space has the following parameters. 
    # 
    #       delta_p = precip_timing_index   in [-1,1]
    #       f_s     = snow_precip_fraction  in [0,1]
    #       phi     = aridity_index         in [0, infinity]
    #-------------------------------------------------------------
    #       Set ORIGINAL=False and USE_B3=True to use expanded
    #       class boundaries and to introduce a new class, B3,
    #       so that classes span the full 3-parameter space,
    #       without any overlap.  The class B3 "stacks above"
    #       the classes B1 and B2, just as A3 does for A1 & A2.
    #       So far, this extra class/region has not been needed.
    #-------------------------------------------------------------
    #       A Mathematica program has been written to visualize
    #       the set of "cuboid" regions in the parameter space.
    #-------------------------------------------------------------   
    d_p = delta_p  # synonym
    swb_class = 'None'

    #-------------------------------------------------------------
    # All A classes have -0.4 as the same max value of delta_p.
    # All C classes have 0.0 as the same min value of delta_p.
    # Classes D2 and D3 have -0.1 as the same min value of
    # delta_p, while the min value for D1 is -0.4.
    #-------------------------------------------------------------
    # This leaves a big void region in the 3D parameter space of
    # (delta_p, f_s, phi) that can be closed if we set these 5
    # delta_p values to be the same value in [-0.4, -0.1].
    #-------------------------------------------------------------
    A_dp_min  = -1.0   # (Berghuijs figure 7 & table 3)
    B_dp_min  = -1.0   # (Berghuijs figure 7 & table 3)
    B_dp_max  = 1.0    # (Berghuijs figure 7; not table 3)
    C_dp_max  = 1.0    # (Berghuijs figure 7 & table 3)
    D1_dp_min = -0.4   # (Berghuijs figure 7 & table 3)
            
    if (ORIGINAL):
        A_dp_max  = -0.4   # (Berghuijs figure 7 & table 3)
        # B_dp_max = 0.0   # (Berghuijs table 3)
        C_dp_min  = 0.0    # (Berghuijs table 3)
        # C_dp_min = 0.3   # (Berghuijs figure 7)
        D2_dp_min = -0.1   # (Berghuijs figure 7 & table 3)
        D3_dp_min = -0.1   # (Berghuijs figure 7 & table 3)
    else:
        #-----------------------------------------------------------
        # Note: For mid_delta_p, we can choose any value between
        #       -0.4 and -0.1 to "close the gap".  At first,
        #       this midpoint was set to -0.2.  However, when
        #       applied to all GAGES2 CONUS basins, the counts
        #          for A1, A2, and A3 were: 444, 188, and 22 &
        #          for D1, D2, and D3 were: 442, 3345, and 1684.
        #       To make the counts more equal, we are now 
        #       setting it to -0.1.  The new counts are:
        #          A1, A2, A3 = 553, 227, and 25 &
        #          D1, D2, D3 = 442, 3266, and 1668.
        #       This also means that C_dp_min is now less negative
        #       when its value was zero in Berghuis et al. (2014).
        #-----------------------------------------------------------
        ###### mid_delta_p = -0.2
        mid_delta_p = -0.1
        A_dp_max    = mid_delta_p 
        C_dp_min    = mid_delta_p
        D2_dp_min   = mid_delta_p
        D3_dp_min   = mid_delta_p
    
    #---------------------------------------------------------
    # The max value of delta_p for all the D classes can be
    # raised to 1.0 without intersecting any other classes.
    #---------------------------------------------------------
    if (ORIGINAL):
        ###############################################
        ## Currently, d_p = 1.0 often when f_s = 0.0,
        ## so D1 is not getting any basins.
        ###############################################
        D1_dp_max = 0.3      # (Berghuijs table 3)
        ## D1_dp_max = 0.25  # (Berghuijs figure 7)
        D2_dp_max = 0.3      # (Berghuijs figure 7 & table 3)
        D3_dp_max = 0.4      # (Berghuijs figure 7 & table 3)
    else:
        #--------------------------------------------
        # Helps: 951 goes down to 333 unclassified.  RECHECK ##########
        #-------------------------------------------- 
        D1_dp_max = 1.0   # (max allowed value)
        D2_dp_max = 1.0
        D3_dp_max = 1.0
    
    #--------------------------------------------------------------
    # The min value of phi for classes A1, B1, D1, D2, & D3 can
    # all be lowered to 0 without intersecting any other classes.
    #--------------------------------------------------------------
    if (ORIGINAL):
        A1_phi_min = 0.0      # (Berghuijs figure 7)
        ## A1_phi_min = 0.35  # (Berghuijs table 3)
        B1_phi_min = 0.0      # (Berghuijs figure 7)  #### BUG FIX 2023-11-29
        ## B1_phi_min = 0.4   # (Berghuijs table 3)
        D1_phi_min = 0.5      # (Berghuijs figure 7 & table 3)
        D2_phi_min = 0.5      # (Berghuijs figure 7 & table 3)
        D3_phi_min = 0.4      # (Berghuijs figure 7 & table 3)
    else:
        #-------------------------------------------
        # Helps: 333 goes down to 18 unclassified.
        #------------------------------------------- 
        A1_phi_min = 0.0
        B1_phi_min = 0.0
        D1_phi_min = 0.0
        D2_phi_min = 0.0
        D3_phi_min = 0.0
    
    #-------------------------------------------------------------
    # In Berghuijs et al. (2014), the min value of phi for class
    # C1 is 0.9, in both figure 7 or table 3, and the max value
    # of phi for all D classes is 0.9.  There is a void in the
    # parameter space for phi values < 0.9 below the C1 box.
    # But if we were to lower C1_phi_min, the box for class C1
    # would overlap all the D-class boxes.
    #-------------------------------------------------------------
    # A better alternative to filling this void seems to be to
    # increase the maximum value of delta_p for all of the D
    # classes to the max possible value of 1.0.   
    #-------------------------------------------------------------
    C1_phi_min = 0.9    # (Berghuijs figure 7 & table 3)
    D1_phi_max = 0.9    # (Berghuijs figure 7 & table 3)
    D2_phi_max = 0.9    # (Berghuijs figure 7 & table 3)
    D3_phi_max = 0.9    # (Berghuijs figure 7 & table 3)
  
    #-----------------------------------------------------------     
    # In Berghuijs, the max value of f_s for class A3 is 0.45.
    # The max value of f_s for classes C1 and C2 is 0.25.
    # All of these max values of f_s could be raised to 1.0
    # without intersecting any other class boxes.
    #-----------------------------------------------------------
    # However, if we introduce a class B3 analogous to A3,
    # as done now, then we cannot raise A3_fs_max.
    # If we USE_B3 and don't raise A3_fs_max, previously
    #   unclassified points will be assigned to class B3.
    # If we don't USE_B3 and raise A3_fs_max, previously
    #   unclassified points wil be assigned to class A3. 
    #-----------------------------------------------------------
    # As classification code is written now, we assume that
    # C1 and C2 both use the same C_fs_max.
    #-----------------------------------------------------------
    A_fs_min  = 0.0      # (Berghuijs figure 7 & table 3)
    A1_fs_max = 0.45     # (Berghuijs figure 7 & table 3)
    A2_fs_max = 0.45     # (Berghuijs figure 7 & table 3)
    #----------------
    B_fs_min  = 0.45     # (Berghuijs figure 7 & table 3)
    B_fs_max  = 1.0      # (Berghuijs figure 7 & table 3)
    #----------------
    C_fs_min  = 0.0      # (Berghuijs figure 7 & table 3)
    #----------------
    D1_fs_val  = 0.0     # (Berghuijs figure 7 & table 3)
    D2_fs_min  = 0.0     # (Berghuijs figure 7 & table 3)
    D2_fs_max  = 0.2     # (Berghuijs figure 7 & table 3)
    D3_fs_min  = D2_fs_max
    D3_fs_max  = 0.45    # (Berghuijs figure 7 & table 3)
     
    if (ORIGINAL or USE_B3):
        A3_fs_max = 0.45   # (Berghuijs figure 7 & table 3)
    else:
        A3_fs_max = 1.0
    #---------------------
    if (ORIGINAL):
        C_fs_max = 0.25    # (Berghuijs figure 7 & table 3)
    else:
        #-----------------------------------------
        # Helps: 
        #--------------------------------------------------
        # Can't set this higher than 0.45 or we intersect
        # the new bounds of B1, B2, etc.
        #--------------------------------------------------
        C_fs_max = 0.45

    #----------------------------------------------------------------    
    # In Berghuijs et al. (2014), the maximum value of phi for the
    # classes A3 and C2 are slightly different, namely, 5 and 5.3.
    # They can be safely raised to any higher value without
    # intersecting any other class box.  The max value of phi for
    # class B2 is 1.75, but it could also be raised to any higher
    # value.  Alternately, we could introduce a class B3 that has
    # phi values ranging from 1.75 to 5.3 (or some other phi_max),
    # and that is done here when USE_B3=True.
    #----------------------------------------------------------------
    A1_phi_max = 0.75        # (Berghuijs figure 7 & table 3)
    A2_phi_min = A1_phi_max
    A2_phi_max = 1.75        # (Berghuijs figure 7 & table 3)
    A3_phi_min = A2_phi_max
    #-----------------------
    B1_phi_max = 0.75        # (Berghuijs figure 7 & table 3)
    B2_phi_min = B1_phi_max
    B2_phi_max = 1.75        # (Berghuijs figure 7 & table 3)
    B3_phi_min = B2_phi_max  # (no B3 in Berghuijs)
    #-----------------------
    C1_phi_max = 1.5         # (Berghuijs figure 7 & table 3)
    C2_phi_min = C1_phi_max

    if (ORIGINAL):
        A3_phi_max = 5.0     # (Berghuijs figure 7 & table 3)
        #------------------------------------------------------
        C2_phi_max = 5.3     # (Berghuijs table 3)
        ## C2_phi_max = 5.0  # (Berghuijs figure 7)
        B3_phi_max = 5.3     # (No B3 in original SWB method)
    else:  
        #----------------------------------------------------------
        # In the GAGES-II CONUS, after converting units from mm
        # to cm, the max value of phi is 5.711.
        # The theoretical range is 0 to Infinity.
        # In the Berghuijs (2014) paper, phi max is 5.3, so here
        # we increase the max allowed to 5.8.
        #---------------------------------------------------------
        # Helps: 1 goes down to 0 unclassified.
        #----------------------------------------
        A3_phi_max = 5.8
        C2_phi_max = 5.8
        B3_phi_max = 5.8

    #--------------------------------
    # Check "A" classes for a match
    # "Precipitation out of phase"
    #--------------------------------
    if ((A_dp_min < d_p) and (d_p <= A_dp_max)):
        if ((A1_phi_min < phi) and (phi <= A1_phi_max)):
            if ((A_fs_min < f_s) and (f_s <= A1_fs_max)):
                swb_class = 'A1'   
        elif ((A2_phi_min < phi) and (phi <= A2_phi_max)):
            if ((A_fs_min < f_s) and (f_s <= A2_fs_max)):
                swb_class = 'A2' 
        elif ((A3_phi_min < phi) and (phi <= A3_phi_max)):
            #-------------------------------------------
            # We can either introduce a class B3, or
            # we can raise A3_fs_max from 0.45 to 1.0
            #-------------------------------------------
            if ((A_fs_min < f_s) and (f_s <= A3_fs_max)):                
                swb_class = 'A3'

    #--------------------------------
    # Check "B" classes for a match
    # These are snow-dominated.
    #--------------------------------
    if ((B_dp_min < d_p) and (d_p <= B_dp_max)):
        if ((B_fs_min < f_s) and (f_s <= B_fs_max)):
           if ((B1_phi_min < phi) and (phi <= B1_phi_max)):
               swb_class = 'B1'
           elif ((B2_phi_min < phi) and (phi <= B2_phi_max)):
               swb_class = 'B2'
           elif ((B3_phi_min < phi) and (phi <= B3_phi_max)):
               #-------------------------------------------
               # B3 is not one of the original 10 classes
               # but it is analogous to A3.
               #-------------------------------------------
               if (USE_B3):
                   swb_class = 'B3'
                
    #--------------------------------
    # Check "C" classes for a match
    # "Precipitation in phase"
    #--------------------------------
    if ((C_dp_min < d_p) and (d_p <= C_dp_max)):
        if ((C_fs_min <= f_s) and (f_s <= C_fs_max)):
            if ((C1_phi_min < phi) and (phi <= C1_phi_max)):
                swb_class = 'C1'
            elif ((C2_phi_min < phi) and (phi <= C2_phi_max)):
                swb_class = 'C2'

    #------------------------------------------------
    # Check "D" classes for a match
    # "Mild seasonality and humid"
    # Higher number indicates higher snow fraction.
    #------------------------------------------------
    if (f_s == D1_fs_val):
        if ((D1_dp_min < d_p) and (d_p <= D1_dp_max)):
            if ((D1_phi_min < phi) and (phi <= D1_phi_max)):
                swb_class = 'D1'
    #--------------------------------------------------
    elif ((D2_fs_min < f_s) and (f_s <= D2_fs_max)):
        if ((D2_dp_min < d_p) and (d_p <= D2_dp_max)):
            if ((D2_phi_min < phi) and (phi <= D2_phi_max)):
                swb_class = 'D2'
    #--------------------------------------------------
    elif ((D3_fs_min < f_s) and (f_s <= D3_fs_max)):
        if ((D3_dp_min < d_p) and (d_p <= D3_dp_max)):
            if ((D3_phi_min < phi) and (phi <= D3_phi_max)):
                swb_class = 'D3'

#     #--------------------------------------------------                
#     if ((D1_dp_min < d_p) and (d_p <= D1_dp_max)):
#         if (f_s <= 0):
#             if ((D1_phi_min < phi) and (phi <= D1_phi_max)):
#                 swb_class = 'D1'
#     #--------------------------------------------------
#     if ((D2_dp_min < d_p) and (d_p <= D2_dp_max)):
#         if ((0 < f_s) and (f_s <= 0.2)):
#             if ((D2_phi_min < phi) and (phi <= D1_phi_max)):
#                 swb_class = 'D2'
#     #--------------------------------------------------
#     if ((D3_dp_min < d_p) and (d_p <= D3_dp_max)):
#         if ((0.2 < f_s) and (f_s <= 0.45)):
#             if ((D3_phi_min < phi) and (phi <= D1_phi_max)):
#                 swb_class = 'D3'
   
    if (swb_class == 'None'):
        if (REPORT or (ORIGINAL == False)):
            print('### SORRY, No matching SWB class for:')
            print('    d_p =', d_p)
            print('    f_s =', f_s)
            print('    phi =', phi)
            print()
    elif (REPORT):
        print('SWB class =', swb_class)
        print()

    return swb_class

#   get_swb_class() 
#---------------------------------------------------------------------
def get_seasonal_water_balance_class( delta_p, f_s, phi):

    return get_swb_class( delta_p, f_s, phi)
  
#   get_seasonal_water_balance_class()
#---------------------------------------------------------------------

