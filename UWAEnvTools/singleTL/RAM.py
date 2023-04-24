# -*- coding: utf-8 -*-
"""

Feed the environmental info in to a function here,
and compute a single frequency response over the domain.

Build 1D first, 2D afterwards.

Other modules can then put this / these function(s) in a loop over frequency
or other parameters that can be passed. But this should work on its own too
from __main__ call

Does not handle storing results to file!

"""

import numpy as np
import matplotlib.pyplot as plt

import pydal._variables as _vars
import pydal._directories_and_files as _dirs

#my modules
from UWAEnvTools.environment import Environment_ARL,Environment_PYAT,Environment_RAM
from UWAEnvTools.bathymetry import Bathymetry_CHS_2
from UWAEnvTools.ssp import SSP_Blouin_2015, SSP_Isovelocity
from UWAEnvTools.seabed import SeaBed
from UWAEnvTools.surface import Surface
from UWAEnvTools.locations import Location
from UWAEnvTools.source import Source

def compute_RAM_perimeter_to_hyd(
        p_freq, #integer or float singleton
        p_tx_latlon_tuple,
        p_dir_RAM,
        p_RAM_delta_r,
        p_TX_depth = 1.7,
        p_BASIS_SIZE_depth = 100,
        p_BASIS_SIZE_distance = 100,
        p_depth_offset = 0, # legazy correction for bathy alignment with ssp, 0 should work but keep as parameter.
        p_N_points_lon = 200, #for 1x1 m resolution should be ~80
        p_N_points_lat = 80, # for 1x1 m resolution should be ~200 
        p_location = 'Patricia Bay'):
    """
    
    Build up the north and south hydrophone environments for RAM processing.
    
    This is for a the perimineter of the range box for each north and south.

    Parameters
    ----------
    p_freq : TYPE
        DESCRIPTION.
    #integer or float singleton        p_tx_latlon_tuple : TYPE
        DESCRIPTION.
    p_dir_RAM : TYPE
        DESCRIPTION.
    p_RAM_delta_r : TYPE
        DESCRIPTION.
    p_TX_depth : TYPE, optional
        DESCRIPTION. The default is 1.7.
    p_BASIS_SIZE_depth : TYPE, optional
        DESCRIPTION. The default is 100.
    p_BASIS_SIZE_distance : TYPE, optional
        DESCRIPTION. The default is 100.
    p_depth_offset : TYPE, optional
        DESCRIPTION. The default is 0.
    # legazy correction for bathy alignment with ssp : TYPE
        DESCRIPTION.
    0 should work but keep as parameter.        p_N_points_lon : TYPE, optional
        DESCRIPTION. The default is 200.
    #for 1x1 m resolution should be ~80        p_N_points_lat : TYPE, optional
        DESCRIPTION. The default is 80.
    # for 1x1 m resolution should be ~200        p_location : TYPE, optional
        DESCRIPTION. The default is 'Patricia Bay'.

    Returns
    -------
    env_RAM_S : TYPE
        DESCRIPTION.
    env_RAM_N : TYPE
        DESCRIPTION.

    """

    the_location = Location(p_location) 

    # RECEIVER LOCATIONS (depths can be retrieved later)
    rx_loc_N = (the_location.hyd_1_lat,the_location.hyd_1_lon)
    rx_z_N   = the_location.hyd_1_z
    rx_loc_S = (the_location.hyd_2_lat,the_location.hyd_2_lon)
    rx_z_S   = the_location.hyd_2_z

    bathy = Bathymetry_CHS_2()
    bathy.get_2d_bathymetry_trimmed( #og variable values
        p_location_as_object = the_location,
        p_num_points_lon = p_N_points_lon,
        p_num_points_lat = p_N_points_lat,
        p_depth_offset = 0
        )

    ssp = SSP_Blouin_2015()

    # # DEFINE THE SOURCES FOR PYRAM
    source_RAM_S = Source()
    source_RAM_S.set_name()
    source_RAM_S.set_depth(p_TX_depth) #Default is 1.7m
    source_RAM_S.set_speed()
    source_RAM_S.generate_course_perim_pat_bay(  
        the_location.LAT_RANGE_CORRIDOR_TUPLE,
        the_location.LON_RANGE_CORRIDOR_TUPLE,
        p_num_lat_points = p_N_points_lat,
        p_num_lon_points = p_N_points_lon,
        p_hyd_side = 'South'
          )

    source_RAM_N = Source()
    source_RAM_N.set_name()
    source_RAM_N.set_depth(p_TX_depth) #Default is 1.7m
    source_RAM_N.set_speed()
    source_RAM_N.generate_course_perim_pat_bay(  
        the_location.LAT_RANGE_CORRIDOR_TUPLE,
        the_location.LON_RANGE_CORRIDOR_TUPLE,
        p_num_lat_points = p_N_points_lat,
        p_num_lon_points = p_N_points_lon,
        p_hyd_side = 'North'
          )


    # # DEFINE THE ENVIRONMENT FOR SOUTH HYD
    env_RAM_S = Environment_RAM(
        the_location,p_N_points_lat,p_N_points_lon
        )  
    env_RAM_S.set_model_save_directory(p_dir_RAM)
    env_RAM_S.set_bathymetry_common(bathy)
    env_RAM_S.set_seabed_common()
    env_RAM_S.set_ssp_common(ssp)
    env_RAM_S.set_calc_params(p_RAM_delta_r)
    env_RAM_S.set_freqs_common([p_freq])
    env_RAM_S.set_rx_location_common(rx_loc_S)
    env_RAM_S.set_rx_depth_common(rx_z_S)

    # # DEFINE THE ENVIRONMENT FOR NORTH HYD    
    env_RAM_N = Environment_RAM(
        the_location,p_N_points_lat,p_N_points_lon
        )  
    env_RAM_N.set_model_save_directory(p_dir_RAM)
    env_RAM_N.set_bathymetry_common(bathy)
    env_RAM_N.set_seabed_common()
    env_RAM_N.set_ssp_common(ssp)
    env_RAM_N.set_calc_params(p_RAM_delta_r)
    env_RAM_N.set_freqs_common([p_freq])
    env_RAM_N.set_rx_location_common(rx_loc_N)
    env_RAM_N.set_rx_depth_common(rx_z_N)

    # SET SOURCE POINTS
    # THEN COMPUTE RESULTS
    env_RAM_S.set_source_common(source_RAM_S)
    env_RAM_N.set_source_common(source_RAM_N)

    env_RAM_S.calculate_exact_TLs(
        BASIS_SIZE_DEPTH = p_BASIS_SIZE_depth, 
        BASIS_SIZE_DISTANCE = p_BASIS_SIZE_distance
        )
    env_RAM_N.calculate_exact_TLs(
        BASIS_SIZE_DEPTH = p_BASIS_SIZE_depth, 
        BASIS_SIZE_DISTANCE = p_BASIS_SIZE_distance
        )

    return env_RAM_S, env_RAM_N


def compute_RAM_single_latlon_to_hyd(
        p_freq, #integer or float singleton
        p_tx_latlon_tuple,
        p_dir_RAM,
        p_RAM_delta_r,
        p_TX_depth = 1.7,
        p_BASIS_SIZE_depth = 100,
        p_BASIS_SIZE_distance = 100,
        p_depth_offset = 0, # legazy correction for bathy alignment with ssp, 0 should work but keep as parameter.
        p_N_points_lon = 200, #for 1x1 m resolution should be ~80
        p_N_points_lat = 80, # for 1x1 m resolution should be ~200 
        p_location = 'Patricia Bay'):
    """
    
    Build up the north and south hydrophone environments for RAM processing.
    
    This is for a single tx_latlon point.

    Parameters
    ----------
    p_freq : TYPE
        DESCRIPTION.
    #integer or float singleton        p_tx_latlon_tuple : TYPE
        DESCRIPTION.
    p_dir_RAM : TYPE
        DESCRIPTION.
    p_RAM_delta_r : TYPE
        DESCRIPTION.
    p_TX_depth : TYPE, optional
        DESCRIPTION. The default is 1.7.
    p_BASIS_SIZE_depth : TYPE, optional
        DESCRIPTION. The default is 100.
    p_BASIS_SIZE_distance : TYPE, optional
        DESCRIPTION. The default is 100.
    p_depth_offset : TYPE, optional
        DESCRIPTION. The default is 0.
    # legazy correction for bathy alignment with ssp : TYPE
        DESCRIPTION.
    0 should work but keep as parameter.        p_N_points_lon : TYPE, optional
        DESCRIPTION. The default is 200.
    #for 1x1 m resolution should be ~80        p_N_points_lat : TYPE, optional
        DESCRIPTION. The default is 80.
    # for 1x1 m resolution should be ~200        p_location : TYPE, optional
        DESCRIPTION. The default is 'Patricia Bay'.

    Returns
    -------
    env_RAM_S : TYPE
        DESCRIPTION.
    env_RAM_N : TYPE
        DESCRIPTION.

    """

    the_location = Location(p_location) 

    # RECEIVER LOCATIONS (depths can be retrieved later)
    rx_loc_N = (the_location.hyd_1_lat,the_location.hyd_1_lon)
    rx_z_N   = the_location.hyd_1_z
    rx_loc_S = (the_location.hyd_2_lat,the_location.hyd_2_lon)
    rx_z_S   = the_location.hyd_2_z

    bathy = Bathymetry_CHS_2()
    bathy.get_2d_bathymetry_trimmed( #og variable values
        p_location_as_object = the_location,
        p_num_points_lon = p_N_points_lon,
        p_num_points_lat = p_N_points_lat,
        p_depth_offset = 0
        )

    ssp = SSP_Blouin_2015()

    # # DEFINE THE SOURCE FOR PYRAM
    source_RAM = Source()
    source_RAM.set_name()
    source_RAM.set_depth(p_TX_depth) #Default is 1.7m
    source_RAM.set_speed()
    source_RAM.set_singleton_course(p_tx_latlon_tuple)

    # # DEFINE THE ENVIRONMENT FOR SOUTH HYD
    env_RAM_S = Environment_RAM(
        the_location,p_N_points_lat,p_N_points_lon
        )  
    env_RAM_S.set_model_save_directory(p_dir_RAM)
    env_RAM_S.set_bathymetry_common(bathy)
    env_RAM_S.set_seabed_common()
    env_RAM_S.set_ssp_common(ssp)
    env_RAM_S.set_calc_params(p_RAM_delta_r)
    env_RAM_S.set_freqs_common([p_freq])
    env_RAM_S.set_rx_location_common(rx_loc_S)
    env_RAM_S.set_rx_depth_common(rx_z_S)

    # # DEFINE THE ENVIRONMENT FOR NORTH HYD    
    env_RAM_N = Environment_RAM(
        the_location,p_N_points_lat,p_N_points_lon
        )  
    env_RAM_N.set_model_save_directory(p_dir_RAM)
    env_RAM_N.set_bathymetry_common(bathy)
    env_RAM_N.set_seabed_common()
    env_RAM_N.set_ssp_common(ssp)
    env_RAM_N.set_calc_params(p_RAM_delta_r)
    env_RAM_N.set_freqs_common([p_freq])
    env_RAM_N.set_rx_location_common(rx_loc_N)
    env_RAM_N.set_rx_depth_common(rx_z_N)

    # SET SOURCE POINTS
    # THEN COMPUTE RESULT
    env_RAM_S.set_source_common(source_RAM)
    env_RAM_N.set_source_common(source_RAM)

    env_RAM_S.calculate_exact_TLs(
        BASIS_SIZE_DEPTH = p_BASIS_SIZE_depth, 
        BASIS_SIZE_DISTANCE = p_BASIS_SIZE_distance
        )
    env_RAM_N.calculate_exact_TLs(
        BASIS_SIZE_DEPTH = p_BASIS_SIZE_depth, 
        BASIS_SIZE_DISTANCE = p_BASIS_SIZE_distance
        )

    return env_RAM_S, env_RAM_N


    # env_RAM_S.interpolate_RAM_data(100, 100)
    # env_RAM_S.plot_TL_interpolation_with_couse(p_r=100,
    #                                            p_x_track=np.linspace(-100,100,100),
    #                                            p_y_track = (-1/np.sqrt(3)) * np.linspace(-100,100,100))
