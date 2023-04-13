# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 14:45:39 2021

@author: Jasper
"""

import numpy as np

import haversine

class Source():
    """
    A source can be a ship, point model, or something towed.
    
    Setting aside the source's spectrum, it will have
    course, in lat-lon tuples.
    depth, in m
    speed, in m/s for each tuple above
    name
    """
    
    def __init__(self):
        self.course = 'not set'
        self.depth = 'not set'
        self.speed = 'not set' #m/s
        self.name = 'point source'
        
    def set_depth(self,p_z = 1.7):
        self.depth = p_z
        
    def set_speed(self,p_v = 3):
        self.speed = p_v
        
    def set_name(self,p_name = 'Point source on grid'):
        self.name = p_name
        
    def generate_course_sailed(self,
                        p_CPA_lat_lon,
                        p_CPA_deviation_m,
                        p_CPA_deviation_heading=haversine.Direction.EAST,
                        p_course_heading=0.45*np.pi, #(0,2pi), mathematical not navigation angles
                        p_distance=200,
                        p_divisions = 25):

        new_center = haversine.inverse_haversine(p_CPA_lat_lon, 
                                                 p_CPA_deviation_m, 
                                                 p_CPA_deviation_heading,
                                                 unit=haversine.Unit.METERS)
        # returns tuple (lat, lon)
        course_start = haversine.inverse_haversine(new_center,
                                                   p_distance//2,
                                                   p_course_heading,
                                                   unit=haversine.Unit.METERS)
        course_end = haversine.inverse_haversine(new_center,
                                                 p_distance//2,
                                                 p_course_heading+np.pi,
                                                 unit=haversine.Unit.METERS)
        
        lats = np.linspace(course_start[0],course_end[0],p_divisions)
        lons = np.linspace(course_start[1],course_end[1],p_divisions)
        course = []
        for la, lo, in zip(lats,lons):
            course.append((la,lo))
        #Leave as list of tuples, tuples are input to other things.
        self.course = course
        return course
    
    def generate_course_from_grid(self,
                        p_lat_minmax_tuple,
                        p_lon_minmax_tuple,
                        p_num_lat_points,
                        p_num_lon_points
                        ):
        """
        Basically takes extend information about a lat-lon bounded box
        and turns it in to a list of paired points.
        """
        self.lat_basis = np.linspace(p_lat_minmax_tuple[0],
                                p_lat_minmax_tuple[1],
                                p_num_lat_points)
        self.lon_basis = np.linspace(p_lon_minmax_tuple[0],
                                p_lon_minmax_tuple[1],
                                p_num_lon_points)
        
        lats = []
        lons = []        
        self.course = []
        for la in self.lat_basis:
            for lo in self.lon_basis:
                lats.append(la)
                lons.append(lo)
                self.course.append((la,lo))
        #Leave as list of tuples, tuples are input to other things.
        
        self.lats = np.array(lats)
        self.lons = np.array(lons)        
        return self.course
    
    
    def generate_course_perim_pat_bay(self,
                                        p_lat_minmax_tuple,
                                        p_lon_minmax_tuple,
                                        p_num_lat_points,
                                        p_num_lon_points,
                                        p_hyd_side = 'South'):
        """
        Pass the parameters defining a range corridor,
        but build a set of points that are the perimeter rather than the
        area of the thing.
        Must also supply a cut-off or something.
        
        FOR NOW: Assume p_hyd_side south, then take only the north perimeter
        (don't care about perimeter on the southern bound)
         
        Would need to be extended for each of the four possibilities here.         
        """
        lats = []
        lons = []        
        course = []

        self.lat_basis = np.linspace(p_lat_minmax_tuple[0],
                                p_lat_minmax_tuple[1],
                                p_num_lat_points
                                )
        self.lon_basis = np.linspace(p_lon_minmax_tuple[0],
                                p_lon_minmax_tuple[1],
                                p_num_lon_points)

        # Handle the southern and northern edges separately.
        # North hyd needs south edge, south hyd needs north edge.
        if p_hyd_side == 'South':
            # Handle the northern edge because it is hyd loc dependent.
            for lo in self.lon_basis:                
                la = np.max(self.lat_basis)
                lats.append(la)
                lons.append(lo)
                course.append((la,lo))
            
        if p_hyd_side =='North':
            # Handle the southern edge because it is hyd loc dependent.
            for lo in self.lon_basis:                
                la = np.min(self.lat_basis)
                lats.append(la)
                lons.append(lo)
                course.append((la,lo))
                        
        for la in self.lat_basis:
            # must append "left" and "right" sides.
            if la == np.max(self.lat_basis) or la == np.min(self.lat_basis): 
                continue # this entry is already in the lon loop.
            lats.append(la)
            lons.append(np.min(self.lon_basis))
            course.append((la,np.min(self.lon_basis)))
            lats.append(la)
            lons.append(np.max(self.lon_basis))
            course.append((la,np.max(self.lon_basis)))
        #Leave as list of tuples, tuples are input to other things.        
        self.course = course
        self.lats = np.array(lats)
        self.lons = np.array(lats)        
        return self.course
            
            
            
            
            
            
    
    