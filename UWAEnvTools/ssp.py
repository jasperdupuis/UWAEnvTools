import numpy as np
import ast
import matplotlib.pyplot as plt
import pandas as pd
import os

class SSP():
    """
    All SSP sources must map to this class.
    """
    
    class SSP_Error(Exception):
        pass
    
    def __init__(self):
        self.depths = r'not set'
        self.dict = dict()
        self.lat = r'not set'
        self.lon = r'not set'
        self.source = r'not set'
        
    def read_profile(self,fname):
        pass
    
    def get_summer(self):
        return self.dict['Summer']
    
    def get_winter(self):
        return self.dict['Winter']


    
class SSP_Blouin_2015(SSP):
    """
    Takes as input Stephane Blouin's Bedford Basin coefficients.
    
    Depth-dependent relation only. It can take an arbitrary depth basis function.
    
    The file in which they are has already been edited for this script.    
    
    Third order approximation only.
    """

    def set_depths(self,p_depths):
        self.depths = p_depths
    
    def read_profile(self,
                     fname,
                     summer = r'August',
                     winter = r'February'):
        if self.depths == r'not set':
            raise SSP.SSP_Error("Depths must be set before calling SSP_Blouin.2015.read_profile.")
        self.source = 'Blouin 2015'
        self.lat = 44.693611
        self.lon = -63.640278
        
        with open(fname) as f:
            data = f.readlines()
        
        coeffs_dict = dict()
        for line in data:
            strs = line.split('[')
            coefs = ast.literal_eval('['+strs[1])
            strs[0].split(' ')[0]
            coeffs_dict[strs[0].split(' ')[0]] = coefs

        summer_c = np.zeros(len(self.depths))
        winter_c = np.zeros(len(self.depths))
        for index in range(len(self.depths)):
            summer_c[index] = self.third_order_estimate(
                coeffs_dict[summer],
                self.depths[index])
            winter_c[index] = self.third_order_estimate(
                coeffs_dict[winter],
                self.depths[index])
        self.dict['Summer'] = summer_c
        self.dict['Winter'] = winter_c
            

    def third_order_estimate(self,p_coeffs,p_depth):
        c = 0
        for index in range(len(p_coeffs)):
            c = c + (p_coeffs[index] * p_depth**(index))
        return c

class SSP_Munk(SSP):
    """
    Implement the standard Munk Profile
    
    Must provide a depth vector.
    """
    
    def set_depths(self,p_depths):
        self.depths = p_depths
        self.epsilon = 0.00737
        self.z_scale = 1300 #m/s , typical around 1500m depth.
    
    def read_profile(self,fname='none_required'):
        """
        Hard coded Munk profile generation.
        
        Sets both Winter and Summer profiles to be the same.
        """
        c = []
        for d in self.depths:
            c.append(self.munk(d))
        c = np.array(c)
        self.dict['Summer'] = c
        self.dict['Winter'] = c
        
            
    def munk(self,z):
        z_tilde = 2 * (z - self.z_scale)/self.z_scale
        c = 1500 * (1 + self.epsilon*(z_tilde - 1 + np.exp(-z_tilde)))
        return c


class SSP_Isovelocity(SSP):
    """
    Implement an isovelocity profile, 1500m/s
 
    Must provide a depth vector.
    """
    
    def set_depths(self,p_depths):
        self.depths = p_depths
    
    def read_profile(self,fname='none_required'):
        """
        Hard coded isovelocity at 1500.
        
        Sets both Winter and Summer profiles to be the same.
        """
        c = np.ones_like(self.depths)
        c = c * 1500.
        self.dict['Summer'] = c
        self.dict['Winter'] = c
        

class SSP_Measured(SSP):
    """
    A set of measured profiles, or a mean of a profile.
    
    Must provide a depth vector.
    """
    
    
    def tolerant_mean(p_list):
        """
        For averaging together different-length arrays in a list.
        e.g. here the CTD is not uniform depth every cast.
        """
        lens = [len(i) for i in p_list]
        arr = np.ma.empty((np.max(lens),len(p_list)))
        arr.mask = True
        for idx, l in enumerate(p_list):
            arr[:len(l),idx] = l
        return arr.mean(axis = -1), arr.std(axis=-1)
    
    def set_depths(self,p_depths):
        self.depths = p_depths
    
    def generate_mean_from_ssps_and_plot(self,
                      p_source_dir,
                      p_target_dir,
                      datetime_filter='2019'):
        """
        For a given querty string on filenames in a source directory,
        calculate the mean SSP.
        Mean is calculated on the longest array length.
        See tolerant_mean for details.
        """

        dirname_source = r'C:\Users\Jasper\Documents\Repo\pyDal\UWAEnvTools\data\Pat Bay CTDs\\'
        dirname_target = r'C:\Users\Jasper\Documents\Repo\pyDal\UWAEnvTools\data\interim\ssp\\'
        files = os.listdir(dirname_source)
        line_skip = 28
        fig, ax = plt.subplots(1,1)
        v_list = []
        z_list = []
        for f in files:
            if 'mean' in f: continue
            if datetime_filter in f:
                strs = f.split('_')
                t = strs[1] + strs[2].split('.')[0]
                df = pd.read_csv( dirname_source + f , skiprows=line_skip)
                x = df[df.columns[6]]
                y = df[df.columns[1]]
                ax.plot(x,y, label = t)
                z_list.append(np.asarray(y))
                v_list.append(np.asarray(x))
        v, v_error = self.tolerant_mean(v_list)
        z, _ = self.tolerant_mean(z_list)
        ax.plot(v,z, label = 'MEAN')
        ax.invert_yaxis()
        plt.legend()
        
        df_res = pd.DataFrame(data={'Depth (m)':z,'Sound speed(m/s)':v})
        df_res.to_csv(dirname_target + 'mean_SSP_' + datetime_filter + '.csv')

        
    
    
    def read_profile(self,
                     fname='explicit_must_be_passed'):
        """
        Read in a real SSP, where the first column is depth and the second is
        sound speed in a CSV.
        
        Sets both Winter and Summer profiles to be the same.
        """

        ## TODO
        df = pd.read_csv(fname)
        
        depths = df[df.columns[0]]
        summer_c = df[df.columns[1]]
        winter_c = df[df.columns[1]]

        self.set_depths(depths)
        self.dict['Summer'] = summer_c
        self.dict['Winter'] = winter_c