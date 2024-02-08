'''
Geometry of the bulk.

Part of GRIN module

Author: Devey Singh Rathore
Created: 8 Feb 2024, 11:13 am
Updated:
'''

import numpy as np
import matplotlib.pyplot as plt

class Geometry(object):
    '''
    Needs to have following features:
    1) Shape function
    2) Volume 
    3) Surface area
    '''
    def __init__(self) -> None:
        pass

    def GetVolume(self) -> float:
        pass
    
    def GetSurfaceArea(self) -> float:
        pass
    
    def MakeMesh(self) -> None:
        pass
    
    def Draw(self):
        pass


class Cuboid(Geometry):
    '''
    Cuboid geometry
    '''
    def __init__(self, height, width, length) -> None:
        super().__init__()
        self.height = height
        self.width = width
        self.length = length
        self.volume = 0
        self.surface_area = 0

    def GetVolume(self) -> float:
        '''
        Gets Volume of the geometry
        '''
        self.volume = self.width*self.height*self.length
        return self.volume
    
    def GetSurfaceArea(self) -> float:
        '''
        Gets surface area of the geometry
        '''
        self.surface_area = 2*(self.height*self.width + self.height*self.length + self.width*self.length)
        return self.surface_area
    

    def Draw(self):
        '''
        A 3D plot showing the geometry.
        '''
        pass