from math import exp
import numpy as np
from .utils.configuration import load_config

config = load_config('input/config.json')

class Cell():
  """Description
    Parameters
    ---------
    self: Cell
    the Cell to initialise

    pveg: float
    discrete value of the vegetation type

    pden: float
    discrete value of the vegetation density

    state: IntEnum
    the state in which the cell is initialised in

    Functionality
    ---------
    it initialises the Cell with given FuelAttributes and state
    

    """
  def __init__(self,pveg,pden,state):
      self.pveg = pveg
      self.pden = pden
      self.state = state

  def getState(self):
    return self.state
  
  def setState(self,state):
     self.state = state

  def getFuelAttributes(self):
     return self.pveg,self.pden
  
  def setFuelAttributes(self,pveg,pden):
     self.pveg= pveg;
     self.pden= pden;
  
  def getAltitude(self):
     return self.altitude
  
  def setAltitude(self, grid, x, y):
    mu_x = grid.width / 2
    mu_y = grid.height / 2
    sigma_x = config['terrain']['sigma_x']
    sigma_y = config['terrain']['sigma_y']
    A = config['terrain']['A']
    altitude = A * exp(-(((x - mu_x) ** 2) / (2 * sigma_x ** 2) + ((y - mu_y) ** 2) / (2 * sigma_y ** 2)))
    self.altitude = altitude
  
  def setMoistureContent(self, grid, x, y):
    if y > grid.width/2:
      C_m = np.random.uniform(10,30)
    else:
      C_m = np.random.uniform(60,100)
    self.C_m = C_m
 
  def getMoistureContent(self):
     return self.C_m
