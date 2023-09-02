"""Performs a simulation using cellular automata to predict the spread of fire for random conditions.

Author: Isabella Papageorgiou isabellapap14@gmail.com
Created: Wed 14 Jun, 2023
"""

# imports
from enum import IntEnum
from math import cos, exp, sqrt, tan, pi, atan
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os

# set path of project for input/output management
output_path = os.path.abspath('ca-sim-model')
# classes

# discrete states of cells
class State(IntEnum):
    NO_FUEL = 1
    NOT_BURNING = 2
    BURNING = 3
    BURNED = 4

class Grid():
    
    """Description
    Parameters
    ---------
    self: Grid
    the Grid to initialise

    height: int
    the height of the grid in cells

    width: int
    the width of the grid in cells

    Functionality
    ---------
    it initialises the grid with given dimensions

    """
    def __init__(self,height,width):
        self.height = height
        self.width = width
        self.grid = np.zeros((self.height,self.width), dtype=object)

    """Description
    Parameters
    ---------
    self: Grid
    the Grid to populate

    p_nofuel: float
    the probability of cells being unable to be burned => class NoFuel()
    In the future this is set by the user based on real-life evidence, for now it is random.

    Functionality
    ---------
    it populates the grid with cells that are either able or unable to be burned (class Nofuel() || Fuel() )

    """
    def PopulateGrid(self,p_nofuel):
      # probability to have a cell with no fuel
      self.p_nofuel = p_nofuel
      grid_rep = np.zeros((self.height,self.width), dtype=int)
      for i in range(len(self.grid)):
          for j in range(len(self.grid[i])):
              # Change how the grid is populated
              grid_rep[i][j] = np.random.choice(3,1,p=[0,p_nofuel,1-p_nofuel])
              
              if grid_rep[i][j] == 1:
                pveg = 0
                pden = 0
                state = State.NO_FUEL

              else:
                pveg = np.random.choice([1,2,3],1)
                pden = np.random.choice([1,2,3],1) 
                state = State.NOT_BURNING

              # initialise Cell
              self.grid[i][j] = Cell(pveg,pden,state)

              Cell.setAltitude(self.grid[i][j],self,x=i,y=j)

      self.ShowGrid('dem')
    """Description
    Parameters
    ---------
    self: Grid
    the Grid to print

    colormap: string
    sim for simulation colours
    dem for elevation colours

    iter_num: int
    iteration number, useful for printing the time step on the plot 
    and keeping the output file names unique

    Functionality
    ---------
    it prints the grid using matplotlib's pyplot and colors libraries. 
    The colormap is sliced depending on the number of unique states present on the grid at the time it is called.

    """
    def ShowGrid(self,colormap,iter_num=0):
      grid_rep_sim = np.zeros((self.height,self.width), dtype=int)
      grid_rep_dem = np.zeros((self.height,self.width), dtype=float)
      normalized_grid = np.zeros((self.height,self.width), dtype=float)

      plt.figure(figsize=(10,10))

      #DEM Colormap
      for i in range(len(self.grid)):
        for j in range(len(self.grid[i])):
          grid_rep_dem[i][j] = (Cell.getAltitude(self.grid[i][j]))
      min_altitude = np.min(grid_rep_dem)
      max_altitude = np.max(grid_rep_dem)
      normalized_grid = (grid_rep_dem - min_altitude) / (max_altitude - min_altitude)
      plt.imshow(normalized_grid,cmap='Greens')

      # SIM Colormap
      unique_states = set()
      for i in range(len(self.grid)):
          for j in range(len(self.grid[i])):
              grid_rep_sim[i][j] = (self.grid[i][j].state.value)
              unique_states.add(self.grid[i][j])
      transparent_green = (0, 1, 0, 0)
      colormap = colors.ListedColormap([transparent_green,"red","black"])
      plt.imshow(grid_rep_sim,cmap=colormap,alpha=0.9)

      
      
      plt.title('Grid Visualization (Time Step = {})'.format(iter_num))
      plt.grid(False)
      plt.savefig('./output/{}.png'.format(10 + iter_num))
      plt.close()
    
    """Description
    Parameters
    ---------
    self: Grid
    the Grid on which to perform the simulation

    center_x: int
    x coordinate of the center cell

    center_y: int
    y coordinate of the center cell

    Functionality
    ---------
    it creates a window based on the location of the cell within the grid

    Returns
    ---------
    the "window" numpy array with 0s and 1s where only if window[i,j] == 1 the algorithm will run

    """
    def window(self, center_x, center_y):
      x = center_x
      y = center_y
      window = np.ones((3, 3), dtype=int)
      if x == 0:
          window[0, :] = [0, 0, 0]
          if y == 0:
              window[:, 0] = [0, 0, 0]
          elif y == self.width - 1:
              window[:, -1] = [0, 0, 0]
      if x == self.height - 1:
          window[-1, :] = [0, 0, 0]
          if y == 0:
              window[:, 0] = [0, 0, 0]
          elif y == self.width - 1:
              window[:, -1] = [0, 0, 0]
      if y == 0:
          window[:, 0] = [0, 0, 0] 
      elif y == self.width - 1:
          window[:, -1] = [0, 0, 0]
    
      window[1, 1] = 0

      return window
         
    def ignite(self,ignition_x,ignition_y):
      x = ignition_x
      y = ignition_y

      if self.grid[x][y].state == State.NO_FUEL:
        sys.exit("You tried to burn a cell with no fuel")
      else:
        self.grid[x][y].setState(State.BURNING)

    
    def check_neighbours(self, center_x, center_y):
      x = center_x
      y = center_y

      window = self.window(x,y)
      for i in range(len(window)):
        for j in range(len(window[i])):
            if window[i][j] == 1:
              if self.grid[x+i-1][y+j-1].state == State.NOT_BURNING:
                if x == x+i-1:
                   position = 0
                else:
                   position = 1
                p_burn = Grid.calculatePburn(self,self.grid[x][y],self.grid[x+i-1][y+j-1], position)
                if np.random.rand() < p_burn:
                    self.ignite(x+i-1,y+j-1)

                 
        """Description
    Parameters
    ---------
    self: Grid
    the Grid on which to perform the simulation

    Functionality
    ---------
    it changes the cell values of the input Grid (specifically the cell status) 
    based on the Cellular Automata Rules.
    """
        
    def ca_simulation(self, count):
      burn_count = count
      burning_cells = []
      # check neighbours
      cells_to_check = []
      for x in range(len(self.grid)):
        for y in range(len(self.grid[x])):
            cells_to_check.append((x, y))
            np.random.shuffle(cells_to_check)
      
      for cell in cells_to_check:
        x,y = cell
        if self.grid[x][y].state == State.BURNING:
          self.check_neighbours(x,y)
          burning_cells.append((x,y))
      
      # enforce rule 2: if a cell is in state 3 - BURN, it will be BURNED down in the next time step
      for cell in burning_cells:
        x,y = cell
        # if x==self.width/2 and y ==self.height/2:
        #   print(burn_count)
        # if x==self.width and y==self.height:
        #   print(burn_count)
        self.grid[x][y].setState(State.BURNED)

      self.ShowGrid('sim',burn_count)

    def calculatePburn(self, cell_from, cell_to, position):
      P_h = 0.58
      p_veg, p_den = Cell.getFuelAttributes(cell_to)

      # c1 = 0.045
      # c2 = 0.131
      # ft = exp(V*c2*(cos(theta)-1))
      # P_w = exp(c1*V)*ft

      a = 0.9
      # cell altitudes
      E1 = Cell.getAltitude(cell_from)
      E2 = Cell.getAltitude(cell_to)
      l = 25
      if E1 - E2 == 0:
         theta_s = 0
      else:  
        if position == 0:
          theta_s = atan((E1 - E2) / l)
        elif position == 1:
          theta_s = atan((E1 - E2) / (l * sqrt(2)))
      P_s = exp(a*theta_s)

      P_burn = P_h*(1+p_veg)*(1+p_den)*P_s


      return P_burn

       


      

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
    sigma_x = grid.width / 6
    sigma_y = grid.height / 6
    A=0
    altitude = A * exp(-(((x - mu_x) ** 2) / (2 * sigma_x ** 2) + ((y - mu_y) ** 2) / (2 * sigma_y ** 2)))
    self.altitude = altitude

    if x==mu_x and y ==mu_y:
      print(self.altitude)




# modules

# run when file is directly executed

if __name__ == "__main__":
  height = 120
  width = 120
  grid = Grid(height,width)
  count = 1

  grid.PopulateGrid(0)
  grid.ignite(1,1)
  for i in range(50):
    grid.ca_simulation(count)
    count += 1