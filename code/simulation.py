"""Performs a simulation using cellular automata to predict the spread of fire for random conditions.

Author: Isabella Papageorgiou isabellapap14@gmail.com
Created: Wed 14 Jun, 2023
"""

# imports
from enum import IntEnum
from math import cos, exp, sqrt, tan, pi
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
                self.grid[i][j] = noFuel()
                noFuel.setAltitude(self.grid[i][j],self,x=i+1,y=j+1)

              else:
                self.grid[i][j] = Fuel(np.random.choice([-0.3,0,0.4],1),np.random.choice([-0.4,0,0.3],1))
                Fuel.setAltitude(self.grid[i][j],self,x=i+1,y=j+1)

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

      if colormap == 'sim':
        unique_states = set()
        for i in range(len(self.grid)):
            for j in range(len(self.grid[i])):
                grid_rep_sim[i][j] = (self.grid[i][j].state.value)
                unique_states.add(self.grid[i][j])
        num_states = len(unique_states)
        colormap_to_use = colors.ListedColormap(["grey","green","red","black"][:num_states])
        plt.figure(figsize=(10,10))
        plt.imshow(grid_rep_sim,cmap=colormap_to_use)
      elif colormap == 'dem':
         for i in range(len(self.grid)):
            for j in range(len(self.grid[i])):
                if isinstance(self.grid[i][j],noFuel):
                  grid_rep_dem[i][j] = (noFuel.getAltitude(self.grid[i][j]))
                else:
                   grid_rep_dem[i][j] = (Fuel.getAltitude(self.grid[i][j]))
         min_altitude = np.min(grid_rep_dem)
         max_altitude = np.max(grid_rep_dem)
         colormap_to_use = 'Greens'
         normalized_grid = (grid_rep_dem - min_altitude) / (max_altitude - min_altitude)

         plt.figure(figsize=(10,10))
         plt.imshow(normalized_grid,cmap=colormap_to_use)
      
      plt.title('Grid Visualization (Time Step = {})'.format(iter_num))
      plt.grid(False)
      plt.savefig('./output/{}.png'.format(10 + iter_num))
    
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
                if p_burn > 0:
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
        self.grid[x][y].setState(State.BURNED)

      self.ShowGrid('sim',burn_count)

    def calculatePburn(self, cell_from, cell_to, position):
      P_h = 0.58
      p_veg, p_den = Fuel.getFuelAttributes(cell_to)

      # c1 = 0.045
      # c2 = 0.131
      # ft = exp(V*c2*(cos(theta)-1))
      # P_w = exp(c1*V)*ft

      a = 0.078
      # cell altitudes
      E1 = Fuel.getAltitude(cell_from)
      E2 = Fuel.getAltitude(cell_to)
      l = 25
      if E1 - E2 == 0:
         theta_s = 0
      else:  
        if position == 0:
          theta_s = pow(tan(E1-E2/l),-1)
        elif position == 1:
          theta_s = pow(tan(E1-E2/sqrt(2*l)),-1)
      P_s = exp(a*theta_s)

      P_burn = P_h*(1+p_veg)*(1+p_den)*P_s

      return P_burn

       


      

class Fuel():
   
  def __init__(self,pveg,pden):
      self.pveg = pveg
      self.pden = pden
      self.state = State.NOT_BURNING
    
  def setState(self,state):
     self.state = state

  def getFuelAttributes(self):
     return self.pveg,self.pden
  
  def getAltitude(self):
     return self.altitude
  def setAltitude(self,grid,x,y):
     center_x = grid.width / 2
     center_y = grid.height / 2
     distance_from_center = sqrt((x - center_x)**2 + (y - center_y)**2)
     sigma = -10
     self.altitude = (1/(2*pi*sigma**2))*exp(-(distance_from_center**2)/2*sigma**2)


class noFuel():
  def __init__(self):
        self.state = State.NO_FUEL

  def setAltitude(self,grid,x,y):
     center_x = grid.width / 2
     center_y = grid.height / 2
     distance_from_center = sqrt((x - center_x)**2 + (y - center_y)**2)
     sigma = -10
     self.altitude = (1/(2*pi*sigma**2))*exp(-(distance_from_center**2)/2*sigma**2)

  def getAltitude(self):
     return self.altitude

# modules

# run when file is directly executed

if __name__ == "__main__":
  height = 120
  width = 100
  grid = Grid(height,width)
  count = 1

  grid.PopulateGrid(0.02)
  grid.ignite(10,20)
  for i in range(25):
    grid.ca_simulation(count)
    count += 1