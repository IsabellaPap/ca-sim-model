"""Performs a simulation using cellular automata to predict the spread of fire for random conditions.

Author: Isabella Papageorgiou isabellapap14@gmail.com
Created: Wed 14 Jun, 2023
"""

# imports
from enum import IntEnum
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
    NO_BURN = 2
    BURNING = 3
    BURNED = 4

class Grid():
    
    def __init__(self,height,width):
        self.height = height
        self.width = width
        self.grid = np.zeros((self.height,self.width), dtype=object)

    # populate the Grid with fuel and no fuel cells
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
              else:
                self.grid[i][j] = Fuel(np.random.choice([-0.3,0,0.4],1),np.random.choice([-0.4,0,0.3],1))

    # visualise the Grid
    def ShowGrid(self,iter_num=0):
      grid_rep = np.zeros((self.height,self.width), dtype=int)
      for i in range(len(self.grid)):
          for j in range(len(self.grid[i])):
              grid_rep[i][j] = (self.grid[i][j].state.value)

      colormap = colors.ListedColormap(["grey","black","green","red",])
      plt.figure(figsize=(10,10))
      plt.imshow(grid_rep,cmap=colormap)
      plt.title('Grid Visualization (Time Step = {})'.format(i))
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
        print("this cell {},{} is now burning. look! {}".format(x,y,grid.grid[x][y].state))

    
    def check_neighbours(self, center_x, center_y):
      x = center_x
      y = center_y

      window = self.window(x,y)
      for i in range(len(window)):
        for j in range(len(window[i])):
            if window[i][j] == 1:
              if self.grid[x+i-1][y+j-1].state == State.NO_BURN:
                p_burn = np.random.normal(0,1)
                print("burning cell {},{} with probability {}".format(i,j,p_burn))
                if p_burn > 0:
                  self.ignite(x+i-1,y+j-1)

                 
        """Description
    Parameters
    ---------
    self: Grid
    the Grid on which to perform the simulation

    ignition_x: int
    x coordinate of the ignition cell

    ignition_y: int
    y coordinate of the ignition cell
    Functionality
    ---------
    it changes the cell values of the input Grid (specifically the cell status) 
    based on the Cellular Automata Rules.
    """
        
    def ca_simulation(self):
      burn_count = 0

      # check neighbours
      for x in range(len(self.grid)):
        for y in range(len(self.grid[x])):
          if self.grid[x][y].state == State.BURNING:
            self.check_neighbours(x,y)
      
      # enforce rule 2: if a cell is in state 3 - BURN, it will be BURNED down in the next time step
      for i in range(len(self.grid)):
        for j in range(len(self.grid[i])):
          if self.grid[i][j].state == State.BURNING:
            print("cell {},{} has burned down".format(i,j))
            burn_count += 1
            self.grid[i][j].state == State.BURNED

      self.ShowGrid(burn_count)


      

class Fuel():
   
  def __init__(self,pveg,pden):
      self.pveg = pveg
      self.pden = pden
      self.state = State.NO_BURN
    
  def setState(self,state):
     self.state = state

class noFuel():
  def __init__(self):
        self.state = State.NO_FUEL

# modules

# run when file is directly executed

if __name__ == "__main__":
  height = 120
  width = 100
  grid = Grid(height,width)

  grid.PopulateGrid(0.2)
  grid.ignite(35,46)
  for i in range(25):
    grid.ca_simulation()