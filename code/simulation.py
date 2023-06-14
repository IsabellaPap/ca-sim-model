"""Performs a simulation using cellular automata to predict the spread of fire for random conditions.

Author: Isabella Papageorgiou isabellapap14@gmail.com
Created: Wed 14 Jun, 2023
"""

# imports
from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os

# set path of project for input/output management
output_path = os.path.abspath('ca-sim-model')
# classes

# discrete states of cells
class State(Enum):
    NO_FUEL = 1
    NO_BURN = 2
    BURN = 3
    BURNED = 4

class Grid():
    
    def __init__(self,height,width):
        self.height = height
        self.width = width
        self.grid = np.zeros((self.height,self.width), dtype=int)

    # populate the Grid with fuel and no fuel cells
    def PopulateGrid(self,p_nofuel):
      # probability to have a cell with no fuel
      self.p_nofuel = p_nofuel  

      for i in range(len(self.grid)):
          for j in range(len(self.grid[i])):
              # Change how the grid is populated
              self.grid[i][j] = np.random.choice(3,1,p=[0,p_nofuel,1-p_nofuel])

    # visualise the Grid

    def ShowGrid(self):
      colormap = colors.ListedColormap(["grey","green"])
      plt.figure(figsize=(10,10))
      plt.imshow(self.grid,cmap=colormap)
      plt.title('Grid Visualization (Time Step = {})')
      plt.grid(False)
      plt.savefig('./output/initial_state.png')
      plt.show()

# modules


# run when file is directly executed

if __name__ == "__main__":
  height = 12
  width = 10
  grid = Grid(height,width)

  grid.PopulateGrid(0.2)
  grid.ShowGrid()