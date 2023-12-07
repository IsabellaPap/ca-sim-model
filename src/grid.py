from enum import IntEnum
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from .cell import Cell

class State(IntEnum):
    # state of cells (data model)
    NO_FUEL = 1
    NOT_BURNING = 2
    BURNING = 3
    BURNED = 4

class Grid():
    # setting up the plain
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
    #initialize
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
                pveg = 0
                pden = 0
                state = State.NOT_BURNING

              # initialise Cell
              self.grid[i][j] = Cell(pveg,pden,state)

              Cell.setAltitude(self.grid[i][j],self,x=i,y=j)
              Cell.setMoistureContent(self.grid[i][j],self,x=i,y=j)

      self.ShowGrid()
      
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
    it prints the grid and the dem using matplotlib's pyplot and colors libraries. 
    """

    def ShowGrid(self,iter_num=0):
      grid_rep_sim = np.zeros((self.height,self.width), dtype=int)
      grid_rep_dem = np.zeros((self.height,self.width), dtype=float)
      normalized_grid = np.zeros((self.height,self.width), dtype=float)

      plt.figure(figsize=(10,10))

      #DEM Colormap
      for i in range(len(self.grid)):
        for j in range(len(self.grid[i])):
          grid_rep_dem[i][j] = (Cell.getMoistureContent(self.grid[i][j]))
      plt.imshow(grid_rep_dem,cmap='Greens')

      # SIM Colormap
      # unique_states = set()
      for i in range(len(self.grid)):
          for j in range(len(self.grid[i])):
              grid_rep_sim[i][j] = (self.grid[i][j].state.value)
              # unique_states.add(self.grid[i][j])
      transparent_green = (0, 1, 0, 0)
      colormap = colors.ListedColormap([transparent_green,"red","black"])
      plt.imshow(grid_rep_sim,cmap=colormap,alpha=0.9)

      
      
      plt.title('Grid Visualization (Time Step = {})'.format(iter_num))
      plt.grid(False)
      plt.savefig('./output/{}.png'.format(10 + iter_num))
      plt.close()