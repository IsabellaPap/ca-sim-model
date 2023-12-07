from math import atan, exp, pi, sqrt
from .cell import Cell
from .grid import Grid, State
import sys
import numpy as np


class Simulation():

  def __init__(self,grid,mode,ignition_point):
        self.mode = mode
        self.grid = grid
        self.ignite(ignition_point)

  """Description
  Parameters
  ---------
  grid: Grid
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
  the "window" numpy array with 0s and 1s where 1 indicates a cell to be checked by the ignition algorithm

  """
  def window(self, center_x, center_y):
    x = center_x
    y = center_y
    window = np.ones((3, 3), dtype=int)
    if x == 0:
        window[0, :] = [0, 0, 0]
        if y == 0:
            window[:, 0] = [0, 0, 0]
        elif y == self.grid.width - 1:
            window[:, -1] = [0, 0, 0]
    if x == self.grid.height - 1:
        window[-1, :] = [0, 0, 0]
        if y == 0:
            window[:, 0] = [0, 0, 0]
        elif y == self.grid.width - 1:
            window[:, -1] = [0, 0, 0]
    if y == 0:
        window[:, 0] = [0, 0, 0] 
    elif y == self.grid.width - 1:
        window[:, -1] = [0, 0, 0]
  
    # the center cell is never to be checked (it's the reference cell)
    window[1, 1] = 0

    return window
        
  def ignite(self,ignition_point):
    x,y = ignition_point

    if self.grid.grid[x][y].state == State.NO_FUEL:
      sys.exit("You tried to burn a cell with no fuel")
    else:
      self.grid.grid[x][y].setState(State.BURNING)

  def calculatePburn(self,cell_from, cell_to, position):
    P_h = 1
    p_veg, p_den = Cell.getFuelAttributes(cell_to)

    # c1 = 0.045
    # c2 = 0.131
    # ft = exp(V*c2*(cos(theta)-1))
    # P_w = exp(c1*V)*ft

    a = 0.1

    # cell altitudes
    E1 = Cell.getAltitude(cell_from)
    E2 = Cell.getAltitude(cell_to)
    height_diff = E2-E1

    l = 25.0

    if position == 0:
      theta_s = atan(height_diff / l)
    elif position == 1:
      theta_s = atan(height_diff / (l * sqrt(2)))

    theta_s = theta_s * 180/pi
    alpha = 3.258
    beta = 0.111
    C_m = Cell.getMoistureContent(cell_to)
    P_m = alpha * exp(-beta*C_m)
    # Clip the values to [0, 1]
    P_m = 1 / (1 + exp(-P_m))
    #P_burn = P_h*(1+p_veg)*(1+p_den)*P_s

    P_burn = P_m

    return P_burn
  
  def checkNeighbours(self, center_x, center_y):
    x = center_x
    y = center_y

    window = self.window(x,y)
    for i in range(len(window)):
      for j in range(len(window[i])):
          if window[i][j] == 1:
            if self.grid.grid[x+i-1][y+j-1].state == State.NOT_BURNING:
              # if i-1 = 0 then the cells are on the same x axis, if y-1 = 0 they are on the same y axis (so adjacent)
              if x == x+i-1 or y == y+j-1:
                  position = 0
              else:
                  # if not they are diagonal
                  position = 1
              p_burn = self.calculatePburn(self.grid.grid[x][y],self.grid.grid[x+i-1][y+j-1], position)
              print("burning cell: {},{} with p_burn = {}".format(x,y,p_burn))
              
              random = np.random.choice(2,1,p=[p_burn,1-p_burn])
              if random == 0:
                  self.ignite((x+i-1,y+j-1))
                  
      
  def startSimulation(self, steps):
    for i in range(steps):
      time_step = i
      burning_cells = []
      burned_cells = []
      # check neighbours
      cells_to_check = []
      for x in range(len(self.grid.grid)):
        for y in range(len(self.grid.grid[x])):
            cells_to_check.append((x, y))
            np.random.shuffle(cells_to_check)
      
      for cell in cells_to_check:
        x,y = cell
        if self.grid.grid[x][y].state == State.BURNING:
          self.checkNeighbours(x,y)
          burning_cells.append((x,y))
      
      # enforce rule 2: if a cell is in state 3 - BURN, it will be BURNED down in the next time step
      for cell in burning_cells:
        x,y = cell
        self.grid.grid[x][y].setState(State.BURNED)
        burned_cells.append(cell)

      self.grid.ShowGrid(time_step)

  