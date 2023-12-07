from enum import Enum
from math import atan, exp, pi, sqrt, cos
from .cell import Cell
from .grid import State
import sys
import numpy as np

class WindDirection(Enum):
    E = (1,0)
    NE = (sqrt(2)/2,sqrt(2)/2)
    N = (0,1)
    NW = (-sqrt(2)/2,sqrt(2)/2)
    W = (-1,0)
    SW = (-sqrt(2)/2,-sqrt(2)/2)
    S = (0,-1)
    SE = (sqrt(2)/2,-sqrt(2)/2)
    
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

  # Calculation of different P function

  # P moisture
  def calculatePm(self,cell_to):
    alpha = 3.258
    beta = 0.111
    C_m = Cell.getMoistureContent(cell_to)
    P_m = alpha * exp(-beta*C_m)
    # Clip the values to [0, 1]
    P_m = 1 / (1 + exp(-P_m))

    return P_m
  
  # P slope (terrain)
  def calculatePs(self,cell_from, cell_to, position):
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

    P_s = exp(a*theta_s)/2

    # Clip the values to [0, 1]
    P_s = 1 / (1 + exp(-a * theta_s))
    
    return P_s
  
  # P wind
  def calcualtePw(self, V, theta):
    c1 = 0.045
    c2 = 0.131
    ft = exp(V*c2*(cos(theta)-1))
    P_w = exp(c1*V)*ft
  
    if P_w>1:
      P_w = 1
    
    return P_w

  # Extra calculations for the parameters
  def getTheta(dir_w,cell_window_position):
    match cell_window_position:
      case (1,2):
          dir_p = WindDirection.E.value
      case (0,2):
          dir_p = WindDirection.NE.value
      case (0,1):
          dir_p = WindDirection.N.value
      case (0,0):
          dir_p = WindDirection.NW.value
      case (1,0):
          dir_p = WindDirection.W.value
      case (2,0):
          dir_p = WindDirection.SW.value
      case (2,1):
          dir_p = WindDirection.S.value
      case (2,2):
          dir_p = WindDirection.SE.value
    # note case (1,1) will never happen because window[1][1] is always 0 as shown in window() method. 

    theta = np.arccos(np.clip(np.dot(dir_w, dir_p) / (np.linalg.norm(dir_w) * np.linalg.norm(dir_p)), -1.0, 1.0))
    theta = theta * (180 / pi)

    return theta

  def checkNeighbours(self, center_x, center_y):
    x = center_x
    y = center_y

    window = self.window(x,y)
    for i in range(len(window)):
      for j in range(len(window[i])):
          if window[i][j] == 1:
            cell_to = self.grid.grid[x+i-1][y+j-1]
            current_cell = self.grid.grid[x][y]

            if cell_to.state == State.NOT_BURNING:
              
              # if i-1 = 0 then the cells are on the same x axis, if y-1 = 0 they are on the same y axis (so adjacent)
              if x == x+i-1 or y == y+j-1:
                  position_flag = 'ADJ'
              else:
                  # if not they are diagonal
                  position_flag = 'DIAG'

              if self.mode == 'humidity':
                 p_burn = self.calculatePm(cell_to)
              elif self.mode == 'terrain':
                p_burn = self.calculatePs(current_cell,cell_to, position_flag)
              elif self.mode == 'wind':
                dir_w = WindDirection.SE.value
                V = 3
                theta = self.getTheta(dir_w,(i,j))
                p_burn = self.calcualtePw(V,theta)

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

  