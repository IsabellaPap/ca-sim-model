"""Performs a simulation using cellular automata to predict the spread of fire for random conditions.

Author: Isabella Papageorgiou isabellapap14@gmail.com
Created: Wed 14 Jun, 2023
"""
from .grid import Grid
from .simulation import Simulation

height = 120
width = 120
grid = Grid(height,width)
count = 1

grid.PopulateGrid(0)
mode = 'humidity'
ignition_point = (1,1)
sim = Simulation(grid, mode, ignition_point)
steps = 60
sim.startSimulation(steps)
