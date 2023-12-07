"""Performs a simulation using cellular automata to predict the spread of fire for random conditions.

Author: Isabella Papageorgiou isabellapap14@gmail.com
Created: Wed 14 Jun, 2023
"""
from .grid import Grid
from .simulation import Simulation
from .utils.configuration import load_config
import numpy

config = load_config('input/config.json')

seed_value = config['random_seed']
numpy.random.seed(seed=seed_value)
grid_config = config['grid']

grid = Grid(grid_config['height'],grid_config['width'])

sim_config = config['simulation']

mode = sim_config['mode']
grid.PopulateGrid(grid_config['p_no_fuel'], mode)
ignition_point = (sim_config['ignition_point']['x'],sim_config['ignition_point']['y'])
sim = Simulation(grid, mode, ignition_point,config[mode])
steps = sim_config['steps']
sim.startSimulation(steps)
