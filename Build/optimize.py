# coding=utf-8
# Copyright 2023 Frank Latos AC8P
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#%%---------------------------------------------------

#
# Yagi optimization code used for performance comparisons of NEC5 executables
#
# -->  pymoo single-objective optimization, NEC5 executes 10,100 times
#
# Based on optim6.py
#



import numpy as np
import subprocess
from pathlib import Path
import time
import re
import matplotlib.pyplot as plt
import functools
import itertools
import io
import os  
import pandas as pd
import json
import pickle
import copy
import logging

import necutil

from numba import njit, prange, types

from dataclasses import dataclass

from necutil import in2m, ft2m, nec5_sim, plot_z, plot_complex_z, plot_elevation, nec5_sim_stdio2, nec5_sim_stdio3
from necutil import vswr, plot_vswr, plot_azimuth2, gen_tapered_el2, gen_tapered_el2_sym, gen_tapered_el2_h, find_de_resfreq
from necutil import plot_azimuth3


#%%------------------------------------------------
from pymoo.algorithms.moo.nsga2 import binary_tournament
from pymoo.algorithms.moo.nsga2 import RankAndCrowdingSurvival
# from pymoo.operators.survival.rank_and_crowding import RankAndCrowding
from pymoo.util.misc import has_feasible
from pymoo.algorithms.base.genetic import GeneticAlgorithm
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.selection.tournament import compare, TournamentSelection
from pymoo.termination.default import DefaultMultiObjectiveTermination
from pymoo.util.display.multi import MultiObjectiveOutput
from pymoo.core.repair import NoRepair
from pymoo.core.duplicate import DefaultDuplicateElimination
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.algorithms.soo.nonconvex.ga import comp_by_cv_and_fitness
from pymoo.core.problem import Problem
from pymoo.termination.xtol import DesignSpaceTermination
from pymoo.termination.robust import RobustTermination




#%%------------------------------------------------

# FB_MIN_CONSTRAINT = 20.0
FB_MIN_CONSTRAINT = 21.8
VSWR_MAX_CONSTRAINT = 3.0
BOOM_LENGTH_CONSTRAINT = 13.5 * 12      # in

# Increasing this pushes max vswr within band from approx 1.45 down to 1.28
# VSWR_WEIGHT = 0.1
VSWR_WEIGHT = 0.25


class MyProblem(Problem):
    
# Design variables:
#       0   len of reflector (in)
#       1   ratio of DE / R
#       2   ratio of D1 / DE
#       3   ratio of D2 / D1
#       4   spacing, R to DE (in)
#       5   spacing, DE to D1 (in)
#       6   spacing, D1 to D2 (in)

    def __init__(self, direct_feed=False, **kwargs):
        super().__init__(n_var=7,
                         n_obj=1,
                         n_ieq_constr=3,
                         **kwargs)

        # New design proposals produced in _evaluate()
        self.design_proposals = None

        # Template for element sections (length,diam), half element (in)
        self.secs = [[None,0.5],[18,0.625],[24,0.75]]
        self.inner_sections_len = sum((s[0] for s in self.secs[1:]))

        self.nsegs = 15                 # Number of NEC segments per element
        self.z = 0                      # z position (height) (m)
        self.z0 = 50                    # TL impedance
        xl_max = 200                    # Max inductive matching reactance considered (ohms)
        self.x_m = np.arange(1,xl_max+1,dtype=float)             # Matching reactance values considered, e.g. 1,2,3...200 ohms
        self.direct_feed = direct_feed
        
        # NEC5 design deck template
        f_min = 28.0                # Band of interest
        f_max = 29.0
        f_num = 9                   # Num of freqs evaluated within band
        f_center = np.mean([f_min,f_max])
        self.necpre = 'CE Yagi\n'
        self.necpost = ''.join(
            ['GX 100 010\n',                                            # Reflect across xz plane
             'GE 0 0\n'                                                 # End of geometry
             'EX 4 {tag} {seg} 2 1.0 0.0\n',                            # Excitation
             f'FR 0 {f_num} 0 0 {f_min} {(f_max-f_min)/(f_num-1)}\n',   # Freqs for XQ (feedpoint Z) simulation
             'XQ 0\n',                                                  # Simulate feedpoint Z
             f'FR 0 1 0 0 {f_center } 0\n'                              # Freqs for radiation pattern sim
             ])
        self.necend = ''.join(                  # For optimization,
            ['RP 0 1 2 0000 90 0 0 180\n',      # theta=90 (horizon), phi=0,180 (azimuth)
             'EN\n'])
        self.necdispend = ''.join(              # For display of results,
            ['RP 0 91 181 0000 0 0 1 1\n',      #  full rad pat with 1 deg grid
             'EN\n'])

        # Pre-calculate 1/Xl for range of matching inductive reactances considered
        freqs = np.linspace(f_min, f_max, num=f_num)        # Freqs used in feedpoint z simulation
        # Array of (1 / matching reactance) at each freq; shape = (len(x_m), len(freqs))
        self.Xlmr = np.reciprocal(((1.0j)*self.x_m)[:,None] @ (freqs[None,:] / f_center ))         # L






    # 
    # Create NEC5 card deck from design variables
    # 
    #  Args:    x       array of design variables for a single design instance
    #                   e.g. [el0_len, el1/el0, el2/el1,...spacing01, spacing12,...] (inches)
    #           display True = calculate rad pat on finer grid for display 
    #
    #   Returns:    NEC5 design as a single string
    #
    def make_nec5_design(self, x, display=False):
        n_el = len(x) // 2 + 1              # Number of elements

        # x positions of elements (meters) (DE is at x=0)
        xpos = in2m( np.cumsum(np.hstack([np.zeros(1),x[n_el:]])) - x[n_el] )

        # Element lengths
        ellens = np.cumprod(x[:n_el])

        necstrs = [self.necpre]
        for i in range(n_el):
            self.secs[0][0] = ellens[i] / 2 - self.inner_sections_len    # Tip length
            s, ex = gen_tapered_el2_sym(self.secs, xpos[i], self.z, i+1, self.nsegs, inches=True)
            necstrs.append(s)
            if i == 1:
                exseg = ex
        necstrs.append(self.necpost.format(tag=2, seg=exseg))
        if display:
            necstrs.append(self.necdispend)
        else:
            necstrs.append(self.necend)
        necstr = ''.join(necstrs)
        return necstr


    # 
    # Find optimum matching reactance and resulting max VSWR within band of interest
    # 
    #  Args:    x        complex z for each frequency in simulation (row vector)
    #
    #  Returns:     [matching inductive reactance, max vswr in band]  (row vector)
    #
    def find_opt_match(self, x):    
        rzs = np.reciprocal(x)                          # 1/z
        zmatch = np.reciprocal(self.Xlmr + rzs)         # 1 / ( 1/z + 1/Xl )
        abs_refl_coef = np.abs((zmatch-self.z0) / (zmatch+self.z0))     # Reflection coefs
        vswr_curves = (1 + abs_refl_coef) / (1 - abs_refl_coef)         # vswr
        max_vswr = np.max(vswr_curves, axis=1)          # Max vswr for any row (reactance value)
        idx = np.argmin(max_vswr)                       # Index of min vswr value
        # Return 
        return np.array([self.x_m[idx], max_vswr[idx]])

    # 
    # Compute max vswr for a vector of feedpoint impedances
    # 
    #  Args:    x       complex z for each frequency in simulation (row vector)
    #
    #  Returns:         max vswr in band
    #
    def direct_feed_vswr(self, x):    
        abs_refl_coef = np.abs((x-self.z0) / (x+self.z0))     # Reflection coefs
        vswr_curve = (1 + abs_refl_coef) / (1 - abs_refl_coef)         # vswr
        return np.max(vswr_curve)          # Max vswr in band



    # Evaluate the designs in X    
    def _evaluate(self, X, out, *args, **kwargs):

        # Make a list of NEC decks, one per row of X
        # *** np.apply_along_axis might truncate strings!
        # designs = list(np.apply_along_axis(self.make_nec5_design, 1, X))
        designs = [self.make_nec5_design(x) for x in X]

        # Run the simulations
        res = nec5_sim_stdio3(designs, timelimit=10000.0)


        # Note: our 'RP' card produces results at theta=90, phi=0,180, e.g.:
        #   array([[ 90.  ,   0.  ,   6.68],
        #          [ 90.  , 180.  ,   2.82]])

        # Extract forward (phi=0) and reverse (phi=180) gains and return as column vectors
        f_gain = np.array([x[1][0][0][1][0,2] for x in res])[:,None]
        r_gain = np.array([x[1][0][0][1][1,2] for x in res])[:,None]

        # Constraint: (F - R) > FB_MIN_CONSTRAINT
        # Calculate forward (phi=0) minus back (phi=180)    (negative values, since constraints are satisfied if <= 0)
        constr_min_fb = -(f_gain - r_gain - FB_MIN_CONSTRAINT)
                    
        # Constraint:  boom length < BOOM_LENGTH_CONSTRAINT
        # Computes the amount boom length exceeds max limit (or 0 if constraint is satisfied)
        n_el = X.shape[1] // 2 + 1              # Number of elements
        constr_boom_len = np.maximum(np.sum(X[:,n_el:], axis=1)[:,None] - BOOM_LENGTH_CONSTRAINT, 0)


        # Our 'XQ' card produces feedpoint impedances that look like:
        #  res[design#][0][0] = 
                    # [[28.0, (18.156-28.716j)],
                    # [28.05, (18.262-26.878j)],
                    # [28.1, (18.338-25.051j)],
                    # [28.15, (18.383-23.23j)],
                    # [28.2, (18.396-21.412j)],
                    # [28.25, (18.377-19.592j)],
                    # [28.3, (18.325-17.767j)],
                    # [28.35, (18.243-15.933j)],
                    # [28.4, (18.129-14.086j)]]

        # Extracts feedpoint complex z for each design --> complex array of shape (#designs, #freqs)
        zs = np.array([[freq[1] for freq in des[0][0]] for des in res])

        if not self.direct_feed:
            # Compute optimum matching reactance, max vswr within band for each design
            xl_mv = np.apply_along_axis(self.find_opt_match, 1, zs)
            vswr = xl_mv[:,1][:,None]                   # Max vswr for each design
            opt_xl = xl_mv[:,0][:,None]                 # Corresponding matching reactance
        else:
            vswr = np.apply_along_axis(self.direct_feed_vswr, 1, zs)[:,None]
            

        # Constraint: max vswr in band < some limit
        constr_max_vswr = vswr - VSWR_MAX_CONSTRAINT

        # Return the target: gain plus a reward for lower vswr
        #  (negative values, since we're minimizing)
        out["F"] = -f_gain + VSWR_WEIGHT * np.minimum(constr_max_vswr, 0.0)
        # print(constr_min_fb.shape, constr_max_vswr.shape, constr_boom_len.shape)
        # Return the constraints
        out["G"] = [constr_min_fb, constr_max_vswr, constr_boom_len]

        # Return the optimum matching reactances, max vswr, fwd and rev gains
        if not self.direct_feed:
            out["XL"] = opt_xl
        out["VSWR"] = vswr
        out["FGAIN"] = f_gain
        out["RGAIN"] = r_gain



problem = MyProblem(xl=np.array([200.0, 0.85, 0.85, 0.85, 40, 10, 40]),
                    xu=np.array([225.0, 1.0, 1.0, 1.0, 80, 30, 100]),
                    direct_feed=True)





#%%------------------------------------------------

from pymoo.termination import get_termination
from pymoo.termination.default import DefaultSingleObjectiveTermination
termination = get_termination("n_gen", 100)


# Called in algorithm._post_advance() (after 'survival')        
from pymoo.core.callback import Callback
class MyCallback(Callback):
    def _update(self, algorithm):
        # print('****** ', len(algorithm.pop))
        pass

from pymoo.optimize import minimize


algorithm = GA(callback=MyCallback(), pop_size=200, n_offsprings=100)

# Run the optimization
timestart = time.monotonic()
res = minimize(problem,
               algorithm,
               termination,
               save_history=False,
               verbose=True)
exectime = time.monotonic() - timestart
print(f'Exec time: {exectime}')


#%%------------------------------------------------

