
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

#
# Utility functions for use with NEC5
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
import pandas as pd
import copy
import logging

import os  
from dataclasses import dataclass

from numba import njit

# Base name for auto-generated files
# Generate 8-digit prefix, verify that no existing files match it
cwd = Path.cwd()
while (True):
    unique_base = str(np.random.randint(1e7,1e8))    
    if list(cwd.glob(unique_base+'*')) == []:
        break
unique_count = 0

logger = logging.getLogger('bayes')

#%%------------------------------------------------

# Match INPUT LINE line (acts as section separator)
#  m[0] entire match, m[1] integer line number, m[2] command, e.g. 'EX', 'FR'
#  m[3]-m[6] integer fields, m[7]-m[12] float fields
input_line_pat = r'^\s*\*\*\*\*\* INPUT LINE\s+([0-9]+)\s+([A-Z]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)'

# Search a file for a line matching a pattern
# Return match object if found, or None if end of section ('INPUT LINE')
#   else raise EOFError
# Args:  f:a file handle,  pat:a regexp
# Returns: (match_obj|None, line_that_produced_match)
def find_string(f, pat):
    for s in f:
        m = re.match(input_line_pat, s)
        if m:
            return None, s
        m = re.match(pat, s)
        if m:
            return m, s
    raise EOFError

# Search a file for one of a list of lines matching a pattern
# Return match object if found, or None if end of section ('INPUT LINE')
#   else raise EOFError
# Args:  f:a file handle,  pats:list of regexps
# Returns: (match_obj|None, line_that_produced_match. index of match in list)
def find_strings(f, pats):
    for s in f:
        m = re.match(input_line_pat, s)
        if m:
            return None, s, 0
        for idx,pat in enumerate(pats):
            m = re.match(pat, s)
            if m:
                return m, s, idx
    raise EOFError


# Skip 'n' lines, while also checking for end of section
# If skipped 'n' lines without encountering INPUT LINE, return (True, last_line_read)
#   if INPUT LINE encountered, (False, last_line_read=INPUT LINE)
#   else raise EOFError
def skip_lines(f, n):
    for count in range(n):
        s = f.readline()
        # print(f'SL {s}')
        if not s:
            raise EOFError
        m = re.match(input_line_pat, s)
        if m:
            return False, s
    return True, s

# Search a file for an INPUT LINE line (acts as section separator)
# Return match object if found, else raise EOFError
def find_input_line(f):
    for s in f:
        m = re.match(input_line_pat, s)
        if m:
            return m, s
    raise EOFError

# Check if next line (only) matches a pattern
# Return match object if found, raise EOFError if EOF, ValueError if no match
def match_string(f, pat):
    s = f.readline()
    # print(f'MS {s}')
    if not s:
        raise EOFError
    m = re.match(pat, s)
    if m:
        return m, s
    raise ValueError(s)



# Process NEC5 output file results
#
# Deck can contain XQ (for computing input Z) and/or RP (radiation pattern) cards
# Single or multiple frequencies can be specified
#
# Args:
#   outf            output file from NEC5 run (as a filename string or a Path object)
#
def nec5_read_output_file(outf):
    with open(outf, 'r', encoding="utf-8") as f:
        res = nec5_read_output(f)
    return res

# Process NEC5 output file results
#
# Deck can contain XQ (for computing input Z) and/or RP (radiation pattern) cards
# Single or multiple frequencies can be specified
#
# Args:
#   st              output from NEC5 as an in-memory string
#
def nec5_read_output_str(st):
    with io.StringIO(st) as f:
        res = nec5_read_output(f)
    return res


# Process NEC5 output file results
#
# Deck can contain XQ (for computing input Z) and/or RP (radiation pattern) cards
# Single or multiple frequencies can be specified
#
# Args:
#   outf            output file from NEC5 run (as a string or a Path object)
#
def nec5_read_output(f):

    XQresults = []
    RPresults = []
    freq = 0.0          # Save last FREQUENCY value
    
    try:

        # Find next 'INPUT LINE' line in file -- 
        #   marks start of output resulting from a command in input deck
        # Loop until it's found, or an EOFError exception is raised
        lastline = None
        while True:
            # Scan for INPUT LINEs in file or in lastline string if it exists
            m = None
            if lastline:
                m = re.match(input_line_pat, lastline)      # Returns match obj or None
                lastline = None
            if not m:                           # If lastline didn't match...
                m, _ = find_input_line(f)

            # Section type is in m[2]: 'XQ', 'RP', etc.
            if m[2] == 'XQ':


                # Read data for (possibly) multiple frequencies
                # Loop terminates when we encounter a new section ('INPUT LINE' line) or reach EOF
                XQSectionResults = []
                while True:
                    # Read the frequency field 'FREQUENCY ='
                    m, lastline = find_string(f, r'^\s+FREQUENCY= ([0-9.E+]+)\s+([A-Z]+)')
                    if not m:       # If end of section, break
                        break
                    if m[2] != 'MHZ':
                        raise ValueError(lastline)
                    freq = float(m[1])

                    m, lastline = find_string(f, r'^\s+- - - ANTENNA INPUT PARAMETERS - - -')
                    if not m:       # If end of section, break
                        break
                    m, lastline = skip_lines(f, 3)
                    if not m:       # If end of section, break
                        break
                    # Read the data for this frequency:
                    #   tag seg end Vr,i Ir,i Zr,i Yr,i pwr
                    m, lastline = match_string(f, r'^\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)')

                    # Segment currents and other stuff found here -- not currently being used

                    # Add freq, Z to results 
                    XQSectionResults = XQSectionResults + [[freq, complex(float(m[8]), float(m[9]))]]

                # Done with this XQ section, so add results to XQresults
                XQresults = XQresults + [XQSectionResults]

                # When we arrive here we have reached a new section (new input command's results)
                # The string 'lastline' contains the new section's 'INPUT LINE' line


            elif m[2] == 'RP':
                # Process the 'RP' radiation pattern output
                ntheta = int(m[4])          # Number of steps for theta, phi
                nphi = int(m[5])
                
                RPSectionResults = []
                # Read data for (possibly) multiple frequencies
                # Loop terminates when we encounter a new section ('INPUT LINE' line) or reach EOF
                while True:
                    # # Read the frequency field 'FREQUENCY ='
                    # m, lastline = find_string(f, r'^\s+FREQUENCY= ([0-9.E+]+)\s+([A-Z]+)')
                    # if not m:       # If end of section, break
                    #     break
                    # if m[2] != 'MHZ':
                    #     raise ValueError(lastline)
                    # freq = float(m[1])

                    # # Check for presence of header stuff
                    # m, lastline = find_string(f, r'^\s+- - - RADIATION PATTERNS - - -')
                    # if not m:       # If end of section, break
                    #     raise ValueError(lastline)
                    # m, lastline = skip_lines(f, 4)
                    # if not m:       # If end of section, break
                    #     raise ValueError(lastline)

                    # *** Careful: FREQUENCY section might be missing, so look for either FREQUENCY or RADIATION PATTERNS
                    m, lastline, patidx = find_strings(f, [r'^\s+FREQUENCY= ([0-9.E+]+)\s+([A-Z]+)', r'^\s+- - - RADIATION PATTERNS - - -'])
                    if not m:       # If end of section, break
                        break
                    if patidx == 0:         # FREQUENCY
                        if m[2] != 'MHZ':
                            raise ValueError(lastline)
                        freq = float(m[1])

                        m, lastline = find_string(f, r'^\s+- - - RADIATION PATTERNS - - -')     # Now wait for RADIATION PATTERN header
                        if not m:       # If end of section, error
                            break

                    # Found RADIATION PATTERNS section
                    m, lastline = skip_lines(f, 4)
                    if not m:       # If end of section, break
                        raise ValueError(lastline)


                    # Read the data table
                    # Rows (e.g. 48) = phi steps x theta steps
                    # Columns:  'THETA DEGREES', 'PHI DEGREES', 'POWER GAINS - TOTAL DBS'
                    rpdat = np.zeros((ntheta * nphi,3))
                    rp_table_row = r'^\s+([0-9.+-]+)\s+([0-9.]+)\s+([0-9.E-]+)\s+([0-9.E-]+)\s+([0-9.E-]+)\s+([0-9.E-]+)\s+([0-9.E-]+)\s+([A-Z]*)\s*([0-9.E+-]+)\s+([0-9.-]+)\s+([0-9.E+-]+)\s+([0-9.-]+)'
                    for r in range(rpdat.shape[0]):
                        m, lastline = match_string(f, rp_table_row)       # May raise EOFError, ValueError
                        rpdat[r,0] = float(m[1])
                        rpdat[r,1] = float(m[2])
                        rpdat[r,2] = float(m[5])


                    # Add freq, RP table to results 
                    RPSectionResults = RPSectionResults + [[freq, rpdat]]

                # Done with this RP section, so add results to XQresults
                RPresults = RPresults + [RPSectionResults]


            # Some other section type we're not currently using....
            else:
                pass




    except EOFError:
        pass
    except ValueError as e:
        # Do something?
        raise

    return XQresults, RPresults




#%%------------------------------------------------
# Run NEC5 simulation(s) in parallel 
#
# Deck can contain XQ (for computing input Z) and/or RP (radiation pattern) cards
# Single or multiple frequencies can be specified
#
# Note: This version does process I/O via on-disk files,
#  which may be desirable for debugging
# See 'nec5_sim_stdio', 'nec5_sim_stdio2' for versions that communicate via
#  process' stdin, stdout instead
#
# Args:
#   designs         list of designs to be simulated, as multi-line strings
#   clean           True = delete input and output files created
#
def nec5_sim(designs, clean=True, timelimit=100.0):
    global unique_count

    # Start parallel NEC5 processes
    for i in range(len(designs)):
        infile = '{fn}_{suffix:d}.dat'.format(fn=unique_base, suffix=unique_count)
        outfile = '{fn}_{suffix:d}.out'.format(fn=unique_base, suffix=unique_count)
        unique_count += 1
        with open(infile, 'w', encoding="utf-8") as f:
            f.write(designs[i])

        p = subprocess.Popen(['../NEC/nec5', infile, outfile], 
                                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # Save a tuple of process obj, design string, in/out filenames
        designs[i] = (p, designs[i], infile, outfile)

    timestart = time.monotonic()
    timelim = timestart + timelimit

    # Wait for procs to complete, or eventually throw TimeoutError
    while not np.all([isinstance(des[0].poll(),int) for des in designs]):
        if time.monotonic() > timelim:
            raise TimeoutError('Child process timed out')

    # Read output files
    outdata = []
    for t in designs:
        XQresults, RPresults = nec5_read_output_file(t[3])
        if clean:
            Path(t[2]).unlink()
            Path(t[3]).unlink()
        outdata = outdata + [[XQresults, RPresults]]

    return outdata

#%%------------------------------------------------


# Try using stdin/stdout instead of temporary files
# Based on ngspice code in ~/Projects/ngspice/test 
# Useful: https://stackoverflow.com/questions/7756609/pass-stdout-as-file-name-for-command-line-util

# Size of stdout buffer? https://stackoverflow.com/questions/10904067/in-c-whats-the-size-of-stdout-buffer
# May be a concern -- have to incrementally decode output for all processes?

#
# Run NEC5 simulation(s) in parallel, using process stdin, stdout for input, output
#
# Args:
#   designs         list of designs to be simulated, as multi-line strings
#
# Returns: list of lists of results, one per design:
#   [ [Popen obj, resultcode, output as bytearray, input_design],... ]
#
# *** see 'nec5_sim_stdio2' below for an improved version
#
def nec5_sim_stdio(designs, timelimit=100.0):

    # Start parallel NEC5 processes
    for i in range(len(designs)):
        p = subprocess.Popen(['../NEC/nec5', '/dev/stdin', '/dev/stdout'], 
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p.stdin.write(bytes(designs[i], 'utf-8'))
        p.stdin.close()

        # Save a tuple of process obj, return code, bytearray for stdout, design string
        designs[i] = [p, None, bytearray(), designs[i]]

    timestart = time.monotonic()
    timelim = timestart + timelimit

    # Wait for procs to complete, or eventually throw TimeoutError
    alldone = False
    while not alldone:   # While some retcodes are still 'None'...
        if time.monotonic() > timelim:
            raise TimeoutError('Child process timed out')
        alldone = True
        for des in designs:
            if des[1] == None:
                des[1] = des[0].poll()      # Before checking streams (so we don't lose anything in buffer)
                if des[1]==None:            # If this proc not done, we're not done with loop
                    alldone = False
                if des[0].stdout.readable():
                    des[2] += bytearray(des[0].stdout.read())
    return designs


#
# Run NEC5 simulation(s) in parallel, using process stdin, stdout for input, output
# Improved version of above - don't create processes all at once
#
# Args:
#   designs         list of designs to be simulated, as multi-line strings
#
# Returns: list of NecProcess objects:
#       index:  position in input list of designs
#       design: input design (str)
#       p:      resultcode (int)
#       outp:   output as bytearray
#

@dataclass
class NecProcess:
    index:      int
    design:     str
    p:          object
    outp:       bytearray       # Apparently initializing this here makes it a class var (?)
    def __lt__(self, other):
        return self.index < other.index

def nec5_sim_stdio2(designs, timelimit=100.0):
    PROCLIMIT = 100                 # Max processes created at a time
    running = []                    # Running processes as NecProcess objs
    complete = []                   # Completed processes as NecProcess objs
    ndesigns = len(designs)         # Number of input designs
    nrunning = 0
    ncomplete = 0
    index = 0                       # Pointer into input design list
    timestart = time.monotonic()
    timelim = timestart + timelimit

    # Loop until all simulations have completed
    while ncomplete < ndesigns:

        # Start some processes if room
        while nrunning < PROCLIMIT and index < ndesigns:
            # p = subprocess.Popen(['../NEC/nec5', '/dev/stdin', '/dev/stdout'], 
            #                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p = subprocess.Popen(['/home/flatos/Projects/VSCodeProjects/necopt/nec5s', '/dev/stdin', '/dev/stdout'], 
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p.stdin.write(bytes(designs[index], 'utf-8'))       # Send input design to process
            p.stdin.close()
            running.append(NecProcess(index, designs[index], p, bytearray()))        # Add to running queue
            index += 1
            nrunning += 1


        # Check running processes for output data or completion
        i = 0
        while i < nrunning:
            result = running[i].p.poll()        # Before checking streams (so we don't lose anything in buffer)
            if running[i].p.stdout.readable():  # If output data available, read and append to buffer
                running[i].outp += bytearray(running[i].p.stdout.read())
            if result != None:                  # If process completed:
                pobj = running[i].p                 # Get the Popen object, save return code,
                running[i].p = pobj.returncode      #  then delete object
                del pobj                            # * this is important to free file handles
                complete.append(running.pop(i))     # Move NecProcess obj to 'complete' list
                nrunning -= 1
                ncomplete += 1
            else:
                i += 1
            
        # Check timeout
        if time.monotonic() > timelim:
            raise TimeoutError('Child process timed out')

    complete.sort()         # Sort NecProcess objs by index number (to correspond to input designs)
    return complete




#
# Same as nec5_sim_stdio2() but allows the MKL_NUM_THREADS env var to be set
# ---> Doesn't seem to be useful -- no discernible difference in performance
#
def nec5_sim_stdio2t(designs, nthreads, timelimit=100.0):
    PROCLIMIT = 100                 # Max processes created at a time
    running = []                    # Running processes as NecProcess objs
    complete = []                   # Completed processes as NecProcess objs
    ndesigns = len(designs)         # Number of input designs
    nrunning = 0
    ncomplete = 0
    index = 0                       # Pointer into input design list
    timestart = time.monotonic()
    timelim = timestart + timelimit

    # Add MKL_NUM_THREADS to environment of subprocesses
    env = os.environ
    env['MKL_NUM_THREADS'] = str(nthreads)

    # Loop until all simulations have completed
    while ncomplete < ndesigns:

        # Start some processes if room
        while nrunning < PROCLIMIT and index < ndesigns:
            # p = subprocess.Popen(['../NEC/nec5', '/dev/stdin', '/dev/stdout'], 
            p = subprocess.Popen(['../NEC/OpenMPI/nec5mp', '/dev/stdin', '/dev/stdout'], 
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                    env=env)
            p.stdin.write(bytes(designs[index], 'utf-8'))       # Send input design to process
            p.stdin.close()
            running.append(NecProcess(index, designs[index], p, bytearray()))        # Add to running queue
            index += 1
            nrunning += 1


        # Check running processes for output data or completion
        i = 0
        while i < nrunning:
            result = running[i].p.poll()        # Before checking streams (so we don't lose anything in buffer)
            if running[i].p.stdout.readable():  # If output data available, read and append to buffer
                running[i].outp += bytearray(running[i].p.stdout.read())
            if result != None:                  # If process completed:
                pobj = running[i].p                 # Get the Popen object, save return code,
                running[i].p = pobj.returncode      #  then delete object
                del pobj                            # * this is important to free file handles
                complete.append(running.pop(i))     # Move NecProcess obj to 'complete' list
                nrunning -= 1
                ncomplete += 1
            else:
                i += 1
            
        # Check timeout
        if time.monotonic() > timelim:
            raise TimeoutError('Child process timed out')

    complete.sort()         # Sort NecProcess objs by index number (to correspond to input designs)
    return complete







#
# Run NEC5 simulation(s) in parallel, using process stdin, stdout for input, output
# Same as nec5_sim_stdio2() but adds parsing of NEC output files 
#
# Args:
#   designs         list of designs to be simulated, as multi-line strings
#   debug           if True, return the list of NecProcess objs (e.g. to access NEC5 output)
#                    
#
# Returns: list of XQ results, RP results for each design:
#       [[XQresults, RPresults]]
#
def nec5_sim_stdio3(designs, timelimit=100.0, debug=False):

    res = nec5_sim_stdio2(designs)      # Returns list of NecProcess objs

    # Process output files
    outdata = []
    for r in res:
        XQresults, RPresults = nec5_read_output_str(r.outp.decode())
        outdata = outdata + [[XQresults, RPresults]]
    if debug:
        return outdata, res
    return outdata


#%%------------------------------------------------

# Display azimuthal radiation pattern at one or more elevations
#
#  arr      ndarray produced by RP command
#  elevs    list of elevations to plot
#  title    title to display
#
def plot_azimuth(rparr, elevs, title='', colors=['black','blue','green','red']):

    # Create the figure
    fig = plt.figure(figsize=(5, 5), facecolor='xkcd:sky blue')
    ax = fig.add_subplot(projection='polar', facecolor="lightgoldenrodyellow")
    ax.tick_params(grid_color="palegoldenrod")
    ax.set_theta_direction('clockwise')
    ax.set_theta_zero_location('N')

    # Create one or more plots
    phi = rparr[rparr[:,0]==rparr[0,0]][:,1]                # 0 ... 350
    phi = np.hstack([phi,phi[0]])                           # Repeat first entry at end to complete loop
    phi_rad = phi * np.pi / 180
    powerdb_max = np.max(rparr[:,2])                        # Max power at any elevation, azimuth

    lines = []                          # Save plots for use with legend
    labels = []
    for elev,color in zip(elevs,colors):
        # Remember: theta in RP array is from pos z axis, not horizon
        arr = rparr[rparr[:,0]==(90-elev)] 
        arr = np.vstack([arr,arr[0]])                       # Repeat first entry at end to complete loop
        powerdb = np.maximum(arr[:,2] - powerdb_max, -40)      # Relative power (clip minimum at -40 dB)        
        line, = ax.plot(phi_rad, powerdb, color=color)
        lines.append(line)
        labels.append(f'{elev}°')


    # dB range = 0 / -40
    ax.grid(True)
    ax.set_rmax(0)
    ax.set_rmin(-40)

    # Radial grid (place your own labels - can't get font resizing to work right)
    ax.set_rgrids([0,-10,-20,-30,-40], labels=['']*5)
    for db,pos in zip([10,20,30,40], [-10,-20,-30,-40]):
        ax.text(np.pi/6,pos, f'-{db}dB', fontsize=7, color='black')

    # Phi (azimuth) labels: adjust color, position
    for label in ax.xaxis.get_ticklabels():
        label.set_color('red')
        label.set_fontsize(6)
        t, r = label.get_position()         # Positions: theta 0-2pi, radius 0 (normal pos) - 1.0 (at origin)
        label.set_position((t,0.05))

    # Legend
    ax.legend(lines, labels, title='Elevation', fontsize=6, title_fontsize=7, loc='upper right')

    ax.set_title(title, va='bottom')
    plt.show()




# Display azimuthal radiation pattern at one or more elevations
#  --> plot_azimuth2() handles multiple input RP arrays
#
#  rparrs   list of ndarrays produced by RP command
#  elevs    list of elevations to plot (of respective arrays in rparrs)
#  title    title to display
#
def plot_azimuth2(rparrs, elevs, tags, title='', legend_title='', colors=['black','blue','green','red']):

    # Create the figure
    fig = plt.figure(figsize=(5, 5), facecolor='xkcd:sky blue')
    ax = fig.add_subplot(projection='polar', facecolor="lightgoldenrodyellow")
    ax.tick_params(grid_color="palegoldenrod")
    ax.set_theta_direction('clockwise')
    ax.set_theta_zero_location('N')

    # Create one or more plots
    powerdb_max = np.max([np.max(rparrs[n][:,2]) for n in range(len(rparrs))])  # Max power at any elevation, azimuth

    lines = []                          # Save plots for use with legend
    for rparr,elev,color in zip(rparrs,elevs,colors):

        phi = rparr[rparr[:,0]==rparr[0,0]][:,1]                # 0 ... 350
        phi = np.hstack([phi,phi[0]])                           # Repeat first entry at end to complete loop
        phi_rad = phi * np.pi / 180

        # Remember: theta in RP array is from pos z axis, not horizon
        arr = rparr[rparr[:,0]==(90-elev)] 
        arr = np.vstack([arr,arr[0]])                       # Repeat first entry at end to complete loop
        powerdb = np.maximum(arr[:,2] - powerdb_max, -40)      # Relative power (clip minimum at -40 dB)        
        line, = ax.plot(phi_rad, powerdb, color=color)
        lines.append(line)


    # dB range = 0 / -40
    ax.grid(True)
    ax.set_rmax(0)
    ax.set_rmin(-40)

    # Radial grid (place your own labels - can't get font resizing to work right)
    ax.set_rgrids([0,-10,-20,-30,-40], labels=['']*5)
    for db,pos in zip([10,20,30,40], [-10,-20,-30,-40]):
        ax.text(np.pi/6,pos, f'-{db}dB', fontsize=7, color='black')

    # Phi (azimuth) labels: adjust color, position
    for label in ax.xaxis.get_ticklabels():
        label.set_color('red')
        label.set_fontsize(6)
        t, r = label.get_position()         # Positions: theta 0-2pi, radius 0 (normal pos) - 1.0 (at origin)
        label.set_position((t,0.05))

    # Legend
    ax.legend(lines, tags, title=legend_title, fontsize=6, title_fontsize=7, loc='upper right')

    ax.set_title(title, va='bottom')
    plt.show()



# Display azimuthal radiation pattern at one or more elevations
#  --> handles multiple input RP arrays
#  
# ** Same as plot_azimuth2(), but adds ability to mirror 'half-hemisphere' data
#  --> if largest azimuth angle is 180, plot will be mirrored to produce full 360deg plot
#
#  rparrs   list of ndarrays produced by RP command
#  elevs    list of elevations to plot (of respective arrays in rparrs)
#  title    title to display
#
def plot_azimuth3(rparrs, elevs, tags, title='', legend_title='', colors=['black','blue','green','red']):

    # Create the figure
    fig = plt.figure(figsize=(5, 5), facecolor='xkcd:sky blue')
    ax = fig.add_subplot(projection='polar', facecolor="lightgoldenrodyellow")
    ax.tick_params(grid_color="palegoldenrod")
    ax.set_theta_direction('clockwise')
    ax.set_theta_zero_location('N')

    # Max power for all plots to be displayed
    max_per_plot = []
    for rparr,elev in zip(rparrs,elevs):
        max_per_plot.append( np.max(rparr[rparr[:,0]==(90-elev)][:,2]) )
    powerdb_max = np.max(max_per_plot)

    lines = []                          # Save plots for use with legend
    for rparr,elev,color in zip(rparrs,elevs,colors):


        # All rows at selected elevation
        arr = rparr[rparr[:,0]==(90-elev)] 
        maxaz = np.max(arr[:,1])
        if maxaz == 180.0:
            marr = np.flipud(arr[:-1].copy())                  # Reversed copy of arr (except for last row)
            marr[:,1] = 360 - marr[:,1]                 # Mirror the azimuths
            arr = np.vstack([arr,marr])                 # Add to end

        elif maxaz != 360:                              # Typically if azimuth = 0 - 359
            arr = np.vstack([arr,arr[0]])               # Repeat first entry at end to complete loop
            
        phi_rad = arr[:,1] * np.pi / 180
        powerdb = np.maximum(arr[:,2] - powerdb_max, -40)      # Relative power (clip minimum at -40 dB)        
        line, = ax.plot(phi_rad, powerdb, color=color)
        lines.append(line)


    # dB range = 0 / -40
    ax.grid(True)
    ax.set_rmax(0)
    ax.set_rmin(-40)

    # Radial grid (place your own labels - can't get font resizing to work right)
    ax.set_rgrids([0,-10,-20,-30,-40], labels=['']*5)
    for db,pos in zip([10,20,30,40], [-10,-20,-30,-40]):
        ax.text(np.pi/6,pos, f'-{db}dB', fontsize=7, color='black')

    # Phi (azimuth) labels: adjust color, position
    for label in ax.xaxis.get_ticklabels():
        label.set_color('red')
        label.set_fontsize(6)
        t, r = label.get_position()         # Positions: theta 0-2pi, radius 0 (normal pos) - 1.0 (at origin)
        label.set_position((t,0.05))

    # Legend
    ax.legend(lines, tags, title=legend_title, fontsize=6, title_fontsize=7, loc='upper right')

    ax.set_title(title, va='bottom')
    plt.show()



# Display thumbnail azimuthal radiation patterns in a grid (one pattern per plot)
#
#  rparrs   list of ndarrays produced by RP command
#  title    list of titles for lots
#  elevs    single elevation value for all plots
#  vh       (rows,cols) of grid
#
# Example:  nzip frequencies, rad pat arrays
#               freqs,arrs = zip(*result[0][1][0])
#               plot_azimuth_thumb(arrs, vh=(1,6),titles=freqs)
#
def plot_azimuth_thumb(rparrs, titles, vh=(1,5), size=2, elev=0, nolabels=False):
    plt.rc('xtick', labelsize=8) 
    plt.rc('ytick', labelsize=8) 
    fig,axs = plt.subplots(*vh, figsize=(size*vh[1],size*vh[0]), subplot_kw={'projection': 'polar'})
    if axs.ndim == 1:
        axs = axs[None,:]
    for row in range(axs.shape[0]):
        for col in range(axs.shape[1]):
            ax = axs[row,col]
            ax.set_rmax(0)
            ax.set_rmin(-40)
            idx = row*vh[1] + col
            if idx < len(rparrs):
                ax.set_theta_direction('clockwise')
                ax.set_theta_zero_location('N')

                arr = rparrs[idx]
                arr = arr[arr[:,0]==(90-elev)]                  # RP for this elevation
                arr = np.vstack([arr,arr[0]])                       # Repeat first entry at end to complete loop
                phi_rad = arr[:,1] * np.pi / 180                    # phi (azimuth) in radians
                powerdb_max = np.max(arr[:,2])                      # Max power at this elevation

                powerdb = np.maximum(arr[:,2] - powerdb_max, -40)      # Relative power (clip minimum at -40 dB)        
                ax.plot(phi_rad, powerdb)
                if nolabels:
                    ax.tick_params(labelleft = False, labelbottom = False)
                ax.set_title(titles[idx], va='bottom', fontsize=8)
    fig.tight_layout()
plt.show()


#%%------------------------------------------------

# Display one or more elevation radiation patterns
#
#  rparrs       list of one or more ndarrays produced by RP command
#  tags         list of display tags corresponding to above
#  maxgain      max ('outer ring') gain to display
#               if not specified, max gain 
#  title        main title
#  legend_title title for legend box
#
def plot_elevation(rparrs, tags, title='', legend_title='', maxgain=None, colors=['black','blue','green','red']):

    # Create the figure
    fig = plt.figure(figsize=(5, 5), facecolor='xkcd:sky blue')
    ax = fig.add_subplot(projection='polar', facecolor="lightgoldenrodyellow")
    ax.tick_params(grid_color="palegoldenrod")
    ax.set_theta_direction('clockwise')
    ax.set_theta_zero_location('N')
    ax.set_thetalim([-np.pi/2-0.001,np.pi/2+0.001])


    # Elevation pattern in phi=0.0 plane (theta) for each input array
    arrs = []
    for arr in rparrs:
        arr0 = arr[arr[:,1]==0]             # 0 deg direction
        arr180 = arr[arr[:,1]==180]         # 180 deg direction
        arr180[:,0] = -arr180[:,0]          # Convert 0...90 to 0...-90
        arr = np.vstack([arr0,arr180[1:]])  # combine (but omit redundant 0 deg entry)
        arr = arr[np.argsort(arr[:,0])]     # Sort -90...90
        arrs.append(arr)

    # Max power over all plots
    theta = arr[:,0]               # -90 ... 90
    theta_rad = theta * np.pi / 180
    if maxgain == None:
        powerdb_max = np.max([np.max(a[:,2]) for a in arrs])
    else:
        powerdb_max = maxgain

    # Draw the plots
    lines = []                          # Save plots for use with legend
    for arr,color in zip(arrs,colors):
        powerdb = np.minimum(np.maximum(arr[:,2] - powerdb_max, -40), powerdb_max)      # Relative power (clip minimum at -40 dB)
        line, = ax.plot(theta_rad, powerdb, color=color)
        lines.append(line)

    # Radial grid
    ax.set_rgrids([0,-10,-20,-30,-40], labels=['']*5)
    for db,pos in zip([10,20,30,40], [-10,-20,-30,-40]):
        ax.text(np.pi/6,pos, f'-{db}dB', fontsize=7, color='black')
    ax.grid(True)
    ax.set_rmax(0)
    ax.set_rmin(-40)

    # Elevation grid and title
    ax.set_title(title, va='bottom')
    nsects = 18
    ax.set_thetagrids(np.arange(nsects+1)*180.0/nsects - 90, 
                    [f'{d+90}°' if d<= 0 else f'{90-d}°' for d in range(-90,91,180//nsects)])
    for label in ax.xaxis.get_ticklabels():
        label.set_color('red')
        label.set_fontsize(6)
        t, r = label.get_position()
        label.set_position((t,0.05))

    # Legend
    ax.legend(lines, tags, title=legend_title, fontsize=6, title_fontsize=7, loc='upper right')

    plt.show()



#%%------------------------------------------------

# Display one or more impedance-vs-frequency plots (as magnitude)
#
#  zlists       list of lists of [freq, cplx_z] pairs
#  tags         list of display tags corresponding to above
#  title        main title
#  legend_title title for legend box
#
def plot_z(zlists, tags, title='', legend_title='', colors=['black','blue','green','red']):
    fig, ax = plt.subplots(figsize=(5, 2.7))
    lines = []
    for zlist,color in zip(zlists, colors):
        fs,zs = list(zip(*zlist))                   # Tuples of freqs, complex zs
        zs = np.abs(zs)                             # Mag of complex z
        line, = ax.plot(fs, zs, color=color)
        lines.append(line)
    ax.set_xlabel('Freq MHz')
    ax.set_ylabel(r'Mag Z $\Omega$')
    ax.set_title(title)
    ax.grid(True)
    ax.legend(lines, tags, title=legend_title)


# Display one or more impedance-vs-frequency plots (as real, complex parts)
#
#  zlists       list of lists of [freq, cplx_z] pairs
#  tags         list of display tags corresponding to above
#  title        main title
#  legend_title title for legend box
#
def plot_complex_z(zlists, tags, title='', legend_title='', colors=['black','blue','green','red']):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 5.4))
    lines = []
    for zlist,color in zip(zlists, colors):
        fs,zs = list(zip(*zlist))                   # Tuples of freqs, complex zs
        line, = ax1.plot(fs, np.real(zs), color=color)
        ax2.plot(fs, np.imag(zs), color=color)
        lines.append(line)
    ax2.set_xlabel('Freq MHz')
    ax1.set_ylabel(r'Re(Z) $\Omega$')
    ax2.set_ylabel(r'Im(Z) $\Omega$')
    ax1.set_title(title)
    ax1.grid(True)
    ax2.grid(True)
    ax1.legend(lines, tags, title=legend_title)



#
# Reflection coefficient for given Zl, Z0
#
def refl_coef(zl, z0=50):
    return (zl-z0) / (zl+z0)

#
# VSWR for given Zl, Z0
#
def vswr(zl, z0=50):
    arc = np.abs(refl_coef(zl, z0))
    return (1 + arc) / (1 - arc)

# Display one or more SWR plots
#
#  zlists       list of lists of [freq, cplx_z] pairs
#  tags         list of display tags corresponding to above
#  title        main title
#  legend_title title for legend box
#
def plot_vswr(zlists, tags, title='', legend_title='', colors=['black','blue','green','red'], z0=50):
    fig, ax = plt.subplots(figsize=(5, 2.7))
    lines = []
    for zlist,color in zip(zlists, colors):
        fs,zs = list(zip(*zlist))                   # Tuples of freqs, complex zs
        zs = list(map(lambda z: vswr(z, z0), zs))
        line, = ax.semilogy(fs, zs, color=color)
        lines.append(line)
    ax.set_xlabel('Freq MHz')
    ax.set_ylabel(r'VSWR')
    ax.set_title(title)
    ax.set_ybound(1.0,20.0)
    yticks = [1.1,1.5,2.0,3.0,5.0,10.0]
    ax.set_yticks(yticks, minor=False)
    ax.set_yticklabels([str(n) for n in yticks])
    # ax.set_yticks([20.0], minor=True)
    ax.yaxis.grid(True, which='major')
    # ax.yaxis.grid(True, which='minor')
    ax.grid(True)
    if len(tags)!=0:
        ax.legend(lines, tags, title=legend_title)

# Display one or more SWR plots
# Version 2: takes VSWR values instead of complex Z's
#
#  zlists       list of lists of [freq, vswr] pairs
#  tags         list of display tags corresponding to above
#  title        main title
#  legend_title title for legend box
#
def plot_vswr_2(zlists, tags, title='', legend_title='', colors=['black','blue','green','red'], z0=50):
    fig, ax = plt.subplots(figsize=(5, 2.7))
    lines = []
    for zlist,color in zip(zlists, colors):
        fs,vs = list(zip(*zlist))                   # Tuples of freqs, complex zs
        line, = ax.semilogy(fs, vs, color=color)
        lines.append(line)
    ax.set_xlabel('Freq MHz')
    ax.set_ylabel(r'VSWR')
    ax.set_title(title)
    ax.set_ybound(1.0,20.0)
    yticks = [1.1,1.5,2.0,3.0,5.0,10.0]
    ax.set_yticks(yticks, minor=False)
    ax.set_yticklabels([str(n) for n in yticks])
    # ax.set_yticks([20.0], minor=True)
    ax.yaxis.grid(True, which='major')
    # ax.yaxis.grid(True, which='minor')
    ax.grid(True)
    if len(tags)!=0:
        ax.legend(lines, tags, title=legend_title)


#%%------------------------------------------------

def in2m(inch):
    return inch / 39.3701

def ft2m(ft):
    return ft / (39.3701 / 12)

#
# Generate NEC cards for a tapered antenna element
# Elements are parallel to the y axis
#
#   seclist     list of (len, dia) pairs (all dimensions in meters)
#   x           x position
#   z           z (height)
#   tag         tag used for all sections of element
#   segs        # of segments within each section
#   inches      True: elem lengths/dias in inches (default = meters)
#
#   Returns:    (multi-line string of NEC5 cards, seg# for center of element)
#    Use far end of returned seg# to drive element at center, e.g.
#    EX 1 <tag> <seg#> 2 ...
#
def gen_tapered_el(seclist, x, z, tag, segs, inches=False):
    if inches:
        # seclist = list(map(lambda t:tuple(map(in2m, t)), seclist))
        seclist = [(in2m(t[0]), in2m(t[1])) for t in seclist]
    totlen = sum(t[0] for t in seclist)        # Sum of section lengths
    y = -totlen / 2.0
    s = ''
    for ln,dia in seclist:
        s = s + f'GW {tag} {segs} {x} {y} {z} {x} {y+ln} {z} {dia/2}\n'
        y = y + ln
    return s, int(len(seclist)*segs) // 2

#
# This version takes total number of segments to use for element (not per-section)
# Segments distributed proportionately among sections,
#  middle section will have an even number of segments
#
def gen_tapered_el2(seclist, x, z, tag, totsegs, inches=False):
    if inches:
        seclist = [(in2m(t[0]), in2m(t[1])) for t in seclist]
    totlen = sum(t[0] for t in seclist)                                 # Sum of section lengths
    nsegs = [np.max((round(t[0]/totlen*totsegs), 1)) for t in seclist]  # List of seg count per section (minimum of 1)
    if len(nsegs)%2 == 1 and nsegs[len(nsegs)//2]%2 == 1:
        nsegs[len(nsegs)//2] += 1
    y = -totlen / 2.0
    s = ''
    for (ln,dia),segs in zip(seclist,nsegs):
        s = s + f'GW {tag} {segs} {x} {y} {z} {x} {y+ln} {z} {dia/2}\n'
        y = y + ln
    return s, sum(nsegs) // 2




#
# Same, but seclist contains dimensions of *half* of element
#   E.g. [ <tip elem>, <next elem>, ... <half of center element>]
#
def gen_tapered_el2_h(seclist, x, z, tag, totsegs, inches=False):
    slst = seclist[:-1] + [(seclist[-1][0]*2, seclist[-1][1])] + seclist[-2::-1]
    return gen_tapered_el2(slst, x, z, tag, totsegs, inches)



#
# Generate NEC cards for a tapered antenna element
# Elements are parallel to the y axis
#
# ** This version accepts a description of half of the element
#  and relies on the deck having a subsequent 'GX nnn 010' card to complete the geometry
#
#   seclist     list of (len, dia) pairs (element from tip to midpoint)
#   x           x position
#   z           z (height)
#   tag         tag used for all sections of element
#   totsegs     segments in _full_ element, distributed proportionately among sections
#   inches      True: elem lengths/dias in inches (default = meters)
#
#   Returns:    (multi-line string of NEC5 cards, seg# for center of element)
#    Use far end of returned seg# to drive element at center, e.g.
#    EX 1 <tag> <seg#> 2 ...
#
def gen_tapered_el2_sym(seclist, x, z, tag, totsegs, inches=False):
    if inches:
        seclist = [(in2m(t[0]), in2m(t[1])) for t in seclist]
    totlen = sum(t[0] for t in seclist) * 2                                 # Full element length
    slst = seclist[:-1] + [(seclist[-1][0]*2, seclist[-1][1])]              # Double the length of the last ('middle') section
    nsegs = [np.max((round(t[0]/totlen*totsegs), 1)) for t in slst]         # List of seg count per section (minimum of 1) 
    if nsegs[-1]%2 == 1:     
        nsegs[-1] += 1
    nsegs[-1] //= 2
    # y = -sum(t[0] for t in seclist)
    y = -totlen / 2.0
    s = ''
    for (ln,dia),segs in zip(seclist,nsegs):
        yend = y + ln
        if np.isclose(yend, 0.0, atol=1e-5):                                   # Fix fp errors if close to plane
            yend = 0.0
        s = s + f'GW {tag} {segs} {x} {y} {z} {x} {yend} {z} {dia/2}\n'
        y = yend
    return s, sum(nsegs)





#%%------------------------------------------------

#
# Convert RP radiation pattern array returned from nec5_sim_stdio3(), etc.:
#       [[theta0, phi0, db], [theta1, phi1, db].....]
#  into numpy array of shape (#_theta_vals, #_phi_vals)
#  Note: theta is angle with z axis (0=zenith, 90=horizon), phi is azimuth
#
# Original non-numba version -- still useful if you need ph, th values returned
def convert_rp_array_0(arr, return_vals=False):
    thvals = np.unique(arr[:,0].astype(int))
    phvals = np.unique(arr[:,1].astype(int))
    rarr = np.reshape(arr[:,2], (len(phvals), len(thvals))).T
    if return_vals:
        return rarr, thvals, phvals
    return rarr


#
# Convert RP radiation pattern array returned from nec5_sim_stdio3(), etc.:
#       [[theta0, phi0, db], [theta1, phi1, db].....]
#  into numpy array of shape (#_theta_vals, #_phi_vals)
#  Note: theta is angle with z axis (0=zenith, 90=horizon), phi is azimuth
#
# New numba version
@njit
def convert_rp_array(arr):
    nthvals = np.unique(arr[:,0]).shape[0]
    nphvals = np.unique(arr[:,1]).shape[0]
    rarr = np.empty((nthvals, nphvals))
    for th in np.arange(nthvals):
        for ph in np.arange(nphvals):
            rarr[th,ph] = arr[ph*nthvals + th, 2]
    return rarr








# Notes on integrating power density over ranges of angles
#
# Area of horizontal ring: width = r*d_theta, circum = 2*pi*r*sin(theta)
#  integrate from th1 to th2:  2pi*r^2*sin(theta)*d_theta  -->  2pi*r^2[cos(th2) - cos(th1)]



#
# Integrate power density produced by NEC5 RP card, over a range of theta, phi values
# Reminder: theta is angle with pos z axis (0-90), phi is azimuth (0-359)
#
# Note: This (original) version uses RP covering entire upper hemisphere
# See 'integrate_power_density_3' for improved half-hemisphere version
# 
# Args:     arr     array of shape (91, 360), power in dB as returned by NEC5
#           [thetal, thetah), [phil, phih)    region of interest (degrees)
#
#           Note: if phil < phih, integrates through 0 degrees
#
# Returns   power summed over specified region (divide by area for mean power density)
#           area of specified region (for a sphere of radius=1)
#
def integrate_power_density(arr, thetal=0, thetah=91, phil=0, phih=360):
    thetas = np.pi / 180.0 * np.arange(91.0)        # Theta vals from 0-90 deg (in radians)
    # Areas of rings centered on thetas from 0-90deg
    weight = 2*np.pi * (np.cos(thetas - 0.5 * np.pi / 180) - np.cos(thetas + 0.5 * np.pi / 180))
    # Strip at horizon: only half is in upper hemisphere
    weight[90] *= 0.5
    # Circle at zenith also needs its area corrected
    weight[0] = 2*np.pi * (1 - np.cos(0.5 * np.pi / 180))
        # Mult by weight (as col vec)
        # pwr = np.sum( (weight[:,None] * np.power(10, arr / 10))[thetal:thetah, phil:phih] ) / 360
        # area = np.sum(weight[thetal:thetah])*(phih-phil) / 360
        # return pwr, area
    # May have to sum over two regions if we cross phi=0
    phi_ranges = [(phil,phih)]
    if phil > phih:
        phi_ranges = [(phil,360), (0,phih)]
    pwr, area = (0,0)
    # Convert dB power array to ratio, mult by weight based on area of ring for theta val
    #  Note: weight[:,None] is a column vec, broadcast along rows
    weighted_pwr = weight[:,None] * np.power(10, arr / 10)      # Mult by weight (as col vec)
    for pl,ph in phi_ranges:
        pwr += np.sum( weighted_pwr[thetal:thetah, pl:ph] ) / 360
        area += np.sum(weight[thetal:thetah])*(ph-pl) / 360
    return pwr, area


#
# Integrate power density produced by NEC5 RP card
#
# Version 2: Improve speed by considering only half of upper hemisphere (theta 0-90, phi 0-180)
#             and accommodating a coarser grid (must be divisor of 90: 2,3,5,6,9,10,15,18,30,45)
#
# Reminder: theta is angle with pos z axis (0-90), phi is azimuth (0-359)
# 
# Args:     arr     array of shape (90/grid+1, 180/grid+1), power in dB as returned by NEC5
#           grid    theta, phi increment, 2,3,5,6,9,10,15,18,30,45 deg
#           [thetal, thetah), [phil, phih)    region of interest (degrees)
#
# Returns   relative power summed over specified region 
#            If upper/lower, left/right symmetry holds, summing entire half-hemisphere shouid return 0.25
#
# def integrate_power_density_2(arr, grid, thetal=0, thetah=91, phil=0, phih=181):
#     thgrid = np.array(range(0,91,grid))
#     phgrid = np.array(range(0,181,grid))
#     assert arr.shape == (len(phgrid), len(thgrid))
#     thetas = np.pi / 180.0 * thgrid             # Theta vals from 0-90 deg (in radians)
#     # Areas of rings centered on thetas from 0-90deg
#     half_grid_rad = grid/2 * np.pi / 180        # Half of grid width in radians
#     # weight = 2*np.pi * (np.cos(thetas - half_grid_rad) - np.cos(thetas + half_grid_rad))
#     weight = (np.cos(thetas - half_grid_rad) - np.cos(thetas + half_grid_rad))
#     # Strip at horizon: only half is in upper hemisphere
#     weight[-1] *= 0.5
#     # Circle at zenith also needs its area corrected
#     weight[0] = 2*np.pi * (1 - np.cos(half_grid_rad))
#     # Sum area and power weighted by area
#     #   Convert dB power array to ratio, mult by weight based on area of ring for theta val
#     #    Note: weight[:,None] is a column vec, broadcast along rows
#     weighted_pwr = weight[:,None] * np.power(10, arr / 10)      # Mult by weight (as col vec)
#     weighted_pwr[:,0] *= 0.5                # Strips at 0, 180 deg azimuth are half-width
#     weighted_pwr[:,-1] *= 0.5
#     # Boolean masks to select desired rows, columns
#     thbool = np.logical_and(thgrid >= thetal, thgrid < thetah )
#     phbool = np.logical_and(phgrid >= phil, phgrid < phih )
#     # Sum over selected area
#     pwr = np.sum(weighted_pwr[phbool][:,thbool])
#     return pwr

#
# Integrate power density produced by NEC5 RP card
#
# Version 2: Improve speed by considering only half of upper hemisphere (theta 0-90, phi 0-180)
#             and accommodating a coarser grid (must be divisor of 90: 2,3,5,6,9,10,15,18,30,45)
#
# Reminder: theta is angle with pos z axis (0-90), phi is azimuth (0-359)
#
# Note: *Deprecated* -- integrate_power_density_3 is an improved version
#
# 
# Args:     arr     array of shape (90/grid+1, 180/grid+1), power in dB as returned by NEC5
#           grid    theta, phi increment, 2,3,5,6,9,10,15,18,30,45 deg
#           [thetal, thetah), [phil, phih)    region of interest (degrees)
#
# Returns   relative power summed over specified region 
#           area of specified region
#
# The area of the quarter-sphere is pi, so the call
#   integrate_power_density_2(arr, grid, thetal=0, thetah=91, phil=0, phih=181)
#       will return (pi, pi).....pi/pi=1.0 will be the relative power density for any pattern
#
def integrate_power_density_2(arr, grid, thetal=0, thetah=91, phil=0, phih=181):
    thgrid = np.array(range(0,91,grid))         # e.g. [0, 10, ..., 90] for grid = 10
    phgrid = np.array(range(0,181,grid))        #      [0, 10, ..., 180]
    assert arr.shape == (len(thgrid), len(phgrid))
    thetas = np.pi / 180.0 * thgrid             # Theta vals from 0-90 deg (in radians)
    # Areas of rings centered on thetas from 0-90deg
    half_grid_rad = grid/2 * np.pi / 180        # Half of grid width in radians
    # Actual areas of patches (grid x grid) degrees in size
    weight = (np.cos(thetas - half_grid_rad) - np.cos(thetas + half_grid_rad))
    # Strip at horizon: only half is in upper hemisphere
    weight[-1] *= 0.5
    # Circle at zenith also needs its area corrected
    weight[0] = 1 - np.cos(half_grid_rad)
    weight = weight * np.pi / (180/grid)
    # Sum area and power weighted by area
    #   Convert dB power array to ratio, mult by weight based on area of ring for theta val
    #    Note: weight[:,None] is a column vec, broadcast along rows
    weighted_pwr = weight[:,None] * np.power(10, arr / 10)      # Mult by weight (as col vec)
    weighted_pwr[:,0] *= 0.5                # Strips at 0, 180 deg azimuth are half-width
    weighted_pwr[:,-1] *= 0.5
    # Boolean masks to select desired rows, columns
    thbool = np.logical_and(thgrid >= thetal, thgrid < thetah )
    phbool = np.logical_and(phgrid >= phil, phgrid < phih )
    # Sum over selected area
    pwr = np.sum(weighted_pwr[thbool][:,phbool])
    # Area of selected region
    nsec = np.sum(phbool)           # Number of sectors (in phi) 
    if phbool[0]:                   # But 0, 180 deg sectors are only 1/2 width
        nsec -= 0.5
    if phbool[-1]:
        nsec -= 0.5
    area = np.sum(weight[thbool]) * nsec
    return pwr, area




#
# Integrate power density produced by NEC5 RP card
#
# Version 3: Improve speed by considering only half of upper hemisphere (theta 0-90, phi 0-180)
#             and accommodating a coarser grid: 
#            theta grid increment must be divisor of 90: 2,3,5,6,9,10,15,18,30,45)
#            phi increment must be divisor of 180: 2,3,4,5,6,9,10,12,15,18,20,30,45,90)
#
# Reminder: theta is angle with pos z axis (0-90), phi is azimuth (0-359)
#
# Note: Improved version: use any feasible grid increments, angle limits are now any float 0.0 - 90.0 or 180.0
#
# 
# Args:     arr     array of shape (90/tgrid+1, 180/pgrid+1), power in dB as returned by NEC5
#           [thetal, thetah], [phil, phih]    region of interest (degrees)
#
# Returns   relative power summed over specified region 
#           area of specified region
#
# The area of the quarter-sphere is pi, so the call
#   integrate_power_density_3(arr, thetal=0.0, thetah=90.0, phil=0.0, phih=180.0)
#       will return (pi, pi).....pi/pi=1.0 will be the relative power density for any pattern
#
def integrate_power_density_3(arr, thetal=0.0, thetah=90.0, phil=0.0, phih=180.0):

    # Compute span of angles associated with each grid point of radiation pattern
    #   grid    angles on the grid calculated by RP card, e.g. [0,10,20,...,180] or [0,10,20,...,90]
    #   startang, endang    range of interest (any float value)
    def get_angle_ranges(grid, startang, endang):
        ngrid = len(grid)
        gridinc = grid[1] - grid[0]
        angs = np.zeros((2,ngrid))       # Start, end angles for each grid point
        for i in range(ngrid):
            angs[0,i] = np.clip(grid[i] - gridinc/2, grid[0], grid[-1])
            angs[1,i] = np.clip(grid[i] + gridinc/2, grid[0], grid[-1])
        angs[angs < startang] = startang
        angs[angs > endang] = endang
        return angs

    assert arr.shape[0] in [91, 46, 31, 19, 16, 11, 10, 7, 6, 4, 3]
    assert arr.shape[1] in [181, 91, 61, 46, 37, 31, 21, 19, 16, 13, 11, 10, 7, 5, 3]
    thinc = 90 // (arr.shape[0] - 1)            # Grid increment (theta)
    phinc = 180 // (arr.shape[1] - 1)           # Grid increment (phi)
    thgrid = np.array(range(0,91,thinc))        # e.g. [0, 10, ..., 90] for grid = 10
    phgrid = np.array(range(0,181,phinc))       #      [0, 10, ..., 180]

    # 'tweight' = areas of horizontal strips associated with each grid point (within thetal,thetah limits)
    thlimits = get_angle_ranges(thgrid, thetal, thetah)         # Angle range to integrate for each grid point
    tweight = np.pi * np.apply_along_axis(lambda col: np.cos(np.deg2rad(col[0])) - np.cos(np.deg2rad(col[1])), 
                                 axis=0, arr=thlimits)

    # 'pweight' = frac of total area represented by vertical strips associated with each grid point (within phil,phih limits)
    phlimits = get_angle_ranges(phgrid, phil, phih)         # Angle range to integrate for each grid point
    pweight = np.apply_along_axis(lambda col: (col[1] - col[0]) / 180, axis=0, arr=phlimits)

    # Sum area and power weighted by area
    #   Convert dB power array to ratio, mult by corresponding area 
    #    -> pweight broadcast along cols, weight[:,None] is a column vec, broadcast along rows
    weighted_pwr = tweight[:,None] * (pweight * np.power(10, arr / 10))      # Mult by weight (as col vec)
    # Total relative power
    pwr = np.sum(weighted_pwr)
    # Area of selected region
    area = np.sum(tweight) * np.sum(pweight)

    return pwr, area



#
# Integrate power density produced by NEC5 RP card
#
# Version 4: numba-compiled equivalent of integrate_power_density_3()  (see above)
#
# Reminder: requires RP data covering half of upper hemisphere (theta 0-90, phi 0-180)
@njit
def integrate_power_density_4(arr, thetal=0.0, thetah=90.0, phil=0.0, phih=180.0):

    # Compute span of angles associated with each grid point of radiation pattern
    #   grid    angles on the grid calculated by RP card, e.g. [0,10,20,...,180] or [0,10,20,...,90]
    #   startang, endang    range of interest (any float value)
    def get_angle_ranges(grid, startang, endang):
        ngrid = len(grid)
        halfgrid = (grid[1] - grid[0]) / 2
        angs = np.zeros((2,ngrid))
        angs[0,:] = grid-halfgrid
        angs[1,:] = grid+halfgrid
        angs = np.minimum(np.maximum(angs,startang), endang)
        return angs

    assert arr.shape[0] in [91, 46, 31, 19, 16, 11, 10, 7, 6, 4, 3]
    assert arr.shape[1] in [181, 91, 61, 46, 37, 31, 21, 19, 16, 13, 11, 10, 7, 5, 3]
    thinc = 90 // (arr.shape[0] - 1)            # Grid increment (theta)
    phinc = 180 // (arr.shape[1] - 1)           # Grid increment (phi)
    thgrid = np.arange(0,91,thinc)        # e.g. [0, 10, ..., 90] for grid = 10
    phgrid = np.arange(0,181,phinc)       #      [0, 10, ..., 180]

    # 'tweight' = areas of horizontal strips associated with each grid point (within thetal,thetah limits)
    thlimits = np.cos(np.deg2rad(get_angle_ranges(thgrid, thetal, thetah)))         # Angle range to integrate for each grid point
    tweight = np.pi * (thlimits[0,:] - thlimits[1,:])

    # 'pweight' = frac of total area represented by vertical strips associated with each grid point (within phil,phih limits)
    phlimits = get_angle_ranges(phgrid, phil, phih)         # Angle range to integrate for each grid point
    pweight = (phlimits[1,:] - phlimits[0,:]) / 180
    
    # Sum area and power weighted by area
    #   Convert dB power array to ratio, mult by corresponding area 
    #    -> pweight broadcast along cols, weight[:,None] is a column vec, broadcast along rows
    weighted_pwr = tweight[:,None] * (pweight * np.power(10, arr / 10))      # Mult by weight (as col vec)
    # Total relative power
    pwr = np.sum(weighted_pwr)
    # Area of selected region
    area = np.sum(tweight) * np.sum(pweight)

    return pwr, area



#%%------------------------------------------------


#
# Find locations of specified NEC5 cards in a list of strings 
#
# Args:     strlst  e.g. ['str1', 'str2', ...]  (possibly multiline strings)
#           cardlist  e.g. ['EX', 'FR', ...]
# Returns list of indexes into list  (multiple cards of a single type not supported)
# Raises ValueError if any not found
#
def find_nec_cards(strlst, cardlst):
    cardidx = [None]*len(cardlst)
    for idx,s in enumerate(strlst):
        for cidx,card in enumerate(cardlst):
            if re.search('^'+card, s, re.MULTILINE):
                cardidx[cidx] = idx
    for idx in range(len(cardlst)):
        if cardidx[idx] == None:
            raise ValueError(f"find_nec_cards(): {cardlst[idx]} card not found")
    return cardidx



#
# Find antenna length for resonance at a specified frequency
#
# Args:
#   neclist     list of strings that will be concatenated to form NEC design deck
#               * must contain placeholder string for DE at index 'de_idx'
#   elemsects   list of lengths, diameters of sections of driven element (e.g. Yagi tubing sections)
#   freq        frequency of interest (MHz)
#   de_idx      pos in neclist of driven element
#   de_nsegs    DE's number of segments (for NEC5 simulation)
#   de_x, de_z  x,z positions of DE
#   de_tag      tag for DE
#   divs        # of points to evaluate within a range of DE lengths
#

def find_de_resonance(neclist, elemsects, freq, de_idx=2, de_nsegs=20, de_x=0, de_z=0, de_tag=2, divs=10, depth=2):
 
    desects = copy.deepcopy(elemsects)              # Make a copy of section lengths (we're going to modify it)
    de_inner_len = sum([t[0] for t in elemsects[1:-1]])          # DE length, excluding tips (inches)

    llim, ulim = 492*12/freq * np.array([0.8, 1.2])     # Initial guess: 1/2wl -10% / + 30%
    for deep in range(depth):

        delens = np.linspace(llim, ulim, divs)          # The DE lengths to evaluate
        designs = []
        for delen in delens:
            tiplen = (delen - de_inner_len) / 2         # DE tip length
            desects[0][0] = tiplen                      # Modify end section lengths            
            desects[-1][0] = tiplen
            destr, exseg = gen_tapered_el2(desects, de_x, de_z, de_tag, de_nsegs, inches=True)
            neclist[de_idx] = destr                     # NEC cards for DE
            # Create single multiline string, filling in params {freq}, etc.
            designs.append(''.join(neclist).format(tag=de_tag, seg=exseg, nfreq=1, fstart=freq, fstep=''))

        # Execute the designs in parallel
        res = nec5_sim_stdio3(designs)

        # Find the range where sign of reactance changes
        if res[0][0] == []:                                 # Check for failed NEC execution
            logger.info('Warning: find_de_resonance(): NEC5 execution failed')
            return (None, None)
        zs = [d[0][0][0][1] for d in res]                   # List of impedances
        xs = np.imag(zs)                                    # Array of reactances
        if xs[0] >= 0 or xs[-1] < 0:                        # Error if no sign change
            logger.info('Warning: find_de_resonance(): Resonant freq out of range')
            logger.info(f'{zs}')
            return (None, None)
        first_pos = np.nonzero(xs >= 0)[0][0]
        llim = delens[first_pos-1]                          # Range where resonance occurs
        ulim = delens[first_pos]


    return (llim, ulim)





#
# Find antenna length, shunt matching reactance for minimum VSWR at a specified frequency
#
# Args:
#   neclist     list of strings that will be concatenated to form NEC design deck
#               * must contain placeholder string for DE at index 'de_idx'
#   elemsects   list of lengths, diameters of sections of driven element (e.g. Yagi tubing sections)
#   freqs       tuple of (base_freq, #_of_freqs, freq_step)
#   de_idx      pos in neclist of driven element
#   de_nsegs    DE's number of segments (for NEC5 simulation)
#   de_x, de_z  x,z positions of DE
#   de_tag      tag for DE
#   divs        # of points to evaluate within a range of DE lengths
#   z_0         feedline Z
#

def find_de_min_vswr(neclist, elemsects, freqs, z_0=50, de_idx=2, de_nsegs=20, de_x=0, de_z=0, de_tag=2, divs=10, depth=2, match_X='L', peak=False):

    desects = copy.deepcopy(elemsects)              # Make a copy of section lengths (we're going to modify it)
    de_inner_len = sum([t[0] for t in elemsects[1:-1]])          # DE length, excluding tips (inches)
    # llim, ulim = 492*12/freqs[0] * np.array([0.9, 1.3])     # Initial guess: 1/2wl -10% / + 30%
    llim, ulim = 492*12/freqs[0] * np.array([0.8, 1.3])     # Initial guess: 1/2wl -20% / + 30%

    # Create array of 1/x_m values for all frequencies and matching reactances
    #   x_m = shunt matching reactance at lowest freq (0 - 199ohms, inductive or capacitive)
    #   shape = (200, n_freqs)
    if match_X =='L':
        inv_xms = np.reciprocal((np.arange(200)[...,None] + 1e-3)*1j)         # Pos reactances (column vector)
        fscale = freqs[0] / (np.arange(0,freqs[1])*freqs[2]+freqs[0])[None,...]     # Reactance scaling by freq (row)
        inv_x_m = inv_xms @ fscale      # shape (#reactances, #freqs)
    if match_X =='C':
        inv_xms = np.reciprocal((np.arange(200)[...,None] + 1e-3)*(-1j))         # Neg reactances (column vector)
        fscale = (np.arange(0,freqs[1])*freqs[2]+freqs[0])[None,...] / freqs[0]           # Reactance scaling by freq (row)
        inv_x_m = inv_xms @ fscale      # shape (#reactances, #freqs)

    retvals = []
    for deep in range(depth):

        delens = np.linspace(llim, ulim, divs)          # The DE lengths to evaluate
        designs = []
        for delen in delens:
            tiplen = (delen - de_inner_len) / 2         # DE tip length
            desects[0][0] = tiplen                      # Modify end section lengths            
            desects[-1][0] = tiplen
            destr, exseg = gen_tapered_el2(desects, de_x, de_z, de_tag, de_nsegs, inches=True)
            neclist[de_idx] = destr                     # NEC cards for DE
            # Create single multiline string, filling in params {freq}, etc.
            designs.append(''.join(neclist).format(tag=de_tag, seg=exseg, nfreq=freqs[1], fstart=freqs[0], fstep=freqs[2]))
        # return designs
        # Execute the designs in parallel
        res = nec5_sim_stdio3(designs)
        # res = nec5_sim_stdio3(designs, debug=True)
        # return res
        if res[0][0] == []:                                 # Check for failed NEC execution
            return None

        # Create an array of complex impedances calculated by NEC5, shape = (# of lengths, # of freqs)
        zarr = np.array( [[f[1] for f in d[0][0]]   for d in res] )

        # Scan for optimum VSWR across freq range, at each element length
        # @njit(parallel=True)      *** numba worsens execution time, for some reason
        def vswr_scan(arr, vswr_peak=True):

            # Results: min VSWR (mean or peak) and corresponding shunt reactance at each element length
            # min_vswr = np.zeros(arr.shape[0])                   # Min vswr
            min_vswr = np.ones(arr.shape[0]) * 1e9                   # Min vswr
            min_vswr_x = np.zeros(arr.shape[0], dtype=np.int64)      # Shunt reactance for min vswr

            # for idx in prange(arr.shape[0]):            # Iterate over all element lengths
            for idx in range(arr.shape[0]):            # Iterate over all element lengths

                inv_zs = np.reciprocal(arr[idx])        # Z at each freq
                matched_zs = np.reciprocal(inv_x_m + inv_zs)    # Z of element + shunt component
                gamma = (matched_zs - z_0) / (matched_zs + z_0) 
                abs_gamma = np.abs(gamma)
                # abs_gamma[abs_gamma==1.0] -= 1e-9
                vswr = (1 + abs_gamma) / (1 - abs_gamma)

                if vswr_peak:
                    vswr_criterion = np.amax(vswr, axis=1)
                else:
                    vswr_criterion = np.mean(vswr, axis=1)
                min_vswr_x[idx] = np.argmin(vswr_criterion)      # Reactance for lowest mean vswr
                min_vswr[idx] = vswr_criterion[min_vswr_x[idx]]

            # Alternative to above that works with numba:
                # for xval,vswr_row in enumerate(vswr):
                #     if vswr_peak:
                #         vswr_criterion = np.max(vswr_row)
                #     else:
                #         vswr_criterion = np.mean(vswr_row)
                #     if vswr_criterion < min_vswr[idx]:
                #         min_vswr[idx] = vswr_criterion
                #         min_vswr_x[idx] = xval

            return min_vswr, min_vswr_x

        # Returns array of vswr and matching reactance for each length
        vswrs, shuntx = vswr_scan(zarr, vswr_peak=peak)
        # Append to return vals:  element lengths, vswrs at those lengths, matching reactances
        retvals.append([delens,vswrs,shuntx])
        
        # Find min vswr and (possibly) iterate
        min_vswr = np.argmin(vswrs)
        clip_vswr = np.clip(min_vswr, 1, vswrs.shape[0]-2)
        if min_vswr != clip_vswr:
            logger.info('Warning: find_de_min_vswr(): VSWR minimum at end of element length range')
        llim = delens[clip_vswr-1]
        ulim = delens[clip_vswr+1]

    return retvals




#
# Find resonant frequency of an antenna
#
# Args:
#   necdesign               NEC design deck (as a single string)
#   nfreq, fstart, fstep    specify frequency range (MHz)
#

def find_de_resfreq(necdesign, nfreq, fstart, fstep, depth=2):
    
    for deep in range(depth):

        freqs = np.arange(nfreq) * fstep + fstart          # freqs to evaluate
        design = necdesign.format(nfreq=nfreq, fstart=fstart, fstep=fstep)
        res = nec5_sim_stdio3([design])           # Simulate
        
        # Find the range where sign of reactance changes
        if res[0][0] == []:                                 # Check for failed NEC execution
            logger.info('find_de_resfreq(): NEC5 execution failed')
            return (None, None)
        
        zs = [d[1] for d in res[0][0][0]]                   # List of impedances
        xs = np.imag(zs)                                    # Array of reactances
        if xs[0] >= 0 or xs[-1] < 0:                        # Error if no sign change
            logger.info('find_de_resfreq(): Resonant freq not within range')
            return (None, None)
        first_pos = np.nonzero(xs >= 0)[0][0]               # Index where x turns positive
        fstart = freqs[first_pos-1]                          # Range where resonance occurs
        fstep = fstep / (nfreq - 1)

    return (freqs[first_pos-1] , freqs[first_pos] )


#%%------------------------------------------------

import plotly.graph_objects as go


# Slightly simplified syntax for specifying our arrays of segments:
#  make_linear_element((0,0,0), (1,1,1), (2,3,4))   -->    [np.array([[0,0,0], [1,1,1], [2,3,4]])]
def make_linear_element(*args):
    return [np.array(args)]


# Total wire length in a design
def total_wire_len(des):
    return np.sum( [np.sum(np.linalg.norm(np.diff(arr,axis=0),axis=1)) for arr in des] )


#
# Translate and rotate antenna elements  (either a single array or list of arrays)
#
def _rot(elem,arr):
    if isinstance(elem, list):
        return [el @ arr for el in elem]
    return elem @ arr
def rot_x(elem,ang):
    return _rot(elem, np.array([[1,0,0],[0,np.cos(ang),np.sin(ang)],[0,-np.sin(ang),np.cos(ang)]]))
def rot_y(elem,ang):
    return _rot(elem, np.array([[np.cos(ang),0,-np.sin(ang)],[0,1,0],[np.sin(ang),0,np.cos(ang)]]))
def rot_z(elem,ang):
    return _rot(elem, np.array([[np.cos(ang),np.sin(ang),0],[-np.sin(ang),np.cos(ang),0],[0,0,1]]))

def translate(elem, dxyz):
    rowarr = np.array([dxyz])
    if isinstance(elem, list):
        return [el + rowarr for el in elem]
    return elem + rowarr



#
# Create a cage dipole element: cylindrical arrangement of 'n' wires with conical sections at ends
# Element will be constructed along y axis, between points y=y0 and y=y1
#       a:      height of conical section
#       r:      cylinder radius
#       n:      number of parallel wires
#
def make_cage_element(y0, y1, a, r, n):
    seg = np.array([[0,y0,0],[0,y0+a,r],[0,y1-a,r],[0,y1,0]])
    segs = [seg]
    for i in range(1,n):
        rads = i*2*np.pi/n
        segs.append(rot_y(seg, rads))
    return segs






#
# Simple wire antenna visualizer
#
#   segs        list of numpy arrays, each containing x,y,z coords of nodes in a segment
#   x,y,z       range of each axis, e.g. x=(-10,10) --> x-axis spans -10meters - +10meters
#   name        displayed when hovering
#
def wire_ant_visualize(segs,x,y,z,name='',width=800, height=700, mirror=True):

    # Create some visual elements: sxy = 'ground level', sxz = mirroring plane (transparent)
    # To change colors, see https://plotly.com/python/builtin-colorscales/
    sxy = go.Surface(x=x, y=y, z=np.full((2,2),0), colorscale='Greens', surfacecolor=np.full((2,2),0.69), showscale=False,cmin=0,cmax=1)
    sxz = go.Mesh3d(x=(x[0],x[0],x[1],x[1]), y=(0,0,0,0), z=(z[0],z[1],z[1],z[0]), i=(0,0),j=(1,2),k=(2,3),color='lightpink', opacity=0.1, hoverinfo='skip')
    
    # Make a (possibly) multi-segment trace segment for display
    #  arr: numpy array of x,y,z coords, shape = (#points, 3)
    def make_trace(arr):
        return go.Scatter3d(
            x=arr[:,0], y=arr[:,1], z=arr[:,2],
            marker=dict(size=2, color='red'),               # Set color, size of nodes
            line=dict(color='darkblue', width=2),           # Set color, thickness of line
            name=name, showlegend=False                     # 'name' is displayed when you hover
        )
    # Mirror list of arrays across xz plane
    def mirror_xz(arrs):
        mxz = np.array([[-1,-1,1]])
        return [ar * mxz for ar in arrs]
    
    fig = go.Figure(data=[sxy, sxz])                        # Create the figure, with fixed visual elements
    for seg in segs:
        fig.add_trace(make_trace(seg))                  # Add the wire segments
    if mirror:
        for seg in mirror_xz(segs):
            fig.add_trace(make_trace(seg))                  # ...and their mirrored versions

    # Set up appearance: size of viewer, range for each axis, etc.
    fig.update_layout(
        width=width, height=height, autosize=False,
        scene=dict(
            # 'eye' is your initial viewing angle
            camera=dict(up=dict(x=0,y=0,z=1),  eye=dict(x=1,y=.3,z=0.3)),
            # Relative display sizes of axes; this produces a square patch in xy plane, with z (elevation) half that length
            # aspectratio = dict( x=1, y=1, z=0.5 ),
            aspectratio = dict( x=1, y=1, z=1.0 ),
            aspectmode = 'manual',
            # Limits of each axis, in meters
            xaxis=dict(range=x, showbackground=False),
            yaxis=dict(range=y, showbackground=False),
            zaxis=dict(range=z, showbackground=False),
        ),
    )

    fig.show()





# 
# gen_nec5_str()
# Convert antenna descriptions (list of numpy arrays) to NEC5 cards (single string)
#
# Accepts a list of numpy arrays; each array describes one or more segments connected in series
# E.g. np.array([[0,0,0], [1,2,3], [2,5,4]]) --> two segments, connecting (0,0,0) to (1,2,3) to (2,5,4)
#
# Args:
#       arrs        a list of arrays as described above (all dimensions in meters)
#       segs_per_m  number of NEC segments per meter of wire length
#                   * these are the NEC segments used in the simulation;
#                      I've been confusingly also using the term for the wire conncting two points (see above)
#       tag         tag used on all NEC cards generated
#       radius      wire radius (meters)
#

# Clunky numba version: typical exec time 50-60us vs 300us for more elegant version that follows
from numba import njit
@njit
def gen_nec5_str_inner(arr, segs_per_m):
    s = []
    for i in range(arr.shape[0]-1):
        lng = np.sqrt(np.sum(np.square(arr[i+1] - arr[i])))
        s.append(( np.maximum(int(np.rint(lng*segs_per_m)), 1), arr[i,0],arr[i,1],arr[i,2],arr[i+1,0],arr[i+1,1],arr[i+1,2] ))
    return s
def gen_nec5_str(arrs, segs_per_m, radius, tag=1):
    ls = []
    gw = 'GW {0} {1} {2} {3} {4} {5} {6} {7} {8}\n'
    # For each array...
    for arr in arrs:
        data = gen_nec5_str_inner(arr, segs_per_m)
        for d in data:
            ls.append( gw.format(tag, *d, radius))
    return ''.join(ls)

            



#%%------------------------------------------------

# These contain various 'assert' statements to assure valid data; catch exceptions like this:
# try:
#     <your code>
# except AssertionError:
#     <deal with problem>

#
# Locate rising or falling edges of a function (used by res freq, BW estimation functions)
#  (rising = crossing from < 0 to >= 0)
#
# Args:
#   func        takes an array of frequencies, returns same-size array of real values
#   fs, xs      the initial arrays of frequencies and values
#   tol         iterates until consecutive freq estimates are less than tol (MHz)
#   reverse     if False, find first rising edge, else find last falling edge
#
def zoom_to_edge(func, fs, xs, tol=0.001, reverse=False, debug=False):
    assert np.any(xs>=0) and ((not reverse and xs[0]<0) or (reverse and xs[-1]<0))        # First (or last) value must be low
    xs = xs.copy()

    freq, niter = (-1, 0)
    while (True):
        niter += 1
        assert niter < 100

        xge0 = np.flatnonzero(xs>=0)                    # Indices of freqs with xs >= 0
        i = xge0[-1] if reverse else xge0[0]-1          # i = index just before edge
        flast = freq
        freq = (fs[i]*xs[i+1] - fs[i+1]*xs[i]) / (xs[i+1] - xs[i])      # Interpolate to get best freq estimate
        if np.abs(freq-flast) < tol:
            break

        # Zoom in
        fs = np.linspace(fs[i], fs[i+1], num=len(fs))   # Zoom in between freqs closest to edge
        xs[[0,-1]] = xs[[i,i+1]]                        # Retain prev values for start,end points
        xs[1:-1] = func(fs[1:-1])                       #  ...and fill in middle values

    if debug:
        return (freq, niter)
    return freq


# #
# # Test rig for above
# #
# for iter in range(10):
#     freq0 = np.random.rand()*40 + 1.01                      # Some random freq
#     freq1 = np.random.rand()*(50-freq0-6) + freq0 + 6   # Some larger freq
#     def func(fs):
#         return ((fs >= freq0) & (fs < freq1)) * 1.234 - 0.4321
#     fs = np.linspace(1, 51, num=11)                     # Initial freqs
#     xs = func(fs)                                       #  ...and values
#     print('f ', freq0, zoom_to_edge(func, fs, xs, tol=0.001, reverse=False, debug=True))
#     print('r ', freq1, zoom_to_edge(func, fs, xs, tol=0.001, reverse=True, debug=True))



#
# Find resonant frequency of a design
#
# Args:
#   necstr      design as a single string; must include placeholders 'flow', 'fstep', 'fnum'
#   flow, fhigh     extremes of frequency range to consider
#   fnum        number of freq points in range [flow, fhigh]
#   tol         freq tolerance to stop iteration (default: 0.001 MHz)
#
def find_res_freq(necstr, flow, fhigh, fnum=11, tol=0.001, debug=False):
    # Func that converts freqs to reactance
    def func(fs):
        res = nec5_sim_stdio3([necstr.format(fnum=len(fs), flow=fs[0], fstep=(fs[-1]-fs[0])/(len(fs)-1))])
        return np.imag([d[1] for d in res[0][0][0]])   
    fs = np.linspace(flow, fhigh, num=fnum)                 # Initial freqs
    xs = func(fs)                                           #  ...and values
    return zoom_to_edge(func, fs, xs, tol=0.001, reverse=False, debug=debug)


# #
# # Test rig for above
# #
# # Simple dipole template
# nec_dipole_start = """CE Dipole
# """
# gw_card = 'GW 1 10 0 0 {z} 0 {y} {z} 0.001\n'
# nec_dipole_end = """GX 100 010
# GE 1 0
# GD 0 0 0 0 13 0.005 0 0
# EX 4 1 1 1 1.0 0.0
# FR 0 {fnum} 0 0 {flow} {fstep}
# XQ 0
# EN
# """
# #
# # Generate test dipoles for freqs of approx 1 - 50 MHz
# for freq in range(1,51):
#     y = 142.646 / freq / 2                  # Approx half-dipole length (m)
#     z = np.random.rand() * 200 + 10         # Some random height (m)
#     str = nec_dipole_start + gw_card.format(z=z,y=y) + nec_dipole_end

#     f, iters = find_res_freq(str, freq*0.5, freq*1.5, tol=0.001, debug=True)
#     print(freq, f, iters)


#
# Find bandwidth of a design (VSWR less than some value)
#
# Args:
#   necstr      design as a single string; must include placeholders 'flow', 'fstep', 'fnum'
#   vswr_th     vswr threshold to use for computing bandwidth
#   flow, fhigh     extremes of frequency range to consider
#   fnum        number of freq points in range [flow, fhigh]
#   tol         freq tolerance to stop iteration (default: 0.001 MHz)
#
def find_vswr_bw(necstr, vswr_th, flow, fhigh, fnum=11, tol=0.001, debug=False):
    # Func that converts freqs to -(vswr - vswr_th)
    def func(fs):
        res = nec5_sim_stdio3([necstr.format(fnum=len(fs), flow=fs[0], fstep=(fs[-1]-fs[0])/(len(fs)-1))])
        xs = np.array([d[1] for d in res[0][0][0]])         # Feedpoint impedances
        return -(vswr(xs) - vswr_th)
    fs = np.linspace(flow, fhigh, num=fnum)                 # Initial freqs
    xs = func(fs)                                           #  ...and values

    bw_low = zoom_to_edge(func, fs, xs, tol=0.001, reverse=False)
    bw_high = zoom_to_edge(func, fs, xs, tol=0.001, reverse=True)

    if debug:
        return (bw_high-bw_low, bw_low, bw_high)
    return (bw_high-bw_low)

  
# #
# # Test rig for above
# #
# # Simple dipole template
# nec_dipole_start = """CE Dipole
# """
# gw_card = 'GW 1 10 0 0 {z} 0 {y} {z} 0.001\n'
# nec_dipole_end = """GX 100 010
# GE 1 0
# GD 0 0 0 0 13 0.005 0 0
# EX 4 1 1 1 1.0 0.0
# FR 0 {fnum} 0 0 {flow} {fstep}
# XQ 0
# EN
# """
# #
# # Generate test dipoles for freqs of approx 1 - 50 MHz
# freq = 10
# y = 142.646 / freq / 2                  # Approx half-dipole length (m)
# # z = np.random.rand() * 200 + 10         # Some random height (m)
# z = 50                                  # height (m)
# str = nec_dipole_start + gw_card.format(z=z,y=y) + nec_dipole_end

# t = find_vswr_bw(str, 2.0, freq*0.8, freq*1.2, fnum=11, tol=0.001, debug=True)
# print(t)

# res = nec5_sim_stdio3([str.format(fnum=100, flow=9, fstep=2/100)])
# zs = res[0][0][0]         # freqs, zs
# plot_vswr([zs], tags=[''])



#
# Trim a design to resonate at a specified frequency
#
# Args:
#   necstr      design as a single string; must include placeholders 'len', 'flow', 'fstep', 'fnum'
#               'len' is some design dimension that scales approximately as 1/f (e.g. dipole length)
#   freq        desired resonant freq
#   init_len    initial estimate of 'len'
#   tol         freq tolerance to stop iteration (default: 0.001 MHz)
#   f_range     range to search for resonance:  (1/f_range)*freq --> (f_range)*freq
#
def trim_res_freq(necstr, freq, init_len, tol=0.001, f_range=1.2, debug=False):

    x, iters = (init_len, 0)
    while True:
        iters += 1
        assert iters < 100

        s = necstr.format(len=x, flow='{flow}', fstep='{fstep}', fnum='{fnum}')
        f = find_res_freq(s, (1/f_range)*freq, (f_range)*freq, fnum=11, tol=tol/5)
        if np.abs(f - freq) < tol:
            break

        x *= f / freq       # Adjust dimension

    if debug:
        return (x, f, s, iters)
    return x




# #
# # Test rig for above
# #
# # Simple dipole template
# necstr = """CE Dipole
# GW 1 10 0 0 50 0 {len} 50 0.001
# GX 100 010
# GE 1 0
# GD 0 0 0 0 13 0.005 0 0
# EX 4 1 1 1 1.0 0.0
# FR 0 {fnum} 0 0 {flow} {fstep}
# XQ 0
# EN
# """
# #
# freq = 1.5
# y = 142.646 / freq / 2                  # Approx half-dipole length (m)

# (x, f, s, iters) = trim_res_freq(necstr, freq, y, tol=0.001, f_range=1.2, debug=True)
# print(freq, y, f, x, iters)

# res = nec5_sim_stdio3([necstr.format(len=x, fnum=100, flow=f*0.8, fstep=f*0.4/100)])
# zs = res[0][0][0]         # freqs, zs
# plot_vswr([zs], tags=[''])




#%%------------------------------------------------

#
# Implements the series-section-matching equations from Regier (1971)
#   zl      feedpoint impedances (numpy array)
#   z0a     z0 of TL nearest antenna (and feedline to radio)
#   z0b     z0 of matching section
#
#   Returns:    electrical lengths of TL sections (degrees):
#       [len(a), len(b), len(a)(second solution), len(b)(second solution)]
#       shape (len(zl), 4)

def ssm_thetas(zl, z0a=50, z0b=75):
    n = z0b/z0a
    r = np.real(zl) / z0a
    x = np.imag(zl) / z0a
    thetab = np.arctan(np.sqrt(((r-1)**2 + x**2) / (r*(n-1/n)**2 - (r-1)**2 - x**2)))
    thetaa = np.arctan(((n-r/n) * np.tan(thetab) + x) / (r + x*n*np.tan(thetab) - 1))
    if not isinstance(zl, np.ndarray):
        if thetaa < 0:
            thetaa += np.pi
    else:
        thetaa[thetaa < 0] += np.pi
    thetabn = np.pi - thetab
    thetaan = np.arctan(((n-r/n) * np.tan(thetabn) + x) / (r + x*n*np.tan(thetabn) - 1))
    if not isinstance(zl, np.ndarray):
        if thetaan < 0:
            thetaan += np.pi
    else:
        thetaan[thetaan < 0] += np.pi
    return np.rad2deg(  np.array((thetaa,thetab,thetaan,thetabn)).transpose() )


#%%------------------------------------------------

# Precompute some arrays
import numpy as np
from numba import njit

# Precompute arrays A,B,C,D to speed impedance transformation calculations (lossless TL)
#   zoa     z of first matching section (connected to antenna)
#   zob     z of second matching section
#   flow, fhigh, nfreq      freq band of interest
#
#   Usage:  precompute A,B,C,D once for specified parameters
#       A,B,C,D = series_match_precompute(zoa=50,zob=75,nfreq=9,flow=3.5,fhigh=4.0)
#
#       z = (zs*A[a,b] + B[a,b]) / ((1j)*zs*C[a,b] + D[a,b])
#           where   zs is a row vector of complex zs, shape (1,nfreq)
#                   a,b   lengths of zoa,zob matching section (degrees)
#
@njit
def series_match_precompute(zoa=50,zob=75,nfreq=9,flow=3.5,fhigh=4.0):
    # Scale phase delay for freqs across band of interest
    phscale = (np.linspace(flow,fhigh,num=nfreq) / ((fhigh+flow)/2))[None,:]
    A = np.empty((181,181,nfreq))
    B = np.empty((181,181,nfreq), dtype=np.complex128)
    C = np.empty((181,181,nfreq))
    D = np.empty((181,181,nfreq))
    for a in range(181):
        tana = np.tan(np.deg2rad(a*phscale))
        for b in range(181):
            tanb = np.tan(np.deg2rad(b*phscale))
            A[a,b] = zoa*zob - zob**2 * tana * tanb
            B[a,b] = zoa*zob*(zoa*tana + zob*tanb)*(1j)
            C[a,b] = zob*tana + zoa*tanb
            D[a,b] = zoa*(zob - zoa * tana * tanb)
    return A,B,C,D
# Exec time 35ms


# Scan combinations of matching-section lengths to find optimum vswr in band of interest
#
# Args:
#   zs      row vector of complex zs, shape (1,nfreq)
#   step    step size (degrees) of matching section lengths,
#           e.g. if step=2 all combinations of a=(0,2,4,...) and b=(0,2,4,...) will be tried
#
# Returns:  optimal lengths of matching sections zoa,zob for minimum vswr in freq band
#
@njit
def series_match_scan(zs, A,B,C,D, step=1, z0=50.0):
    aopt, bopt, vswr_max_opt, vswr_curve_opt = (0, 0, 99999.0, None)
    for a in range(0,A.shape[0],step):
        for b in range(0,A.shape[1],step):

            z = (zs*A[a,b] + B[a,b]) / ((1j)*zs*C[a,b] + D[a,b])
            arf = np.abs((z - z0) / (z + z0))         # Reflection coefs
            vswr_curve = (1 + arf) / (1 - arf)
            vswr_max = np.max(vswr_curve)
            if vswr_max < vswr_max_opt:
                vswr_max_opt = vswr_max
                aopt = a
                bopt = b
                vswr_curve_opt = vswr_curve
    return aopt,bopt,vswr_curve_opt,vswr_max_opt





#%%------------------------------------------------
#%%------------------------------------------------
#%%------------------------------------------------

