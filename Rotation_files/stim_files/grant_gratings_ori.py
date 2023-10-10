"""

gratingsexp.py

@author: derricw

updaded: 2014/6/8

This is an example of a grating stimulus that varies is screen location and
    orientation.  It also features a sync square that flashes black to white
    every 60 frames.

"""

from aibs.SweepStim import SweepStim
from psychopy import visual
from aibs.Projector.ProjectorWindow import *
from psychopy.visual.windowwarp import Warper


params = {
    'runs': 40,  # number of runs
    'shuffle': True,  # shuffle sweep tables
    'preexpsec': 2,  # seconds at the start of the experiment
    'postexpsec': 2,  # seconds at the end of the experiment
    'sweeplength': 1.5,  # length of sweeps
    'postsweepsec': 0.5,  # black period after sweeps (foreground remains)
    'syncsqr': False,  # display a flashing square for synchronization
    'syncsqrloc': (30, 200),
    'syncsqrfreq': 60,  # square flips from black to white every x frames
    'syncsqrcolorsequence': [-1, 1],
    'nidevice': 'Dev1',
    'syncpulse':True,
    'doport':0,
    'diport':2,
    'syncpulselines':[1,2],
    'eyetracker':False,
    'eyetrackerip': "W7DTMJ82MX5", 
    'eyetrackerport': 10000,
    'mouseid': 'CHANGE_THIS_GRANT',
    }

#INITIALIZE WINDOWS
#WINDOW
settings = WindowSettingsFromStimCfg()
mon = monitors.Monitor('uvprojector',width=115.2,distance=30,gamma=None,notes=None,useBits=None,verbose=True)
mon.setSizePix([1024,768])
win =  visual.Window(monitor = 'testMonitor',color=(0,0,0),screen=1,fullscr=True,useFBO=True)
warper = Warper(win, warpfile=r'C:\Users\denma\Desktop\dan\warps\uvProjector2_Warp2.data',
                                warp='warpfile',flipHorizontal=True)

#CREATE BACKGROUND STIMULUS
grating = visual.GratingStim(win, tex="sqr", mask='raisedCos', texRes=256,
                             size=[1024, 768], sf=1, ori=0, units='deg',
                             )

#CREATE BACKGROUND SWEEP PARAMETERS ([possible_values], column_in_table)
bgSweep = {
    'Ori': ([0,22.5,45,67.5,90,112.5,135,157.5], 6),
    'SF': ([0.08,0.16], 3),
    'Contrast': ([1], 0),
    'TF': ([2], 2),
    'Phase': ([0], 4),
    'PosX': ([0], 5),
    'PosY': ([0], 1),
    }
#CREATE FORAGING CLASS INSTANCE
g = SweepStim(window=win, params=params, bgStim=grating, bgSweep=bgSweep)
#RUN IT
g.run()
