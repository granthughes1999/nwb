# -*- coding: utf-8 -*-
"""
Created on Mon Jan 04 12:43:51 2016

@author: stim1
"""
import numpy as np
from psychopy import visual, monitors
# from aibs.Projector.ProjectorWindow import *
from psychopy.visual.windowwarp import Warper
import _pickle as pkl
import datetime, os
try:
    from toolbox.toolbox.IO.nidaq import DigitalOutput
except:
    print("Could not import iodaq.")

from aibs.Core import ImageStimNumpyuByte
from aibs.SweepStim import SweepStim
import time


matrix = 'ternary'
frameDuration=3 # in frames
stimDuration=600 #in seconds
fps=60

#set up the window and some variables
mon = monitors.Monitor('uvprojector',width=115.2,distance=30,gamma=None,notes=None,
useBits=None,verbose=True)
mon.setSizePix([1024,768]) 
# win =  ProjectorWindow(monitor = 'testMonitor',color=(-1,-1,-10),screen=1,fullscr=True,warpfile=r'C:\Users\stim1\Desktop\dan\warps\uvProjector2_Warp2.data',warp=Warp.Warpfile,flipHorizontal=True)
win = visual.Window(monitor = 'testMonitor',screen=1,fullscr=True,useFBO=True)

warper = Warper(win,warpfile=r'C:\Users\denma\Desktop\dan\warps\uvProjector2_Warp2.data',warp='warpfile',flipHorizontal=True)
warper.changeProjection(warpfile=r'C:\Users\denma\Desktop\dan\warps\uvProjector2_Warp2.data',warp='warpfile')
# win =  ProjectorWindow(monitor = 'testMonitor',color=(0,0,0),screen=0,fullscr=False,warp=Warp.Disabled,flipHorizontal=True)

# win = visual.Window(monitor = 'testMonitor',screen=0,fullscr=False)
# warper = Warper(win,warpfile=r'C:\Users\stim1\Desktop\dan\warps\uvProjector2_Warp2.data',
#                     warp='warpfile',
#                     flipHorizontal=True)



numFrames = int(stimDuration * fps / frameDuration)

#generate the matrices
dim=64
movie= (np.round(np.random.uniform(0,1,(numFrames,dim,dim)))*-2).astype(np.ubyte)
print(movie.shape)
print(np.max(movie))
print(np.min(movie))
print('saving matrix to: '+os.path.join(r'C:\AibStim\arbMatrix',''.join((str(datetime.datetime.now()).replace(':','').replace('.','').replace(' ','-'),'.pkl'))))
fl = open(os.path.join(r'C:\AibStim\arbMatrix',''.join((str(datetime.datetime.now()).replace(':','').replace('.','').replace(' ','-'),'.pkl'))),'wb')
pkl.dump(movie,fl)
fl.close()
print('saved.')
# matrices = {
#     "ternary": np.floor(np.random.uniform(-1,2,(numFrames,dim,dim))),
#     "binary" : br,
#     "white" : np.round(np.random.uniform(0,1,(numFrames,dim,dim))),
#     "black" : np.round(np.random.uniform(1,1,(numFrames,dim,dim))),
#     "pink" : np.round(np.random.uniform(0,-1,(numFrames,dim,dim))),
# }    
# if isinstance(matrix,str):
#     movie = matrices[matrix].astype(np.ubyte)
#     print('saving matrix to: '+os.path.join(r'C:\AibStim\arbMatrix',''.join((str(datetime.datetime.now()).replace(':','').replace('.','').replace(' ','-'),'.pkl'))))
#     fl = open(os.path.join(r'C:\AibStim\arbMatrix',''.join((str(datetime.datetime.now()).replace(':','').replace('.','').replace(' ','-'),'.pkl'))),'wb')
#     pkl.dump(movie,fl)
#     fl.close()
#     print('saved.')
# else:
#     movie=matrix.astype(np.ubyte) #duh


# do = DigitalOutput('Dev1',0)#device 1, port 0
# do.StartTask()
# do.WriteBit(7,0)

params = {
    'mouseid': 'CHANGE_THIS_GRANT',
    "sweeplength": float(frameDuration)/60,  
    "postsweepsec": 0,
    "savesweeptable": False,
    'syncsqr': False,  # display a flashing square for synchronization
    'syncsqrloc': (0, 200),
    'syncsqrfreq': 60,  # square flips from black to white every x frames
    'syncsqrcolorsequence': [-1, 1],
    'nidevice': 'Dev1',
    'syncpulse':True,
    'doport':0,
    'diport':2,
    'invertdo':False,
    'syncpulselines':[1,2],
    'eyetracker':False,
    'eyetrackerip': "W7DTMJ82MX5", 
    'eyetrackerport': 10000,
    
    'runs':1,
}
movie_res = movie.shape[1:]
print('the movie is: '+str(np.shape(movie)[0])+' frames')
win.flip()
#CREATE BACKGROUND STIMULUS (Change size to stretch movie.)
bgStim = ImageStimNumpyuByte(win, image=movie[0], mask="None",
                             size=[800,600], ori=0, pos=(0, 0),
                             units='pix', flipVert=False)  

bgSweep = {
    'ReplaceImage': (movie, 0),
}

t0 = time.time()
g = SweepStim(win, params=params, bgStim=bgStim, bgSweep=bgSweep)
g.run()
win.flip()
print('that took: '+str(time.time()-t0)+' sec')
#return stack

