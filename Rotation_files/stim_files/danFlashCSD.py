# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 20:38:56 2014

@author: djd
"""
import numpy as np
from psychopy import visual, monitors
from aibs.Projector.ProjectorWindow import *
import datetime
from aibs.iodaq import DigitalInput, DigitalOutput, AnalogInput, AnalogOutput


def danFlash(number,contrast=1):
      
    mon = monitors.Monitor('uvprojector',width=115.2,distance=30,gamma=None,notes=None,useBits=None,verbose=True)
    mon.setSizePix([1024,768])
    win =  visual.Window(monitor = 'testMonitor',color=(0,0,0),screen=1,fullscr=True,useFBO=True)
    isi=2.95
    fps=60
    stimDuration = 0.05#in seconds
    stimDuration = int(fps*stimDuration);
    
    totalFrames = int(round(number*isi*fps));print(totalFrames)
    c=1
    
    
    #stim = visual.GratingStim(win,color=(0,0,0),mask=None, units="pix",size = (912,1140))
    white = np.ones((4,4))*contrast
    black = np.ones((4,4))*contrast*-1
    grey = np.zeros((4,4))
    stim = visual.GratingStim(win, tex=grey,mask=None, units="pix",size = (1024,768))
    col = 'grey'

    #initialize digital I/O and bit states
    do = DigitalOutput('Dev1',0)#device 1, port 0
    do.StartTask()
    do.WriteBit(1,0)
    do.WriteBit(2,0)
    do.WriteBit(3,0)
    do.WriteBit(5,0)
    
    #start eye tracking 
##    try:
##        eyetrackerip = "W7DTMJ82MX5"
##        eyetrackerport = 10000
##        trackeyepos = False
##        from aibs.Eyetracking.EyetrackerClient import Client
##        eyetracker = Client(outgoing_ip=eyetrackerip,
##                                 outgoing_port=eyetrackerport,
##                                 output_filename=str(datetime.datetime.now()).replace(':','').replace('.','').replace(' ','-'))
##        eyetracker.setup()
##        eyedatalog = []
##        if trackeyepos:
##            eyeinitpos = None
##    except Exception, e:
##        print "Could not initialize eyetracker:", e
##        eyetracker = None    
##    eyetracker.recordStart()
    
    #stim.setColor((0,0,0),colorSpace='rbg')
    
    do.WriteBit(1,1)
    do.WriteBit(4,0)
    for i in range(number):
        opto = False;#np.round(np.random.uniform(0,1))
        if col is 'white':
            c=black
            col = 'black'
            c2=-1
        else:
            c=white
            col = 'white'
            c2=1
        print(str(i)+':  '+col)
        if opto:
            do.WriteBit(4,1)
        for j in range(stimDuration):
            stim.setTex(c)                
            #stim.setColor([c2,c2,c2])#,colorSpace='rbg')        
            stim.draw()           
            win.flip()
            do.WriteBit(2,1)

        for k in range(int(round(isi*fps))):
            stim.setTex(grey) 
            #stim.setColor([0,0,0])#,colorSpace='rbg')  
            stim.draw()
            win.flip()
            do.WriteBit(2,0)
            if opto:
                if k>9: #number of frames after offset to keep opto on for; 9 = 150msec
                    do.WriteBit(4,0)
  
    #cleanup   
    do.WriteBit(1,0)     
    win.close()
##    eyetracker.recordStop()
    do.StopTask()
    do.ClearTask()
                
if __name__ == "__main__":
    d=danFlash(100,1.0)
