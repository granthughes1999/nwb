from psychopy import visual
import os
from aibs.Foraging import Foraging
from aibs.Projector.ProjectorWindow import *
from psychopy.visual.windowwarp import Warper


params = {
    'runs':50,
    'blanksweeps':0,
    'shuffle':False,
    'preexpsec':2, #seconds at the start of the experiment
    'postexpsec':0, #seconds at the end of the experiment
    'sweeplength':0.5, #length of sweeps
    'postsweepsec':0.5, #black period after sweeps (foreground remains)
    'mouseid': "CHANGE_THIS_GRANT",                                #name of the mouse
    'userid': "danield",                        #name of the user    
    'nidevice': 'Dev1',               #NI device name
    'blanksweeps': 0,             #blank sweep every x sweeps
    'bgcolor': 'black',                #background color
    'syncsqr': False,                  #display a flashing square for synchronization
    'syncsqrloc': (30,200),
    'syncsqrfreq': 60,                   #from black to white every X frames
    'syncpulse':True,
    'doport':0,
    'diport':2,
    'syncpulselines':[0,1,2],
    "savesweeptable": False,
    'eyetracker':False,
    'eyetrackerip': "W7DTMJ82MX5", 
    'eyetrackerport': 10000,
}

settings = WindowSettingsFromStimCfg()
mon = monitors.Monitor('uvprojector',width=115.2,distance=30,gamma=None,notes=None,useBits=None,verbose=True)
mon.setSizePix([1024,768])
win =  visual.Window(monitor = 'testMonitor',color=(0,0,0),screen=1,fullscr=True,useFBO=True)
warper = Warper(win, warpfile=r'C:\Users\denma\Desktop\dan\warps\uvProjector2_Warp2.data',
                                warp='warpfile',flipHorizontal=True)

path = r"C:\Users\denma\Desktop\stimuli_static\CAM_stimuli\selected_images"
imagefiles = [os.path.join(path,f) for f in os.listdir(path) if len(f) > 4 and f[-4:] in ['.jpg','.png','.tif','tiff']]
# print imagefiles[0][imagefiles[0].find('NaturalImages\\')+len('NaturalImages\\'):imagefiles[0].rfind('.')]
sequence = []
for i in imagefiles:
    filename = i[i.find('selected_images\\')+len('selected_images\\'):i.rfind('.')]
    firstcharacterinimagestring = filename[0]
    if firstcharacterinimagestring.lower() == 'm' or firstcharacterinimagestring.lower() == 'p':
        sequence.append(visual.ImageStim(win,image = i,flipVert=True,flipHoriz=True,size=(228*1.2,170*1.2),pos=(-50,-50),units='pix'))
    else:
        sequence.append(visual.ImageStim(win,image = i,flipVert=True,flipHoriz=True,size=(256*1.2,170*1.2),pos=(-50,-50),units='pix'))


g = Foraging(win, params, bgStim=sequence)
g.imagefiles = imagefiles


g.run()


