# flashes
#rf mapping
#ori grating
#natural scenes
import subprocess

scripts = [r'C:\Users\denma\Desktop\dan\danFlashCSD.py',
           r'C:\Users\denma\Desktop\dan\arbMatrix_highspeed.py',
           r'C:\Users\denma\Desktop\dan\grant_gratings_ori.py',
           r'C:\Users\denma\Desktop\dan\grant_scene_flicker.py',
           ]


#go through all of the stimuli
for script in scripts:
    subprocess.call(['python',script])


        
# start the opto stimulus
# on for 1 seconds every 2 seconds for 30 seconds every 5 mins
import datetime, time

from aibs.iodaq import DigitalInput, DigitalOutput, AnalogInput, AnalogOutput



#initialize digital I/O and bit states
do = DigitalOutput('Dev1',0)#device 1, port 0
do.StartTask()
do.WriteBit(4,0) #cyclops is plugged in to ch 7, port 0, dev1
for i in range(15):
    do.WriteBit(4,1)
    time.sleep(1)
    do.WriteBit(4,0)
    time.sleep(1)
do.WriteBit(4,0)
del(do)

scripts = [r'C:\Users\denma\Desktop\dan\danFlashCSD.py',
           r'C:\Users\denma\Desktop\dan\arbMatrix_highspeed.py',
           r'C:\Users\denma\Desktop\dan\grant_gratings_ori.py',
           r'C:\Users\denma\Desktop\dan\grant_scene_flicker.py',
           ]
#go through all of the stimuli
for script in scripts:
    subprocess.call(['python',script])

    do = DigitalOutput('Dev1',0)#device 1, port 0
    do.StartTask()
    do.WriteBit(4,0) #cyclops is plugged in to ch 7, port 0, dev1
    for i in range(15):
        do.WriteBit(4,1)
        time.sleep(1)
        do.WriteBit(4,0)
        time.sleep(1)
    do.WriteBit(4,0)
    del(do)
        