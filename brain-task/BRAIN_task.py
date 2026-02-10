
from datetime import datetime
import pygame
import os
import pygame, sys  # NOTE: After installing pygame, visual studio code needs to be restarted
import pygame.locals
import pandas
import time

import numpy as np

# VERY IMPORTANT: If they start tapping too early (before the go), 
# so that the key press is registered in the code window, instead, you need 
# to click on the square window, so that the key presses are registered (or it 
# too late, abort it, pressing Escape)

# ===========================
# Instruction:
# TRY TO ACCURATELY PRESS THE TWO KEYS IN ALTERNATION AS MANY TIMES AS YOU CAN
# WITHIN THE GIVEN TIME
#
# We'll first do 15 sec training with each hand to get you used to the go and stop sounds, 
# and then 30 sec of actual recording twice with each hand, which is what we will analyse.
#
# HAVE YOUR HAND READY, and start tapping as soon as you hear "GO", and 
# STOP TAPPING when you hear a trumpet

# ===========================
# Experimenter Input 
ID  =  'P1' 
TRAINING_ON = 0
hand = 'right'
#hand = 'left'

exerc = 'preExerc'
#exerc = 'postExerc'

if TRAINING_ON: 
     task_duration = 15
else:
     task_duration = 30
# ===========================

data_folder = 'C:/Users/ti21392/OneDrive - University of Bristol/Bristol/Experiments/PD_warrior/Pre_Post_BrainTask_Protocol/BRAIN_task/data/'
EEG_REC = False

# Create filename
today = datetime.now()
today_str = today.strftime("%d_%m_%Y_%H_%M_%S_")
filename  = today_str + ID + '_' + exerc + '_' + hand

WHITE = (255, 255, 255)
GREEN = (50, 255, 50)
GREY = (100, 100, 100)
BLACK = (0, 0, 0)

os.environ['SDL_VIDEO_WINDOW_POS'] = "50,50"# this will set a environment default window 
# set environment variable to specify the position of window
#'SDL_VIDEO_WINDOW_POS' is a specific environment variable used 
# by certain graphics libraries (such as SDL, Simple DirectMedia Layer) 
# to determine the initial window position.
#" %d,%d" % (x, y) is a string formatting operation 
# that replaces %d placeholders with the values of x and y. 
# In this case, it generates a string based on the values of x and `y.
# https://www.pygame.org/wiki/SettingWindowPosition


# Initiate the window settings
pygame.mixer.init()
pygame.init()
WIDTH = 200
HEIGHT = 200
surface = pygame.display.set_mode((WIDTH,HEIGHT),0,32)
pygame.display.set_caption('BRAIN TASK')
# surface color
surface.fill(WHITE)

pygame.display.flip() # displaying the window color

start_time = pygame.time.get_ticks() # starter tick

# Load sound files
start_sound = pygame.mixer.Sound("Go.wav")
end_sound   = pygame.mixer.Sound("Trumpet.wav")
start_sound.play() 
start_sound_duration = 1 # will be added to the task_duration

# Initialize data frame and values
data = pandas.DataFrame(columns=['strokeOnset','dwellTime','letters'])
seconds = 0
 
while seconds < (task_duration + start_sound_duration):
     seconds = (pygame.time.get_ticks() - start_time) / 1000
     
     for event in pygame.event.get():
          
          if event.type == pygame.KEYDOWN:           
               if event.key == pygame.K_ESCAPE: # quit program if you press "Escape"
                    pygame.quit()
                    sys.exit()
               elif event.key == pygame.K_s: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 's'
               elif event.key == pygame.K_SEMICOLON: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'semicolon'
               elif event.key == pygame.K_a: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'a'
               elif event.key == pygame.K_w: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'w'
               elif event.key == pygame.K_e: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'e'
               elif event.key == pygame.K_d: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'd'
               elif event.key == pygame.K_x: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'x'
               elif event.key == pygame.K_z: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'z'
               elif event.key == pygame.K_l: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'l'
               elif event.key == pygame.K_p: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'p'
               elif event.key == pygame.K_LEFTBRACKET: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'leftbracket'
               elif event.key == pygame.K_QUOTE: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'quote'
               elif event.key == pygame.K_SLASH: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'slash'
               elif event.key == pygame.K_PERIOD: # replace the 'p' to whatever key you wanted to be pressed
                    currKey = 'period'
               else:
                    currKey = 'other'

               strokeOnset_time = (pygame.time.get_ticks() - start_time)/1000
               
               if EEG_REC:
                    surface.fill(WHITE)
                    pygame.display.flip()
                             
          elif event.type == pygame.KEYUP:                             
               key_release_time = (pygame.time.get_ticks() - start_time)/1000
               dwell_time = key_release_time - strokeOnset_time
               # surface.fill(GREY) 
               # pygame.display.flip()   
               
               data.loc[len(data)] = [strokeOnset_time, dwell_time, currKey]

               # Save each button press
               #t = time.time()
               data.to_csv(data_folder + filename +'.csv', header=True, index=False, encoding='utf-8')
               #elapsed = time.time() - t
               #print(elapsed)
               
          elif event.type == pygame.locals.QUIT:
               pygame.quit()
               sys.exit()

end_sound.play()
time.sleep(2) # ensure enough time for playing trumpet sound
total_taps = len(data)

if total_taps > 0:
     total_dwell_time = np.sum(data['dwellTime'])
     distance = 15 # cm
     tap_duration  = np.max(data['strokeOnset']) - np.min(data['strokeOnset'])
     #tap_duration = tap_onsets[-1] - tap_onsets[0]
     mean_velocity = ((total_taps-1)*distance) / tap_duration 

     print('\nTotal Taps: '+ str(total_taps))
     print(f'Mean Velocity: {mean_velocity:.1f} cm/sec')
     print(f'Total dwell time: {total_dwell_time:.2f} sec')  
          
pygame.quit()
sys.exit() 