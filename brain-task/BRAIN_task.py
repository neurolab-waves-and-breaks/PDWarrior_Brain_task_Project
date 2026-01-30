# input 
ID  =  'P5' 
hand = '_right'
#hand = “_left”
exerc = '_preExerc'
#exerc = '_postExerc'
 
 
from datetime import date
today = date.today()
today_str = today.strftime("%d_%m_%Y")
file_name = ''
# patient_id = input("Enter the patient's identifier: ")
filename = today_str + ID  +  exerc + hand

import pygame
import os
import pygame, sys  # NOTE: After installing pygame, visual studio code needs to be restarted
import pygame.locals
import pandas
import time

import numpy as np


WHITE = (255, 255, 255)
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


# initiating the window settings
pygame.mixer.init()
pygame.init()
WIDTH = 300
HEIGHT = 300
surface = pygame.display.set_mode((WIDTH,HEIGHT),0,32)
pygame.display.set_caption('BRAIN TASK')
# surface color
surface.fill(WHITE)
pygame.display.flip()
#Move the window to the top left corner of the screen


# Define the position for the fixation point (center of the window)
fixation_color = (BLACK)
fixation_position = (WIDTH // 2, HEIGHT // 2)
line_length = 10  # Adjust the length of the cross as needed
line_thickness = 2  # Adjust the line thickness as needed

    # Vertical line
pygame.draw.line(surface, fixation_color, (fixation_position[0], fixation_position[1] - line_length),(fixation_position[0], fixation_position[1] + line_length), line_thickness)

    # Horizontal line
pygame.draw.line(surface, fixation_color, (fixation_position[0] - line_length, fixation_position[1]),(fixation_position[0] + line_length, fixation_position[1]), line_thickness)


pygame.display.flip() # displaying the window color




start_time = pygame.time.get_ticks() # starter tick


#windowSurface = pygame.display.set_mode((WIDTH, HEIGHT), 0, 32)  #  TODO: This needs to be moved into top left corner
#surface = pygame.display.set_mode((100, 100))


#load sound files
start_sound = pygame.mixer.Sound("beep.wav")
end_sound = pygame.mixer.Sound("Trumpet.wav")
start_sound.play() 
start_sound_time = pygame.time.get_ticks() + 1000  # 1000 milliseconds piror to the start of task
task_duration = 6



data = pandas.DataFrame(columns=['strokeOnset','dwellTime','letters'])
seconds = 0
taps = 0
total_dwell_time = 0
# at the begining play beep sound
#trumpet_sound_played = False # trumpet sound played only once after the task 
while seconds < task_duration:
     seconds = (pygame.time.get_ticks() - start_time)/1000
     for event in pygame.event.get():
          if event.type == pygame.KEYDOWN:           
               if event.key == pygame.K_ESCAPE: # quit program if you press "Escape"
                    #data.to_csv('file_name.csv', header=True, index=False, encoding='utf-8')
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
               #surface.fill(WHITE)
               #pygame.display.flip()
               print(currKey + ' pressed')   
               strokeOnset_time = (pygame.time.get_ticks() - start_time)/1000
               
          elif event.type == pygame.KEYUP:  
                            
               key_release_time = (pygame.time.get_ticks() - start_time)/1000
               dwell_time = key_release_time - strokeOnset_time
               print(dwell_time)
               # surface.fill(GREY) 
               # pygame.display.flip()   
               
               data.loc[len(data)] = [strokeOnset_time, dwell_time, currKey]

               # Save each button press
               # Save the data with the patient's identifier
               data.to_csv(f'D:\PD_warrior_data\PD_warrior_tests\{filename}.csv', header=True, index=False, encoding='utf-8')

          elif event.type == pygame.locals.QUIT:
               pygame.quit()
               sys.exit()

if seconds >= task_duration: # when task is finished, play the trumpet sound
     end_sound.play()
     time.sleep(2) # ensure enough time for playing trumpet sound
     #pygame.quit()
     #sys.exit()
     total_taps = len(data)
     if total_taps > 0:
          total_dwell_time = np.sum(data['dwellTime'])
     mean_velocity = total_taps / total_dwell_time if total_dwell_time > 0 else 0

     #print(f"Total Taps: {total_taps}")
     #print(f"Mean Velocity: {mean_velocity:.2f} taps/second")
     print('Total Taps: '+ str(total_taps))
     print('Mean Velocity: ' + str(mean_velocity )+ 'taps/second')
    
     pygame.quit()
     sys.exit() 
         
       #if not trumpet_sound_played: 
            #end_sound.play()
            #time.sleep(2) # ensure enough time for the trumpet sound 
            #trumpet_sound_played = True 
# Quit the task
#pygame.quit()
#sys.exit()
    # Vertical line
    # pygame.draw.line(surface, fixation_color, (fixation_position[0], fixation_position[1] - line_length),(fixation_position[0], fixation_position[1] + line_length), line_thickness)

    # Horizontal line
   #  pygame.draw.line(surface, fixation_color, (fixation_position[0] - line_length, fixation_position[1]),(fixation_position[0] + line_length, fixation_position[1]), line_thickness)
    # pygame.display.update()

# Quit the task     
#pygame.quit()
#ys.exit()