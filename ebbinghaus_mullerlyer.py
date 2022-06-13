#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2022.1.4),
    on Mon Jun 13 16:09:31 2022
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

from psychopy.hardware import keyboard



# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)
# Store info about the experiment session
psychopyVersion = '2022.1.4'
expName = 'ebbinghaus_mullerlyer'  # from the Builder filename that created this script
expInfo = {
    'participant': '',
}
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='/Users/workingman/Documents/ASD_and_Synesthesia/Python/psychopy_shared_github/illusions_ebbinghaus_mullerlyer/ebbinghaus_mullerlyer.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.DEBUG)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# Setup the Window
win = visual.Window(
    size=[1280, 800], fullscr=True, screen=0, 
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True, 
    units='deg')
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess
# Setup ioHub
ioConfig = {}
ioSession = ioServer = eyetracker = None

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard(backend='ptb')

# Initialize components for Routine "welcome"
welcomeClock = core.Clock()
welcome_txt = visual.TextStim(win=win, name='welcome_txt',
    text="In this test, you will do tasks which each involves a pair of objects.\n\nIn each task, you are to make sure that the objects' sizes are as similar as possible. To do this, you adjust the size of one of the objects.\n\nTo increase the size of an object, use the up arrow key. To reduce its size, use the down arrow key. When you think the objects are the same size, press space to confirm.\n\nBefore each task you will see a small green box. It shows where the object that you can resize will be located.\n\nWhen you have finished reading, press space to begin.",
    font='Arial',
    units='deg', pos=(0, 0), height=1, wrapWidth=30, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
welcome_keyb = keyboard.Keyboard()

# Initialize components for Routine "target_indicator_1000ms"
target_indicator_1000msClock = core.Clock()
target_rectangle = visual.Rect(
    win=win, name='target_rectangle',units='deg', 
    width=(1, 1)[0], height=(1, 1)[1],
    ori=0, pos=(0, 0), anchor='center',
    lineWidth=1,     colorSpace='rgb',  lineColor=[0,1,0], fillColor=[0,1,0],
    opacity=1, depth=0.0, interpolate=True)

# Initialize components for Routine "illusions_trials"
illusions_trialsClock = core.Clock()
from random import sample
from math import asin, acos
from numpy import sign
from copy import deepcopy

# information fetched from
# Burghoorn, F., Dingemanse, M., Van Lier, R., & Van Leeuwen,
# T. M. (2019). The relation between the degree of synaesthesia,
# autistic traits, and local/global visual perception.
# Journal of Autism and Developmental Disorders.
# Advance Online Publication.

# all size/coordinate specifications in the experiment use deg units

# set positions of left/right central circles
left_center = (-7, 0)
right_center = (7, 0)

# set size (diameter) of "reference" circle, i. e. the circle whose
# size is to be matched
ref_circ_diam = 1.25
# set length of "reference" line
ref_line_len = 3
# set diameter of "reference" star
ref_star_diam = 2
# set how long the smaller arrow components should be in reference to the 'reference' line
# (basing this on manual measurements of image in Burghoorn)
small_line_len = 0.1851 * ref_line_len
# set angle (in rad) of smaller line from tip of horizontal line (is 45 deg according to Burghoorn)
small_line_angle = pi/4
# set orientations for small components in inward/outward arrow
inward_oris =  [x*small_line_angle/pi*180 for x in [1, -1, 1, -1]]
outward_oris = [x * small_line_angle/pi*180 for x in [-1, 1, -1, 1]]

# set minimum initial size of adjustable circle (.68, as per Burghoorn article)
adjust_circle_init_min = .68
# set maximum initial size of adjustable circle
adjust_circle_init_max = 1.82
# set minimum initial size of adjustable line (2.43, as per Burghoorn article)
adjust_line_init_min = 2.43
# set maximum initial size of adjustable line
adjust_line_init_max = 3.86
# set minimum initial size of adjustable star (not specified in Burghoorn)
adjust_star_init_min = 0.7
# set maximum initial size of adjustable star
adjust_star_init_max = 2.5

# set minimum/maximum sizes that stimuli can be adjusted to (this isn't
# specified in the Burghoorn article, instead they were chosen based
# on what seems reasonable
adjust_circle_min = .59 * adjust_circle_init_min
adjust_circle_max = adjust_circle_init_max
adjust_line_min = 0.62 * adjust_line_init_min
adjust_line_max = 1.63 * adjust_line_init_max
adjust_star_min = 0.71 * adjust_star_init_min
adjust_star_max = 1.20 * adjust_star_init_max

# set distance of small/large context circles from
# central circles
small_circ_dist = 1.25
large_circ_dist = 2.08

# set distance of small arrow components' centres from centre of
# long/middle component
# (basing these distances on manual measurements of image in Burghoorn,
# and using the rule of cosines)
# - when arrows are 'outward'
outward_dist = ((ref_line_len/2)**2 + (small_line_len/2)**2 - \
               2*(ref_line_len/2)*small_line_len/2*cos(pi-small_line_angle))**(1/2)
# - when arrows are 'inward'
inward_dist = ((ref_line_len/2)**2 + (small_line_len/2)**2 - \
               2*(ref_line_len/2)*small_line_len/2*cos(small_line_angle))**(1/2)

# set angle (in rad) going from center of long component of arrow to
# center of top-right outward arrow component
outward_angle = asin(small_line_len/2*sin(pi/4)/
                     outward_dist)

# set angle (in rad) going from center of long component of arrow to
# center of top-right inward arrow component
inward_angle = asin(small_line_len/2*sin(pi/4)/
                     inward_dist)


# set size (diameter) of smaller/larger context circles
small_circ_diam = .42
large_circ_diam = 1.67

# generate positions of the 8 smaller context circles,
# by using the fact that
# cos(angle) = (distance from center)/(xcoord)
# and
# sin(angle) = (distance from center)/(ycoord)
# when appearing on the left side
small_circ_pos_l = []
for c in range(8): # 8 circles are to be generated
    x_displ = cos(c*(2*pi)/8) * small_circ_dist
    y_displ = sin(c*(2*pi)/8) * small_circ_dist
    new_pos = (left_center[0]+x_displ, left_center[1]+y_displ)
    small_circ_pos_l.append(new_pos)
# when appearing on the right side
small_circ_pos_r = [(x + right_center[0] - left_center[0], y + right_center[1] - left_center[1]) for x, y in small_circ_pos_l]

# use an analogous process to generate positions of
# the 4 larger context circles
# when appearing on the left side
large_circ_pos_l = []
for c in range(4):
    x_displ = cos(c*(2*pi)/4) * large_circ_dist
    y_displ = sin(c*(2*pi)/4) * large_circ_dist
    new_pos = (left_center[0]+x_displ, left_center[1]+y_displ)
    large_circ_pos_l.append(new_pos)
# when appearing on the right side
large_circ_pos_r = [(x + right_center[0] - left_center[0], y + right_center[1] - left_center[1]) for x, y in large_circ_pos_l]

# generate positions of inward arrow short/small components
inward_line_pos_l = []
inward_angles = [inward_angle, pi-inward_angle, # list of all angles from long component center to short components' centers
                 pi+inward_angle, -inward_angle]
for alpha in inward_angles: # 4 lines are to be generated
    x_displ = cos(alpha) * inward_dist
    y_displ = sin(alpha) * inward_dist
    new_pos = [left_center[0]+x_displ, left_center[1]+y_displ]
    inward_line_pos_l.append(new_pos)
# when appearing on the right side
inward_line_pos_r = [[x + right_center[0] - left_center[0], y + right_center[1] - left_center[1]] for x, y in inward_line_pos_l]

# generate positions of outward arrow short/small components
outward_line_pos_l = []
outward_angles = [outward_angle, pi-outward_angle, # list of all angles from long component center to short components' centers
                 pi+outward_angle, -outward_angle]
for alpha in outward_angles: # 4 lines are to be generated
    x_displ = cos(alpha) * outward_dist
    y_displ = sin(alpha) * outward_dist
    new_pos = [left_center[0]+x_displ, left_center[1]+y_displ]
    outward_line_pos_l.append(new_pos)
# when appearing on the right side
outward_line_pos_r = [[x + right_center[0] - left_center[0], y + right_center[1] - left_center[1]] for x, y in outward_line_pos_l]


# create adjustable circle polygon
adjust_circle = visual.Polygon(
    win=win, name='adjust_circle',units='deg',
    edges=120, size=[1.0, 1.0],
    ori=0, pos=(0, 0),
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

# create a circle polygon that will be used for drawing all static
# circles (i. e. all circles except the adjustable one)
static_circle = visual.Polygon(
    win=win, name='static_circle',units='deg',
    edges=120, size=[1.0, 1.0],
    ori=0, pos=(0, 0),
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

# create rectangle polygon that will be used as the adjustable component
# of muller-lyer arrows. the height value is about
# 2.36% of that of the reference line's width, since that is what manual measurement
# of image in Burghoorn indicates.
adjust_line = visual.Rect(
    win=win, name='adjust_line',
    width=ref_line_len, height=0.02364 * ref_line_len,
    ori=0, pos=(0, 0),
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-2.0, interpolate=True,
    units='deg')

# create a context rectangle polygon
# that will be used for static lines in muller-lyer
# arrows, (with length that is 18.51%
# of that of the long reference muller-lyer arrow component when drawing
# shorter arrow components, since that is what
# manual measurement of stimuli images in Burghoorn indicates)
static_line = visual.Rect(
    win=win, name='static_line',
    width=0.1851 * ref_line_len, height=0.02364 * ref_line_len,
    ori=-45, pos=(0, 0),
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-2.0, interpolate=True,
    units='deg')

# create an adjustable star stimulus that will be used
# in "star"/practice trials
adjust_star = visual.ShapeStim(
    win=win, name='adjust_star', vertices='star7',
    size=(0.5, 0.5),
    ori=0, pos=(0, 0),
    lineWidth=1, lineColor=[1,1,-1], lineColorSpace='rgb',
    fillColor=[1,1,-1], fillColorSpace='rgb',
    opacity=1, depth=-2.0, interpolate=True,
    units='deg')

# create a reference star stimulus that will be used
# in "star"/practice trials
static_star = visual.ShapeStim(
    win=win, name='static_star', vertices='star7',
    size=(ref_star_diam, ref_star_diam),
    ori=0, pos=(0, 0),
    lineWidth=1, lineColor=[1,1,-1], lineColorSpace='rgb',
    fillColor=[1,1,-1], fillColorSpace='rgb',
    opacity=1, depth=-2.0, interpolate=True,
    units='deg')

# generate tuples with specifications of if context/distractor
# stimuli should or shouldn't be drawn, based on
# "four task sessions: two sessions of the Ebbinghaus task
# and two sessions of the Müller- Lyer task, in counterbalanced
# order [...] Of the 20 experimental trials per session 10 trials
# were context-free and 10 trials had context. The context- and
# context-free trials were blocked, and across the two sessions per
# task, it was counterbalanced whether participants started with
# 10 context-trials or 10 context-free trials. Within each 10
# trials trial order was randomised."
# NOTE: 'con' means context. 'noc' means no context. 'ebbing'
# means Ebbinghaus illusion, 'muller' Muller-Lyer illusion.
# this generates a list where trial_spec[0][0] is the type of
# illusion for first block, [0][1][0] is what happens for first
# 10 trials, [0][1][0][0] describes if context should/should not be
# included for 10 trials, [0][1][0][1] describes where small circles
# or inward arrows should be (on left or right side)
trial_spec = []

# randomly select order of illusion type blocks - also making a separate list of
# context positions for context trials (for use when forming list of to-be-adjusted
# stimulus positions)
con_positions = [] # separate list of context positions for context trials
con_noc_order = [] # separate list for storing order of 'context'/'nocontext' groupings

# choose one of the two possible orderings of ebbinghaus/muller-lyer tasks
ill_type_order = sample([['ebbing', 'muller', 'ebbing', 'muller'], ['muller', 'ebbing', 'muller', 'ebbing']], 1)[0]

# randomly deciding for each task type if it its first or second block should be 
# the one where the first 10 trials have 'no context'
cont_spec = sample(([['noc', 'con'], ['con', 'noc']], [['con', 'noc'], ['noc', 'con']]), 2)
cont_spec = cont_spec[0][0], cont_spec[1][0], cont_spec[0][1], cont_spec[1][1]

# for each illusion type value, generate specifications for a block of trials
for i, ill_type in enumerate(ill_type_order): 
    # initialize a list that will hold two tuples - one for 10 'context' trials,
    # another for 10 'nocontext' trials
    cont_ls = []
    for cont_type in cont_spec[i]: # for 'no context'/'noc' and 'context'/'con' trials
        loop_pos = sample(['l']*5+['r']*5, 10) # create a list of 5 'r' and 5 'l' vals, randomly shuffled
        cont_ls.append( (cont_type, loop_pos)) # append a tuple with the content type spec and the 10 'l'/'r' vals
        if cont_type == "con": # if it's a set of trials where context is included
            con_positions += loop_pos # add the list of 'r' and 'l' values to total list of context pos vals
        con_noc_order.append(cont_type)
    # append, to the trial specifications list a tuple with the illusion type and
    # the list of 'noc'/'con' trials tuples
    trial_spec.append( (ill_type, cont_ls) )

# initialize counters that will be used for stepping through the
# trial_spec lists
l_r_count = 0 # for stepping through left/right values
first_last_10_count = 0 # for stepping through first 10/last 10 trials of each block
block_count = 0 # for stepping through blocks/illusion type values

# forming a list of to-be-adjusted stimulus positions, where it’s
# counter-balanced on a “per-illusion” level whether
# A the to-be-adjusted stimulus is embedded in a ‘large’/‘small’ context
# (‘outward’/‘inward’ arrow)
# B the to-be-adjusted stimulus is on the left or right hand side of the screen
# 1. extract 1st block+3rd block and 2nd block+4th block context trial, ‘l’/‘r’
# context position lists
# 2. get indices of all ‘l’ values in 1st+3rd block list
# 3. form an empty list of length 20 (block_1_3_ls, later block_2_4)
# 4. randomly sample from (‘l’, ‘r’)*5
# 5. put the values from step 4 into the list from step 3, at the indices
# fetched in step 2
# 6. repeat steps 4-5, only this time using the indices of ‘r’ values in
# 1st+3rd block list
# (i. e. the complement to the indices from step 2)
# 7. repeat steps 2-6, only this time doing it with the 2nd block+4th block
# ‘l’/‘r’ list
# 8. Join the two lists produced by steps 1-7, first putting them all together,
# then separating them using random samplings of 10 'l'/'r' values for groupings
# of 'no context' trials
con_pos_1_3 = con_positions[0:10] + con_positions[20:30]
con_pos_2_4 = con_positions[10:20] + con_positions[30:40]
l_1_3 = [ind for ind, val in enumerate(con_pos_1_3) if val=='l']
r_1_3 = [ind for ind, val in enumerate(con_pos_1_3) if val=='r']
l_2_4 = [ind for ind, val in enumerate(con_pos_2_4) if val=='l']
r_2_4 = [ind for ind, val in enumerate(con_pos_2_4) if val=='r']
adj_stim_lr_con_1_3 = ['']*20
adj_stim_lr_con_2_4 = ['']*20
for side_indices_ls in [l_1_3, r_1_3]:
    rand_10_pos = sample(['l', 'r'] * 5, 10)
    for rand_pos_ind, ind in enumerate(side_indices_ls):
        adj_stim_lr_con_1_3[ind] = rand_10_pos[rand_pos_ind]
for side_indices_ls in [l_1_3, r_1_3]:
    rand_10_pos = sample(['l', 'r'] * 5, 10)
    for rand_pos_ind, ind in enumerate(side_indices_ls):
        adj_stim_lr_con_2_4[ind] = rand_10_pos[rand_pos_ind]
adj_stim_lr = [] # initialize list for storing all to-be-adjusted stimulus positions
# form list that holds to-be-adjusted stimulus positions for all 'context' trials
adj_stim_lr_con = adj_stim_lr_con_1_3[0:10] + adj_stim_lr_con_2_4[0:10] + \
adj_stim_lr_con_1_3[10:20] + adj_stim_lr_con_2_4[10:20]
lower_index_loop = 0 # initialize a counter to be used for stepping through adj_stim_lr_con values
for context_or_no in con_noc_order:
    if context_or_no=="con": # if it's a grouping of 10 context trials
        adj_stim_lr += adj_stim_lr_con[lower_index_loop:lower_index_loop+10] # add 10 values from adj_stim_lr_con
        lower_index_loop += 10 # increase the counter for stepping through adj_stim_con_lr values
    else: # if it's a grouping of 10 no-context trials
        adj_stim_lr += sample(['l', 'r']*5, 10) # randomly sample from 5 'l' values and 5 'r' values

# initialize counter for stepping through to-be-adjusted left/right values
adj_lr_count = 0

# ensure that the first block of trials starts with a 'star'/practice trial
do_star = True

# initialize keyboard that is to be used for trial routines
adjust_keyboard = keyboard.Keyboard()

time_holder_text = visual.TextStim(win=win, name='time_holder_text',
    text=None,
    font='Arial',
    pos=(0, 0), height=0.000001, wrapWidth=None, ori=0, 
    color='grey', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);

# Initialize components for Routine "goodbye"
goodbyeClock = core.Clock()
end_text = visual.TextStim(win=win, name='end_text',
    text='Thank you for your participation.',
    font='Arial',
    units='deg', pos=(0, 0), height=1.5, wrapWidth=26, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "welcome"-------
continueRoutine = True
# update component parameters for each repeat
welcome_keyb.keys = []
welcome_keyb.rt = []
_welcome_keyb_allKeys = []
# keep track of which components have finished
welcomeComponents = [welcome_txt, welcome_keyb]
for thisComponent in welcomeComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
welcomeClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "welcome"-------
while continueRoutine:
    # get current time
    t = welcomeClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=welcomeClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *welcome_txt* updates
    if welcome_txt.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        welcome_txt.frameNStart = frameN  # exact frame index
        welcome_txt.tStart = t  # local t and not account for scr refresh
        welcome_txt.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(welcome_txt, 'tStartRefresh')  # time at next scr refresh
        welcome_txt.setAutoDraw(True)
    
    # *welcome_keyb* updates
    waitOnFlip = False
    if welcome_keyb.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        welcome_keyb.frameNStart = frameN  # exact frame index
        welcome_keyb.tStart = t  # local t and not account for scr refresh
        welcome_keyb.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(welcome_keyb, 'tStartRefresh')  # time at next scr refresh
        welcome_keyb.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(welcome_keyb.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(welcome_keyb.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if welcome_keyb.status == STARTED and not waitOnFlip:
        theseKeys = welcome_keyb.getKeys(keyList=['space'], waitRelease=False)
        _welcome_keyb_allKeys.extend(theseKeys)
        if len(_welcome_keyb_allKeys):
            welcome_keyb.keys = _welcome_keyb_allKeys[-1].name  # just the last key pressed
            welcome_keyb.rt = _welcome_keyb_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in welcomeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "welcome"-------
for thisComponent in welcomeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# the Routine "welcome" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
illusions_trials_loop = data.TrialHandler(nReps=len(ill_type_order) * 21, method='random', 
    extraInfo=expInfo, originPath=-1,
    trialList=[None],
    seed=None, name='illusions_trials_loop')
thisExp.addLoop(illusions_trials_loop)  # add the loop to the experiment
thisIllusions_trials_loop = illusions_trials_loop.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisIllusions_trials_loop.rgb)
if thisIllusions_trials_loop != None:
    for paramName in thisIllusions_trials_loop:
        exec('{} = thisIllusions_trials_loop[paramName]'.format(paramName))

for thisIllusions_trials_loop in illusions_trials_loop:
    currentLoop = illusions_trials_loop
    # abbreviate parameter names if possible (e.g. rgb = thisIllusions_trials_loop.rgb)
    if thisIllusions_trials_loop != None:
        for paramName in thisIllusions_trials_loop:
            exec('{} = thisIllusions_trials_loop[paramName]'.format(paramName))
    
    # ------Prepare to start Routine "target_indicator_1000ms"-------
    continueRoutine = True
    routineTimer.add(1.000000)
    # update component parameters for each repeat
    # randomly select if the adjustable stimulus (and hence, the green 
    # rectangle indicating where the adjustable stimulus will be shown) 
    # goes on the right or left side if it’s a ‘star’/practice trial, otherwise
    # fetch value from adjustable stimulus l/r values list
    if do_star:
        adjust_pos_side = adjust_pos_side = sample(['l', 'r'], 1)[0] # randomly chooses 'l' or 'r'
    else:
        adjust_pos_side = adj_stim_lr[adj_lr_count]
    if adjust_pos_side=='l':
        adjust_pos = left_center
    else:
        adjust_pos = right_center
    target_rectangle.setPos(adjust_pos)
    # keep track of which components have finished
    target_indicator_1000msComponents = [target_rectangle]
    for thisComponent in target_indicator_1000msComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    target_indicator_1000msClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "target_indicator_1000ms"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = target_indicator_1000msClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=target_indicator_1000msClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *target_rectangle* updates
        if target_rectangle.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            target_rectangle.frameNStart = frameN  # exact frame index
            target_rectangle.tStart = t  # local t and not account for scr refresh
            target_rectangle.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(target_rectangle, 'tStartRefresh')  # time at next scr refresh
            target_rectangle.setAutoDraw(True)
        if target_rectangle.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > target_rectangle.tStartRefresh + 1-frameTolerance:
                # keep track of stop time/frame for later
                target_rectangle.tStop = t  # not accounting for scr refresh
                target_rectangle.frameNStop = frameN  # exact frame index
                win.timeOnFlip(target_rectangle, 'tStopRefresh')  # time at next scr refresh
                target_rectangle.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in target_indicator_1000msComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "target_indicator_1000ms"-------
    for thisComponent in target_indicator_1000msComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    
    # ------Prepare to start Routine "illusions_trials"-------
    continueRoutine = True
    # update component parameters for each repeat
    # get value specifying if smaller circles/inward arrows go on
    # right or left side
    l_r_val = trial_spec[block_count][1][first_last_10_count][1][l_r_count]
    # get value specifying if context or no context is to be used
    noc_con = trial_spec[block_count][1][first_last_10_count][0]
    
    # set the reference line/circle's position to be the side where the
    # adjustable stimulus isn't placed (the 'adjust_pos' variable is assigned
    # its value in the code snippet for the target indicator routine)
    if adjust_pos_side=='r':
        reference_pos = left_center
    else:
        reference_pos = right_center
    
    # initialize/reset has_changed variable, which indicates if participant has given a response
    has_changed = False
    
    # initialize/reset counter for the number of adjustments the participant makes
    adjustment_counter = 0
    
    # store routine start time
    routine_start_t = t
    
    # flush keyboard input, to prevent pre-presentation keyboard presses from being registered
    _ = adjust_keyboard.getKeys(keyList=['up', 'down', 'space'], waitRelease=False)
    
    # if a 'star'/practice trial is to be done:
    if do_star:
        adjust_size = adjust_star_init_min + (adjust_star_init_max - adjust_star_init_min) * random()
        ref_size = ref_star_diam
        if adjust_pos_side == 'l':
            adjust_star.pos = left_center
            static_star.pos = right_center
        elif adjust_pos_side == 'r':
            adjust_star.pos = right_center
            static_star.pos = left_center
    # if ebbinghaus block:
    elif trial_spec[block_count][0]=="ebbing":
        # set adjustable circle diameter variable initial value
        adjust_size = adjust_circle_init_min + (adjust_circle_init_max - adjust_circle_init_min) * random()
        adjust_circle.pos = adjust_pos
        ref_size = ref_circ_diam # reference size is reference circle diameter
        if noc_con == "con": # if context is to be used
            if l_r_val == "l": # if small circles should be on the left
                small_circ_pos = small_circ_pos_l
                large_circ_pos = large_circ_pos_r
            elif l_r_val == "r":
                small_circ_pos = small_circ_pos_r
                large_circ_pos = large_circ_pos_l
    # if muller-lyer block:
    elif trial_spec[block_count][0]=="muller":
        # set adjustable arrow length to random initial value (using values specified in begin_exp code)
        adjust_size = adjust_line_init_min + (adjust_line_init_max - adjust_line_init_min) * random()
        # get difference between reference line's length and adjustable line's length
        adjust_diff_len = adjust_size - ref_line_len
        ref_size = ref_line_len # reference size is reference line length
        adjust_line.pos = adjust_pos
        if noc_con == "con": # if context is to be used
            if l_r_val == "l": # if inward arrows should be on the left
                if adjust_pos_side == "l": # if the adjustable stimulus goes on the left
                    # in this case, inward lines should go on the left, and be adjustable
                    adjust_con_pos = deepcopy(inward_line_pos_l)  # adjustable context positions list
                    static_con_pos = outward_line_pos_r # static context positions list
                    adjust_con_oris = inward_oris
                    static_con_oris = outward_oris
                elif adjust_pos_side == "r": # if the adjustable stimulus goes on the right
                    adjust_con_pos = deepcopy(outward_line_pos_r)
                    static_con_pos = inward_line_pos_l
                    adjust_con_oris = outward_oris
                    static_con_oris = inward_oris
            elif l_r_val == "r": # if inward arrows should be on the right
                if adjust_pos_side == "l": # if the adjustable stimulus goes on the left
                    # in this case, inward lines should go on the right, and be static
                    adjust_con_pos = deepcopy(outward_line_pos_l)
                    static_con_pos = inward_line_pos_r
                    adjust_con_oris = outward_oris
                    static_con_oris = inward_oris
                elif adjust_pos_side == "r":
                    adjust_con_pos = deepcopy(inward_line_pos_r)
                    static_con_pos = outward_line_pos_l
                    adjust_con_oris = inward_oris
                    static_con_oris = outward_oris
            # shift initial positions of adjustable small arrow components based on adjustable arrow length's initial value
            for pos in adjust_con_pos:
                # if position of small arrow component has a greater x-coordinate than the position of the
                # long component, increase the small arrow component's x-coordinate by half the adjust_diff_len
                # value (bringing the small component further away if adjust_diff_len is positive,
                # closer if the value is negative)
                if pos[0] > adjust_pos[0]:
                    pos[0] += adjust_diff_len / 2
                else:
                    pos[0] -= adjust_diff_len / 2
            has_changed = False
    
    # keep track of which components have finished
    illusions_trialsComponents = [time_holder_text]
    for thisComponent in illusions_trialsComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    illusions_trialsClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "illusions_trials"-------
    while continueRoutine:
        # get current time
        t = illusions_trialsClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=illusions_trialsClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # check if participant has pressed up/down/space
        waitOnFlip = False
        theseKeys = adjust_keyboard.getKeys(keyList=['up', 'down', 'space'], waitRelease=False)
        if len(theseKeys):
            if 'up' in theseKeys and ((trial_spec[block_count][0]=="ebbing" and adjust_size < adjust_circle_max) or (trial_spec[block_count][0]=="muller" and adjust_size < adjust_line_max) or (do_star and adjust_size < adjust_star_max)):
                adjust_change = 0.05 * adjust_size
                adjust_size += adjust_change # increase diameter/arrow length by 10%
                has_changed = True
                adjustment_counter += 1
            elif 'down' in theseKeys and ((trial_spec[block_count][0] == "ebbing" and adjust_size > adjust_circle_min) or (trial_spec[block_count][0] == "muller" and adjust_size > adjust_line_min) or (do_star and adjust_size > adjust_star_min)):
                adjust_change = -0.05 * adjust_size
                adjust_size += adjust_change # decrease diameter/arrow length by 10%
                has_changed = True
                adjustment_counter += 1
            elif 'space' in theseKeys:
                continueRoutine = False # end routine
        
        # if a 'star'/practice trial is to be done
        if do_star:
            # draw reference star
            static_star.draw()
            # draw adjustable star
            adjust_star.size = (adjust_size, adjust_size)
            adjust_star.draw()
        # if in a muller-lyer block at the moment (and no 'star' trial is to be done)
        elif trial_spec[block_count][0]=="muller":
            # adjust positions of adjustable small arrow components, if participant has given input
            if noc_con == "con" and has_changed:
                for pos in adjust_con_pos:
                    if pos[0] > adjust_pos[0]:
                        pos[0] += adjust_change/2 # change position by half as much that the long component's width changes
                    else:
                        pos[0] -= adjust_change/2
                has_changed = False
            # draw context components
            if noc_con == "con":
                adjust_line.width, static_line.width = [0.1851 * ref_line_len] * 2
                for pos_index in range(len(adjust_con_pos)):
                    adjust_line.pos = adjust_con_pos[pos_index]
                    adjust_line.ori = adjust_con_oris[pos_index]
                    adjust_line.draw()
                for pos_index in range(len(static_con_pos)):
                    static_line.pos = static_con_pos[pos_index]
                    static_line.ori = static_con_oris[pos_index]
                    static_line.draw()
            # draw adjustable long arrow component
            adjust_line.width = adjust_size
            adjust_line.pos = adjust_pos
            adjust_line.ori = 0
            adjust_line.draw()
            # draw static long arrow component
            static_line.width = ref_size
            static_line.pos = reference_pos
            static_line.ori = 0
            static_line.draw()
        # if in an ebbinghaus block at the moment
        elif trial_spec[block_count][0]=="ebbing":
            if noc_con == "con": # if context is to be used
                # draw large context circles
                static_circle.size = [large_circ_diam, large_circ_diam]
                for pos in large_circ_pos:
                   static_circle.pos = pos
                   static_circle.draw()
                # draw small context circles
                static_circle.size = [small_circ_diam, small_circ_diam]
                for pos in small_circ_pos:
                   static_circle.pos = pos
                   static_circle.draw()
            # draw reference circle
            static_circle.size = [ref_size, ref_size]
            static_circle.pos = reference_pos
            static_circle.draw()
            # draw adjustable circle
            adjust_circle.size = (adjust_size, adjust_size)
            adjust_circle.draw()
        
        
        # *time_holder_text* updates
        if time_holder_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            time_holder_text.frameNStart = frameN  # exact frame index
            time_holder_text.tStart = t  # local t and not account for scr refresh
            time_holder_text.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(time_holder_text, 'tStartRefresh')  # time at next scr refresh
            time_holder_text.setAutoDraw(True)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in illusions_trialsComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "illusions_trials"-------
    for thisComponent in illusions_trialsComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    illusions_trials_loop.addData('trial.duration', t - routine_start_t)
    illusions_trials_loop.addData('trial.adjustable_side', adjust_pos_side)
    illusions_trials_loop.addData('trial.adjustments_counts', adjustment_counter)
    if do_star: # if a 'star'/practice trial just finished
        do_star = False # tell PsychoPy that the next trial is not to be a 'star' trial
        illusions_trials_loop.addData('trial.illusion_type', 'star')
        illusions_trials_loop.addData('trial.context_used', 'no_context')
        illusions_trials_loop.addData('trial.small_circle_inward_arrow_side', 'no_context')
    else: # if a 'non-star'/test trial just finished
        l_r_count += 1 # bump context left/right values counter, i. e. proceed to next trial
        adj_lr_count += 1 # bump to-be-adjusted-stimulus left/right values
        illusions_trials_loop.addData('trial.illusion_type', trial_spec[block_count][0])
        illusions_trials_loop.addData('reference_minus_adjustable_size', ref_size - adjust_size)
        if noc_con == 'con':
            illusions_trials_loop.addData('trial.context_used', 'context')
            illusions_trials_loop.addData('trial.small_circle_inward_arrow_side', l_r_val)
        else:
            illusions_trials_loop.addData('trial.context_used', 'no_context')
            illusions_trials_loop.addData('trial.small_circle_inward_arrow_side', 'no_context')
    
    if l_r_count >= 10: # if gone through all l/r values
        first_last_10_count += 1 # go from context trials to no context trials or vice versa
        l_r_count = 0 # reset l/r values counter
        if first_last_10_count >= 2: # if already done both context and no context trials
            block_count += 1 # proceed to next block of trials
            do_star = True # ensure that the next block of trials starts with a 'star'/practice trial
            first_last_10_count = 0 # reset no context/context specification trial
    
    # this is outside of the above if/else statements because the data to be saved is the same for
    # 'star'/'non-star' trials
    illusions_trials_loop.addData('reference_minus_adjustable_size', ref_size - adjust_size)
    
    # the Routine "illusions_trials" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()
    
# completed len(ill_type_order) * 21 repeats of 'illusions_trials_loop'


# ------Prepare to start Routine "goodbye"-------
continueRoutine = True
routineTimer.add(10.000000)
# update component parameters for each repeat
# keep track of which components have finished
goodbyeComponents = [end_text]
for thisComponent in goodbyeComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
goodbyeClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "goodbye"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = goodbyeClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=goodbyeClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *end_text* updates
    if end_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        end_text.frameNStart = frameN  # exact frame index
        end_text.tStart = t  # local t and not account for scr refresh
        end_text.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(end_text, 'tStartRefresh')  # time at next scr refresh
        end_text.setAutoDraw(True)
    if end_text.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > end_text.tStartRefresh + 10-frameTolerance:
            # keep track of stop time/frame for later
            end_text.tStop = t  # not accounting for scr refresh
            end_text.frameNStop = frameN  # exact frame index
            win.timeOnFlip(end_text, 'tStopRefresh')  # time at next scr refresh
            end_text.setAutoDraw(False)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in goodbyeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "goodbye"-------
for thisComponent in goodbyeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
if eyetracker:
    eyetracker.setConnectionState(False)
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
