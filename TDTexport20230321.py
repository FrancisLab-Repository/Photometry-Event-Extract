# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 15:14:03 2021

@author: Chase
"""

"____ User Inputs ___________________________________________________________"
photometryOffset = 1000
offset = 0 + photometryOffset # how much the video and photometry recording is offset in msec
BehStartTime = 20000 # time behavior starts relative to photometry rec in msec
LengthOfBaseline = 10 #set length of baseline for Z-Scoring in sec
#SampleRate = 1017.2526245117188 #found in data.streams._405G

"____ Package Import ________________________________________________________"

# Imports for file management
import os
import datetime
import time
from openpyxl import load_workbook
from tdt import read_block, epoc_filter
from datetime import datetime

# Imports for computation
import numpy as np
import pandas as pd
import scipy as sp
import scipy.stats as stats
from scipy import signal
from pandas import Series, DataFrame
from scipy.signal import find_peaks

# Imports for user interface
import tkinter as tk
from tkinter import simpledialog as dial
from tkinter import filedialog
from tkinter.filedialog import askdirectory
from tkinter import *
from pathlib import Path


"____ Import File ___________________________________________________________"

def GetBlockPath():
    BLOCKPATH = askdirectory()
    return([BLOCKPATH])

root = Tk()
BLOCKPATH = askdirectory(title='Select tank folder')
root.destroy()

SampleRate = data.streams._405G.fs

# extract and define signals
data = read_block(os.path.abspath(BLOCKPATH))
isosb = DataFrame(data.streams._405G.data) #, columns = {'405'})

combined = {'time':isosb.index / SampleRate, 
            '405':data.streams._405G.data,
            '470':data.streams._470G.data}
combined = pd.DataFrame(combined)

print('Extracted ' + BLOCKPATH)

"____ Pre-processing ________________________________________________________"

# calc dF/F and align to behavioral start time
combined['dF'] = combined['470'] - combined['405']
combined['dFF'] = combined['dF'].div(combined['405'].values)

StartingIter = BehStartTime + offset
dFF = combined.drop(combined.index[0:StartingIter], axis = 0, inplace = False)
dFF.reset_index(inplace = True, drop = True)

# Find slope and normalize
dFFdetrend = signal.detrend(dFF['dFF'])
dFFNorm = pd.DataFrame(dFFdetrend, columns = ['dFF'])

# extract only dF/F and time and normalize baseline to zero
dFFNorm.rename(columns={'dFF':'dFFnorm'}, inplace=True)
dFFNormComb = pd.concat([dFF, dFFNorm], axis = 1, join='inner')
dFFtimeNorm = dFFNormComb.drop(['405', '470', 'dF', 'dFF'], 
                               axis = 1, inplace = False)

#set datetime for file path
timestr = time.strftime("%Y%m%d")

#set file path
Path(BLOCKPATH + '/EventFiles' + '_' + timestr).mkdir(parents=True, 
                                                      exist_ok=True)
SaveDir = BLOCKPATH + '/EventFiles' + '_' + timestr + '/'

# save raw as csv file
dFFtimeNorm.to_csv(SaveDir + 'dFF.csv', index=False)
print('...dFF file saved')

"____ Event Extraction ______________________________________________________"
# Select file with event start and end
BehavEvents = filedialog.askopenfilename(
    title="Select behavior events .csv file")
Events = pd.read_csv(BehavEvents) #, names = ['Start', 'End'])

dFFtimeNorm['time']= dFFtimeNorm['time']
Events2 = Events.multiply(SampleRate) 
dFFtimeNorm = dFFtimeNorm.reset_index(drop = True)

# iterate through behavioral event ranges
EventExtract2 = []
for i in range(len(Events2)):
    TimeEpochs = dFFtimeNorm.iloc[np.r_[Events2.iloc[i,0]:Events2.iloc[i,1]]]
    dFFevents = TimeEpochs['dFFnorm']
    EventExtract2.append(dFFevents.reset_index(drop = True))

dFFeventExtract = pd.concat(EventExtract2, axis=1) # convert list to dataframe
#.fillna('') 

"____ Z-score Events ________________________________________________________"
#extract baseline for z-scoring
BaselineForZscore = LengthOfBaseline * 1000
ZscoreBaseline = dFFeventExtract.drop( 
    dFFeventExtract.index[BaselineForZscore:len(dFFeventExtract)],
                                      axis = 0, inplace = False)
# z-score
MeanBaselineForZ = ZscoreBaseline.mean(axis = 0)
MeanBaselineForZdf = pd.DataFrame(MeanBaselineForZ).T

StdBaselineForZ = np.std(ZscoreBaseline, axis = 0)
StdBaselineForZdf = pd.DataFrame(StdBaselineForZ).T

ZscoreEventsSub = dFFeventExtract.subtract(MeanBaselineForZ, axis = 1)
ZscoreEvents = (ZscoreEventsSub / StdBaselineForZdf.loc[0])
#ZscoreEvents = ZscoreEvents #.fillna('')

# rejoin time column
timeColumn = pd.DataFrame(range(len(ZscoreEvents)), columns = ['time'])
timeColumn = timeColumn.divide(SampleRate)

# z-score
dFFeventsZ = pd.concat([timeColumn, ZscoreEvents], axis = 1, join='inner')

# dFF
dFFevents = pd.concat([timeColumn, dFFeventExtract], axis = 1, join='inner')

# Save events

dFFeventsZ.to_csv(SaveDir + 'Events_Zscore.csv', index=False)
dFFevents.to_csv(SaveDir + 'Events.csv', index=False)

dFFeventsZavg = dFFeventsZ.groupby(np.arange(len(dFFeventsZ))//200).mean()
dFFeventsZavg.to_csv(SaveDir + 'Events_Zscore200msec.csv')

print('...events saved')

"____ Peak Detection ________________________________________________________"
# set parameters
width_event = 50 #width of peak in msec
dist_event = 100 #minimum time between each event in msec

#median and mean absolute deviation
FullSignal = dFFtimeNorm.iloc[:,1]
Median = np.median(FullSignal)
FullSignal_series = pd.Series(FullSignal)
MedDev = stats.median_abs_deviation(FullSignal)
#MedDev = FullSignal_series.mad() #this is actually mean deviation
threshold = Median + ((MedDev) * (StdBaselineForZ[1]*1.5))
                      #*1.28151) #20% both tails, 10% top tails #25% of events is 1.15034

#finding peaks from entire trace
Peaks , _ = find_peaks(FullSignal, distance = dist_event,
                            prominence = threshold,
                            height = -1, width = width_event)
Peaks_df = pd.DataFrame(Peaks)

#event extraction
Peaks_df_events = []
for i in range(len(Events2)):
    Epochs = (Peaks[(Peaks[:len(Peaks)] >= Events2.iloc[i,0]) &
                  (Peaks[:len(Peaks)] < Events2.iloc[i,1])]) - Events2.iloc[i,0]
    Epochs2 = pd.Series(Epochs) #making series so can concatenate later
    Peaks_df_events.append(Epochs2)

Peaks_df_events = pd.concat(Peaks_df_events, axis=1) #make dataframe instead of list 
Peaks_df_events.to_csv(SaveDir + 'EventPeaks.csv')

print('...peaks saved')