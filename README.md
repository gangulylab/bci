# bci

<br>

## Installation
To download or clone the full repository follow use this [link](https://github.com/gangulylab/bci.git).

There are a few dependencies that you must have installed to run bci:

* [Matlab](https://www.mathworks.com/)
* [Matlab Arduino Support](https://www.mathworks.com/hardware-support/arduino-matlab.html)
* [Psychtoolbox](http://psychtoolbox.org/)

<br>

## BCI coding environment

---

### General Principles:

1. Blackrock neural acquisition system and API (cbmex) reads in neural data.
2. Psychtoolbox controls the drawing / timing through cbmex files
3. Matlab code is used to control the task flow, signal processing, data saving, etc...

To run an experiment, type 

```matlab
ExperimentStart(task_name, subject, control_mode, blackrock, debug)
```

where

* *task_name* is a string containing the name of a valid task,
* *subject* is a string containing the id of a subject (use 'test' or 'Test' to avoid saving a ton of useless data), 
* *control_mode* is an integer {1 - mouse position control, 2 - mouse joystick control, 3 - full kalman filter, 4 - velocity kalman filter}, 
* *blackrock* is a flag, if true it attempts to use the Blackrock API to acquire neural data.
* *debug* is a flag, when true it invokes a debugging environment in which the screen is a bit smaller, etc...

An Experiment can paused by pressing *'p'* during the task flow. To exit the task gracefully, press the *Escape* key at the pause screen.

There are two important folders that this code base relies on:

1. **persistence**
2. **Data**

The paths to the folders are defined in ExperimentStart.

**persistence** is a folder that contains the following files:

* ch_stats: mean and variance of the signal on each channel / trial
* feature_stats:  mean and variance of the signal of each feature / trial
* kf_params:  reduced dimensionality 2D kalman filter parameters / trial
* full_kf_params:  full dimensionality 2D kalman filter parameters / experiment
* kf_params_1D:  reduced dimensionality 1D kalman filter parameters / trial
* full_kf_params_1D:  full dimensionality 1D kalman filter parameters / experiment

Each of these files are saved after each trial (or experiment). This way, if these parameters are changed during an experiment, the most recent parameters are present for future experiments.

**Data** is a folder that contains the data from each trial in a separate file.

---

### Task Flow
This is an example of how the functions call each other in an example task. Starred functions, are specific to each task module. Also note, the function names listed here might be different depending on whether a 1D or 2D kalman filter is used, whether a clicker is being used, whether typing is being used, etc...

---> ExperimentStart.............................# Initializes Experiment

------> *GetParams..................................# Loads Task Specific Parameters

------> GetNeuroParams.........................# Loads Neural Processing Parameters

------> LoadFeatureMask........................# Loads FeatureMask if using one

------> RunBaseline..................................# Runs at Params.UpdateRate, no other task flow

---------> NeuroPipeline............................# Loads and processes neural data

------> *RunTask........................................# Loads Clicker, Fits dimensionality reduction and kalman filter if needed

---------> *RunLoop....................................# Sets up block, saves data and persistent params after each trial

------------> *RunTrial...................................# Main Loop (Updates at Params.UpdateRate)

---------------> NeuroPipeline.......................# Loads and processes neural data

---------------> UpdateClicker .......................# Uses neural data to discretely decode clicking

---------------> UpdateCursor........................# Uses neural data to continuously decode velocity and apply dynamics

---> ExperimentStop................................# Ends an experiment. Can be called early by user.

<br>

## Matlab Modules

1. task_modules - control task flow in different environments
	* CenterOut
	* CenterOut1D
	* GridTask
	* RadialTask
	* RadialTyping
	* RadialTypingMultiClick
2. task_helpers - general tools to help with task flow across task modules
3. neuro - controls neural processing
4. arduino - controls arduino for time-sync and exo system
5. kalman_filter - continuous decoder, used to control cursor pos/vel 
6. clicker - discrete decoder, used to discretely control task flow
7. typing_env - controls typing parameters: keyboard layout, character/word selection, etc...
8. exo_control - controls exo system

Each module will be explained in more detail below.

---

### task_modules

Each task module contains at least four files:

* GetParams.m
* RunTask.m
* RunLoop.m
* RunTrial.m

*GetParams.m* is used to set the parameters for a given task / experimental session. It is meant to be updated before each experimental run depending on what the experimenter wants to test. It returns a struct, *Params* that is meant to be passed to a large range of functions so that those functions are modular, but can remain static (no task coding) once a task is set up. Parameters fall into rough sections. These sections are heavily commented. I will expand on a few important ones here:

* Control - this section contains a few highly relevant parameters that change often. Some parameters pertain to cursor control (e.g., degree of assistance, type of cursor control adaptation, whether the cursor is "reset" after each trial) while others pertain to collection of neural features (e.g., spatial filtering, baseline time, feature mask)

*RunTask.m*, *RunLoop.m*, and *RunTrial.m* control the task flow for each task. Importantly, *RunTrial.m* contains a main for loop that occasionally (at Params.UpdateRate Hz) calls the NeuroPipeline (see neuro module).

---> *RunTask........................................# Loads Clicker, Fits dimensionality reduction and kalman filter if needed

------> *RunLoop....................................# Sets up block, saves data and persistent params after each trial

---------> *RunTrial...................................# Main Loop (Updates at Params.UpdateRate)

---

### task_helpers

This module contains helper functions to control task flow (e.g., *ExperimentPause.m*) and to control data structures (e.g., *UpdateCursor.m*) across many tasks. A few of these are outlined below.

---> CheckPause ..............................# Polls keyboard, checks for *'p'*

------> ExperimentPause .................# Pause Screen

---------> ExperimentStop ................# Gracefully shuts down experiment

<br>

---> UpdateClicker .......................# Uses neural data to discretely decode clicking

------> Cicker.Func .......................# sets function in RunTask

<br>

---> UpdateCursor........................# Uses neural data to continuously decode velocity and apply dynamics

------> OptimalCursorUpdate.......# Defines optimal velocity vector

------> UpdateRmlKF......................# updates kalman filter weights in a recursive maximum likelihood framework

---

### neuro

This module acts on the *Neuro* struct. It contains a neural processing pipeline that can be modified by changing parameters in the *GetParams.m* file. The main function here is *NeuroPipeline.m*, which collects neural data (through the Blackrock API), performs signal preprocessing, computes neural features, etc...

*GetNeuroParams.m* contains parameters that change how neural data is processed, including details of neural filter bank, flags for processing/saving steps (whether to save raw, processed, reduced, etc... data) and defines the edges of filters, which filtered signals will get averaged together, which are computed through a hilbert vs. traditional power estimate, etc...

---> GetNeuroParams...................# Loads Neuro processing parameters

<br>

---> NeuroPipeline........................# collects and processes neural data / note that these steps are turned on/off in GetNeuroParams

------> ReadBR................................# collects neural data through BlackRock API

------> RefNeuralData....................# references neural data

------> UpdateChStats....................# updates mean and variance of signal on each channel

------> ZscoreChannels...................# z-scores signal using channel statistics

------> ApplyFilterBank...................# filters signal from each channel across multiple freq bands

------> UpdateNeuroBuf.................# buffers low freqs

------> UpdateFeatureStats............# updates mean and variance of each feature

------> ZscoreFeatures....................# z-scores features using feature statistics

------> VelToNeuralFeatures...........# generates fake neural features for testing / debugging

------> ApplyFeatureMask...............# not an independent function - reduces dimensionality by applying binary mask to features

------> Neuro.DimRed.F...................# reduces dimensionality by applying function (loaded in RunTask)

---

### arduino

This module contains a few helper functions for interacting through Matlab with an Arduino device. This is useful for time synchronization (sending arduino pulses), and for control of the planar exo system.

---

### kalman_filter

This module contains code for loading, fitting, and updating parameters of a Kalman Filter. It contains separate functions for 1D and 2D Kalman Filters. It acts on the *KF* struct.

#### Loading/Initializing Kalman Filter
---> GetParams..............................# Sets experimental parameters

------> LoadKF2dDynamics............# Sets kalman Filter A and W matrices

<br>

---> FitKF........................................# Either loads data and fits new kalman C and Q matrices or loads existing matrices from persistence folder

#### Running Kalman Filter
---> UpdateCursor........................# Applies current kalman filter to cursor state  

------> UpdateRmlKF......................# updates kalman filter weights in a recursive maximum likelihood framework  
*In reality, sufficient statistics for generating a kalman filter are tracked and updated at each time bin. Then the C and Q matrices are computed from these sufficient statistics.*

#### Saving Kalman Filter
---> SavePersistence....................# Saves current current kalman filter after each trial (kf_params or kf_params_1D)

---> ExperimentStop....................# Updates adapted parameters within a full kalman filter after each experiment ( full_kf_params or full_kf_params_1D). Leaves unadapted parameters alone.

---

### clicker

This module contains classification models and functions that apply classification to neural data and return discrete states. These functions act on the *Clicker* struct. Currently, *click_classifier.m* returns

---

### typing_env

This module contains functions to control the flow of a typing environment that is overlaid on the traditional tasks. For example, character selections, word selections, keyboard displays, undo functionality, etc... These functions primarily act on the *Keyboard* struct in the *Params* struct.

---

### exo_control

This module contains functions that support interfacing with and controling a custom planar exoskeleton device. These functions primarily act on the *Arduino* struct in the *Params* struct.

---


