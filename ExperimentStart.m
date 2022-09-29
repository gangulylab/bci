function ExperimentStart(Task,Subject,ControlMode,BLACKROCK,DEBUG)
% function ExperimentStart(Subject,ControlMode)
% Task - string for task
% Subject - string for the subject id
% ControlMode - [1,2,3,4] for mouse pos, mouse vel, pos/vel kalman, vel
%   kalman
% BLACKROCK - [0,1] if 1, collects, processes, and saves neural data
% DEBUG - [0,1] if 1, enters DEBUG mode in which screen is small and cursor
%   remains unhidden

%% Clear All and Close All
clearvars -except Task Subject ControlMode BLACKROCK DEBUG
clearvars -global
clc
warning off

AssertOpenGL;
KbName('UnifyKeyNames');

% Deal with inputs
valid_tasks = {...
    'CenterOut',...
    'CenterOut1D',...  
    'GridTask',...
    'GridTyping',...
    'RadialTask',...
    'RadialTyping',...
    'RadialTypingMultiClick',...
    'ExoControl1D',...
    'MultiStateDiscrete',...
    'RadialTaskMultiStateDiscrete',...
    'RadialTaskMultiStateDiscreteArrow',...
    'Robot',...
    'Robot3D',...
    'Robot3DArrow',...
    'RobotDistance',...
    'HandImagined',...
    'Gripper',...
    'RobotRR',...
    'RobotIntention',...
    'RobotR2GAuto',...
    'RealRobot3D',...
    'RobotReachStop',...
    'RobotStop',...
    'RobotR2GModeSwitch',...
    'RealRobot3D',...
    'RealRobotBatch',...
    'Robot3DPath',...
    'RobotEngagement',...
    'RobotLateralR2G',...
    'Robot3DWaypoint',...
    'RobotEtchASketch',...
    'RobotBeta',...
    'Robot3DRot',...
    'HandOnline',...
    'RealRobotBetaStop',...
    'RealRobotPath',...
    'RealRobotBatch_RL'};
assert(any(strcmp(Task,valid_tasks)), 'Unknown task')
if ~exist('Subject','var'), Subject = 'Test'; DEBUG = 1; end
if ~exist('ControlMode','var'), ControlMode = 2; end
if ~exist('BLACKROCK','var'), BLACKROCK = 0; end
if ~exist('DEBUG','var'), DEBUG = 0; end

if strcmpi(Subject,'Test'), Subject = 'Test'; end % ignore case

%% Set up path for running task dependent code
if IsWin,
    %homedir = 'C:\Users\ganguly-lab2\Documents\bci_data_test';
    %projectdir = fullfile('C:\Users\ganguly-lab2\Documents\MATLAB\bci');
    homedir='C:\Users\Nikhlesh\Documents\bci_data_test';
    projectdir = fullfile('C:\Users\Nikhlesh\Documents\GitHub\bci');
elseif IsOSX,
    homedir = '/Users/daniel/';
    projectdir = '/Users/daniel/Projects/bci/';
else,
    homedir = '~';
    projectdir = '~/Projects/bci/';
%     projectdir = '/home/sarah/Documents/bci';

end

% reset path
rmpath(genpath(fullfile(projectdir)));
% add non-task modules
addpath(genpath(fullfile(projectdir,'arduino')));
addpath(genpath(fullfile(projectdir,'audio')));
addpath(genpath(fullfile(projectdir,'clicker')));
addpath(genpath(fullfile(projectdir,'kalman_filter')));
addpath(genpath(fullfile(projectdir,'neuro')));
addpath(genpath(fullfile(projectdir,'task_helpers')));
addpath(genpath(fullfile(projectdir,'typing_env')));
addpath(genpath(fullfile(projectdir,'exo_control')));
% add task module
addpath(genpath(fullfile(projectdir,'task_modules',Task)))

%% Retrieve Parameters from Params File
Params.Task = Task;
Params.Subject = Subject;
Params.ControlMode = ControlMode;
Params.BLACKROCK = BLACKROCK;
Params.DEBUG = DEBUG;
Params.Homedir = homedir;
Params = GetParams(Params);
Params = GetNeuroParams(Params);
Params = LoadFeatureMask(Params);

%% Initialize Sound
% 
% InitializePsychSound(1);
% Params.sound_nrchannels = 2;                                                         % Number of channels and Frequency of the sound
% Params.sound_freq = 48000;
% Params.sound_pahandle = PsychPortAudio('Open', [], 1, 1, Params.sound_freq, Params.sound_nrchannels);
% PsychPortAudio('Volume', Params.sound_pahandle, 0.5);                                % Set the volume to half for this demo
% 
% Params.beepHigh    = MakeBeep(1000, 0.1, Params.sound_freq);
% % Params.beepLow     = MakeBeep(500, 0.25, Params.sound_freq);
% freq = 500;
% duration = .1;
% samplingRate = 48000;
% 
% s = sin(2*pi*freq*(0:duration*samplingRate-1)/samplingRate);
% 
% Params.beepLow  = [s, zeros(1,1000),s];


%% Verbose
if Params.DEBUG, Params.Verbose = true;
else, Params.Verbose = false;
end

%% Data Saving
now = datetime; % get today's date
Params.YYYYMMDD = sprintf('%i',yyyymmdd(now));

H = now.Hour;
M = now.Minute;
S = round(now.Second);

Params.HHMMSS = sprintf('%02i%02i%02i',H,M,S);

if  strcmp(Task(1:4), 'Real')
    Params.udp = udpport("LocalPort", 43210);
    Params.pythonPort = 5006;
    write(Params.udp, [0,6,H,M,S,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);   
end


% if Subject is 'Test' or 'test' then can write over previous test
if strcmpi(Params.Subject,'Test'),
    Params.YYYYMMDD = 'YYYYMMDD';
    Params.HHMMSS = 'HHMMSS';
end

% create folders for saving
Params.Datadir = fullfile(Params.Homedir,'Data',...
    Params.Subject,Params.YYYYMMDD,Params.Task,Params.HHMMSS);
mkdir(Params.Datadir);

% persistence folder
Params.Persistencedir = fullfile(Params.Homedir,'persistence');

%% Initialize Blackrock System
if BLACKROCK,
    if IsWin,
        addpath('C:\Program Files (x86)\Blackrock Microsystems\Cerebus Windows Suite')
    elseif IsLinux,
        addpath('/usr/local/CereLink');
    end
    cbmex('close'); % always close
    cbmex('open'); % open library
    cbmex('trialconfig', 1); % empty the buffer
else,
    Params.GenNeuralFeaturesFlag = true;
end

% load channel layout file and store in params
f = load(fullfile('gridmaps',Params.ChMapFile));
Params.ChMap = f.ecog_grid;

%% Initialize Sync to Blackrock
if Params.ArduinoSync,
    Params.ArduinoPtr = arduino;
    Params.ArduinoPin = 'D13';
    writeDigitalPin(Params.ArduinoPtr, Params.ArduinoPin, 0); % make sure the pin is at 0
    PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,20);
end

%% Neural Signal Processing
% create neuro structure for keeping track of all neuro updates/state
% changes
Neuro.Blackrock         = Params.BLACKROCK;
Neuro.GenNeuralFeaturesFlag = Params.GenNeuralFeaturesFlag;
Neuro.ZscoreRawFlag     = Params.ZscoreRawFlag;
Neuro.ZscoreFeaturesFlag= Params.ZscoreFeaturesFlag;
Neuro.NumFeatureBins    = Params.NumFeatureBins;
Neuro.DimRed            = Params.DimRed;
Neuro.CLDA              = Params.CLDA;
Neuro.SaveProcessed     = Params.SaveProcessed;
Neuro.SaveRaw           = Params.SaveRaw;
Neuro.FilterBank        = Params.FilterBank;
Neuro.NumChannels       = Params.NumChannels;
Neuro.BufferSamps       = Params.BufferSamps;
Neuro.BadChannels       = Params.BadChannels;
Neuro.ReferenceMode     = Params.ReferenceMode;
Neuro.NumPhase          = Params.NumPhase;
Neuro.NumPower          = Params.NumPower;
Neuro.NumBuffer         = Params.NumBuffer;
Neuro.NumHilbert        = Params.NumHilbert;
Neuro.NumFeatures       = Params.NumFeatures;
Neuro.LastUpdateTime    = GetSecs;
Neuro.UpdateChStatsFlag = Params.UpdateChStatsFlag;
Neuro.UpdateFeatureStatsFlag = Params.UpdateFeatureStatsFlag;
Neuro.ChMap             = Params.ChMap;
Neuro.SpatialFiltering  = Params.SpatialFiltering;
Neuro.FeatureMask       = Params.FeatureMask;
Neuro.FeatureBufferSize = Params.FeatureBufferSize;
Neuro.SmoothDataFlag = Params.SmoothDataFlag;
Neuro.biLSTMFlag = Params.biLSTMFlag;
Neuro.LSTMBufferSize = Params.LSTMBufferSize;
Neuro.SaveLSTMFeatures = Params.SaveLSTMFeatures;
Neuro.lpFilt = Params.lpFilt;

% initialize filter bank state
for i=1:length(Params.FilterBank)
    Neuro.FilterBank(i).state = [];
end

% initialize feature buffer
Neuro.FeatureDataBuffer = zeros(Neuro.FeatureBufferSize,Neuro.NumFeatures*Neuro.NumChannels);
Neuro.FilteredFeatures = zeros(Neuro.NumFeatures*Neuro.NumChannels,1);

% initialize the number of data bins for LSTM model
Neuro.LSTMBuffer = 1e-5*randn(128,Neuro.LSTMBufferSize);
Neuro.LSTMFeatures = 1e-5*randn(Neuro.LSTMBufferSize/10,256);

% initialize stats for each channel for z-scoring
Neuro.ChStats.mean      = zeros(1,Params.NumChannels); % estimate of mean for each channel
Neuro.ChStats.var       = zeros(1,Params.NumChannels); % estimate of variance for each channel
Neuro.ChStats.Idx       = 0;
Neuro.ChStats.BufSize   = Params.ZBufSize * Params.UpdateRate;
Neuro.ChStats.Buf       = cell(1,Neuro.ChStats.BufSize);

% initialize stats for each feature for z-scoring
Neuro.FeatureStats.mean     = zeros(1,Params.NumFeatures*Params.NumChannels); % estimate of mean for each channel
Neuro.FeatureStats.var      = zeros(1,Params.NumFeatures*Params.NumChannels); % estimate of variance for each channel
Neuro.FeatureStats.Idx      = 0;
Neuro.FeatureStats.BufSize  = Params.ZBufSize * Params.UpdateRate;
Neuro.FeatureStats.Buf      = cell(1,Neuro.FeatureStats.BufSize);

% create low freq buffers
Neuro.FilterDataBuf = zeros(Neuro.BufferSamps,Neuro.NumChannels,Neuro.NumBuffer);
if Neuro.NumFeatureBins>1,
    Neuro.NeuralFeaturesBuf = zeros(Neuro.NumFeatures*Neuro.NumChannels,...
        Neuro.NumFeatureBins);
end

% params  setting for LSTM
Params.BaselineRunningFlag=false;



%% Kalman Filter
if Params.ControlMode>=3,
    KF = Params.KF;
    KF.CLDA = Params.CLDA;
    KF.Lambda = Params.CLDA.Lambda;
else,
    KF = [];
end

%% Check Important Params with User
LogicalStr = {'off', 'on'};
IMStr = {'imagined mvmts', 'shuffled imagined mvmts', 'prev mvmts', 'prev adapted'};
DimRedStr = {'PCA', 'FA'};

fprintf('\n\nImportant Experimental Parameters:')
fprintf('\n\n  Task Parameters:')
fprintf('\n    - task: %s', Params.Task)
fprintf('\n    - subject: %s', Params.Subject)
fprintf('\n    - control mode: %s', Params.ControlModeStr)
fprintf('\n    - blackrock: %s', LogicalStr{Params.BLACKROCK+1})
fprintf('\n    - debug mode: %s', LogicalStr{Params.DEBUG+1})
fprintf('\n    - arduino sync: %s', LogicalStr{Params.ArduinoSync+1})

fprintf('\n\n  Neuro Processing Pipeline:')
if Params.GenNeuralFeaturesFlag,
    fprintf('\n    - generating neural features!')
else,
    fprintf('\n    - reference mode: %s', Params.ReferenceModeStr)
    fprintf('\n    - zscore raw: %s', LogicalStr{Params.ZscoreRawFlag+1})
    fprintf('\n    - zscore features: %s', LogicalStr{Params.ZscoreFeaturesFlag+1})
    fprintf('\n    - save raw data: %s', LogicalStr{Params.SaveRaw+1})
    fprintf('\n    - save filtered data: %s', LogicalStr{Params.SaveProcessed+1})
end
fprintf('\n    - dimensionality reduction: %s', LogicalStr{Params.DimRed.Flag+1})
if Params.DimRed.Flag,
    fprintf('\n      - method: %s', DimRedStr{Params.DimRed.Method})
    fprintf('\n      - before clda: %s', LogicalStr{Params.DimRed.InitAdapt+1})
    fprintf('\n      - before fixed: %s', LogicalStr{Params.DimRed.InitFixed+1})
end

fprintf('\n\n  BCI Parameters:')
fprintf('\n    - Imagined Movements: %s', LogicalStr{double(Params.NumImaginedBlocks>0) +1})
fprintf('\n      - initialization mode: %s', IMStr{Params.InitializationMode})
fprintf('\n    - Adaptation Decoding: %s', LogicalStr{double(Params.NumAdaptBlocks>0) +1})
if Params.NumAdaptBlocks>0,
    fprintf('\n      - adapt type: %s', Params.CLDA.TypeStr)
    fprintf('\n      - adapt change type: %s', Params.CLDA.AdaptType)
end
fprintf('\n    - Fixed Decoding: %s', LogicalStr{double(Params.NumFixedBlocks>0) +1})


str = input('\n\nContinue? (''n'' to quit, otherwise continue)\n' ,'s');
if strcmpi(str,'n'),
    fprintf('\n\nExperiment Ended\n\n')
    return
end

%% Initialize Window
Screen('Preference', 'SkipSyncTests', 0);

if strcmp(Task(1:5), 'Robot') || strcmp(Task(1:4), 'Hand') || strcmp(Task(1:4), 'Grip') || strcmp(Task(1:4), 'Real') 
    Params.Center = [0,0,0];
else
    if DEBUG
        [Params.WPTR, Params.ScreenRectangle] = Screen('OpenWindow', 0, 0, [50 0 1750 1000]);
    else
        [Params.WPTR, Params.ScreenRectangle] = Screen('OpenWindow', max(Screen('Screens')), 0);
    end
    Params.Center = [mean(Params.ScreenRectangle([1,3])),mean(Params.ScreenRectangle([2,4]))];

    % Font
    Screen('TextFont',Params.WPTR, 'Arial');
    Screen('TextSize',Params.WPTR, 28);
end

%% Initialze keyboard
typing_tasks = {'RadialTyping','RadialTypingMultiClick', 'GridTyping'};
if any(strcmp(Task,typing_tasks)),
    Params = SetKeyboardParams(Params);
end

%% Start
try
    % Baseline 
    if Params.BaselineTime>0,
        % turn on update stats flags
        Neuro.UpdateChStatsFlag = true;
        Neuro.UpdateFeatureStatsFlag = true;
        Neuro.DimRed.Flag = false;

        % collect data during baseline period
        %Neuro = RunBaseline(Params,Neuro);        
        Params.BaselineRunningFlag=true;
        Neuro = RunBaseline(Params,Neuro);
        Params.BaselineRunningFlag=false;
        
        % set flags back to original vals
        Neuro.UpdateChStatsFlag = Params.UpdateChStatsFlag;
        Neuro.UpdateFeatureStatsFlag = Params.UpdateFeatureStatsFlag;
        Neuro.DimRed.Flag = Params.DimRed.Flag;

        % save of useful stats and params
        SavePersistence(Params,Neuro,KF,0)
        
        pause(5.0)
        
        
    else, % if baseline is set to 0, just load stats
        f=load(fullfile(Params.Persistencedir, 'ch_stats.mat'));
        Neuro.ChStats = f.ch_stats;
        f=load(fullfile(Params.Persistencedir, 'feature_stats.mat'));
        Neuro.FeatureStats = f.feature_stats;
        clear('f');
    end
    
    % Imagined Cursor Movements Loop
    if Params.NumImaginedBlocks>0,
        [Neuro,KF,Params] = RunTask(Params,Neuro,1,KF);
    end
    
    % Adaptation Loop
    if Params.NumAdaptBlocks>0,
        [Neuro,KF,Params] = RunTask(Params,Neuro,2,KF);
    end
    
    % Fixed Decoder Loop
    if Params.NumFixedBlocks>0,
        [Neuro,KF,Params] = RunTask(Params,Neuro,3,KF);
    end
    
    % Pause and Finish!
%     PsychPortAudio('Close', [,  Params.sound_pahandle]);
    ExperimentStop(Params,0);

catch ME % handle errors gracefully
    Screen('CloseAll')
    for i=length(ME.stack):-1:1
        if i==1
            errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s\n\n', ...
                ME.stack(1).name, ME.stack(1).line, ME.message);
        else
            errorMessage = sprintf('Error in function %s() at line %d.\n\n', ...
                ME.stack(i).name, ME.stack(i).line);
        end
        fprintf(1,'\n%s\n', errorMessage);
    end
    
    keyboard;
end

end % ExperimentStart
