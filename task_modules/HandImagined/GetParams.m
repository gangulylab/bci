function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Experiment
Params.Task = 'HandImagined';
switch Params.ControlMode,
    case 1, Params.ControlModeStr = 'MousePosition';
    case 2, Params.ControlModeStr = 'MouseVelocity';
    case 3, Params.ControlModeStr = 'KalmanPosVel';
    case 4, Params.ControlModeStr = 'KalmanVelocity';
end

%% Control
Params.CenterReset      = true;
Params.Assistance       = 0; %0.05; % value btw 0 and 1, 1 full assist
Params.DaggerAssist 	= true;

Params.CLDA.Type        = 3; % 0-none, 1-refit, 2-smooth batch, 3-RML
Params.CLDA.AdaptType   = 'linear'; % {'none','linear'}, affects assistance & lambda for rml

Params.InitializationMode   = 4; % 1-imagined mvmts, 2-shuffled imagined mvmts, 3-choose dir, 4-most recent KF
Params.BaselineTime         = 0; % secs
Params.BadChannels          = [];
Params.SpatialFiltering     = false;
Params.UseFeatureMask       = true;
Params.GenNeuralFeaturesFlag= false; % if blackrock is off, automatically sets to true

%% Cursor Velocity
Params.Gain                     = 7;
Params.OptimalVeloctityMode     = 1; % 1-vector to target, 2-LQR
Params.VelocityTransformFlag    = false;
Params.MaxVelocityFlag          = false;
Params.MaxVelocity              = 200;

%% Cursor Click
Params.ClickerBins = -1; % set to -1 to use target hold time instead of click
Params.DecisionBoundary= -0.5;
Params.ClickerDataCollection = true; % if true, does not use clicker, freezes cursor when in target
if Params.ClickerDataCollection,
    Params.ClickerBins = -1; % must override to not use clicker
end

%% Sync to Blackrock
Params.ArduinoSync = true;

%% Update rate in pixels if decoded correctly 
% expressed as a percentage of the overall target distance
Params.PixelLength = 0.05;

%% Neural feature smoothing
Params.SmoothDataFlag = true;
Params.FeatureBufferSize = 5;

%% Timing
Params.ScreenRefreshRate = 5; % Hz
Params.UpdateRate = 5; % Hz

%% Discrete Decoder name
Params.DiscreteDecoder = 'clicker_svm_mdl_6Dir_hG.mat';

%% Multi State Decision Boundary
% set this to negative values. I would say -0.3 to -0.6 would be okay
Params.MultiDecisionBoundary = 0; 



%% Neural network classifier option
% set this to true to use neural network
% also set the softmax option
Params.NeuralNetFlag = true;
if Params.NeuralNetFlag
    Params.NeuralNetSoftMaxThresh = 0.7;    
    Params.Use3Features = true;
    Params.NeuralNetFunction = 'smallerMLP_6DoF_PlusOK_Trained4mOnlineData_20210201';
    %Params.NeuralNetFunction = 'MLP_6DoF__PlusOK_Trained4mOnlineData_20210201';    
    %Params.NeuralNetFunction = 'MLP_6DoF_PlusOK_Trained4mAllData_20210201';
        
    % LAST USED VERSION OF NN CLASSIFIER
    %Params.NeuralNetFunction = 'MLP_6DoF_Trained4mOnlineData_PlusStop';

else
    Params.NeuralNetSoftMaxThresh = 0;
end
%% biLSTM classifier option
Params.biLSTMFlag = false;
if Params.biLSTMFlag
    Params.biLSTMSoftMaxThresh = 0.45;
end

%% biLSTM classifier option
Params.biLSTMFlag = false;
if Params.biLSTMFlag
    Params.biLSTMSoftMaxThresh = 0.45;
end

Params.LSTMFunctionName = 'net_bilstm_20220824';%'net_bilstm_20220929_update';% or use 'net_bilstm_20220824';
Params.LSTM = load(fullfile('clicker',Params.LSTMFunctionName));
Params.LSTM = Params.LSTM.net_bilstm_20220824; %net_bilstm_20220929_update; % or use net_bilstm_20220824
Params.LSTMBufferSize = 1000;
Params.SaveLSTMFeatures = false;

Params.LSTM_Output_Method = false;
if Params.LSTM_Output_Method
    f = load(fullfile('clicker','lstm_output_pattern.mat'));
    Params.lstm_output_pattern = f.lstm_output_pattern;
    Params.LSTM_Output_Method_Thresh = 0.85;
end
%% Targets: radial layout
Params.NumReachTargets   = 6;
Params.TargetSpacing     = 10; % px
Params.OuterCircleRadius = 350; % defines outer edge of target
Params.InnerCircleRadius = 150; % defines inner edge of target
% Params.ReachTargetRadius = .5*(Params.InnerCircleRadius + Params.OuterCircleRadius);

Params.ReachTargetRadius = 250;

d2 = sqrt(1/2);
d3 = sqrt(1/3);

Params.ReachTargetPositions = [Params.ReachTargetRadius, 0, 0;...
    0, Params.ReachTargetRadius, 0; ...
    -Params.ReachTargetRadius, 0, 0;...
    0, -Params.ReachTargetRadius, 0; ...
    0,0,Params.ReachTargetRadius;...
    0, 0,-Params.ReachTargetRadius;...
    d2*Params.ReachTargetRadius, d2*Params.ReachTargetRadius, 0;...
    -d2*Params.ReachTargetRadius, d2*Params.ReachTargetRadius, 0;...
    -d2*Params.ReachTargetRadius, -d2*Params.ReachTargetRadius, 0;...
    d2*Params.ReachTargetRadius, -d2*Params.ReachTargetRadius, 0;...
    0,    d2*Params.ReachTargetRadius, d2*Params.ReachTargetRadius;...
    0, -d2*Params.ReachTargetRadius, d2*Params.ReachTargetRadius;...
    0, -d2*Params.ReachTargetRadius, -d2*Params.ReachTargetRadius;...
    0, d2*Params.ReachTargetRadius, -d2*Params.ReachTargetRadius;...
        d2*Params.ReachTargetRadius, 0, d2*Params.ReachTargetRadius;...
    -d2*Params.ReachTargetRadius, 0, d2*Params.ReachTargetRadius;...
    -d2*Params.ReachTargetRadius, 0, -d2*Params.ReachTargetRadius;...
    d2*Params.ReachTargetRadius, 0, -d2*Params.ReachTargetRadius;...
        d3*Params.ReachTargetRadius, d3*Params.ReachTargetRadius, -d3*Params.ReachTargetRadius;...
    -d3*Params.ReachTargetRadius, d3*Params.ReachTargetRadius, -d3*Params.ReachTargetRadius;...
    -d3*Params.ReachTargetRadius, -d3*Params.ReachTargetRadius, -d3*Params.ReachTargetRadius;...
    d3*Params.ReachTargetRadius, -d3*Params.ReachTargetRadius, -d3*Params.ReachTargetRadius;...
        d3*Params.ReachTargetRadius, d3*Params.ReachTargetRadius, -d3*Params.ReachTargetRadius;...
    -d3*Params.ReachTargetRadius, d3*Params.ReachTargetRadius, d3*Params.ReachTargetRadius;...
    -d3*Params.ReachTargetRadius, -d3*Params.ReachTargetRadius, d3*Params.ReachTargetRadius;...
    d3*Params.ReachTargetRadius, -d3*Params.ReachTargetRadius, d3*Params.ReachTargetRadius];
    
    


%% Kalman Filter Properties
Params.SaveKalmanFlag = false; % if true, saves kf at each time bin, if false, saves kf 1x per trial
G = Params.Gain;
t = 1/Params.UpdateRate;
a = 0.91;%.825;
w = 120;
if Params.ControlMode>=3,
    Params = LoadKF2dDynamics(Params, G, t, a, w);
end

%% LQR Optimal Velocity Controller
if Params.OptimalVeloctityMode==2,
    Params = LoadLQR2dDynamics(Params, G, t, a);
end

%% Velocity Command Online Feedback
Params.DrawVelCommand.Flag = true;
Params.DrawVelCommand.Rect = [-425,-425,-350,-350];

%% Trial and Block Types
Params.NumImaginedBlocks    = 1;
Params.NumAdaptBlocks       = 0;
Params.NumFixedBlocks       = 0;

% Cardinal Directions
% Params.NumTrialsPerBlock    = 6;
% Params.TargetOrder          = [1:6];
% 
% % Diagonals in the Horizontal Plane
Params.NumTrialsPerBlock    = 20;
% % Params.TargetOrder          = [7:10];
% 
% 
% Params.TargetOrder = Params.TargetOrder(randperm(length(Params.TargetOrder)));  % randomize order
% Params.TargetOrder          = [Params.TargetOrder, 1];

%% CLDA Parameters
TypeStrs                = {'none','refit','smooth_batch','rml'};
Params.CLDA.TypeStr     = TypeStrs{Params.CLDA.Type+1};

Params.CLDA.UpdateTime = 80; % secs, for smooth batch
Params.CLDA.Alpha = exp(log(.5) / (120/Params.CLDA.UpdateTime)); % for smooth batch
 
% Lambda
Params.CLDA.Lambda = 5000; % for RML
FinalLambda = 5000; % for RML
DeltaLambda = (FinalLambda - Params.CLDA.Lambda) ...
    / ((Params.NumAdaptBlocks-3)...
    *Params.NumTrialsPerBlock...
    *Params.UpdateRate * 3); % bins/trial;

Params.CLDA.DeltaLambda = DeltaLambda; % for RML
Params.CLDA.FinalLambda = FinalLambda; % for RML

Params.CLDA.FixedRmlFlag = false; % for RML during fixed
Params.CLDA.FixedLambda = FinalLambda; % for RML during fixed

switch Params.CLDA.AdaptType,
    case 'none',
        Params.CLDA.DeltaLambda = 0;
        Params.CLDA.DeltaAssistance = 0;
    case 'linear',
        switch Params.CLDA.Type,
            case 2, % smooth batch
                Params.CLDA.DeltaAssistance = ... % linearly decrease assistance
                    Params.Assistance...
                    /(Params.NumAdaptBlocks*Params.NumTrialsPerBlock*5/Params.CLDA.UpdateTime);
            case 3, % RML
                Params.CLDA.DeltaAssistance = ... % linearly decrease assistance
                    Params.Assistance...
                    /((Params.NumAdaptBlocks-1)*Params.NumTrialsPerBlock);
            otherwise, % none or refit
                Params.CLDA.DeltaAssistance = 0;
        end
end

%% Hold Times
Params.TargetHoldTime = 1;
Params.InterTrialInterval = 2;
Params.InstructedDelayTime = 1;
Params.CueTime = 1.5;
Params.MaxStartTime = 25;
Params.MaxReachTime = 25 ;
Params.InterBlockInterval = 10; % 0-10s, if set to 10 use instruction screen
Params.ImaginedMvmtTime = 3;

%% Feedback
Params.FeedbackSound = false;
Params.ErrorWaitTime = 0;
Params.ErrorSound = 1000*audioread('buzz.wav');
Params.ErrorSoundFs = 8192;
[Params.RewardSound,Params.RewardSoundFs] = audioread('reward1.wav');
% play sounds silently once so Matlab gets used to it
sound(0*Params.ErrorSound,Params.ErrorSoundFs)

%% Robotics 
Params.RobotMode            = 3;  % 0: Horizontal, 1: Vertical+Gripper, 3: 3D robot 
Params.RobotDirectionLines  = 1;  % 0: No lines, 1: Lines
Params.RunningModeBinNum    = 3;  % 1: No filtering, 3+: running mode filter of last n bins: Try 4 bins?


%% Hand
%Target order: 1:thumb, 2:index, 3:middle, 5:ring, 5:pinky, 6:power, ...
%7:pinch, 8:tripod, 9:wrist add, 10:wrist abd, 11: wrist flex, 12 : wrist
%extend, 13: wrist pronate, 14: wrist supinate

Params.NumTrialsPerBlock    = 12;
Params.TargetOrder          = [1:12];
Params.TargetOrder          = [Params.TargetOrder, 1];

Params.handVis = 1;
% Params.ActionDuration = 4;    % seconds
% 
% delta = 1/(0.5*Params.ActionDuration*Params.UpdateRate);
% 
% Params.angles = [0.1:-delta:-1,-1:delta:0.1, 0.1, 0.1, 0.1];
% Params.angles1d = [0.1:-0.5*delta:-1, -1, -1, -1];


Params.CycleDuration    = 2.6;    % seconds
Params.NumCycles        = 3;
Params.NumBinsPause     = 0;

% delta = 1/(0.5*Params.CycleDuration*Params.UpdateRate);
numBins = (Params.CycleDuration + 1/Params.ScreenRefreshRate)*Params.ScreenRefreshRate;

st = 0.1;
en = -1;
d1 = abs(st-en);
d2 = 2*d1;

delta1 = d1/(numBins - 1);
delta2 = d2/(numBins - 1);

Params.angles = [];
Params.angles1d = [];

% full range 0.1: 1

for i = 1:Params.NumCycles

Params.angles = [Params.angles, st:-delta2:en,(en+ delta2):delta2:(st-delta2), st,st,st,st];
Params.angles1d = [Params.angles1d, st:delta1:en];

end

Params.angles = [Params.angles, st,st,st,st];

Params.d1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0];
Params.showDecodeLines = 1;
end % GetParams
