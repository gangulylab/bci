function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Experiment
Params.Task = 'Robot';
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
Params.FeatureBufferSize = 4;

%% Timing
Params.ScreenRefreshRate = 8; % Hz
Params.UpdateRate = 8; % Hz

%% Discrete Decoder name
Params.DiscreteDecoder = 'clicker_svm_mdl_6Dir_3Feat_462021.mat';

%% Multi State Decision Boundary
% set this to negative values. I would say -0.3 to -0.6 would be okay
Params.MultiDecisionBoundary = -2; 


%% Neural network classifier option
% set this to true to use neural network
% also set the softmax option
Params.NeuralNetFlag = true;
if Params.NeuralNetFlag
    Params.NeuralNetSoftMaxThresh = 0.5;       
    Params.Use3Features = true;
    Params.NeuralNetFunction = 'MLP_7DoF_PnP_2022Feb_2norm';
%  Params.NeuralNetFunction = 'MLP_PreTrained_7DoF_1006_AM2';
else
    Params.NeuralNetSoftMaxThresh = 0;
end


%% Neural network 2 classifier option
% Trained in a different way using different optimizer

Params.NeuralNet2Flag = false;
if Params.NeuralNet2Flag
    Params.NeuralNet2SoftMaxThresh = 0.4;       
    Params.Use3Features = true;
    Params.NeuralNet2 = load(fullfile('clicker','net_mlp_7DoF_Feb2022')); % 7DoF classifier trained in a different way
    
else
    Params.NeuralNet2SoftMaxThresh = 0;
end

%% NORMALIZING THE NEURAL FEATURES
Params.Norm2 = true;


%% POOLING CHANNELS FOR CONTROL
% set this 1 only during online control
Params.ChPooling = true; 


%% Targets: radial layout
Params.NumReachTargets   = 4;
Params.TargetSpacing     = 10; % px
Params.OuterCircleRadius = 350; % defines outer edge of target
Params.InnerCircleRadius = 150; % defines inner edge of target
% Params.ReachTargetRadius = .5*(Params.InnerCircleRadius + Params.OuterCircleRadius);

Params.ReachTargetRadius = 250;

Params.TargetsColor        = [100,100,100]; % all targets
Params.CuedTargetColor   = [55,255,0]; % cued target

% triangles
da = 360/Params.NumReachTargets; % degrees btw target centers
b = da/2; % degrees btw targets edges
Params.ReachTargetAngles = 0:da:(360-da);
for i=1:Params.NumReachTargets,
    a = Params.ReachTargetAngles(i);
    Params.ReachTargetPositions(i,:) = ...
        Params.ReachTargetRadius * [cosd(a) sind(a)]; % center of targets
    
    Params.ReachTargetVerts{i} = ...
        Params.OuterCircleRadius .* ...
        [0, 0
        cosd(a-b) sind(a-b)
        cosd(a+b) sind(a+b)];
end

% inner circle
sz = Params.InnerCircleRadius;
Params.InnerCircleColor = 0;
Params.InnerCircleRect = [-sz, -sz, sz, sz];

% useful for optimal cursor updates, not used for targets here
Params.TargetSize = Params.OuterCircleRadius - Params.InnerCircleRadius;

%% Cursor
Params.CursorColor = [0,102,255];
Params.InTargetColor    = [255,55,0];
Params.CursorSize = 10;
Params.CursorRect = [-Params.CursorSize -Params.CursorSize ...
    +Params.CursorSize +Params.CursorSize];

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
Params.NumImaginedBlocks    = 0;
Params.NumAdaptBlocks       = 0;
Params.NumFixedBlocks       = 1;
Params.NumTrialsPerBlock    = 10;

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
Params.InterTrialInterval = 1;
Params.InstructedDelayTime = 1;
Params.CueTime = 0.75;
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

Params.RobotTargetRadius    = 40;
Params.RobotMode            = 0;  % 0: Horizontal, 1: Vertical+Gripper
Params.RobotDirectionLines  = 1;  % 0: No lines, 1: Lines
Params.RunningModeBinNum    = 3;  % 1: No filltering, 3+: running mode filter of last n bins: Try 4 bins?
Params.RunningModeZero      = 1;  % 1: No motion if no winner, 0: maintain prior decision if no winner


if Params.RobotMode == 0
    Params.RobotTargetDim = 2;
elseif Params.RobotMode == 1
    Params.RobotTargetDim = 1;
end
end % GetParams
