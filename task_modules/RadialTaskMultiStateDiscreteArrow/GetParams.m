function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Experiment
Params.Task = 'DiscreteArrow';
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

%% Neural feature smoothing
Params.SmoothDataFlag = true;
Params.FeatureBufferSize = 5;

%% Bins for successful target selection
% The number of bins of successful decodes to hit the target
% Set this to 2/3 bins if enforcing a null class i.e.
% Params.MultiDecisionBoundary <0
Params.ClickCounter=5;

%% Timing
Params.ScreenRefreshRate = 8; % Hz
Params.UpdateRate = 8; % Hz   

%% Discrete Decoder name
Params.DiscreteDecoder = 'clicker_svm_mdl_Hand.mat';
%clicker_svm_mdl_OnlineDays123_OnlyArrowForTraining
%clicker_svm_mdl_Imagined_Smoothing_RtLeg

%% Multi State Decision Boundary for Null Class

% set this to negative values between -2 and -4 IMPORTANT: if this is set
% to large -ve values (enforcing a null class), don't set
% Params.ClickCounter to be > 4 bins, else it might get too difficult for
% subject.

Params.MultiDecisionBoundary =0; 


%% Neural network classifier option
% set this to true to use neural network
% also set the softmax option
Params.NeuralNetFlag = true;
if Params.NeuralNetFlag
   Params.NeuralNetSoftMaxThresh = 0.4;
   Params.NeuralNetFunction = 'MLP_Tong_Lips_PinchHand_Day2A';
   %Params.NeuralNetFunction = 'MLP_4Dir_Imagined_20210217_Day3_AllFeat';
   %Params.NeuralNetFunction = 'multilayer_perceptron_4Dir_MimeUpTongueIn_OnlineData';
else
    Params.NeuralNetSoftMaxThresh = 0;
end

%% To get ERPs
% freezes updates and points arrow towards the target for trial duration
Params.ERPs=false;

%% Targets: radial layout
Params.NumReachTargets   = 4;
Params.TargetSpacing     = 10; % px
Params.OuterCircleRadius = 600; % defines outer edge of target
Params.InnerCircleRadius = 300; % defines inner edge of target
Params.ReachTargetRadius = .5*(Params.InnerCircleRadius + Params.OuterCircleRadius);

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
Params.CursorSize = 50;
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
Params.NumFixedBlocks       = 2;
Params.NumTrialsPerBlock    = 12;

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
Params.InterTrialInterval = 1.0;
Params.InstructedDelayTime = 1.0;
Params.CueTime = 1;
Params.MaxStartTime = 25;
Params.MaxReachTime = 5 ;
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

end % GetParams
