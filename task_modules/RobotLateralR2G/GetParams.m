function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Experiment
Params.Task = 'RobotLateralR2G';
switch Params.ControlMode
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
if Params.ClickerDataCollection
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
Params.UseSVM = false;
Params.DiscreteDecoder = 'clicker_svm_mdl_6Dir_3Feat_462021.mat';

%% Multi State Decision Boundary
% set this to negative values. I would say -0.3 to -0.6 would be okay
Params.MultiDecisionBoundary = 0; 

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

%% BIAS CORRECTION FOR LEFT LEG
% scales the probabilities of the decoder towards a specific action by a
% prespecific amount

Params.NeuralBias = false;
Params.NeuralNetBiasDirection = 2; % class o/p that has the bias. 
Params.NeuralNetBiasCorrection = 0.7; % pulls decision probabilities by this amount



%% ADAPTIVE BASELINE FLAG 
% data is baseline to state 1 data
Params.AdaptiveBaseline = false;

%% CONVOLUTIONAL NEURAL NET OPTION
% set this to true to use neural network
% also set the softmax option
Params.ConvNeuralNetFlag = false;
if Params.ConvNeuralNetFlag
    Params.ConvNeuralNetSoftMaxThresh = 0.6;       
    Params.ConvUse3Features = true;
%     Params.ConvNeuralNetFunctionName = 'CNN_classifier_B1_OnlyLastBins';    
    Params.ConvNeuralNetFunctionName = 'CNN_classifier_B1_OnlyLastBins_AndState2';    
    Params.ConvNeuralNet = load(fullfile('clicker','CNN_classifier'));
else
    Params.ConvNeuralNetSoftMaxThresh = 0;
end


%% POOLING CHANNELS FOR CONTROL
% set this 1 only during online control
Params.ChPooling = true; 

%% Targets: radial layout
Params.NumReachTargets   = 6;
Params.TargetSpacing     = 10; % px
Params.OuterCircleRadius = 350; % defines outer edge of target
Params.InnerCircleRadius = 150; % defines inner edge of target
% Params.ReachTargetRadius = .5*(Params.InnerCircleRadius + Params.OuterCircleRadius);

Params.ReachTargetRadius = 150;

d2 = sqrt(1/2);
d3 = sqrt(1/3);

Params.ReachTargetPositions{1} = [0,0, -200; -250,0, -50];    
Params.LongStartPos{1} =  [200, 0, 0];


Params.ReachTargetPositions{2} = [0,150, -200; -250,0, -50];    
Params.LongStartPos{2} =  [200, 0, 0];


Params.ReachTargetPositions{3} = [0,-150, -200; -250,0, -50];    
Params.LongStartPos{3} =  [200, 0, 0];


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

% Cardinal Directions
Params.NumTrialsPerBlock    = 1;        
Params.TargetOrder          = [1];  %  Set set 1 or 3 

Params.TargetOrder = Params.TargetOrder(randperm(length(Params.TargetOrder)));  % randomize order
Params.TargetOrder          = [Params.TargetOrder, 1];

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

switch Params.CLDA.AdaptType
    case 'none'
        Params.CLDA.DeltaLambda = 0;
        Params.CLDA.DeltaAssistance = 0;
    case 'linear'
        switch Params.CLDA.Type
            case 2 % smooth batch
                Params.CLDA.DeltaAssistance = ... % linearly decrease assistance
                    Params.Assistance...
                    /(Params.NumAdaptBlocks*Params.NumTrialsPerBlock*5/Params.CLDA.UpdateTime);
            case 3 % RML
                Params.CLDA.DeltaAssistance = ... % linearly decrease assistance
                    Params.Assistance...
                    /((Params.NumAdaptBlocks-1)*Params.NumTrialsPerBlock);
            otherwise % none or refit
                Params.CLDA.DeltaAssistance = 0;
        end
end

%% Hold Times
Params.TargetHoldTime = 1;
Params.InterTrialInterval = 2;
Params.InstructedDelayTime = 1;
Params.CueTime = 0.75;
Params.MaxStartTime = 25;
Params.MaxReachTime = 180;
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

Params.flipView = 1;  % 0: standard simulation view, 1: real robot view (across table)

if Params.flipView
    Params.RobotMode = 12;
else
    Params.RobotMode = 10;  
end

Params.boundaryDist     = .01;
Params.boundaryVel      = 0;
Params.RboundaryDist    = .01;
Params.RboundaryVel     = 0;
Params.AssistAlpha      = 0.0;
Params.LongTrial        = 1;

Params.limit = [-256, 256; -256 256; -200 256];

Params.RobotDirectionLines  = 0;  % 0: No lines, 1: Lines
Params.RunningModeBinNum    = 4;  % 1: No filtering, 3+: running mode filter of last n bins: Try 4 bins?
Params.RunningModeZero      = 1;  % 1: Null output if no winner, 0: maintain prior decision if no winner

Params.RobotTargetDim = 1;

Params.ReachTargets      = [1,2,3,4,5,6];
Params.ValidDir          = [1:6,7];

Params.deltaT = 1/Params.UpdateRate;
Params.k_v = 0.9;
Params.k_i = 10.0;       % Can change this to make it more/less responsive  (increase to make faster, normal is 10)

Params.dA = [1 0 0  Params.deltaT 0 0;...
                    0 1 0 0 Params.deltaT 0;...
                    0 0 1 0 0 Params.deltaT;...
                    0 0 0 Params.k_v 0 0;...
                    0 0 0 0 Params.k_v 0;...
                    0 0 0 0 0 Params.k_v];
                
Params.dB = [zeros(3);...
                    eye(3)];
Params.dB = Params.dB*Params.k_i;

% Rotation Dynamics
Params.ra = 0.8;
Params.rb = 0.03;
Params.opMode = 0;
Params.Rlimit = [-0.6, 0.6];

% Target
Params.RobotTargetRadius = [40, 70];  % Radius of targets, can be set separately

% Clicker
Params.RobotClicker     = 1;     % 0: trial ends with hold time, 1: trial ends with click
Params.ClickerBinNum    = 4;
Params.ClickerBinThresh = 0.6;

Params.switchLockCnt    = 10;
Params.GraspBinNum      = 5;
Params.GraspBinThresh   = 0.6;

Params.autoCenterOverTarget = 0;
Params.GripperOrnInit = [0, pi/2,0];

Params.printDecode = 1;

Params.ClickerBreak = 1;
Params.BreakGain = 0.75;  %(0-1; 0: Full stop, 1: No break) 
Params.zeroVelOnModeSwitch = 0;

Params.SubTaskStart = 1;  % 1: full task, 2: start at subtask 2

if Params.SubTaskStart == 2
    Params.LongStartPos{1} =  [0, 0, -200];
    Params.opMode = 1;
end
    
%% For Debug

Params.index = 1;
% Params.clickOrder = [7*(ones(1,7)),2*(ones(1,20)), 4*(ones(1,20)) ];
% 
% Params.clickOrder = [1*(ones(1,20)), 3*(ones(1,20)),2*(ones(1,20)), ...
%     4*(ones(1,20)), 7*(ones(1,7)),  6*(ones(1,20)),5*(ones(1,20)),...
%     2*(ones(1,10)),6*(ones(1,20)),5*(ones(1,20)), 4*(ones(1,10)),...
%     6*(ones(1,20)),5*(ones(1,20))];

Params.clickOrder = [5*(ones(1,40)), 6*(ones(1,40)) ...
    1*(ones(1,8)), 3*(ones(1,8))];



end % GetParams
