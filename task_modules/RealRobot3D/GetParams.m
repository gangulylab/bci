function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Experiment
Params.Task = 'RealRobot3D';
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
Params.DiscreteDecoder = 'clicker_svm_mdl_6Dir_3Feat_462021.mat';

%% Multi State Decision Boundary
% set this to negative values. I would say -0.3 to -0.6 would be okay
Params.MultiDecisionBoundary = 0; 


%% Neural network classifier option
% set this to true to use neural network
% also set the softmax option
Params.NeuralNetFlag = false;
if Params.NeuralNetFlag
    Params.NeuralNetSoftMaxThresh = 0.45;       
    Params.Use3Features = false;
    Params.Use4Features = true; % for low gamma
    Params.NeuralNetFunction = 'MLP_7DoF_PnP_2022July_lg_0817_PM1';%'MLP_7DoF_PnP_2022July_lg_0722_PM3'; %'MLP_7DoF_PnP_2022July_lg_0720_PM';%'MLP_7DoF_PnP_2022July_lg';
%     Params.NeuralNetFunction = 'MLP_FlipView3D_20210817_PM1';
%     Params.NeuralNetFunction = 'MLP_PreTrained_7DoF_PnP4';%'MLP_PreTrained_7DoF_PnP';

%    Params.NeuralNetFunction =  'MLP_7DoF_ZWrist_07012022B '; %'MLP_7DoF_PnP_2022Feb_2norm'; 
%    Params.NeuralNetFunctionName = load(fullfile('clicker','net_new_7DoF_ZWrist_07012022B'));
%    Params.NeuralNet = Params.NeuralNetFunctionName.net_new_7DoF_ZWrist_07012022B;
    
%     Params.NeuralNetFunction = 'multilayer_perceptron_6DoF_Online_Apr16_2021';
    %Params.NeuralNetFunction = 'MLP_6DoF_PlusOK_Trained4mAllData_20210212';    

else
    Params.NeuralNetSoftMaxThresh = 0;
end

%% Use ensemble neural network
 Params.NeuralNetEnsemble = false;
% Params.NeuralNetSoftMaxThresh = 0.45;   
% Params.NeuralNetName = 'net_7DoF_PnP4_ensemble_batch_0520A';%'net_7DoF_PnP4_ensemble_batch';
% Params.NeuralNetFunction = load(fullfile('clicker',Params.NeuralNetName)); 

%% Neural network 2 classifier option
% Trained in a different way using different optimizer

Params.NeuralNet2Flag = false;
if Params.NeuralNet2Flag
    Params.NeuralNet2SoftMaxThresh = 0.45    ;       
    Params.Use3Features = true;
    Params.NeuralNet2 = load(fullfile('clicker','net_new_7DoF_ZWrist_05252022')); % 7DoF classifier trained in a different way
    Params.NeuralNet2.net = Params.NeuralNet2.net_mlp_7DoF_Feb2022;
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


%% CONVOLUTIONAL NEURAL NET OPTION
% set this to true to use neural network
% also set the softmax option
Params.ConvNeuralNetFlag =false;
if Params.ConvNeuralNetFlag
    Params.ConvNeuralNetSoftMaxThresh = 0.6;       
    Params.ConvUse3Features = true;
    Params.ConvNeuralNetFunctionName = 'CNN_classifier_Online_Apr16_2021_B';%'CNN_classifier_B1_16thApr';%'CNN_classifier_B1_OnlyLastBins';    
%     Params.ConvNeuralNetFunctionName = 'CNN_classifier_B1_OnlyLastBins_AndState2';    
    Params.ConvNeuralNet = load(fullfile('clicker','CNN_classifier'));
else
    Params.ConvNeuralNetSoftMaxThresh = 0;
end

%% biLSTM classifier option

Params.biLSTMFlag = true;
if Params.biLSTMFlag
    Params.biLSTMSoftMaxThresh = 0.4;
end

%Params.LSTMFunctionName = 'net_bilstm_robot_20220824B'; %net_bilstm_robot_20220928B %'net_bilstm_robot_20220824B';%net_bilstm_robot_20220824_early_stop
Params.LSTMFunctionName = 'net_bilstm_robot_20220929';% OR TRY 'net_bilstm_robot_20220824B';
Params.LSTM = load(fullfile('clicker',Params.LSTMFunctionName));
Params.LSTM = Params.LSTM.net_bilstm_robot_20220929;%OR TRY net_bilstm_robot_20220824B
Params.LSTMBufferSize = 1000;
Params.SaveLSTMFeatures = false;


%% ADAPTIVE BASELINE FLAG 
% data is baseline to state 1 data
Params.AdaptiveBaseline = false;


%% POOLING CHANNELS FOR CONTROL
% set this 1 only during online control
Params.ChPooling = true; 

%% IMPORT PC WEIGHTS AND MEAN FOR BETA BAND ANALYSIS

Params.BetaWts = load(fullfile('clicker','betawts_stop'));
Params.BetaMean = load(fullfile('clicker','betamean'));
%% Targets: radial layout

Params.ReachTargetRadius = 180;
d2 = sqrt(1/2);
d3 = sqrt(1/3);

h = 260;
% 
% Params.ReachTargetPositions = [Params.ReachTargetRadius, 0, h;...
%     0, Params.ReachTargetRadius, h; ...
%     -Params.ReachTargetRadius, 0, h;...
%     0, -Params.ReachTargetRadius, h; ...
%     0,0,450;...
%     0, 0, 100;...
%     d2*Params.ReachTargetRadius, d2*Params.ReachTargetRadius, 0;...
%     -d2*Params.ReachTargetRadius, d2*Params.ReachTargetRadius, 0;...
%     -d2*Params.ReachTargetRadius, -d2*Params.ReachTargetRadius, 0;...
%     d2*Params.ReachTargetRadius, -d2*Params.ReachTargetRadius, 0];
% 
% Params.ReachTargetPositions = [Params.ReachTargetPositions;...
%     Params.ReachTargetPositions(1,1:2),-150;...
%     Params.ReachTargetPositions(2,1:2),-150;...
%     Params.ReachTargetPositions(3,1:2),-150;...
%     Params.ReachTargetPositions(4,1:2),-150;
%     Params.ReachTargetPositions(3,1:2),-150;...
%     Params.ReachTargetPositions(4,1:2),-150];

Params.ReachTargetPositions = [240, -70, 410;...
240, -220, 410;...
240, 80, 410];


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
Params.NumTrialsPerBlock    = 4;
Params.TargetOrder          = [1:4];
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
            case 3 % RML
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
Params.MaxStartTime = 50;
Params.MaxReachTime = 150;
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

Params.RobotMode    = 2; 
Params.wl           = [-50, -67, 10];
Params.wu           = [5, -15 55];

if Params.RobotMode == 1 % lateral R2G boxes
    Params.ValidDir             = [1:9];
    Params.StartPos             = [-200, 200,220];
    Params.StartPos             = [0, 30, 300];
    Params.NumTrialsPerBlock    = 1;
    Params.TargetOrder          = [1];   
    Params.OperationModeReset   = 1;
    Params.wristStartX          = 3.1415/2*10; 
    Params.wristStartZ          = 0; 
    Params.autoCenterOverTarget = 0;
    Params.autoCenterDist       = 5;
    Params.graspOrientation     = 1;

elseif Params.RobotMode == 2 % vertical R2G
    Params.ValidDir          = [1:9];
    Params.StartPos          = [-150, -150,300];
    Params.NumTrialsPerBlock    = 1;
    Params.TargetOrder          = [1];   
    Params.OperationModeReset = 0;
    Params.wristStartX = 3.1415*10; 
    Params.wristStartZ = 0; 
    Params.autoCenterOverTarget = 0;
    Params.autoCenterDist = 5;
    Params.zlim = 5;
    Params.graspOrientation = 0;
    
    Params.StartWristX = [3.1415]*10;                    
    Params.StartWristZ = 10*[0];
    Params.StartWristY = 10*[0];    
    
Params.wl           = [-50, -67, 16];
Params.wu           = [5, -15 55];
    
    
elseif Params.RobotMode == 3  % lateral R2G wall
    Params.ValidDir             = [1:9];
%     Params.StartPos             = [-200, 200, 220]; % vertically alligned start
    Params.StartPos             = [100, -90, 320]; %
    Params.NumTrialsPerBlock    = 1;
    Params.TargetOrder          = [1];   
    Params.OperationModeReset   = 0;
    Params.wristStartX          = 3.1415/2*10; 
    Params.wristStartZ          = 0; 
    Params.graspOrientation     = 1;
    
    
    Params.StartWristX = [3.1415/2]*10;                    
    Params.StartWristZ = 10*[0];
    Params.StartWristY = 10*[0];        

end

Params.index        = 1;
Params.clickOrder   = [ones(5,1)];
Params.ReachTargets      = [1,2,3,4,5,6];
Params.TargetOrder  = [Params.TargetOrder, 1];

Params.RunningModeBinNum    = 5;  % 1: No filtering, 3+: running mode filter of last n bins: Try 4 bins?
Params.RunningModeZero      = 1;  % 1: No motion if no winner, 0: maintain prior decision if no winner

Params.RobotTargetRadius    = 50;
Params.RobotTargetDim       = 1;

% Robot Dynamics
Params.deltaT   = 1/Params.UpdateRate;
% Params.k_v      = 0.7;
% Params.k_i      = 40;
Params.k_v      = 0.8;
Params.k_i      = 18;    % standard: 18

Params.r_v      = 0.8;
Params.r_i      = 100;

Params.r_v      = 0.7;
Params.r_i      = 90;
% Params.r_v      = 0.7;
% Params.r_i      = 150;

Params.LongTrial = 0;
Params.RobotClicker     = 1;
Params.TargetHoldTime   = 0.25;

Params.boundaryDist     = 0;
Params.boundaryVel      = 0;
Params.AssistAlpha      = 0.0;
Params.AutoGrasp        = 0;
Params.GraspTask        = 1;
Params.lowGainMode      = 0;
Params.autoCenterOverTarget    = 0;
Params.autoCenterDist = 0;

Params.SwitchBinNum     = 8;
Params.SwitchBinThresh  = 0.7;
Params.GraspBinNum      = 8;
Params.GraspBinThresh   = 0.7;

% Beta
Params.UseBetaStop      = 0;
Params.BetaThreshold    = 0.5;

end % GetParams
