function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Experiment
Params.Task = 'Robot3DArrow';
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
Params.ArduinoSync =false;

%% Update rate in pixels if decoded correctly 
% expressed as a percentage of the overall target distance
Params.PixelLength = 0.05;

%% Neural feature smoothing
Params.SmoothDataFlag = true;
Params.FeatureBufferSize = 5;

%% Bins for successful target selection
% The number of bins of successful decodes to hit the target
% Set this to 2/3 bins if enforcing a null class i.e.
% Params.MultiDecisionBoundary <0
Params.ClickCounter = 4;

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
Params.NeuralNetFlag = false;
if Params.NeuralNetFlag
    Params.NeuralNetSoftMaxThresh = 0.45;       
    Params.Use3Features = false; % only delta beta hG
    Params.Use4Features = true; % inlcuding lg
    Params.NeuralNetFunction = 'MLP_9DoF_New_July2023';
    
    Params.NeuralNetFunction = 'MLP_9DoF_07282023Update2';
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
Params.Norm2 = false;

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
    Params.biLSTMSoftMaxThresh = 0.45;
end

% Params.LSTMFunctionName = 'net_bilstm_7DoF_Feb2024_RtWrist_Act4_20240411_2';% or use 'net_bilstm_20220824';
Params.LSTMFunctionName =  'net_bilstm_7DoF_Feb2024_RtWrist_Act4';
Params.LSTM = load(fullfile('clicker',Params.LSTMFunctionName));
% Params.LSTM = Params.LSTM.net_bilstm_7DoF_Feb2024_RtWrist_Act4_20240411_2; %net_bilstm_20220929_update; % or use net_bilstm_20220824
Params.LSTM = Params.LSTM.net_bilstm_7DoF_Feb2024_RtWrist_Act4;

Params.LSTMBufferSize = 1000;
Params.SaveLSTMFeatures = false;

Params.LSTM_Output_Method = false;
if Params.LSTM_Output_Method
    f = load(fullfile('clicker','lstm_output_pattern.mat'));
    Params.lstm_output_pattern = f.lstm_output_pattern;
    Params.LSTM_Output_Method_Thresh = 0.85;
end
%% BAD CHANNELS
% not needed as LSTM weights during training downweights bad channels
%Params.BadChannels=[99];

%% ADAPTIVE BASELINE FLAG 
% data is baseline to state 1 data
Params.AdaptiveBaseline = false;

%% POOLING CHANNELS FOR CONTROL
% set this 1 only during online control
Params.ChPooling = false; 

%% Targets: radial layout
Params.NumReachTargets   = 9;
Params.TargetSpacing     = 10; % px
Params.OuterCircleRadius = 350; % defines outer edge of target
Params.InnerCircleRadius = 150; % defines inner edge of target
% Params.ReachTargetRadius = .5*(Params.InnerCircleRadius + Params.OuterCircleRadius);

Params.ReachTargetRadius = 200;
dl    = Params.ReachTargetRadius/sqrt(2);


Params.ReachTargetPositions = [Params.ReachTargetRadius, 0, 0;...
    0, Params.ReachTargetRadius, 0; ...
    -Params.ReachTargetRadius, 0, 0;...
    0, -Params.ReachTargetRadius, 0; ...
    0,0,Params.ReachTargetRadius;...
    0, 0,-Params.ReachTargetRadius;...
    0,0,0;...
    0.5*Params.ReachTargetRadius,0,40;... % wrist rotation
    -0.5*Params.ReachTargetRadius,0,40;...; % wrist rotation
    dl,-dl,0;...  %  diagonals
    dl, dl,0;
    -dl, dl, 0;...
    -dl, -dl,0;...
    dl,0,dl;...
    0, dl, dl;...
    -dl,0,dl;...
    0, -dl, dl;
    dl,0,-dl;...
    0, dl,-dl;...
    -dl,0,-dl;...
    0, -dl,-dl;...
    0,0,0;...
    0,0,0;...
    0,0,0];
    



%% Kalman Filter Properties
Params.SaveKalmanFlag = false; % if true, saves kf at each time bin, if false, saves kf 1x per trial
G = Params.Gain;
t = 1/Params.UpdateRate;
a = 0.91;%.825;
w = 120;
if Params.ControlMode>=3
    Params = LoadKF2dDynamics(Params, G, t, a, w);
end

%% LQR Optimal Velocity Controller
if Params.OptimalVeloctityMode==2
    Params = LoadLQR2dDynamics(Params, G, t, a);
end

%% Velocity Command Online Feedback
Params.DrawVelCommand.Flag = true;
Params.DrawVelCommand.Rect = [-425,-425,-350,-350];

%% Trial and Block Types
Params.NumImaginedBlocks    = 0;
Params.NumAdaptBlocks       = 0;
Params.NumFixedBlocks       = 1;

 Params.NumTrialsPerBlock    = 27;              % 9 target (wrist rotation)
 Params.TargetOrder          = [1:9,1:9,1:9];
% 
Params.NumTrialsPerBlock    = 21;              % standard 7
Params.TargetOrder          = [1:7,1:7,1:7];

% Params.NumTrialsPerBlock    = 20;                % co-activation - center plane
% Params.TargetOrder          = [10:13, 10:13, 10:13, 10:13, 10:13];

% Params.NumTrialsPerBlock    = 20;                % co-activation - top plane
% Params.TargetOrder          = [14:17, 14:17, 14:17, 14:17, 14:17];
% 
% Params.NumTrialsPerBlock    = 12;                % co-activation - bottom plane
% Params.TargetOrder          = [18:21, 18:21, 18:21];
% 
% Params.NumTrialsPerBlock    = 9;                % co-activation - opposing
% Params.TargetOrder          = [22:24, 22:24, 22:24];



Params.TargetOrder = Params.TargetOrder(randperm(length(Params.TargetOrder)));  % rand order
% Params.TargetOrder          = [Params.TargetOrder, 2];

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
Params.InterTrialInterval = 1;
Params.InstructedDelayTime = 1;
Params.CueTime = 1.0;
Params.MaxStartTime = 25;
Params.MaxReachTime = 5 ;
Params.InterBlockInterval = 10; % 0-10s, if set to 10 use instruction screen
% Params.ImaginedMvmtTime = 3;

%% Feedback
Params.FeedbackSound = false;
Params.ErrorWaitTime = 0;
Params.ErrorSound = 1000*audioread('buzz.wav');
Params.ErrorSoundFs = 8192;
[Params.RewardSound,Params.RewardSoundFs] = audioread('reward1.wav');
% play sounds silently once so Matlab gets used to it
sound(0*Params.ErrorSound,Params.ErrorSoundFs)

%% Robotics 
Params.flipView = 0;
if Params.flipView
    Params.RobotMode        = 14; 
else
    Params.RobotMode            = 4; 
end
Params.LetterMode           = 0;  % 1: letter cues, 0: box cues
Params.RobotTargetRadius    = 40;
Params.RobotDirectionLines  = 1;  % 0: No lines, 1: Lines

Params.RunningModeBinNum    = 1;  % 1: No filtering, 3+: running mode filter of last n bins

if Params.RobotMode == 0
    Params.RobotTargetDim = 2;
elseif Params.RobotMode == 1
    Params.RobotTargetDim = 1;
end

Params.RobotTargetDim = 1;

Params.ReachTargets      = [1,2,3,4,5,6,7];
Params.ValidDir          = [1:7];

%% For Debug
Params.index = 1;
Params.clickOrder = [8*ones(1,5)];
Params.clickOrder = [8*ones(1,5),2*ones(1,5),3*ones(1,5),4*ones(1,5),5*ones(1,5),6*ones(1,5),7*ones(1,5),8*ones(1,5)]';
end % GetParams
