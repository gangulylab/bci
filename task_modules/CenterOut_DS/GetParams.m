  function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Experiment
Params.Task = 'CenterOut_DS';
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
Params.ArduinoSync = false;

%% Update rate in pixels if decoded correctly 
% expressed as a percentage of the overall target distance
Params.PixelLength = 0.05;

%% Neural feature smoothing
Params.SmoothDataFlag       = true;
Params.FeatureBufferSize    = 5;

%% Timing
Params.ScreenRefreshRate    = 5; % Hz
Params.UpdateRate           = 5; % Hz

%% Kalman Filter Properties
Params.SaveKalmanFlag = false; % if true, saves kf at each time bin, if false, saves kf 1x per trial
G = Params.Gain;
t = 1/Params.UpdateRate;
a = 0.1; %0.91;%.825;
w = 120;
if Params.ControlMode>=3
    Params = LoadKF2dDynamics(Params, G, t, a, w);
end
%% Discrete Decoder name
Params.UseSVM = false;
Params.DiscreteDecoder = 'clicker_svm_mdl_6Dir_3Feat_462021.mat';

%% Multi State Decision Boundary
% set this to negative values. I would say -0.3 to -0.6 would be okay
Params.MultiDecisionBoundary = -2; 

%% Neural network classifier option
% set this to true to use neural network
% also set the softmax option
Params.NeuralNetFlag = true;
if Params.NeuralNetFlag
   Params.NeuralNetSoftMaxThresh = 0.4;
%    Params.NeuralNetFunction = 'MLP_9Dir_B1_20240517_CL3_NoPooling';
   Params.NeuralNetFunction = 'MLP_7Dir_B1_20240621_CL3_NoPooling';
   %Params.NeuralNetFunction = 'multilayer_perceptron_4Dir_MimeUpTongueIn_OnlineData';
else
    Params.NeuralNetSoftMaxThresh = 0;
end

%% NEURAL NET2 (USED FOR PNP)

Params.NeuralNet2Flag = false;
Params.NeuralNet2UseAllFeat=true;
if Params.NeuralNet2Flag
    Params.NeuralNet2SoftMaxThresh = 0.4;
    if Params.NeuralNet2UseAllFeat
        Params.NeuralNet2FileName = 'net_B1_253_TrfLearn';
        Params.NeuralNet2 = load(fullfile('clicker',Params.NeuralNet2FileName));
        Params.NeuralNet2 = Params.NeuralNet2.net_B1_253_TrfLearn;
    else
        Params.NeuralNet2FileName = 'net_PnP_hG';
        Params.NeuralNet2 = load(fullfile('clicker',Params.NeuralNet2FileName));
        Params.NeuralNet2  = Params.NeuralNet2.net_PnP_hG;
    end
end

%% BAD CHANNELS

Params.SetBadChannels = [];

%% CONVOLUTIONAL NEURAL NET OPTION
% set this to true to use neural network
% also set the softmax option
Params.ConvNeuralNetFlag = false;
if Params.ConvNeuralNetFlag
    Params.ConvNeuralNetSoftMaxThresh = 0.4;       
    Params.ConvUse3Features = true;
    Params.ConvNeuralNetFunctionName = 'CNN_classifier_B2_new_C';    
    Params.ConvNeuralNet = load(fullfile('clicker',...
         Params.ConvNeuralNetFunctionName));
else
    Params.NeuralNetSoftMaxThresh = 0;
end

%% biLSTM classifier option
Params.biLSTMFlag = false;
if Params.biLSTMFlag
    Params.biLSTMSoftMaxThresh = 0.45;
end

Params.LSTMFunctionName = 'net_bilstm_20230521_lstm';%'net_bilstm_20220929_update';% or use 'net_bilstm_20220824';
Params.LSTM = load(fullfile('clicker',Params.LSTMFunctionName));
Params.LSTM = Params.LSTM.net_bilstm_20230521_lstm; %net_bilstm_20220929_update; % or use net_bilstm_20220824
Params.LSTMBufferSize = 1000;
Params.SaveLSTMFeatures = false;

Params.LSTM_Output_Method = false;
if Params.LSTM_Output_Method
    f = load(fullfile('clicker','lstm_output_pattern.mat'));
    Params.lstm_output_pattern = f.lstm_output_pattern;
    Params.LSTM_Output_Method_Thresh = 0.85;
end

%% LOAD THE CHMAP FILE
tmp = load(fullfile('clicker','ECOG_Grid_8596_000063_B3.mat'));
Params.ChMapB2 = tmp.ecog_grid;


%% 2-norm
Params.Norm2 = false;

%% ADAPTIVE BASELINE FLAG 
% data is baseline to state 1 data
Params.AdaptiveBaseline = false;

%% POOLING CHANNELS FOR CONTROL
% set this 1 only during online control
Params.ChPooling = false; 

%% Targets

Params.TargetSize       = 75;
Params.OutTargetColor   = [55,255,0];
Params.InTargetColor    = [255,55,0];

Params.StartPositions  = [0, 0];
Params.TargetRect = ...
    [-Params.TargetSize -Params.TargetSize +Params.TargetSize +Params.TargetSize];

Params.NumTargets   = 8;
Params.TargetRadius = 450;

Params.ReachTargets = zeros(Params.NumTargets, 2);
p_tmp   = [0,Params.TargetRadius];
for i = 1:Params.NumTargets
    theta = (i-1)*2*pi/Params.NumTargets;
    
    R = [cos(theta), -sin(theta);sin(theta), cos(theta)];    
    Params.ReachTargets(i,:) = R*p_tmp';
end

Params.NumStartPos  = length(Params.StartPositions);


%% Cursor
Params.CursorColor  = [0,102,255];
Params.CursorSize   = 15;  %30
Params.CursorRect   = [-Params.CursorSize -Params.CursorSize ...
    +Params.CursorSize +Params.CursorSize];

%% Trial and Block Types
Params.NumImaginedBlocks    = 0;
Params.NumAdaptBlocks       = 0;
Params.NumFixedBlocks       = 1;

% Cardinal Directions
Params.NumTrialsPerBlock    = 4;
Params.TargetOrder          = [1,3,5,7];
% Params.TargetOrder          = Params.TargetOrder(randperm(length(Params.TargetOrder)));  % randomize order
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
Params.TargetHoldTime       = 1;
Params.InterTrialInterval   = 2;
Params.InstructedDelayTime  = 1;
Params.CueTime              = 0.75;
Params.MaxStartTime         = 25;
Params.MaxReachTime         = 20;
Params.InterBlockInterval   = 10; % 0-10s, if set to 10 use instruction screen
Params.ImaginedMvmtTime     = 5;

%% Feedback
Params.FeedbackSound = false;
Params.ErrorWaitTime = 0;
Params.ErrorSound = 1000*audioread('buzz.wav');
Params.ErrorSoundFs = 8192;
[Params.RewardSound,Params.RewardSoundFs] = audioread('reward1.wav');
% play sounds silently once so Matlab gets used to it
sound(0*Params.ErrorSound,Params.ErrorSoundFs)

%% Task

Params.RunningModeBinNum    = 4;  % 1: No filtering, 3+: running mode filter of last n bins: Try 4 bins?
Params.RunningModeZero      = 0;  % 1: No motion if no winner, 0: maintain prior decision if no winner

Params.ValidDir          = [1:6,7];

% dynamics
Params.deltaT   = 1/Params.UpdateRate;
Params.k_v      = 0.85;
Params.k_i      = 24;   
Params.dA       = [1 0 0  Params.deltaT 0 0;...
                    0 1 0 0 Params.deltaT 0;...
                    0 0 1 0 0 Params.deltaT;...
                    0 0 0 Params.k_v 0 0;...
                    0 0 0 0 Params.k_v 0;...
                    0 0 0 0 0 Params.k_v];               
Params.dB       = [zeros(3); eye(3)];
Params.dB       = Params.dB*Params.k_i;

% % Clicker

Params.UseClicker       = 0;     % 0: trial ends with hold time, 1: trial ends with click
Params.ClickAction      = 7;     % action for click
Params.ClickerBinNum    = 3;
Params.ClickerBinThresh = 0.7;

Params.HoldTime = 1;

end % GetParams