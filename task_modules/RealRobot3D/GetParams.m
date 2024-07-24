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
Params.BadChannels          = [14,20,21,22];
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
Params.SmoothDataFlag = true;
Params.FeatureBufferSize = 5;

%% Timing
Params.ScreenRefreshRate = 5; % Hz
Params.UpdateRate = 5; % Hz

%% Discrete Decoder name
Params.UseSVM = false;
Params.DiscreteDecoder = 'clicker_svm_mdl_6Dir_3Feat_462021.mat';

%% Multi State Decision Boundary for Null Class

% set this to negative values between -2 and -4 IMPORTANT: if this is set
% to large -ve values (enforcing a null class), don't set
% Params.ClickCounter to be > 4 bins, else it might get too difficult for
% subject.

Params.MultiDecisionBoundary =-2; 

%% Neural network classifier option
% set this to true to use neural network
% also set the softmax option
Params.NeuralNetFlag = true;
if Params.NeuralNetFlag
   Params.NeuralNetSoftMaxThresh = 0.4;
    Params.NeuralNetFunction = 'MLP_7Dir_B1_20240621_CL4_NoPooling';
   %Params.NeuralNetFunction = 'MLP_4Dir_Imagined_20210217_Day3_AllFeat';
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
%         Params.NeuralNet2FileName = 'net_PnP';
        Params.NeuralNet2FileName = 'net_PnP_RobotBatchUpdate';
        Params.NeuralNet2 = load(fullfile('clicker',Params.NeuralNet2FileName));
%         Params.NeuralNet2 = Params.NeuralNet2.net_PnP;
        Params.NeuralNet2 = Params.NeuralNet2.net_PnP_RobotBatchUpdate;
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

% %% IMPORT PC WEIGHTS AND MEAN FOR BETA BAND ANALYSIS
% 
% Params.BetaWts = load(fullfile('clicker','betawts_stop'));
% Params.BetaMean = load(fullfile('clicker','betamean'));
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
Params.CueTime = 0.75;
Params.MaxStartTime = 50;
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
Params.RobotMode    = 6; 
Params.wl           = [-50, -67, 10];
Params.wu           = [5, -15 55];
Params.view = 1;   % 1  = far-side of table, 2 = near-side of table % 3 = drawer

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
    Params.UseHeightDist  = 0;
    Params.UseNewAutoGrasp = 0;
    Params.WaitForGraspSignal = 0;
elseif Params.RobotMode == 2 % vertical R2G
    Params.ValidDir          = [1:9];
%     Params.StartPos          = [-250, 0,300];
        Params.StartPos          = [-250, 0,350];
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
    
    Params.wl           = [-50, -67, 20];
    Params.wu           = [5, -15 55];
    Params.UseNewAutoGrasp = 0;
    Params.WaitForGraspSignal = 0;
    Params.UseHeightDist  = 0;
    Params.AutoGraspHorzDist = 10;
    Params.AutoGraspVertDist = 15;
    
    
    Params.AssistMode = 0;
    
elseif Params.RobotMode == 3  % lateral R2G wall
    Params.ValidDir             = [1:9];
%     Params.StartPos             = [-200, 200, 220]; % vertically alligned start
    Params.StartPos             = [100, -70, 380]; %
    Params.NumTrialsPerBlock    = 1;
    Params.TargetOrder          = [1];   
    Params.OperationModeReset   = 0;
    Params.wristStartX          = 3.1415/2*10; 
    Params.wristStartZ          = 0; 
    Params.graspOrientation     = 1;   
    Params.StartWristX = [3.1415/2]*10;                    
    Params.StartWristZ = 10*[0];
    Params.StartWristY = 10*[0]; 
    Params.UseNewAutoGrasp = 0;
    Params.WaitForGraspSignal = 0;
    Params.UseNewAutoGrasp = 0;
    Params.WaitForGraspSignal = 0;
    Params.UseHeightDist  = 0;
    Params.AutoGraspHorzDist = 10;
    Params.AutoGraspVertDist = 15;
Params.autoCenterDist       = 5;
    
elseif Params.RobotMode == 4 % TopDown with assist
    Params.ValidDir          = [1:9];
    Params.StartPos          = [-250, 0,300];  % standard
%     Params.StartPos          = [0, 250,300];  % 
    Params.NumTrialsPerBlock    = 1;
    Params.TargetOrder          = [1];   
    Params.OperationModeReset = 0;
    Params.wristStartX = 3.1415*10; 
    Params.wristStartZ = 0; 
    Params.autoCenterOverTarget = 0;
    Params.autoCenterDist = 10;
    Params.zlim = 5;
    Params.graspOrientation = 0;
    
    Params.StartWristX = [3.1415]*10;                    
    Params.StartWristZ = 10*[0];
    Params.StartWristY = 10*[0];    
    Params.UseNewAutoGrasp = 0;
    Params.UseHeightDist  = 1;
    Params.AutoGraspHorzDist = 10;
    Params.AutoGraspVertDist = 15;
    Params.WaitForGraspSignal   = 1;
    Params.wl           = [-50, -67, 16];
    Params.wu           = [5, -15 50];
    
elseif Params.RobotMode == 5 % Three object assist
    Params.ValidDir          = [1:9];
    Params.StartPos          = [-250, 0,350];
    Params.NumTrialsPerBlock    = 1;
    Params.TargetOrder          = [1];   
    Params.OperationModeReset = 0;
%     Params.wristStartX = 3.1415*10; 
    Params.wristStartZ = 0; 
    Params.autoCenterOverTarget = 0;
    Params.autoCenterDist = 10;
    Params.zlim = 5;
    Params.graspOrientation = 1;
    
%     Params.StartWristX = [3.1415]*10;    
    Params.StartWristX = [3.1415]*10;
    Params.wristStartX = 3.1415*10;
    Params.StartWristZ = 10*[0];
    Params.StartWristY = 10*[0];    
    Params.UseNewAutoGrasp = 2;
    Params.UseHeightDist  = 1;
    Params.AutoGraspHorzDist = 10;
    Params.AutoGraspVertDist = 15;
    Params.WaitForGraspSignal   = 1;
    Params.wl           = [-50, -67, 16];
    Params.wu           = [5, -15 50];
    
    Params.AssistMode = 2;  %1: blended,  2: flexible autocomplete
    Params.view = 2;   % 1  = far-side of table, 2 = near-side of table % 3 = drawer
    
elseif Params.RobotMode == 6 % Drawer Task
    Params.ValidDir          = [1:9];
    Params.StartPos          = [100, -250, 400];
%     Params.StartPos          = [100, -200, 200];
    Params.NumTrialsPerBlock    = 1;
    Params.TargetOrder          = [1];   
    Params.OperationModeReset = 0;
    Params.wristStartX = 3.1415*10; 
    Params.wristStartZ = 0; 
    Params.autoCenterOverTarget = 0;
    Params.autoCenterDist = 10;
    Params.zlim = 5;
    Params.graspOrientation = 0;
    
%     Params.StartWristX = [3.1415]*10;    
    Params.StartWristX = [3.1415]/2*10;
    Params.wristStartX = 3.1415/2*10;
    Params.StartWristZ = 10*[0];
    Params.StartWristY = 10*[0];    
    Params.UseNewAutoGrasp = 2;
    Params.UseHeightDist  = 1;
    Params.AutoGraspHorzDist = 10;
    Params.AutoGraspVertDist = 15;
    Params.WaitForGraspSignal   = 1;
    Params.wl           = [-25, -75, 10];
    Params.wu           = [25, -35,  70];
    
    Params.AssistMode = 2;  %1: blended,  2: flexible autocomplete

    % belief
    Params.BeliefThresh = 0.25;
    Params.distB = .4;
    Params.distK = 3;
    Params.velB = 0.4;
    Params.velK = 0.4;
    Params.pDiag = 0.8;
    Params.slowThresh = 1;
    
    Params.view =3;   % 1  = far-side of table, 2 = near-side of table % 3 = drawer
end

% Params.GoalPos1     = [0, 50, 220];  %red short
Params.GoalPos1     = [-100, -150, 220];  %green short
Params.GoalPos2     = [100, -180, 170]; %red diag
Params.GoalAng      = [pi/2];
Params.index        = 1;


Params.clickOrder = [4*ones(20,1);1*ones(60,1); 7*ones(10,1); 3*ones(5,1);7*ones(10,1); 1*ones(100,1)];
Params.clickOrder = [0*ones(5000,1)];
Params.ReachTargets      = [1,2,3,4,5,6];
Params.TargetOrder  = [Params.TargetOrder, 1];

Params.RunningModeBinNum    = 5;  % 1: No filtering, 3+: running mode filter of last n bins: Try 4 bins?
Params.RunningModeZero      = 1;  % 1: 0 output if not winner, 0: maintain prior decision if no winner

Params.RobotTargetRadius    = 50;
Params.RobotTargetDim       = 1;

% Robot Dynamics
Params.deltaT   = 1/Params.UpdateRate;
% Params.k_v      = 0.7;
% Params.k_i      = 40;
Params.k_v      = 0.8;
Params.k_i      = 18;    % standard: 18

Params.r_v      = 0.75;
Params.r_i      = 90;

% Params.r_v      = 0.7;
% Params.r_i      = 90;
% Params.r_v      = 0.7;
% Params.r_i      = 150;

Params.LongTrial = 0;
Params.RobotClicker     = 1;
Params.TargetHoldTime   = 0.25;

Params.boundaryDist     = 0;
Params.boundaryVel      = 0;
Params.AssistAlpha      = 0.3;
Params.AutoGrasp        = 0;
Params.GraspTask        = 1;
Params.lowGainMode      = 0;
Params.autoCenterOverTarget    = 0;

Params.SwitchBinNum     = 8;
Params.SwitchBinThresh  = 0.74;
Params.GraspBinNum      = 8;
Params.GraspBinThresh   = 0.74;
Params.InputBinNum      = 5;
Params.InputBinThresh   = 0.7;

% Beta
Params.UseBetaStop      = 0;
Params.BetaThreshold    = 0.5;

%Flip
Params.FlipStop         = 0;
Params.FlipBinNum       = 5*8; 
Params.FlipBinThresh    = 3;

% Mode switch sound
Params.UseSoundModeSwitch = 0;

% Lock gripper closed
Params.LockGripperClosed = 1;  % 0: do not lock, 1: lock 

% Slow dynamics when outstretched
Params.UseSlowMode = 1;
Params.SlowDistanceThrehold = 65; % cm
Params.k_v_slow      = 0.7;
Params.k_i_slow      = 12;  

Params.UseModeSwitch = 0;
Params.PreviewTarget = 1;

% Display settings for Eye-gaze Sync
Params.ShowFlash        = 1;
Params.FlashDuration    = 50; %ms

end % GetParams
