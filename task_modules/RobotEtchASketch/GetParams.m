function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Experiment
Params.Task = 'EtchASketch';
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

% %% Neural network classifier option
% % set this to true to use neural network
% % also set the softmax option
% Params.NeuralNetFlag = true;
% if Params.NeuralNetFlag
%     Params.NeuralNetSoftMaxThresh = 0.5;       
%     Params.Use3Features = true;
%     Params.NeuralNetFunction = 'MLP_PreTrained_7DoF_Days1to11';
%     %Params.NeuralNetFunction = 'MLP_6DoF_PlusOK_Trained4mAllData_20210212';    
% 
% else
%     Params.NeuralNetSoftMaxThresh = 0;
% end

%% Neural network classifier option
% set this to true to use neural network
% also set the softmax option
Params.NeuralNetFlag = false;
if Params.NeuralNetFlag
    Params.NeuralNetSoftMaxThresh = 0.55;       
    Params.Use3Features = true;
%     Params.NeuralNetFunction = 'MLP_FlipView3D_20210817_PM1';
    Params.NeuralNetFunction = 'MLP_PreTrained_7DoF_PnP4';%'MLP_PreTrained_7DoF_PnP';
    
%     Params.NeuralNetFunction = 'MLP_PreTrained_7DoF_1006_AM2';
%     Params.NeuralNetFunction = 'multilayer_perceptron_6DoF_Online_Apr16_2021';
    %Params.NeuralNetFunction = 'MLP_6DoF_PlusOK_Trained4mAllData_20210212';    

else
    Params.NeuralNetSoftMaxThresh = 0;
end

%% Neural network 2 classifier option
% Trained in a different way using different optimizer

Params.NeuralNet2Flag = true;
if Params.NeuralNet2Flag
    Params.NeuralNet2SoftMaxThresh = 0.45    ;       
    Params.Use3Features = true;
    Params.NeuralNet2 = load(fullfile('clicker','net_mlp')); % 7DoF classifier trained in a different way
    
else
    Params.NeuralNet2SoftMaxThresh = 0;
end

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
    Params.ConvNeuralNetSoftMaxThresh = 0.4;       
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

Params.LongStartPos =  [200, 0, 0; 200, 0, 0; 0,0,-200; 0,0,-200; ...
    -200,0,0;-200,0,0;0,-200,-200; 0,200,-200];

% Params.Paths{1} = [[0.2, 0, 0]; [0.2, 0.0, -0.2]; [0.0,0,-0.2];  [0.0, 0.2,-0.2]]*1000;

% Params.ReachTargetPositions = [];

% Params.CorrectDecode{1} = [6*ones(1,10), 7*ones(1,10),3*ones(1,10), 2*ones(1,10)];
% pp = Params.Paths{1}(1,:) + [zeros(10,1), zeros(10,1),-[20:20:200]'];
% pp = [pp; pp(end,:) + [zeros(10,1),zeros(10,1), zeros(10,1)]];
% pp = [pp; pp(end,:) + [-[20:20:200]',zeros(10,1), zeros(10,1)]];
% pp = [pp; pp(end,:) + [zeros(10,1),[20:20:200]' zeros(10,1)]];
% 
% Params.PathPoints{1} = [Params.Paths{1}(1,:);pp];

Params.DirOrder = [[6,7,3,2]; [6,3,7,4];[5,4,7,3];[5,1,2,7];...
    [1,7,6,4];[2,6,7,1];[2,3,5,7]; [4,7,5,1]];

stepL = 20;
segL = 200;
numSteps = segL/stepL;
numSegs = 4;

numClickSteps = 5;

for i = 1:8
   pp = Params.LongStartPos(i,:);
   p = Params.LongStartPos(i,:);

        correctD = Params.DirOrder(i,:);
%         Params.CorrectDecode{i} = reshape([correctD'*ones(1,numSteps)]',numSteps*numSegs,1);
    Params.CorrectDecode{i} = [];
    for d = 1:numSegs
        dir = correctD(d);

        if dir < 7
            Params.CorrectDecode{i} = [Params.CorrectDecode{i}; dir*ones(numSteps,1)];
        else
            Params.CorrectDecode{i} = [Params.CorrectDecode{i}; dir*ones(numClickSteps,1)];
        end

        if dir == 1
            pp = [pp; pp(end,:) + [[20:20:200]',zeros(10,1), zeros(10,1)]];
        elseif dir == 3
            pp = [pp; pp(end,:) + [-[20:20:200]',zeros(10,1), zeros(10,1)]];
        elseif dir == 2
            pp = [pp; pp(end,:) + [zeros(10,1),[20:20:200]' zeros(10,1)]];
        elseif dir == 4
            pp = [pp; pp(end,:) + [zeros(10,1),-[20:20:200]' zeros(10,1)]];
        elseif dir == 5
            pp =[pp; pp(end,:) + [zeros(10,1), zeros(10,1),[20:20:200]']];
        elseif dir == 6
           pp = [pp; pp(end,:) + [zeros(10,1), zeros(10,1),-[20:20:200]']];
        elseif dir == 7
            pp = [pp; pp(end,:) + [zeros(numClickSteps,1),zeros(numClickSteps,1), zeros(numClickSteps,1)]];
            
            Params.ReachTargetPositions(i,:) = pp(end,:);
        end
        
        if dir ~= 7
            p = [p; pp(end,:)];
        end
    end
   

        Params.PathPoints{i} = [pp];
        Params.Paths{i} = p;
    
    
end



%% Kal,an Filter Properties
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

% Cardinal Directions
Params.NumTrialsPerBlock    = 8;
Params.TargetOrder          = [1:8];

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
Params.CueTime = 0.75;
Params.MaxStartTime = 25;
Params.MaxReachTime = 120 ;
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

Params.limit = [-256, 256; -256 256; -250 256];
Params.RobotMode            = 15;  
Params.RobotDirectionLines  = 1;  % 0: No lines, 1: Lines
Params.RunningModeBinNum    = 4;  % 1: No filtering, 3+: running mode filter of last n bins: Try 4 bins?
Params.RunningModeZero      = 1;  % 1: Ouput null state if no winner, 0: maintain prior decision if no winner
Params.RobotTargetDim = 2;


Params.RobotTargetDim = 1;

Params.ReachTargets      = [1,2,3,4,5,6];
Params.ValidDir          = [1:6,7];

Params.deltaT = 1/Params.UpdateRate;

Params.LongTrial = 1;

% Target
Params.RobotTargetRadius = 15;  % increase radius if task too hard
Params.TargetHoldTime = 10;

% Clicker
Params.RobotClicker = 1;     % 0: trial ends with hold time, 1: trial ends with click
Params.ClickerBinNum = 10;
Params.ClickerBinThresh = 0.7;
Params.RobotClickerStop = 0;  % 1: decode of 7 will set velocity to zero

Params.boundaryDist = 10;
Params.boundaryVel = 0;
Params.AssistAlpha = 0.0;

Params.PathInd = 0;


Params.GripperOrnInit = [0, pi/2,0];


%% For Debug

Params.index = 1;

d = Params.DirOrder(Params.TargetOrder(1),:);
Params.clickOrder = reshape([d'*ones(1,15)]',60,1);

% Params.clickOrder = [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7];
end % GetParams