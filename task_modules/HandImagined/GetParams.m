function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Experiment
Params.Task = 'Hand';
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
Params.NumImaginedBlocks    = 0;
Params.NumAdaptBlocks       = 0;
Params.NumFixedBlocks       = 1;

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
Params.limit = [-400, 400; -400 400; -350 350];
Params.RobotTargetRadius    = 40;
Params.RobotMode            = 3;  % 0: Horizontal, 1: Vertical+Gripper, 3: 3D robot 
Params.RobotDirectionLines  = 1;  % 0: No lines, 1: Lines
Params.RunningModeBinNum    = 3;  % 1: No filtering, 3+: running mode filter of last n bins: Try 4 bins?
Params.RunningModeZero      = 3;  % 1: No motion if no winner, 0: maintain prior decision if no winner

if Params.RobotMode == 0
    Params.RobotTargetDim = 2;
elseif Params.RobotMode == 1
    Params.RobotTargetDim = 1;
end

Params.RobotTargetRadius = 10;
Params.RobotTargetDim = 1;

Params.ReachTargets      = [1,2,3,4,5,6];
Params.ValidDir          = [1:6,7];

Params.deltaT = 0.1;
Params.k_v = 0.9;
Params.k_i = 10.0;

Params.dA = [1 0 0  Params.deltaT 0 0;...
                    0 1 0 0 Params.deltaT 0;...
                    0 0 1 0 0 Params.deltaT;...
                    0 0 0 Params.k_v 0 0;...
                    0 0 0 0 Params.k_v 0;...
                    0 0 0 0 0 Params.k_v];
                
Params.dB = [zeros(3);...
                    eye(3)];
Params.dB = Params.dB*Params.k_i;

Params.LongTrial = 0;
Params.LongStartPos =  [Params.ReachTargetPositions(3,:);...
    Params.ReachTargetPositions(4,:);...
    Params.ReachTargetPositions(1,:);...
    Params.ReachTargetPositions(2,:);...
    Params.ReachTargetPositions(6,:);...
    Params.ReachTargetPositions(5,:);...
    Params.ReachTargetPositions(9,:);...
    Params.ReachTargetPositions(10,:);...
    Params.ReachTargetPositions(7,:);...
    Params.ReachTargetPositions(8,:)];

Params.RobotClicker = 1;
Params.TargetHoldTime = 3;


%% Hand
Params.NumTrialsPerBlock    = 20;
Params.TargetOrder          = [7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8];
% Params.TargetOrder          = [3,6,3,6,3,6,3,6,3,6,3,6,3,6,3,6,3,6,3,6];

% Params.TargetOrder = Params.TargetOrder(randperm(length(Params.TargetOrder)));  % randomize order
Params.TargetOrder          = [Params.TargetOrder, 1];

Params.wristAxisLim = [-0.5, 0.5;...
                        -0.5, 1.2;...
                        -0.5, 3.15;...
                        0.5, -0.5;...
                         1.2, -0.5;...
                        3.15, -0.5;
                        -1.0, 0.1;...
                        0.1, -1.0];
                    
Params.TrialDur = 1.5;                        
mult = abs(Params.wristAxisLim(:,1) - Params.wristAxisLim(:,2))/Params.UpdateRate/Params.TrialDur;                 

Params.axes         = [1,2,3,1,2,3,1,1];
Params.rotInc       = [1,1,1,-1,-1,-1,1,-1].*mult';
Params.rotDir       = 1;
Params.trialRepeat  = 1;
Params.trialPause   = 0.5;
Params.handVis = 0;

for i = 1:8
    Params.angles{i} = [];
    Params.dir{i} = [];
    for j = 1:Params.trialRepeat
       if Params.rotDir == 1 
            Params.angles{i} = [Params.angles{i} , Params.wristAxisLim(i,1):Params.rotInc(i):Params.wristAxisLim(i,2)];
            Params.dir{i} = [Params.dir{i}, ones(1,length(Params.wristAxisLim(i,1):Params.rotInc(i):Params.wristAxisLim(i,2)))*1];
            Params.angles{i} = [Params.angles{i}, ones(1, round(Params.trialPause* Params.UpdateRate))*Params.angles{i}(end)];
            Params.dir{i} = [Params.dir{i}, zeros(1, round(Params.trialPause* Params.UpdateRate))];
       else if Params.rotDir == -1
            Params.angles{i} = [Params.angles{i} , Params.wristAxisLim(i,2):-Params.rotInc(i):Params.wristAxisLim(i,1)];
            Params.dir{i} = [Params.dir{i}, ones(1,length(Params.wristAxisLim(i,2):-Params.rotInc(i):Params.wristAxisLim(i,1)))*-1];
            Params.angles{i} = [Params.angles{i}, ones(1, round(Params.trialPause* Params.UpdateRate))*Params.angles{i}(end)];
            Params.dir{i} = [Params.dir{i}, zeros(1, round(Params.trialPause* Params.UpdateRate))];
       else
            Params.angles{i} = [Params.angles{i} , Params.wristAxisLim(i,1):Params.rotInc(i):Params.wristAxisLim(i,2)];
            Params.dir{i} = [Params.dir{i}, ones(1,length(Params.wristAxisLim(i,1):Params.rotInc(i):Params.wristAxisLim(i,2)))*1];
            Params.angles{i} = [Params.angles{i}, ones(1, round(Params.trialPause* Params.UpdateRate))*Params.angles{i}(end)];
            Params.dir{i} = [Params.dir{i}, zeros(1, round(Params.trialPause* Params.UpdateRate))];
            Params.angles{i} = [Params.angles{i} , Params.wristAxisLim(i,2):-Params.rotInc(i):Params.wristAxisLim(i,1)];
            Params.dir{i} = [Params.dir{i}, ones(1,length(Params.wristAxisLim(i,2):-Params.rotInc(i):Params.wristAxisLim(i,1)))*-1];
            Params.angles{i} = [Params.angles{i}, ones(1, round(Params.trialPause* Params.UpdateRate))*Params.angles{i}(end)];
            Params.dir{i} = [Params.dir{i}, zeros(1, round(Params.trialPause* Params.UpdateRate))];
       end
   end
    end
end

a = 1;
end % GetParams
