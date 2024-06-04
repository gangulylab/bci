function [Neuro,KF,Params] = RunLoop(Params,Neuro,TaskFlag,DataDir,KF)
% Defines the structure of collected data on each trial
% Loops through blocks and trials within blocks

global Cursor

%% Start Experiment
DataFields = struct(...
    'Params',Params,...
    'Block',NaN,...
    'Trial',NaN,...
    'TrialStartTime',NaN,...
    'TrialEndTime',NaN,...
    'TargetID',NaN,...
    'TargetAngle',NaN,...
    'TargetPosition',NaN,...
    'Time',[],...
    'ChStats',[],...
    'FeatureStats',[],...
    'CursorAssist',[],...
    'CursorState',[],...
    'IntendedCursorState',[],...
    'TaskState',[],...
    'NeuralTime',[],...
    'NeuralTimeBR',[],...
    'NeuralSamps',[],...
    'NeuralFeatures',{{}},...
    'SmoothedNeuralFeatures',{{}},...
    'NeuralFactors',{{}},...
    'BroadbandData',{{}},...
    'Reference',{{}},...
    'ProcessedData',{{}},...
    'KalmanFilter',{{}},...
    'KalmanGain',{{}},...
    'ErrorID',0,...
    'ErrorStr','',...
    'Events',[]...
    );

switch TaskFlag,
    case 1, NumBlocks = Params.NumImaginedBlocks;
    case 2, NumBlocks = Params.NumAdaptBlocks;
    case 3, NumBlocks = Params.NumFixedBlocks;
end

%%  Loop Through Blocks of Trials
Trial = 0;
TrialBatch = {};
tlast = GetSecs;
Cursor.LastPredictTime = tlast;
Cursor.LastUpdateTime = tlast;
for Block=1:NumBlocks, % Block Loop

    % random order of reach targets for each block
    switch Params.TargetSelectionFlag,
        case {1,2,4},
            TargetOrder = Params.TargetFunc(Params.NumTrialsPerBlock);
        case 3,
            TargetOrder = Params.TargetFunc(Block,Params.NumTrialsPerBlock);
    end
    Cursor.State = [0,0,0,0,1]';
    Cursor.IntendedState = [0,0,0,0,1]';
    Cursor.Vcommand = [0,0]';

    for TrialPerBlock=1:Params.NumTrialsPerBlock, % Trial Loop
        % if smooth batch on & enough time has passed, update KF btw trials
        if TaskFlag==2 && Neuro.CLDA.Type==2,
            TrialBatch{end+1} = sprintf('Data%04i.mat', Trial);
            if (GetSecs-tlast)>Neuro.CLDA.UpdateTime,
                Neuro.KF.CLDA = Params.CLDA;
                if Neuro.DimRed.Flag,
                    KF = FitKF(Params,fullfile(Params.Datadir,'BCI_CLDA'),2,...
                        KF,TrialBatch,Neuro.DimRed.F);
                else,
                    KF = FitKF(Params,fullfile(Params.Datadir,'BCI_CLDA'),2,...
                        KF,TrialBatch);
                end
                tlast = GetSecs;
                TrialBatch = {};
                % decrease assistance after batch update
                if Cursor.Assistance>0,
                    Cursor.Assistance = Cursor.Assistance - Cursor.DeltaAssistance;
                    Cursor.Assistance = max([Cursor.Assistance,0]);
                end
            end
        elseif TaskFlag==2 && Neuro.CLDA.Type==3,
            % decrease assistance after batch update
            if Cursor.Assistance>0,
                Cursor.Assistance = Cursor.Assistance - Cursor.DeltaAssistance;
                Cursor.Assistance = max([Cursor.Assistance,0]);
            end
        end
        
        % update trial
        Trial = Trial + 1;
        TrialIdx = TargetOrder(TrialPerBlock);
        
        % set up trial
        TrialData = DataFields;
        TrialData.Block = Block;
        TrialData.Trial = Trial;
        TrialData.TargetID = TrialIdx;
        TrialData.TargetAngle = Params.ReachTargetAngles(TrialIdx);
        TrialData.TargetPosition = Params.ReachTargetPositions(TrialIdx,:);
        % save kalman filter
        if Params.ControlMode>=3 && TaskFlag>=2 && ~Params.SaveKalmanFlag,
            TrialData.KalmanFilter{1}.A = KF.A;
            TrialData.KalmanFilter{1}.W = KF.W;
            TrialData.KalmanFilter{1}.C = KF.C;
            TrialData.KalmanFilter{1}.Q = KF.Q;
            TrialData.KalmanFilter{1}.P = KF.P;
            TrialData.KalmanFilter{1}.Lambda = KF.Lambda;
        end
        % save ch stats and feature stats in each trial
        TrialData.ChStats.Mean = Neuro.ChStats.mean;
        TrialData.ChStats.Var = Neuro.ChStats.var;
        TrialData.FeatureStats.Mean = Neuro.FeatureStats.mean;
        TrialData.FeatureStats.Var = Neuro.FeatureStats.var;
        
        % Run Trial
        TrialData.TrialStartTime  = GetSecs;
        [TrialData,Neuro,KF,Params] = RunTrial(TrialData,Params,Neuro,TaskFlag,KF);
        TrialData.TrialEndTime    = GetSecs;
                
        % Save Data from Single Trial
        save(...
            fullfile(DataDir,sprintf('Data%04i.mat',Trial)),...
            'TrialData',...
            '-v7.3','-nocompression');
        
        % keep track of useful stats and params
        SavePersistence(Params,Neuro,KF,TaskFlag)
        
    end % Trial Loop
    
    % Give Feedback for Block
    if Params.InterBlockInterval >= 10,
        Instructions = [...
            sprintf('\n\nFinished block %i of %i\n\n',Block,NumBlocks),...
            '\nPress the ''Space Bar'' to resume task.' ];
        InstructionScreen(Params,Instructions)
    else,
        WaitSecs(Params.InterBlockInterval);
    end
    
end % Block Loop
%#ok<*NASGU>

end % RunLoop



