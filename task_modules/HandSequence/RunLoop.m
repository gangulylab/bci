function [Neuro,KF,Params,Clicker] = RunLoop(Params,Neuro,TaskFlag,DataDir,KF,Clicker)
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
    'TargetPosition',NaN,...
    'NextTargetID',NaN,...
    'NextTargetPosition',NaN,...
    'SelectedTargetID',NaN,...
    'SelectedTargetPosition',NaN,...
    'Time',[],...
    'ChStats',[],...
    'FeatureStats',[],...
    'HandState',[],...
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
    'ErrorID',0,...
    'ErrorStr','',...    
    'Events',[],...
    'Action',[]...
    );

switch TaskFlag
    case 1, NumBlocks = Params.NumImaginedBlocks;
    case 2, NumBlocks = Params.NumAdaptBlocks;
    case 3, NumBlocks = Params.NumFixedBlocks;
end

%% Open UDP
Params.udp = udp("127.0.0.1", 5006);
fopen(Params.udp)
fwrite(Params.udp, [0,1,0])                  % reset robot
fwrite(Params.udp, [0, 3, Params.showDecodeLines])
% fwrite(Params.udp, [0,2,Params.UpdateRate])  % set update rate
% fwrite(Params.udp, [0,3,Params.RobotMode])   % set robot mode
% fwrite(Params.udp, [0,4,Params.RobotDirectionLines])   % set debug lines

%%  Loop Through Blocks of Trials
Trial = 0;
TrialBatch = {};
tlast = GetSecs;
Cursor.LastPredictTime = tlast;
Cursor.LastUpdateTime = tlast;
for Block=1:NumBlocks % Block Loop

    for TrialPerBlock=1:Params.NumTrialsPerBlock % Trial Loop
        % if smooth batch on & enough time has passed, update KF btw trials
        if TaskFlag==2 && Neuro.CLDA.Type==2
            TrialBatch{end+1} = sprintf('Data%04i.mat', Trial);
            if (GetSecs-tlast)>Neuro.CLDA.UpdateTime
                Neuro.KF.CLDA = Params.CLDA;
                if Neuro.DimRed.Flag
                    KF = FitKF(Params,fullfile(Params.Datadir,'BCI_CLDA'),2,...
                        KF,TrialBatch,Neuro.DimRed.F);
                else
                    KF = FitKF(Params,fullfile(Params.Datadir,'BCI_CLDA'),2,...
                        KF,TrialBatch);
                end
                tlast = GetSecs;
                TrialBatch = {};
                % decrease assistance after batch update
                if Cursor.Assistance>0
                    Cursor.Assistance = Cursor.Assistance - Cursor.DeltaAssistance;
                    Cursor.Assistance = max([Cursor.Assistance,0]);
                end
            end
        elseif TaskFlag==2 && Neuro.CLDA.Type==3
            % decrease assistance after batch update
            if Cursor.Assistance>0
                Cursor.Assistance = Cursor.Assistance - Cursor.DeltaAssistance;
                Cursor.Assistance = max([Cursor.Assistance,0]);
            end
        end
        
        % update trial
        Trial = Trial + 1;
        
        % set up trial
        TrialData       = DataFields;
        TrialData.Block = Block;
        TrialData.Trial = Trial;
        TrialData.TrialActions = Params.BlockTargetOrder{Trial};
        
        % save ch stats and feature stats in each trial
        TrialData.ChStats.Mean = Neuro.ChStats.mean;
        TrialData.ChStats.Var = Neuro.ChStats.var;
        TrialData.FeatureStats.Mean = Neuro.FeatureStats.mean;
        TrialData.FeatureStats.Var = Neuro.FeatureStats.var;
        
        % Run Trial
        TrialData.TrialStartTime  = GetSecs;
        [TrialData,Neuro,KF,Params,Clicker] = ...
            RunTrial(TrialData,Params,Neuro,TaskFlag,KF,Clicker);
        TrialData.TrialEndTime    = GetSecs;
                
        % Save Data from Single Trial
        save(...
            fullfile(DataDir,sprintf('Data%04i.mat',Trial)),...
            'TrialData',...
            '-v7.3','-nocompression');
        
        % keep track of useful stats and params
%         SavePersistence(Params,Neuro,KF,TaskFlag)
        
    end % Trial Loop
    
    % Give Feedback for Block
    if Params.InterBlockInterval >= 10
        Instructions = [...
            sprintf('\n\nFinished block %i of %i\n\n',Block,NumBlocks),...
            '\nPress the ''Space Bar'' to resume task.' ];
%         InstructionScreen(Params,Instructions)
    else
        WaitSecs(Params.InterBlockInterval);
    end
    
end % Block Loop
%#ok<*NASGU>

end % RunLoop


