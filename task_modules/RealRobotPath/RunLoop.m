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
    'CursorAssist',[],...
    'CursorState',[],...
    'IntendedCursorState',[],...
    'ClickerState',[],...
    'ClickerDistance',[],...
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
    'Events',[],...
    'FilteredClickerState',[],...
    'StopState',[],...
    'AssistVel',[],...
    'PathSegInd',[],...
    'CorrectDecode',[],...
    'ClickToSend',[],...
    'BetaScalar', [],...
    'BetaClickerState',[],...
    'LSTMFeatures',{{}});

switch TaskFlag,
    case 1, NumBlocks = Params.NumImaginedBlocks;
    case 2, NumBlocks = Params.NumAdaptBlocks;
    case 3, NumBlocks = Params.NumFixedBlocks;
end

%% Open UDP
% Params.udp = udpport("LocalPort", 43210);
% Params.pythonPort = 5006;
write(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);                  % reset robot

[xa,xb,xc] = doubleToUDP(Params.StartPos(1,1));
[ya,yb,yc] = doubleToUDP(Params.StartPos(1,2)); 
[za,zb,zc] = doubleToUDP(Params.StartPos(1,3) - 256) ;

% write(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, 0], "127.0.0.1", Params.pythonPort) ; % send pos
write(Params.udp, [0,2,Params.RobotMode,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,3,Params.RobotTargetRadius,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,4,Params.TargetHoldTime,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,7,Params.AssistAlpha*10,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,8,Params.AutoGrasp,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
% write(Params.udp, [0,9,Params.ClickerBinNum,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,10,Params.autoCenterOverTarget,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,11,Params.autoCenterDist,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,12,Params.wristStartX,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,13,Params.wristStartZ,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
% write(Params.udp, [0,14,Params.OperationModeReset,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
% write(Params.udp, [0,15,Params.zlim,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,16,Params.lowGainMode,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,17,Params.graspOrientation,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,18,Params.SwitchBinNum,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,19,Params.SwitchBinThresh*10,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,20,Params.GraspBinNum,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,21,Params.GraspBinThresh*10,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 

% set workspace limits
[xa,xb,xc] = doubleToUDP(Params.wl(1));
[ya,yb,yc] = doubleToUDP(Params.wl(2)); 
[za,zb,zc] = doubleToUDP(Params.wl(3)) ;
write(Params.udp, [0,22, xa,xb,xc,ya,yb,yc, za,zb,zc, 0], "127.0.0.1", Params.pythonPort) ; % send pos
[xa,xb,xc] = doubleToUDP(Params.wu(1));
[ya,yb,yc] = doubleToUDP(Params.wu(2)); 
[za,zb,zc] = doubleToUDP(Params.wu(3)) ;
write(Params.udp, [0,23, xa,xb,xc,ya,yb,yc, za,zb,zc, 0], "127.0.0.1", Params.pythonPort) ; % send pos

write(Params.udp, [0,26,Params.k_v*10,Params.k_i,Params.r_v*10,Params.r_i,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 


% pause(2.0)

%%  Loop Through Blocks of Trials
Trial = 0;
TrialBatch = {};
tlast = GetSecs;
Cursor.LastPredictTime = tlast;
Cursor.LastUpdateTime = tlast;
for Block=1:NumBlocks, % Block Loop
    NextTargetID =  Params.TargetOrder(Trial+1);
    % initialize cursor state(s)
    if Params.LongTrial
        Cursor.State = [0,0,0,0,0,0]';
        Cursor.State(1:3) = Params.LongStartPos(NextTargetID,:);
    else
        Cursor.State = [0,0,0,0,0,0]';
        Cursor.State(1:3) = Params.StartPos(Trial+1,:);
        Cursor.IntendedState = [0,0,0,0,1]';
    end
    Cursor.Vcommand = [0,0]';
    Cursor.ClickState = 0;
      
%     NextTargetID = 0;
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
         
        % update target and next target
        TargetID = NextTargetID;
%         while NextTargetID==TargetID,
%             NextTargetID = Params.ReachTargets(randperm(numel(Params.ReachTargets),1))
             NextTargetID =  Params.TargetOrder(Trial+1);
%         end
        
        % set up trial
        TrialData = DataFields;
        TrialData.Block = Block;
        TrialData.Trial = Trial;
        TrialData.TargetID = TargetID;
%         TrialData.TargetPosition = Params.ReachTargetPositions(TargetID,:);
        TrialData.NextTargetID = NextTargetID;
        TrialData.NextTargetPosition = Params.ReachTargetPositions(NextTargetID,:);
        
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
        
        [xa,xb,xc] = doubleToUDP(Params.StartPos(TargetID,1));
        [ya,yb,yc] = doubleToUDP(Params.StartPos(TargetID,2)); 
        [za,zb,zc] = doubleToUDP(Params.StartPos(TargetID,3) - 256) ;
        
        Params.StartPos(TargetID,:)

        write(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, 0], "127.0.0.1", Params.pythonPort) ; % send pos
        write(Params.udp, [0,2,Params.RobotMode,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
        write(Params.udp, [0,10,Params.autoCenterOverTarget,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
        
        [xa,xb,xc] = doubleToUDP(Params.StartWristX(TargetID));
        [ya,yb,yc] = doubleToUDP(Params.StartWristY(TargetID)); 
        [za,zb,zc] = doubleToUDP(Params.StartWristZ(TargetID)) ;
        
        write(Params.udp, [0,14,Params.OperationModeReset(TargetID),0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
        write(Params.udp, [0,31, xa,xb,xc,ya,yb,yc, za,zb,zc, 0], "127.0.0.1", Params.pythonPort) ; % send pos

        write(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);   
        [TrialData,Neuro,KF,Params,Clicker] = ...
            RunTrial(TrialData,Params,Neuro,TaskFlag,KF,Clicker);
        TrialData.TrialEndTime    = GetSecs;
                
        % Save Data from Single Trial
        
        save(...
            fullfile(DataDir,sprintf('Data%04i.mat',Trial)),...
            'TrialData',...
            '-v7.3','-nocompression');
        
    end % Trial Loop
    
    % Give Feedback for Block
    if Params.InterBlockInterval >= 10,
        Instructions = [...
            sprintf('\n\nFinished block %i of %i\n\n',Block,NumBlocks),...
            '\nPress the ''Space Bar'' to resume task.' ];
%         InstructionScreen(Params,Instructions)
    else,
        WaitSecs(Params.InterBlockInterval);
    end
    
end % Block Loop
%#ok<*NASGU>

end % RunLoop



