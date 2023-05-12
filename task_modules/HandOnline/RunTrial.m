function [Data, Neuro, KF, Params, Clicker] = RunTrial(Data,Params,Neuro,TaskFlag,KF,Clicker)
% Runs a trial, saves useful data along the way
% Each trial contains the following pieces
% 1) Get the cursor to the reach target (different on each trial)
% 2) Feedback

global Cursor

%% Set up trial
% ReachTargetPos = Data.TargetPosition;
TargetID = 0; % Target that cursor is in, 0 for no targets

TargetNames = {'Thumb', 'Index', 'Middle', 'Ring', 'Pinky', 'Power Grasp',...
    'Pinch Grasp', 'Tripod Grasp', 'Wrist Out', 'Wrist In','Wrist Flex',...
    'Wrist Extend', 'Wrist Pronate', 'Wrist Supinate'};

% Output to Command Line
fprintf('\nTrial: %i\n',Data.Trial)
fprintf('  Target: %s\n',TargetNames{Data.TargetID})
if Params.Verbose
    if TaskFlag==2
        fprintf('    Cursor Assistance: %.2f%%\n',100*Cursor.Assistance)
        if Params.CLDA.Type==3
            fprintf('    Lambda 1/2 life: %.2fsecs\n',KF.Lambda)
        end
    end
end

% keep track of update times
dt_vec = [];
dT_vec = [];

% grab blackrock data and run through processing pipeline
if Params.BLACKROCK
    Cursor.LastPredictTime = GetSecs;
    Cursor.LastUpdateTime = Cursor.LastPredictTime;
    Neuro = NeuroPipeline(Neuro,[],Params);
end

axis = 1;

if Params.d1(Data.TargetID)
    Params.angles = Params.angles1d;
end

Cursor.State = zeros(8,1);
TargetProgress = zeros(14,1);
Cursor.State(axis) = Params.angles(1);

Cursor.ClickState = 0;
Cursor.ClickDistance = 0;
inTargetOld = 0;

fwrite(Params.udp, [0,16,Data.TargetID]);

%% Instructed Delay
if ~Data.ErrorID && Params.InstructedDelayTime>0
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Instructed Delay';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
    
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;
    
%     Cursor.State(axis) = Params.angles(step);

    [xa,xb,xc] = doubleToUDP(Cursor.State(1)*80);
    [ya,yb,yc] = doubleToUDP(Cursor.State(2)*80); 
    [za,zb,zc] = doubleToUDP(Cursor.State(3)*80) ;
    if Params.handVis
        fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, Data.TargetID]);
        fwrite(Params.udp, [0,2,0])
    end
    
    while ~done
        % Update Time & Position
        tim = GetSecs;
        
        % for pausing and quitting expt
        if CheckPause, [Neuro,Data,Params] = ExperimentPause(Params,Neuro,Data); end
        
        % Update Screen
        if (tim-Cursor.LastPredictTime) > 1/Params.ScreenRefreshRate
            % time
            dt = tim - Cursor.LastPredictTime;
            TotalTime = TotalTime + dt;
            dt_vec(end+1) = dt;
            Cursor.LastPredictTime = tim;
            Data.Time(1,end+1) = tim;
            
            % grab and process neural data
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate)
                dT = tim-Cursor.LastUpdateTime;
                dT_vec(end+1) = dT;
                Cursor.LastUpdateTime = tim;
                
                Data.NeuralTime(1,end+1) = tim;
                [Neuro,Data] = NeuroPipeline(Neuro,Data,Params); 
            end
            
            Data.CursorState(:,end+1) = Cursor.State;
            Data.IntendedCursorState(:,end+1) = Cursor.IntendedState;
            Data.CursorAssist(1,end+1) = Cursor.Assistance;
            
            Cursor.TaskState = 1;
            Data.TaskState(1,end+1)=Cursor.TaskState;
           
            Data.StopState(1,end+1)=0;
            
            % start counting time            
            InTargetTotalTime = InTargetTotalTime + dt;
        end
        
        % end if in start target for hold time
        if InTargetTotalTime > Params.InstructedDelayTime
            done = 1;
        end
    end % Instructed Delay Loop
end % only complete if no errors


%% Cue time
if ~Data.ErrorID && Params.CueTime>0
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Cue';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
       
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;
    ClickDec_Buffer     = zeros(Params.RunningModeBinNum, 1);
    fwrite(Params.udp, [5, Data.TargetID, 0])
    
    while ~done
        % Update Time & Position
        tim = GetSecs;
        
        % for pausing and quitting expt
        if CheckPause, [Neuro,Data,Params] = ExperimentPause(Params,Neuro,Data); end
        
        % Update Screen
        if (tim-Cursor.LastPredictTime) > 1/Params.ScreenRefreshRate,
            % time
            dt = tim - Cursor.LastPredictTime;
            TotalTime = TotalTime + dt;
            dt_vec(end+1) = dt;
            Cursor.LastPredictTime = tim;
            Data.Time(1,end+1) = tim;
            
            % grab and process neural data
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate)
                dT = tim-Cursor.LastUpdateTime;
                dT_vec(end+1) = dT;
                Cursor.LastUpdateTime = tim;
                
                Data.NeuralTime(1,end+1) = tim;
                [Neuro,Data] = NeuroPipeline(Neuro,Data,Params);           
            end

            Data.CursorState(:,end+1) = Cursor.State;
            Data.IntendedCursorState(:,end+1) = Cursor.IntendedState;
            Data.CursorAssist(1,end+1) = Cursor.Assistance;
            
            Cursor.TaskState = 2;
            Data.TaskState(1,end+1)=Cursor.TaskState;      
            Data.StopState(1,end+1)=0;
            
            % start counting time            
            InTargetTotalTime = InTargetTotalTime + dt;
           
        end
        
        % end if in start target for hold time
        if InTargetTotalTime > Params.CueTime
            done = 1;
        end
    end % Instructed Delay Loop
end % only complete if no errors

%% Go to reach target

if ~Data.ErrorID
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Reach Target';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
    
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;
    
    ClickDec_Buffer = zeros(Params.RunningModeBinNum, 1);
    temp_dir = [0,0,0];
    ClickToSend = 0;
    
    step = 0;
    fwrite(Params.udp, [5, Data.TargetID, 1]) % DISPLAY GO CUE
    
    while ~done
        
        % Update Time & Position
        tim = GetSecs;
        
        % for pausing and quitting expt
        if CheckPause, [Neuro,Data,Params] = ExperimentPause(Params,Neuro,Data); end
        
        % Update Screen
        if (tim-Cursor.LastPredictTime) > 1/Params.ScreenRefreshRate
            step = step + 1;
            % time
            dt = tim - Cursor.LastPredictTime;
            TotalTime = TotalTime + dt;
            dt_vec(end+1) = dt;
            Cursor.LastPredictTime = tim;
            Data.Time(1,end+1) = tim;
            
            % grab and process neural data
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate)
                dT = tim-Cursor.LastUpdateTime;
                dT_vec(end+1) = dT;
                Cursor.LastUpdateTime = tim;
                
                Data.NeuralTime(1,end+1) = tim;
                [Neuro,Data] = NeuroPipeline(Neuro,Data,Params);              
                
                % save kalman filter
                if Params.ControlMode>=3 && TaskFlag>1 && Params.SaveKalmanFlag,
                    Data.KalmanGain{end+1} = [];
                    Data.KalmanGain{end}.K = KF.K;
                    Data.KalmanFilter{end+1} = [];
                    Data.KalmanFilter{end}.C = KF.C;
                    Data.KalmanFilter{end}.Q = KF.Q;
                    Data.KalmanFilter{end}.Lambda = KF.Lambda;
                end
            end

            Params.TargetID =  Data.TargetID;
            [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params,Neuro,Clicker);
            %Click_Decision = randperm(12,1);
            targetInd = [1:14];
                
            % Mode filter
            
            Cursor.ClickState = Click_Decision;
            Cursor.ClickDistance = Click_Distance;
            Data.ClickerState(1,end+1) = Cursor.ClickState;
            Data.ClickerDistance(1,end+1) = Cursor.ClickDistance;

            ClickDec_Buffer(1:end-1) = ClickDec_Buffer(2:end);
            ClickDec_Buffer(end) = Click_Decision;
            RunningMode_ClickDec = RunningMode(ClickDec_Buffer);  
            %RunningMode_ClickDec = randperm(12,1);

            ClickToSend = RunningMode_ClickDec;
            Data.FilteredClickerState(1,end+1) = RunningMode_ClickDec;

            
            if RunningMode_ClickDec == targetInd(Data.TargetID)
                step = Params.correctStep;
            else
                step = Params.incorrectStep;
            end    
            
            % Target Progress
            if RunningMode_ClickDec > 0
                TargetProgress(RunningMode_ClickDec) = TargetProgress(RunningMode_ClickDec) + 1;
            end
            Data.TargetProgress(:,end+1) = TargetProgress;
            
            % Finger State
            if RunningMode_ClickDec < 6 &&  RunningMode_ClickDec > 0
                Cursor.State(RunningMode_ClickDec) = Cursor.State(RunningMode_ClickDec) + step;
            elseif RunningMode_ClickDec == 6  % power grasp
                for f = 1:5
                    Cursor.State(f) = Cursor.State(f) + step;
                end
            elseif RunningMode_ClickDec == 7  % pinch grasp
                for f = 1:2
                    Cursor.State(f) = Cursor.State(f) + step;
                end
            elseif RunningMode_ClickDec == 8  % tripod grasp
                for f = 1:3
                    Cursor.State(f) = Cursor.State(f) + step;
                end
            elseif RunningMode_ClickDec == 9  % wrist out
                Cursor.State(6) = Cursor.State(6) - step;
            elseif RunningMode_ClickDec == 10 % wrist in
                Cursor.State(6) = Cursor.State(6) + step;
            elseif RunningMode_ClickDec == 11  % wrist flex
                Cursor.State(7) = Cursor.State(7) - step;
            elseif RunningMode_ClickDec == 12  % wrist extend
                Cursor.State(7) = Cursor.State(7) + step;
            elseif RunningMode_ClickDec == 13 % wrist pronate
                Cursor.State(8) = Cursor.State(8) - step;
            elseif RunningMode_ClickDec == 14 % wrist supinate
                Cursor.State(8) = Cursor.State(8) + step;
            end

            % Boundaries
%             Cursor.State = max(Cursor.State, -1*ones(8,1));
%             Cursor.State = min(Cursor.State, [0,0,0,0,0,1,1,1]');
                
            %%%%% UPDATE CURSOR STATE OR POSITION BASED ON DECODED
            %%%%% DIRECTION
            [ta,tb,tc] = doubleToUDP(Cursor.State(1)*80);
            [va,vb,vc] = doubleToUDP(Cursor.State(2)*80);
            [wa,wb,wc] = doubleToUDP(Cursor.State(3)*80); 
            [xa,xb,xc] = doubleToUDP(Cursor.State(4)*80);
            [ya,yb,yc] = doubleToUDP(Cursor.State(5)*80); 
            [za,zb,zc] = doubleToUDP(Cursor.State(6)*80);

            [wfa,wfb,wfc] = doubleToUDP(Cursor.State(7)*80);
            [wra,wrb,wrc] = doubleToUDP(Cursor.State(8)*80);


            fwrite(Params.udp, [14, ta,tb,tc,va,vb,vc, wa,wb,wc, Data.TargetID]);
            fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, Data.TargetID]);
            fwrite(Params.udp, [8, wfa,wfb,wfc,wra,wrb,wrc,0,0,0 Data.TargetID]);
            fwrite(Params.udp, [15, RunningMode_ClickDec,0,0,0,0,0,0,0,0,0,0]);

            Data.CursorState(:,end+1) = Cursor.State;
            Data.IntendedCursorState(:,end+1) = Cursor.IntendedState;
            Data.CursorAssist(1,end+1) = Cursor.Assistance;
                       
            Cursor.TaskState = 3;
            Data.TaskState(1,end+1)=Cursor.TaskState;

        % end if takes too long
        if TotalTime > Params.MaxReachTime,
            done = 1;
            Data.ErrorID = 3;
            Data.ErrorStr = 'ReachTarget';
            Data.SelectedTargetID = 0;
            Data.SelectedTargetPosition = NaN;
            fprintf('ERROR: %s\n',Data.ErrorStr)
        end
        
       
        targetInd = [1:14];
        
        % end if clicks in a target
%         if (Cursor.State(TargetAxis) - Params.goalState(TargetID)) < 0.1
        if TargetProgress(targetInd(Data.TargetID)) > Params.numCorrectDecode 
            done = 1;
            Data.SelectedTargetID = TargetID;
%             Data.SelectedTargetPosition = Params.ReachTargetPositions(TargetID,:);
        end
        
        % end if in target for hold time (not using clicker)
        if (InTargetTotalTime>=Params.TargetHoldTime) && (Params.ClickerBins==-1),
            done = 1;
            Data.SelectedTargetID = TargetID;
%             Data.SelectedTargetPosition = Params.ReachTargetPositions(TargetID,:); 
        end
        
    end % Reach Target Loop
    
    if step >= length(Params.angles)
        done = 1;
    end
end % only complete if no errors
end
%% Inter trial interval
% blank screen at end of trial but continue collecting data

if Params.InterTrialInterval>0,
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Inter Trial Interval';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
    
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;

    fwrite(Params.udp, [0,1,0])
    while ~done,
        % Update Time & Position
        tim = GetSecs;
        
        % for pausing and quitting expt
        if CheckPause, [Neuro,Data,Params] = ExperimentPause(Params,Neuro,Data); end
        
        % Update Screen
        if (tim-Cursor.LastPredictTime) > 1/Params.ScreenRefreshRate,
            % time
            dt = tim - Cursor.LastPredictTime;
            TotalTime = TotalTime + dt;
            dt_vec(end+1) = dt;
            Cursor.LastPredictTime = tim;
            Data.Time(1,end+1) = tim;
            
            % grab and process neural data
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate),
                dT = tim-Cursor.LastUpdateTime;
                dT_vec(end+1) = dT;
                Cursor.LastUpdateTime = tim;
                
                Data.NeuralTime(1,end+1) = tim;
                [Neuro,Data] = NeuroPipeline(Neuro,Data,Params);   
            end
           
            Data.CursorState(:,end+1) = Cursor.State;
            Data.IntendedCursorState(:,end+1) = Cursor.IntendedState;
            Data.CursorAssist(1,end+1) = Cursor.Assistance;
            
            Cursor.TaskState = 4;
            Data.TaskState(1,end+1)=Cursor.TaskState;
            Data.StopState(1,end+1)=0;
            
            
            % start counting time            
            InTargetTotalTime = InTargetTotalTime + dt;
           
        end
        
        % end if in start target for hold time
        if InTargetTotalTime > Params.InterTrialInterval,
            done = 1;
        end
    end % Instructed Delay Loop
end % only complete if no errors


%% Completed Trial - Give Feedback

% output update times
if Params.Verbose,
    fprintf('      Screen Update: Goal=%iHz, Actual=%.2fHz (+/-%.2fHz)\n',...
        Params.ScreenRefreshRate,mean(1./dt_vec),std(1./dt_vec))
    fprintf('      System Update: Goal=%iHz, Actual=%.2fHz (+/-%.2fHz)\n',...
        Params.UpdateRate,mean(1./dT_vec),std(1./dT_vec))
end

% output feedback
if Data.ErrorID==0,
    fprintf('SUCCESS\n')
    if Params.FeedbackSound,
        sound(Params.RewardSound,Params.RewardSoundFs)
    end
else
    % reset cursor
    Cursor.ClickState = 0;
    Cursor.IntendedState = [0,0,0,0,1]';
    
    fprintf('ERROR: %s\n', Data.ErrorStr)
    
    if Params.FeedbackSound
        sound(Params.ErrorSound,Params.ErrorSoundFs)
    end
    WaitSecs(Params.ErrorWaitTime);
end

end % RunTrial



