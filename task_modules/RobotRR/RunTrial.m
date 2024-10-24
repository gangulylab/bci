function [Data, Neuro, KF, Params, Clicker] = RunTrial(Data,Params,Neuro,TaskFlag,KF,Clicker)
% Runs a trial, saves useful data along the way
% Each trial contains the following pieces
% 1) Get the cursor to the reach target (different on each trial)
% 2) Feedback

global Cursor

set(gcf,'CurrentCharacter','0')

val = double(get(gcf,'CurrentCharacter'));
pause(0.01)

%% Set up trial
ReachTargetPos = Data.TargetPosition;
TargetID = 0; % Target that cursor is in, 0 for no targets

% Output to Command Line
fprintf('\nTrial: %i\n',Data.Trial)
fprintf('  Target: %i\n',Data.TargetPosition)
if Params.Verbose,
    if TaskFlag==2,
        fprintf('    Cursor Assistance: %.2f%%\n',100*Cursor.Assistance)
        if Params.CLDA.Type==3,
            fprintf('    Lambda 1/2 life: %.2fsecs\n',KF.Lambda)
        end
    end
end

% keep track of update times
dt_vec = [];
dT_vec = [];

% grab blackrock data and run through processing pipeline
if Params.BLACKROCK,
    Cursor.LastPredictTime = GetSecs;
    Cursor.LastUpdateTime = Cursor.LastPredictTime;
    Neuro = NeuroPipeline(Neuro,[],Params);
end


Cursor.State = [0,0,0,0,0,0]';
Cursor.State = [Params.StartPos, 0,0,0]';

Cursor.ClickState = 0;
Cursor.ClickDistance = 0;
inTargetOld = 0;
targetNum = 1;

%% Instructed Delay
if ~Data.ErrorID && Params.InstructedDelayTime>0,
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Instructed Delay';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
    

    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;
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
            
            Cursor.TaskState = 1;
            Data.TaskState(1,end+1)=Cursor.TaskState;
           
            Data.StopState(1,end+1)=0;
            
            % start counting time            
            InTargetTotalTime = InTargetTotalTime + dt;
        end
        
        % end if in start target for hold time
        if InTargetTotalTime > Params.InstructedDelayTime,
            done = 1;
        end
    end % Instructed Delay Loop
end % only complete if no errors


%% Cue time
if ~Data.ErrorID && Params.CueTime>0,
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Cue';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end

    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;

    [xa,xb,xc] = doubleToUDP(ReachTargetPos(1));
    [ya,yb,yc] = doubleToUDP(ReachTargetPos(2)); 
    [za,zb,zc] = doubleToUDP(ReachTargetPos(3)) ;

    fwrite(Params.udp, [11, xa,xb,xc,ya,yb,yc,za,zb,zc, 0]);
    
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
            
            Cursor.TaskState = 2;
            Data.TaskState(1,end+1)=Cursor.TaskState;  
            
            Data.StopState(1,end+1)=0;
            
            
            % start counting time            
            InTargetTotalTime = InTargetTotalTime + dt;
           
        end
        
        % end if in start target for hold time
        if InTargetTotalTime > Params.CueTime,
            done = 1;
        end
    end % Instructed Delay Loop
end % only complete if no errors

%% Go to reach target
if ~Data.ErrorID,
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
            
            Cursor.Center = Params.Center;
    
            TargetID = InTargetRobot3D(Cursor,Params.ReachTargetPositions,Params.RobotTargetRadius, Params.RobotTargetDim, Data.TargetID);
      

                Params.TargetID =  Data.TargetID;
                [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params,Neuro,Clicker);
                
                Cursor.ClickState = Click_Decision;
                Cursor.ClickDistance = Click_Distance;
                Data.ClickerState(1,end+1) = Cursor.ClickState;
                Data.ClickerDistance(1,end+1) = Cursor.ClickDistance;
                                
                ClickDec_Buffer(1:end-1) = ClickDec_Buffer(2:end);
                ClickDec_Buffer(end) = Click_Decision;
                RunningMode_ClickDec = RunningMode(ClickDec_Buffer);
                
                ClickToSend = RunningMode_ClickDec;

                Data.FilteredClickerState(1,end+1) = RunningMode_ClickDec;
                
                if RunningMode_ClickDec == 7
                    cs = 'STOP';
                    Cursor.State(4:6) = [0;0;0];
                    
                else 

                A = Params.dA;
                B = Params.dB;
                
                U = zeros(3,1);
                U(1) = int8(RunningMode_ClickDec == 1) - int8(RunningMode_ClickDec == 3);
                U(2) = int8(RunningMode_ClickDec == 2) - int8(RunningMode_ClickDec == 4);
                U(3) = int8(RunningMode_ClickDec == 5) - int8(RunningMode_ClickDec== 6);
                
                vTarget = (ReachTargetPos' - Cursor.State(1:3));
                norm_vTarget = vTarget/norm(vTarget);
                
                AssistVel = Params.AssistAlpha*B*norm_vTarget;
                Data.AssistVel(:,end+1) = AssistVel;
                
                Cursor.State = A*Cursor.State + (1-Params.AssistAlpha)*B*U + AssistVel;
                Cursor.IntendedState = [0 0 0 0 0]';  
                
                % Stop robot at boundaries
                
                if Cursor.State(1) <= Params.limit(1,1)
                   Cursor.State(1) = Params.limit(1,1) + Params.boundaryDist;
                   if Cursor.State(4) < 0
                        Cursor.State(4) = Params.boundaryVel; 
                   end
                elseif Cursor.State(1) >= Params.limit(1,2)
                   Cursor.State(1) = Params.limit(1,2) - Params.boundaryDist;
                   if Cursor.State(4) > 0
                        Cursor.State(4) = -Params.boundaryVel;
                   end   
                elseif Cursor.State(2) <= Params.limit(2,1)
                   Cursor.State(2) = Params.limit(2,1) + Params.boundaryDist;
                    if Cursor.State(5) < 0
                       Cursor.State(5) =  Params.boundaryVel;
                    end
                elseif Cursor.State(2) >= Params.limit(2,2)
                   Cursor.State(2) = Params.limit(2,2) - Params.boundaryDist;
                   if Cursor.State(5) > 044
                        Cursor.State(5) =  -Params.boundaryVel; 
                   end
                elseif Cursor.State(3) <= Params.limit(3,1)
                   Cursor.State(3) = Params.limit(3,1) + Params.boundaryDist;
                   if Cursor.State(6) < 0
                        Cursor.State(6) =  Params.boundaryVel;
                   end
                elseif Cursor.State(3) >= Params.limit(3,2)
                   Cursor.State(3) = Params.limit(3,2) - Params.boundaryDist;
                   if Cursor.State(6) > 0
                        Cursor.State(6) = -Params.boundaryVel;
                   end
                end                               
                end
                
            %%%%% UPDATE CURSOR STATE OR POSITION BASED ON DECODED
            %%%%% DIRECTION

            [xa,xb,xc] = doubleToUDP(Cursor.State(1));
            [ya,yb,yc] = doubleToUDP(Cursor.State(2)); 
            [za,zb,zc] = doubleToUDP(Cursor.State(3)) ;
            
            fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, ClickToSend]);

            Data.CursorState(:,end+1) = Cursor.State;
            Data.IntendedCursorState(:,end+1) = Cursor.IntendedState;
            Data.CursorAssist(1,end+1) = Cursor.Assistance;
                       
            Cursor.TaskState = 3;
            Data.TaskState(1,end+1)=Cursor.TaskState;
  
            % start counting time if cursor is in target
            
             if TargetID==Data.TargetID,
                if inTargetOld == 0
                    if targetNum == 1
                        fwrite(Params.udp, [0, 15, 1])
                    else
                        fwrite(Params.udp, [0, 14, 1])
                    end
                    inTargetOld = 1;
                end
                    InTargetTotalTime = InTargetTotalTime + dt;
                if Params.RobotClicker
                    if RunningMode_ClickDec == 7
                        done = 1;
                        Data.SelectedTargetID = TargetID;
                        Data.SelectedTargetPosition = Params.ReachTargetPositions(TargetID,:); 
%                         fwrite(Params.udp, [0, 6, 0])
                    end
                end
            else
                InTargetTotalTime = 0;
                if inTargetOld
%                     fwrite(Params.udp, [0, 7, 0])
                end
                inTargetOld = 0;
            end
        % end if takes too long
        if TotalTime > Params.MaxReachTime,
            done = 1;
            Data.ErrorID = 3;
            Data.ErrorStr = 'ReachTarget';
            Data.SelectedTargetID = 0;
            Data.SelectedTargetPosition = NaN;
            fprintf('ERROR: %s\n',Data.ErrorStr)
        end
        
        % end if clicks in a target
        if Cursor.ClickState==Params.ClickerBins && TargetID~=0,
            done = 1;
            Data.SelectedTargetID = TargetID;
            Data.SelectedTargetPosition = Params.ReachTargetPositions(TargetID,:);
            if TargetID~=Data.TargetID,
                Data.ErrorID = 4;
                Data.ErrorStr = 'WrongTarget';
            end
        end
        
        % end if in target for hold time (not using clicker)
        if (InTargetTotalTime>=Params.TargetHoldTime) && (Params.ClickerBins==-1),
            
             if       targetNum == 1;
                  targetNum = 2;
                Data.TargetID = 1;
                 fwrite(Params.udp, [0, 15, 2])
                
                 ReachTargetPos = Params.ReachTargetPositions(1,:);

                [xa,xb,xc] = doubleToUDP(ReachTargetPos(1));
                [ya,yb,yc] = doubleToUDP(ReachTargetPos(2)); 
                [za,zb,zc] = doubleToUDP(ReachTargetPos(3)) ;

                fwrite(Params.udp, [12, xa,xb,xc,ya,yb,yc,za,zb,zc, 0]);
             else
                fwrite(Params.udp, [0, 14, 2])
                done = 1;
                Data.SelectedTargetID = TargetID;
                Data.SelectedTargetPosition = Params.ReachTargetPositions(TargetID,:); 
             end
        end
        
    end % Reach Target Loop
end % only complete if no errors

%% Inter trial interval
% blank screen at end of trial but continue collecting data

if Params.InterTrialInterval>0,
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Inter Trial Interval';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
    
%     if TaskFlag==1,
%         OptimalCursorTraj = ...
%             GenerateCursorTraj(StartTargetPos,StartTargetPos,Params.InstructedDelayTime,Params);
%         ct = 1;
%     end
%     
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
            
            % cursor
%             if TaskFlag==1, % imagined movements
%                 Cursor.State(3:4) = (OptimalCursorTraj(ct,:)'-Cursor.State(1:2))/dt;
%                 Cursor.State(1:2) = OptimalCursorTraj(ct,:);
%                 Cursor.Vcommand = Cursor.State(3:4);
%                 ct = ct + 1;
%             end
            
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
else,
    % reset cursor
    Cursor.ClickState = 0;
    % reset cursor
%     Cursor.State = [0,0,0,0,0,0]';
    Cursor.IntendedState = [0,0,0,0,1]';
    
    fprintf('ERROR: %s\n', Data.ErrorStr)
    
    if Params.FeedbackSound,
        sound(Params.ErrorSound,Params.ErrorSoundFs)
    end
    WaitSecs(Params.ErrorWaitTime);
end

end % RunTrial



