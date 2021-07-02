function [Data, Neuro, KF, Params, Clicker] = RunTrial(Data,Params,Neuro,TaskFlag,KF,Clicker)
% Runs a trial, saves useful data along the way
% Each trial contains the following pieces
% 1) Get the cursor to the reach target (different on each trial)
% 2) Feedback

global Cursor

set(gcf,'CurrentCharacter','0')

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

% reset cursor
% if Params.LongTrial
%         Cursor.State = [0,0,0,0,0,0]';
%         Cursor.State(1:3) = Params.LongStartPos(Data.TargetID,:);
% else
%         Cursor.State = [0,0,0,0,0,0]';
%         
% end
Cursor.State = [0,0,0,0,0,0]';
Cursor.IntendedState    = [0,0,0,0,1]';

Cursor.RState = [3.14/2;0]';


Cursor.ClickState = 0;
Cursor.ClickDistance = 0;
inTargetOld = 0;

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
    
    % Send target position
    [xa,xb,xc] = doubleToUDP(ReachTargetPos(1));
    [ya,yb,yc] = doubleToUDP(ReachTargetPos(2)); 
    [za,zb,zc] = doubleToUDP(ReachTargetPos(3)) ;

    fwrite(Params.udp, [11, xa,xb,xc,ya,yb,yc,za,zb,zc, 0]);  
    
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
if ~Data.ErrorID
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Reach Target';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end

    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;
    
    ClickDec_Buffer = zeros(Params.RunningModeBinNum, 1);
    StopClicker_Buffer = zeros(Params.ClickerBinNum, 1);
    Grasp_Buffer = zeros(Params.GraspBinNum, 1);
    
    
    graspnow = 0;
    switchLocked = 0;
    temp_dir = [0,0,0];
    ClickToSend = 0;
    
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
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate),
                dT = tim-Cursor.LastUpdateTime;
                dT_vec(end+1) = dT;
                Cursor.LastUpdateTime = tim;
                
                Data.NeuralTime(1,end+1) = tim;
                [Neuro,Data] = NeuroPipeline(Neuro,Data,Params);              
                
            end
            
            Cursor.Center = Params.Center;
    
            TargetID = InTargetRobot3D(Cursor,Params.ReachTargetPositions,Params.RobotTargetRadius, Params.RobotTargetDim, Data.TargetID);
            
            dist2D = norm([[ReachTargetPos(1); ReachTargetPos(2)] - Cursor.State(1:2)])
            
            if dist2D < Params.RobotTargetRadius
                inTarget = 1;
            end

            Params.TargetID =  Data.TargetID;
            Params.index = Params.index + 1;
            
            [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params,Neuro,Clicker);
            Cursor.ClickState = Click_Decision;
            Cursor.ClickDistance = Click_Distance;
            Data.ClickerState(1,end+1) = Cursor.ClickState;
            Data.ClickerDistance(1,end+1) = Cursor.ClickDistance;
                                
            ClickDec_Buffer(1:end-1) = ClickDec_Buffer(2:end);
            ClickDec_Buffer(end) = Click_Decision;
            RunningMode_ClickDec = RunningMode(ClickDec_Buffer);  
            Data.FilteredClickerState(1,end+1) = RunningMode_ClickDec;
            ClickToSend = RunningMode_ClickDec;
            
            % Open/close 
            if Params.opMode == 1
                Grasp_Buffer(1:end-1) = Grasp_Buffer(2:end);
                Grasp_Buffer(end) = Click_Decision;

                openCnt = sum(Grasp_Buffer == 3);            
                closeCnt = sum(Grasp_Buffer == 1);
                
                if openCnt/Params.GraspBinNum > Params.GraspBinThresh
                    graspnow = 1;
                elseif closeCnt/Params.GraspBinNum > Params.GraspBinThresh
                    graspnow = 2;
                else
                    graspnow = 0;
                end   
            end
            
           
            
            % Mode switching
            
            if switchLocked == 0
                StopClicker_Buffer(1:end-1) = StopClicker_Buffer(2:end);
                StopClicker_Buffer(end) = Click_Decision == 7;
            else
                switchLockCnt = switchLockCnt + 1;
                if switchLockCnt > Params.switchLockCnt
                    switchLocked = 0;
                end
            end
            
            
            if mean(StopClicker_Buffer) > Params.ClickerBinThresh
               if Params.opMode == 0
                   Params.opMode = 1;
                   switchLocked = 1;
                   switchLockCnt = 0;
                   StopClicker_Buffer = zeros(Params.ClickerBinNum, 1);
                   fprintf('Mode: %s\n','Grasp')
                   if inTarget && Params.autoCenterOverTarget
                       Cursor.State(1) = ReachTargetPos(1);
                       Cursor.STate(2) = ReachTargetPos(2);
                       [xa,xb,xc] = doubleToUDP(Cursor.State(1));
                        [ya,yb,yc] = doubleToUDP(Cursor.State(2)); 
                        [za,zb,zc] = doubleToUDP(Cursor.State(3)) ;

                        fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, 0]);
                   end
                   
               elseif Params.opMode == 1
                   Params.opMode = 0;
                   switchLocked = 1;
                   switchLockCnt = 0;
                   StopClicker_Buffer = zeros(Params.ClickerBinNum, 1);
                    fprintf('Mode: %s\n','Translation')

                   
               end
            end
                        
              
            % Click to user input
            
            A = Params.dA;
            B = Params.dB;
            
            U = zeros(3,1);
            U(1) = int8(RunningMode_ClickDec == 1) - int8(RunningMode_ClickDec == 3);
            U(2) = int8(RunningMode_ClickDec == 2) - int8(RunningMode_ClickDec == 4);
            U(3) = int8(RunningMode_ClickDec == 5) - int8(RunningMode_ClickDec== 6);

            
            % Mode 0: Translation
            if Params.opMode == 0

                vTarget = (Data.TargetPosition'- Cursor.State(1:3));
                norm_vTarget = vTarget/norm(vTarget);

                AssistVel = Params.AssistAlpha*B*norm_vTarget;
                Data.AssistVel(:,end+1) = AssistVel;

                Cursor.State = A*Cursor.State + (1-Params.AssistAlpha)*B*U + AssistVel;
                Cursor.IntendedState = [0 0 0 0 0]';  

            % Mode 1: Grasp
            elseif Params.opMode == 1
               
                % Zero horizontal motion
                Cursor.State(4) = 0;
                Cursor.State(5) = 0;

                % Update vertical
                Cursor.State(3) = Cursor.State(3) + A(3,6)*Cursor.State(6);
                Cursor.State(6) = A(6,6)*Cursor.State(6) + B(6,3)*U(3);

                % Update Rotation
                Cursor.RState(1)  = Cursor.RState(1) + A(3,6)*Cursor.RState(2);
                Cursor.RState(2) =  Params.ra*Cursor.RState(2)  + Params.rb*U(2);
                 
            end
            

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
               if Cursor.State(5) > 0
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
                
            %%%%% UPDATE CURSOR STATE OR POSITION BASED ON DECODED
            %%%%% DIRECTION
            
            if Params.opMode == 0
                [xa,xb,xc] = doubleToUDP(Cursor.State(1));
                [ya,yb,yc] = doubleToUDP(Cursor.State(2)); 
                [za,zb,zc] = doubleToUDP(Cursor.State(3)) ;
            
                fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, ClickToSend]);
            elseif Params.opMode == 1
                [xa,xb,xc] = doubleToUDP(Cursor.RState(1)*10);
                [ya,yb,yc] = doubleToUDP(graspnow); 
                [za,zb,zc] = doubleToUDP(Cursor.State(3)) ;
            
                fwrite(Params.udp, [7, xa,xb,xc,ya,yb,yc, za,zb,zc, ClickToSend]);
            end

            Data.CursorState(:,end+1) = Cursor.State;
            Data.IntendedCursorState(:,end+1) = Cursor.IntendedState;
            Data.CursorAssist(1,end+1) = Cursor.Assistance;
                    
            Data.CursorRState(:,end+1) = Cursor.RState;
            Cursor.TaskState = 3;
            Data.TaskState(1,end+1)=Cursor.TaskState;
            
            Data.opMode = Params.opMode;
             
            
        % end if takes too long
        if TotalTime > Params.MaxReachTime,
            done = 1;
            Data.ErrorID = 3;
            Data.ErrorStr = 'ReachTarget';
            Data.SelectedTargetID = 0;
            Data.SelectedTargetPosition = NaN;
            fprintf('ERROR: %s\n',Data.ErrorStr)
        end

        end
        
    end % Reach Target Loop
end % only complete if no errors

pause(0.5)

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






