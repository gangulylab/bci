function [Data, Neuro, KF, Params, Clicker] = RunTrial(Data,Params,Neuro,TaskFlag,KF,Clicker)
% Runs a trial, saves useful data along the way
% Each trial contains the following pieces
% 1) Get the cursor to the reach target (different on each trial)
% 2) Feedback

global Cursor

%% Set up trial
write(Params.udp, [0,40,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); % open file  %set screen to black
% if Params.OperationModeReset == 0

if Data.TargetID == 8
    write(Params.udp, [0,12,3.1415*10,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
elseif Data.TargetID == 9
    write(Params.udp, [0,12,3.1415/2*10,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
end

% Params.StartPos
% [xa,xb,xc] = doubleToUDP(Params.StartPos(1));
% [ya,yb,yc] = doubleToUDP(Params.StartPos(2)); 
% [za,zb,zc] = doubleToUDP(Params.StartPos(3) - 256) ;
% write(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, 0], "127.0.0.1", Params.pythonPort) ; % send pos

% write(Params.udp, [0,2,Params.RobotMode,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
% write(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);                  % reset robot

if Data.TargetID == 3 && Params.OperationModeReset == 1
    l = "LEFT"
    write(Params.udp, [0,27,1,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
end



write(Params.udp, [0,5,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); % open file     
ReachTargetPos = Data.TargetPosition;
TargetID = 0; % Target that cursor is in, 0 for no targets

% Output to Command Line
fprintf('\nTrial: %i\n',Data.Trial)

TargetName = {'Right Thumb', 'Left Leg', 'Left Thumb', 'Right Wrist', 'Lips', 'Tongue', 'Both Middle Fingers', 'Right Wrist', 'Left Wrist'};

fprintf('TARGET: %s\n',TargetName{Data.TargetID})
fprintf('Press any key if robot at start position')

pause()  % pause to make sure robot is at start position


% keep track of update times
dt_vec = [];
dT_vec = [];

% grab blackrock data and run through processing pipeline
Cursor.LastPredictTime = GetSecs;
Cursor.LastUpdateTime = Cursor.LastPredictTime;
if Params.BLACKROCK
    Neuro = NeuroPipeline(Neuro,[],Params);
end

% reset cursor
if Params.LongTrial
        Cursor.State = [0,0,0,0,0,0]';
        Cursor.State(1:3) = Params.LongStartPos(Data.TargetID,:);
else
        Cursor.State = [0,0,0,0,0,0]';
        Cursor.State(1:3) = Params.StartPos;
end
    Cursor.IntendedState = [0,0,0,0,1]';

Cursor.ClickState = 0;
Cursor.ClickDistance = 0;
inTargetOld = 0;


%% Instructed Delay
if ~Data.ErrorID && Params.InstructedDelayTime>0
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Instructed Delay';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
    
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;
    
    if Params.ShowFlash
        % send command for green screen
        write(Params.udp, [0,40,2,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); % open file   
    end     
    
    while ~done
        % Update Time & Position
        tim = GetSecs;
        
        % for pausing and quitting expt
        if CheckPause, [Neuro,Data,Params] = ExperimentPause(Params,Neuro,Data); end
        
        if (tim - tstart) > Params.FlashDuration/1000 
            % send
            write(Params.udp, [0,40,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); % open file  
        end
        
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
    
    write(Params.udp, [0, 29, Data.TargetID,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 

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
    write(Params.udp, [0,40,1,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); % open file  
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
            Params.index = Params.index+1;
            [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params,Neuro,Clicker);
                
            if TaskFlag==1 % imagined movements
%                 if TargetID == Data.TargetID
%                     Click_Decision = 0;
%                     Cursor.State(4:6) = [0;0;0];
%                     Data.StopState(1,end+1)=1;
%                 else
                Click_Decision = Params.TargetID;
                Data.StopState(1,end+1)=0;
%                 end 
            end              
                
                Cursor.ClickState = Click_Decision;
                Cursor.ClickDistance = Click_Distance;
                Data.ClickerState(1,end+1) = Cursor.ClickState;
                Data.ClickerDistance(1,end+1) = Cursor.ClickDistance;
                                
                ClickDec_Buffer(1:end-1) = ClickDec_Buffer(2:end);
                ClickDec_Buffer(end) = Click_Decision;
                RunningMode_ClickDec = RunningMode(ClickDec_Buffer);  
                
                RunningMode_ClickDec = any(Params.ValidDir == RunningMode_ClickDec)*RunningMode_ClickDec; % Filter by allowable directions

                ClickToSend = RunningMode_ClickDec;   
%                 ClickToSend
                
                % Clamp to correct action
                if Params.ClampCorrect
                    if ClickToSend ~= Data.TargetID
                        ClickToSend = 0;
                    end
                end
                
                % Beta stopping
                beta_scalar = betaband_output(Params,Neuro); 
                Data.BetaScalar(1,end+1)    = beta_scalar;
                if Params.UseBetaStop
                     if (beta_scalar >= Params.BetaThreshold)
                         ClickToSend = 0;
                     end
                    Data.BetaClickerState(1,end+1) = ClickToSend;
                end
                
                fprintf('Decode: %i \n',RunningMode_ClickDec)

                Data.FilteredClickerState(1,end+1) = RunningMode_ClickDec;

                A = Params.dA;
                B = Params.dB;
                
                U = zeros(3,1);
                U(1) = int8(RunningMode_ClickDec == 1) - int8(RunningMode_ClickDec == 3);
                U(2) = int8(RunningMode_ClickDec == 2) - int8(RunningMode_ClickDec == 4);
                U(3) = int8(RunningMode_ClickDec == 5) - int8(RunningMode_ClickDec == 6);

                vTarget = (Data.TargetPosition'- Cursor.State(1:3));

                if norm(vTarget) ~= 0
                    norm_vTarget = vTarget/norm(vTarget);
                else
                    norm_vTarget = [0;0;0];
                end
                
                AssistVel = Params.AssistAlpha*B*norm_vTarget;
                Data.AssistVel(:,end+1) = AssistVel;
                
                Cursor.State = A*Cursor.State + (1-Params.AssistAlpha)*B*U + AssistVel;
                Cursor.IntendedState = [0 0 0 0 0]';  
  
            %%%%% UPDATE CURSOR STATE OR POSITION BASED ON DECODED
            %%%%% DIRECTION
%             
%             [xa,xb,xc] = doubleToUDP(Cursor.State(1));
%             [ya,yb,yc] = doubleToUDP(Cursor.State(2)); 
%             [za,zb,zc] = doubleToUDP(Cursor.State(3)-256) ;
%             
%             write(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, ClickToSend], "127.0.0.1", Params.pythonPort); ; % send pos
% 
%             [xa,xb,xc] = doubleToUDP(Cursor.State(4));
%             [ya,yb,yc] = doubleToUDP(Cursor.State(5)); 
%             [za,zb,zc] = doubleToUDP(Cursor.State(6)) ;
            
            write(Params.udp, [5, 0,0,0,0,0,0,0,0,0, ClickToSend], "127.0.0.1", Params.pythonPort); % send vel

            Data.CursorState(:,end+1) = Cursor.State;
            Data.IntendedCursorState(:,end+1) = Cursor.IntendedState;
            Data.CursorAssist(1,end+1) = Cursor.Assistance;
                       
            Cursor.TaskState = 3;
            Data.TaskState(1,end+1)=Cursor.TaskState;
             
        % end if takes too long
        if TotalTime > Params.MaxReachTime
            done = 1;
            Data.ErrorID = 3;
            Data.ErrorStr = 'ReachTarget';
            Data.SelectedTargetID = 0;
            Data.SelectedTargetPosition = NaN;
            fprintf('ERROR: %s\n',Data.ErrorStr)
        end
        
        if Params.udp.NumBytesAvailable > 0
            u = read(Params.udp,Params.udp.NumBytesAvailable);
            u = u(end);
            
            if u == 1
                fprintf("In Target\n")
            elseif u == 2
                done = 1;
                fprintf("SUCCESS\n")
            end
            u = 0;
        end
        
        end 
    end % Reach Target Loop
end % only complete if no errors

%% Inter trial interval
% blank screen at end of trial but continue collecting data

if Params.InterTrialInterval>0
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Inter Trial Interval';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
  
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;

    write(Params.udp, [0,1,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 

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
            
            Cursor.TaskState = 4;
            Data.TaskState(1,end+1)=Cursor.TaskState;
            Data.StopState(1,end+1)=0;
                     
            % start counting time            
            InTargetTotalTime = InTargetTotalTime + dt;     
        end
        
        % end if in start target for hold time
        if InTargetTotalTime > Params.InterTrialInterval
            done = 1;
        end
    end % Instructed Delay Loop
end % only complete if no errors


%% Completed Trial - Give Feedback

% output feedback
if Data.ErrorID==0
%     fprintf('SUCCESS\n')
    if Params.FeedbackSound,
        sound(Params.RewardSound,Params.RewardSoundFs)
    end
else
    % reset cursor
    Cursor.ClickState = 0;
    % reset cursor
%     Cursor.State = [0,0,0,0,0,0]';
    Cursor.IntendedState = [0,0,0,0,1]';
    
    fprintf('ERROR: %s\n', Data.ErrorStr)
    
    if Params.FeedbackSound
        sound(Params.ErrorSound,Params.ErrorSoundFs)
    end
    WaitSecs(Params.ErrorWaitTime);
end

% pause()

end % RunTrial



