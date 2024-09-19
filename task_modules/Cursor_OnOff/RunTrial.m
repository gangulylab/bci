function [Data, Neuro, KF, Params, Clicker] = RunTrial(Data,Params,Neuro,TaskFlag,KF,Clicker)
% Runs a trial, saves useful data along the way
% Each trial contains the following pieces
% 1) Get the cursor to the reach target (different on each trial)
% 2) Feedback

global Cursor

%% Set up trial
ReachTargetPos  = Params.ReachTargets(Data.TargetID,:);
Data.TargetPosition = ReachTargetPos;
targetnum = 1;
TargetID        = Data.TargetID;

% Output to Command Line
fprintf('\nTrial: %i\n',Data.Trial)

% keep track of update times
dt_vec = [];
dT_vec = [];

% grab blackrock data and run through processing 
if Params.BLACKROCK
    Cursor.LastPredictTime = GetSecs;
    Cursor.LastUpdateTime = Cursor.LastPredictTime;
    Neuro = NeuroPipeline(Neuro,[],Params);
end

Cursor.State = [0,0,0,0,0,0]';
Cursor.State(1:2) = Params.StartPositions(1,:);

Cursor.ClickState       = 0;
Cursor.ClickDistance    = 0;
InTargetTotalTime       = 0;
inTargetOld             = 0;
oldTarget               = 0;
%% Instructed Delay
if ~Data.ErrorID && Params.InstructedDelayTime>0
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Instructed Delay';
    
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
        done = 0;
        TotalTime = 0;
        InTargetTotalTime = 0;
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
           
            Cursor.TaskState = 1;
            Data.TaskState(1,end+1) = Cursor.TaskState;
            Data.ClickerState(1,end+1)      = 0;
            Data.ClickerDistance(1,end+1)   = 0;
            Data.FilteredClickerState(1,end+1)   = 0;
            Data.InTarget(:,end+1)  = 0;
            Data.Click(end + 1)     = 0; 

            % start counting time            
            InTargetTotalTime = InTargetTotalTime + dt;
        end
         
        if InTargetTotalTime >= Params.InstructedDelayTime - 0.15
         green = [0 255 0];  % Full green color
         black = BlackIndex(Screen('WindowScreenNumber', Params.WPTR)); % Get black for this window
         % [screenXpixels, screenYpixels] = Screen('WindowSize', window);
         [xCenter, yCenter] = RectCenter(Params.ScreenRectangle);
         baseRect = [0 0 2000 1000];
         centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);
         Screen('FillRect', Params.WPTR, green, centeredRect);
         Screen('Flip', Params.WPTR);
         % WaitSecs(0.05);  % Wait for 50 milliseconds
        end
        % end if in start target for hold time
        if InTargetTotalTime > Params.InstructedDelayTime
            done = 1;
        end
    end % Instructed Delay Loop
     % Return to black screen
         Screen('FillRect', Params.WPTR, black);
         Screen('Flip', Params.WPTR);
   
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

    % draw all targets in gray    
    for k = 1:Params.NumTargets
         TargetPos = Params.ReachTargets(k,:);
         Rect = Params.TargetRect; % centered at (0,0)
         Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
         Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
         Screen('FillOval', Params.WPTR, [200,200,200],  Rect)
    end

   % draw goal target in red
    TargetPos = ReachTargetPos;
    Rect = Params.TargetRect; % centered at (0,0)
    Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
    Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
    Screen('FillOval', Params.WPTR, [255,0,255],  Rect);

    % show drawing
    Screen('DrawingFinished', Params.WPTR);
    Screen('Flip', Params.WPTR);
    
    while ~done
        % Update Time & Position
        tim = GetSecs;
        
        % for pausing and quitting expt
        if CheckPause, [Neuro,Data,Params] = ExperimentPause(Params,Neuro,Data); end
        
        % Update Screen
        if (tim-Cursor.LastPredictTime) > 1/Params.ScreenRefreshRate
            % timet
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
            Cursor.TaskState = 2;
            Data.TaskState(1,end+1)=Cursor.TaskState;
            Data.ClickerState(1,end+1)      = 0;
            Data.ClickerDistance(1,end+1)   = 0;
            Data.FilteredClickerState(1,end+1)   = 0;
            Data.InTarget(:,end+1)  = 0;
            Data.Click(end + 1)     = 0; 

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

    % Initialize 
    done                = 0;
    TotalTime           = 0;
    InTargetTotalTime   = 0;
    lastonofftime       = 0;
    OnState             = 0;
    ClickDec_Buffer     = zeros(Params.RunningModeBinNum, 1);
    StopClicker_Buffer  = zeros(Params.ClickerBinNum, 1);

    temp_dir            = [0,0,0];
    ClickToSend         = 0;

    red                 = [255,0,0];
    green               = [0,255,0];    
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
              
                Cursor.Center = Params.Center;
            
                % check if in any targets
                inTarget = InTargetMulti2D(Cursor,Params.ReachTargets,Params.TargetSize);

                % Update decoder
                [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params,Neuro,Clicker);
                % update the on/off state
                if OnState
                    if (tim-lastonofftime) >= Params.Onwindow
                        OnState = ~OnState;
                        lastonofftime = tim;
                    end
                else
                    if (tim-lastonofftime) >= Params.Offwindow
                        OnState = ~OnState;
                        lastonofftime = tim;
                    end
                end

                % Set Output to 0 if off state
                if ~OnState
                    Click_Decision = 0;
                end

                % Click_Decision
                Cursor.ClickState               = Click_Decision;
                Cursor.ClickDistance            = Click_Distance;
                Data.ClickerState(1,end+1)      = Cursor.ClickState;
                Data.ClickerDistance(1,end+1)   = Cursor.ClickDistance;
                Data.OnState(1,end+1)           = OnState;

                ClickDec_Buffer(1:end-1)        = ClickDec_Buffer(2:end);
                ClickDec_Buffer(end)            = Click_Decision;
                RunningMode_ClickDec            = RunningMode(ClickDec_Buffer);              

                StopClicker_Buffer(1:end-1)     = StopClicker_Buffer(2:end);
                StopClicker_Buffer(end)         = RunningMode_ClickDec == Params.ClickAction;

                ClickToSend                     = RunningMode_ClickDec;
                Data.FilteredClickerState(1,end+1) = RunningMode_ClickDec;

                U    = zeros(3,1);
                U(1) = int8(RunningMode_ClickDec == 1) - int8(RunningMode_ClickDec == 3);
                U(2) = int8(RunningMode_ClickDec == 6) - int8(RunningMode_ClickDec == 5);
                U(3) = 0;

                A = Params.dA;
                B = Params.dB;
                Cursor.State = A*Cursor.State + B*U;
        
                if any(find(inTarget))
                    Data.InTarget(:,end+1)  = find(inTarget);
                else
                    Data.InTarget(:,end+1)  = 0;
                end

                % draw all targets in gray    
                for k = 1:Params.NumTargets
                     TargetPos = Params.ReachTargets(k,:);
                     Rect = Params.TargetRect; % centered at (0,0)
                     Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
                     Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
                     Screen('FillOval', Params.WPTR, [200,200,200],  Rect)
                end

                % draw goal target in Purple
                TargetPos = ReachTargetPos;
                Rect = Params.TargetRect; % centered at (0,0)
                Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
                Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
                Screen('FillOval', Params.WPTR, [255,0,255],  Rect);

                % draw cursor
                CursorRect                  = Params.CursorRect;
                CursorRect([1,3])           = CursorRect([1,3]) + Cursor.State(1) + Params.Center(1); % add x-pos
                CursorRect([2,4])           = CursorRect([2,4]) + Cursor.State(2) + Params.Center(2); % add y-pos
                if OnState
                    Screen('FillOval', Params.WPTR, green, CursorRect);
                else
                    Screen('FillOval', Params.WPTR,red, CursorRect);
                end
                Screen('DrawingFinished', Params.WPTR);
                Screen('Flip', Params.WPTR);
        
                Data.CursorState(:,end+1) = Cursor.State;     
                Cursor.TaskState = 3;
                Data.TaskState(1,end+1) = Cursor.TaskState;
                click = 0;
        
        
                % trial end conditions 

                if find(inTarget)==TargetID  % if in any target
                    if Params.UseClicker % use clicker
                        if (mean(StopClicker_Buffer) > Params.ClickerBinThresh) % if use clicker 
                            click = 1;
                            % done = 1;
                            targetnum = targetnum+1;
                            TargetID = Params.TargetOrder(targetnum);
                            ReachTargetPos = Params.ReachTargets(TargetID,:);

                        end
                    else % use hold time
                        if find(inTarget) == oldTarget
                            InTargetTotalTime = InTargetTotalTime + dt;
                        else
                            oldTarget = find(inTarget);
                            InTargetTotalTime = 0;
                        end

                        if InTargetTotalTime > Params.HoldTime
                            % done = 1;
                            targetnum = targetnum+1;
                            TargetID = Params.TargetOrder(targetnum);
                            ReachTargetPos = Params.ReachTargets(TargetID,:);

                        end
                    end
                end
             

                Data.Click(end + 1) = click; 
           
                % end if takes too long
                if TotalTime > Params.MaxReachTime
                    done = 1;
                    % Data.ErrorID = 3;
                    % Data.ErrorStr = 'Time-out';
                    Data.SelectedTargetID = 0;
                    Data.SelectedTargetPosition = NaN;

                end
        
                % end if clicks in a target or hold
                % if done == 1
                %     Data.SelectedTargetID = find(inTarget);
                %     Data.SelectedTargetPosition = Params.ReachTargets(Data.SelectedTargetID,:);
                % 
                %     if any(Data.SelectedTargetID ~= Data.TargetID)
                %         Data.ErrorID = 4;
                %         Data.ErrorStr = 'WrongTarget';
                % 
                %     end
                % end
        end
    end % Reach Target Loop
end % only complete if no errors

%% Inter trial interval
% blank screen at end of trial but continue collecting data
taskStage = 0;

Screen('DrawingFinished', Params.WPTR);
Screen('Flip', Params.WPTR);

if Params.InterTrialInterval>0
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Inter Trial Interval';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
  
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;

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

            Cursor.TaskState = 4;
            Data.TaskState(1,end+1)=Cursor.TaskState;
                     
            Data.ClickerState(1,end+1)      = 0;
            Data.ClickerDistance(1,end+1)   = 0;
            Data.FilteredClickerState(1,end+1)   = 0;
            Data.InTarget(:,end+1)  = 0;
            Data.Click(end + 1)     = 0; 

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

% output update times
if Params.Verbose
    fprintf('      Screen Update: Goal=%iHz, Actual=%.2fHz (+/-%.2fHz)\n',...
        Params.ScreenRefreshRate,mean(1./dt_vec),std(1./dt_vec))
    fprintf('      System Update: Goal=%iHz, Actual=%.2fHz (+/-%.2fHz)\n',...
        Params.UpdateRate,mean(1./dT_vec),std(1./dT_vec))
end

% output feedback
if Data.ErrorID==0
    fprintf('SUCCESS\n')
    if Params.FeedbackSound
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

end % RunTrial



