function [Data, Neuro, KF, Params, Clicker] = RunTrial(Data,Params,Neuro,TaskFlag,KF,Clicker)
% Runs a trial, saves useful data along the way
% Each trial contains the following pieces
% 1) Get the cursor to the reach target (different on each trial)
% 2) Feedback

global Cursor

%% Set up trial
ReachTargetPos  = Data.TargetPosition;
TargetID        = Data.TargetID; % Target that cursor is in, 0 for no targets
TargetNum1      = Params.TargetList(TargetID,1);
TargetNum2      = Params.TargetList(TargetID,2);
StartPosNum = 1;
% StartPosNum     = Params.TargetList(TargetID,3);

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
Cursor.State(1:2) = Params.StartPositions(StartPosNum,:);

Cursor.ClickState       = 0;
Cursor.ClickDistance    = 0;
inTargetOld             = 0;

% init inference
g = Params.ReachTargets1;
% cosine
p_theta_c = @(theta)exp(Params.velB* Params.velk*(theta+1));

% distance
p_theta_d = @(theta)exp(-Params.distB*Params.distk*theta);

% cos_sim = @(a,b)dot(a,b)/(norm(a)*norm(b));
cos_sim = @(a,b)(1 - acos(dot(a,b)/(norm(a)*norm(b))) / pi);
c = 0.05;
P_g_g = ones(length(g), length(g)) * c/(length(g)-1) + eye(length(g))*((1-c) -  c/(length(g)-1) );

b1(1,:) = ones(1, length(g));
b1(1,:) = [b1(1,:)]/sum(b1(1,:) );  % cosine

b2(1,:) = ones(1, length(g));
b2(1,:) = [b2(1,:)]/sum(b2(1,:) );  % euclidean

b3(1,:) = ones(1, length(g));
b3(1,:) = [b3(1,:)]/sum(b3(1,:) );  %combination

theta_c(1,:) = [0,0,0];
theta_d(2,:) = [0,0,0];

assist  = 0;
userVel = [0,0,0,0,0,0]';

binInd = 0;

taskStage = 0;
Params.RobotClicker = 0;
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
            binInd = binInd + 1;
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
            Data.TaskState(1,end+1)= Cursor.TaskState;
            Data.ClickerState(1,end+1)      = 0;
            Data.ClickerDistance(1,end+1)   = 0;
            Data.FilteredClickerState(1,end+1)   = 0;
            Data.UserVel(:,end+1)   = [0;0];
            Data.AssistVel(:,end+1)   = [0;0];
            Data.Belief(:,end+1)   = b3;
            Data.InTarget(:,end+1)  = 0;
            Data.Assist(:,end+1)    = 0;
            Data.Dyn(:,end+1)       = 0;
            Data.Click(end + 1)     = 0; 
            Data.TaskStage(end + 1)          = taskStage;

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

    % display targets
    
    % goal
    
    for k = 1:Params.NumTargets1
         TargetPos = Params.ReachTargets1(k,:);

         Rect = Params.TargetRect; % centered at (0,0)
         Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
         Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos

         Screen('FillOval', Params.WPTR, [200,200,200],  Rect)
         TargetPos = [0, -200];
    end
    
    for k = 1:Params.NumTargets2
         TargetPos = Params.ReachTargets2(k,:);

         Rect = Params.TargetRect; % centered at (0,0)
         Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
         Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos

         Screen('FillOval', Params.WPTR, [200,200,200],  Rect)
         TargetPos = [0, -200];
    end
 
     % draw goal 1
     TargetPos = Params.ReachTargets1(TargetNum1,:);
     Rect = Params.TargetRect; % centered at (0,0)
     Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
     Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
     Screen('FillOval', Params.WPTR, [255,0,0],  Rect)
    
     % draw goal 2
     TargetPos = Params.ReachTargets2(TargetNum2,:);
     Rect = Params.TargetRect; % centered at (0,0)
     Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
     Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
     Screen('FillOval', Params.WPTR, [255,0,0],  Rect)
      
    Screen('DrawingFinished', Params.WPTR);
    Screen('Flip', Params.WPTR);
    
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
           
            binInd = binInd + 1;
            
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
            Data.UserVel(:,end+1)   = [0;0];
            Data.AssistVel(:,end+1)   = [0;0];
            Data.Belief(:,end+1)   = b3;
            Data.InTarget(:,end+1)  = 0;
            Data.Assist(:,end+1)    = 0;
            Data.Dyn(:,end+1)      = 0;
            Data.Click(end + 1)     = 0; 
            Data.TaskStage(end + 1)          = taskStage;
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
    ClickDec_Buffer     = zeros(Params.RunningModeBinNum, 1);
    StopClicker_Buffer  = zeros(Params.ClickerBinNum, 1);
    ChangeClicker_Buffer  = zeros(Params.ChangeBinNum, 1);
    temp_dir            = [0,0,0];
    ClickToSend         = 0;
    maxV = 0;
    outside = 1;
    assist_target = 0;
    last_ctim = GetSecs;
    taskStage = 1;
    
    t_targetStart = tstart;
    
    while ~done
        % Update Time & Position
        tim = GetSecs;
                
        dt_targetStart = tim - t_targetStart;  %time since target stage start
        
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
            binInd = binInd + 1;
            % grab and process neural data
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate)
                dT = tim-Cursor.LastUpdateTime;
                dT_vec(end+1) = dT;
                Cursor.LastUpdateTime = tim;
                
                Data.NeuralTime(1,end+1) = tim;
                [Neuro,Data] = NeuroPipeline(Neuro,Data,Params);              
                
              
            Cursor.Center = Params.Center;
    
            inTarget = InTargetMulti2D(Cursor,g);

            x = [Cursor.State(1), Cursor.State(2)];
            v = [userVel(4), userVel(5)];
            
            maxV = max([maxV, norm(v)]);
            
            if maxV > 5
            
            % inference
            for cnt_g = 1:length(g)
                dist = norm([x - g(cnt_g,:)]);
                theta_d = norm([x - g(cnt_g,:)]);  %Euclidean Distance

                a1_d = p_theta_d(theta_d);
                a2_d = 0;
                for j = 1:length(g)
                    a2_d = a2_d + P_g_g(cnt_g, j)*b2(j);
                end
                b2_tmp(cnt_g) = a1_d*a2_d;

            end
            b2 = b2_tmp/sum(b2_tmp);
            
            
        %   inference - velocity
            v = [userVel(4), userVel(5)];
            if norm(v) < 10
                v = [0,0];
            end
            for cnt_g = 1:length(g)
                target_d = g(cnt_g,:) - x;

                %Velocity
                if norm(v) == 0
                    theta_c = 0;
                else
                    theta_c = cos_sim(target_d, v); % cosine similarity
                end
                a1_c = p_theta_c(theta_c);
                a2_c = 0;
                for j = 1:length(g)
                    a2_c = a2_c + P_g_g(cnt_g, j)*b1(j);
                end
                
                b1_tmp(cnt_g) = a1_c*a2_c;
            end
                b1 = b1_tmp/sum(b1_tmp);
           
        % inference - combination

        for cnt_g = 1:length(g)
        theta_d = norm([x- g(cnt_g,:)]);  %Euclidean Distance
        target_d =  g(cnt_g,:) -x ;

        a1_d = p_theta_d(theta_d);

        a2_d = 0;
        for j = 1:length(g)
            a2_d = a2_d + P_g_g(cnt_g, j)*b3(j);
        end

        %Velocity
        if norm(v) == 0  || theta_d < 25
            theta_c = 0;
        else
        theta_c = cos_sim(target_d, v); % cosine similarity
        end
        a1_c = p_theta_c(theta_c);

        a2_c = 0;
        for j = 1:length(g)
            a2_c = a2_c + P_g_g(cnt_g, j)*b3(j);
        end
        b3_tmp(cnt_g) = a1_d*a2_d*a1_c;
            
        end
        b3 = b3_tmp/sum(b3_tmp);
            end
            
        if Params.ObservationFunc == 1
            b = b1;
        elseif Params.ObservationFunc == 2
            b = b2;
        elseif Params.ObservationFunc == 3
            b = b3;
        end 
        
        belief_sort = sort(b, 'descend');       
        
        if Params.AssistLock
            if ~assist
                if ((belief_sort(1) - belief_sort(2)) > Params.AssistThresh) & (dt_targetStart > Params.AssistDelay)
                    assist = 1;
                    [~,assist_target] = max(b);
                    ChangeClicker_Buffer  = zeros(Params.ChangeBinNum, 1);
                    last_ctim =GetSecs;
                else
                    assist = 0;
                end
            else % assist on
%                 a = ChangeClicker_Buffer;
%                 ctim = GetSecs;
%                if (sum(a == mode(a))> Params.ChangeBinThresh) && mode(a)>0 && (ctim -last_ctim) > 4 % change target
%                    ChangeClicker_Buffer  = zeros(Params.ChangeBinNum, 1);
%                    last_ctim = ctim;
%                        if mode(a) == 1 % shift right
%                            if assist_target < 3
%                                assist_target = assist_target + 1;
%                            end
%                        elseif mode(a) == 3
%                            if assist_target > 1
%                                assist_target = assist_target - 1;
%                            end
%                        end
%                end
            end

        else
            if (belief_sort(1) - belief_sort(2) > Params.AssistThresh) && (dt_targetStart > Params.AssistDelay)
                assist = 1;
                [~,assist_target] = max(b);
            else
                assist = 0;
            end
        end
            
        Params.TargetID =  Data.TargetID;
        
        % Update decoder
        [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params,Neuro,Clicker);
                            
        Cursor.ClickState               = Click_Decision;
        Cursor.ClickDistance            = Click_Distance;
        Data.ClickerState(1,end+1)      = Cursor.ClickState;
        Data.ClickerDistance(1,end+1)   = Cursor.ClickDistance;

        ClickDec_Buffer(1:end-1)        = ClickDec_Buffer(2:end);
        ClickDec_Buffer(end)            = Click_Decision;
        RunningMode_ClickDec            = RunningMode(ClickDec_Buffer);              

        StopClicker_Buffer(1:end-1)     = StopClicker_Buffer(2:end);
        StopClicker_Buffer(end)         = RunningMode_ClickDec == Params.ClickAction;

        ChangeClicker_Buffer(1:end-1)     = ChangeClicker_Buffer(2:end);
        ChangeClicker_Buffer(end)         = RunningMode_ClickDec;

        ClickToSend                     = RunningMode_ClickDec;
        Data.FilteredClickerState(1,end+1) = RunningMode_ClickDec;

        % Cursor Dynamics
        if Params.AssistLock
            if Params.Assist && Params.SlowAtTarget && any(find(inTarget,1)==assist_target) && any(inTarget)
                dyn = 2;
                A = Params.dA2;
                B = Params.dB2; 
            else
                dyn = 1;
                A = Params.dA;
                B = Params.dB;
            end  
        else   
            if Params.Assist && Params.SlowAtTarget && any(inTarget)
                dyn = 2;
                A = Params.dA2;
                B = Params.dB2; 
            else
                dyn = 1;
                A = Params.dA;
                B = Params.dB;
            end
        end

        U    = zeros(3,1);
        U(1) = int8(RunningMode_ClickDec == 1) - int8(RunningMode_ClickDec == 3);
        U(2) = int8(RunningMode_ClickDec == 2) - int8(RunningMode_ClickDec == 4);
        U(3) = 0;
        
        userVel = A*userVel + B*U;
        if assist && Params.Assist
            vTarget         = (g(assist_target,:) - Cursor.State(1:2)');
            norm_vTarget    = vTarget/norm(vTarget);
            AssistVel       = Params.AssistGain*Params.AssistAlpha*B*[norm_vTarget,0]';
            Cursor.State    = A*Cursor.State + (1-Params.AssistAlpha)*B*U + AssistVel;
        else
            AssistVel = zeros(6,1);
             Cursor.State = A*Cursor.State + B*U;
        end
        
        Data.Belief(:,end+1)    = b;
        Data.UserVel(:,end+1)   = userVel([4,5]);
        Data.AssistVel(:,end+1) = AssistVel([4,5]);
        if any(find(inTarget))
            Data.InTarget(:,end+1)  = find(inTarget);
        else
            Data.InTarget(:,end+1)  = 0;
        end
        Data.Assist(:,end+1)    = assist & Params.Assist;
        Data.Dyn(:,end+1)      = dyn;

        % draw targets
        
        if taskStage == 1
            for k = 1:Params.NumTargets1
             TargetPos = Params.ReachTargets1(k,:);
             Rect           = Params.TargetRect; % centered at (0,0)
             Rect([1,3])    =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
             Rect([2,4])    =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
             Screen('FillOval', Params.WPTR, [200,200,200],  Rect)
             TargetPos = [0, -200];
            end

 % draw trial target color
         TargetPos      = Params.ReachTargets1(TargetNum1,:);
         Rect           = Params.TargetRect; % centered at (0,0)
         Rect([1,3])    =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
         Rect([2,4])    =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
         Screen('FillOval', Params.WPTR, [255,0,0],  Rect)

        end

        for k = 1:Params.NumTargets2
         TargetPos = Params.ReachTargets2(k,:);
         Rect           = Params.TargetRect; % centered at (0,0)
         Rect([1,3])    =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
         Rect([2,4])    =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
         Screen('FillOval', Params.WPTR, [200,200,200],  Rect)
         TargetPos = [0, -200];
        end

        
         % draw trial target color
         TargetPos      = Params.ReachTargets2(TargetNum2,:);
         Rect           = Params.TargetRect; % centered at (0,0)
         Rect([1,3])    =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
         Rect([2,4])    =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
         Screen('FillOval', Params.WPTR, [255,0,0],  Rect)
          
        % color target if cursor inside
        if any(inTarget)
            StartTargetPos      = g(find(inTarget),:);
            StartRect           = Params.TargetRect; % centered at (0,0)
            StartRect([1,3])    = StartRect([1,3]) + StartTargetPos(1) + Params.Center(1); % add x-pos
            StartRect([2,4])    = StartRect([2,4]) + StartTargetPos(2) + Params.Center(2); % add y-pos

            Screen('FillOval', Params.WPTR, [0,200,0], StartRect)
        end
                 
        CursorRect                  = Params.CursorRect;
        CursorRect([1,3])           = CursorRect([1,3]) + Cursor.State(1) + Params.Center(1); % add x-pos
        CursorRect([2,4])           = CursorRect([2,4]) + Cursor.State(2) + Params.Center(2); % add y-pos
        
        % save               
        if assist && Params.ChangeAssistColor && Params.Assist
            Screen('FillOval', Params.WPTR, [255,0,255], CursorRect)
        else
            Screen('FillOval', Params.WPTR, [0,0,255], CursorRect)
        end
                
        % Stop robot at boundaries
        Data.CursorState(:,end+1) = Cursor.State;     
        Cursor.TaskState = 3;
        Data.TaskState(1,end+1) = Cursor.TaskState;
        Data.TaskStage(end + 1)          = taskStage;
        click = 0;
        
        d = 0;
        if any(inTarget)
            if Params.RobotClicker
                if mean(StopClicker_Buffer) > Params.ClickerBinThresh
                    click = 1;
                    d = 1;
                end
            else
                d = 1;
            end
        end
        
        if d == 1
                
            StartTargetPos = g(find(inTarget),:);
            StartRect = Params.TargetRect; % centered at (0,0)
            StartRect([1,3]) = StartRect([1,3]) + StartTargetPos(1) + Params.Center(1); % add x-pos
            StartRect([2,4]) = StartRect([2,4]) + StartTargetPos(2) + Params.Center(2); % add y-pos
            Screen('FillOval', Params.WPTR, [0,200,200], StartRect)
            if taskStage == 1
                Data.SelectedTargetID(1) = find(inTarget);
               taskStage = 2;
               g = Params.ReachTargets2;
                b3(1,:) = ones(1, length(g));
                b3(1,:) = [b3(1,:)]/sum(b3(1,:) );  %combination

                assist  = 0;
                userVel = [0,0,0,0,0,0]';
                Params.RobotClicker = 1;
                t_targetStart = tim;
               % reset inference
            elseif taskStage == 2
                Data.SelectedTargetID(2) = find(inTarget);
                done = 1
                end
            end
        
           Data.Click(end + 1) = click; 
           
        % end if takes too long
        if TotalTime > Params.MaxReachTime
            done = 1;
            Data.ErrorID = 3;
            Data.ErrorStr = 'ReachTarget';
            Data.SelectedTargetID = 0;
            Data.SelectedTargetPosition = NaN;
            fprintf('ERROR: %s\n',Data.ErrorStr)
        end
        
        % end if clicks in a target
        if done == 1
            
            if any(Data.SelectedTargetID ~= [TargetNum1, TargetNum2])
                Data.ErrorID = 4;
                Data.ErrorStr = 'WrongTaget';
            end
        end
 
        Screen('DrawingFinished', Params.WPTR);
        Screen('Flip', Params.WPTR);

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
                binInd = binInd + 1;
                Data.NeuralTime(1,end+1) = tim;
                [Neuro,Data] = NeuroPipeline(Neuro,Data,Params);   
            end
            
            
            Data.CursorState(:,end+1) = Cursor.State;

            Cursor.TaskState = 4;
            Data.TaskState(1,end+1)=Cursor.TaskState;
                     
            Data.ClickerState(1,end+1)      = 0;
            Data.ClickerDistance(1,end+1)   = 0;
            Data.FilteredClickerState(1,end+1)   = 0;
            Data.UserVel(:,end+1)   = [0;0];
            Data.AssistVel(:,end+1)   = [0;0];
            Data.Belief(:,end+1)   = b3;
            Data.InTarget(:,end+1)  = 0;
            Data.Assist(:,end+1)    = 0;
            Data.Dyn(:,end+1)      = 0;
            Data.Click(end + 1)     = 0; 
            Data.TaskStage(end + 1)          = taskStage;
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



