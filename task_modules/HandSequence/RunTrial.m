function [Data, Neuro, KF, Params, Clicker] = RunTrial(Data,Params,Neuro,TaskFlag,KF,Clicker)
% Runs a trial, saves useful data along the way
% Each trial contains the following pieces
% 1) Get the cursor to the reach target (different on each trial)
% 2) Feedback

global Cursor

%% Set up trial
% ReachTargetPos = Data.TargetPosition;
TargetID = 0; % Target that cursor is in, 0 for no targets

% Output to Command Line
fprintf('\nTrial: %i\n',Data.Trial)

% keep track of update times
dt_vec = [];
dT_vec = [];

% grab blackrock data and run through processing pipeline
if Params.BLACKROCK
    Cursor.LastPredictTime = GetSecs;
    Cursor.LastUpdateTime = Cursor.LastPredictTime;
    Neuro = NeuroPipeline(Neuro,[],Params);
end

Cursor.State = zeros(8,1);
Cursor.ClickState       = 0;
Cursor.ClickDistance    = 0;

fwrite(Params.udp, [0,16,Data.TargetID]);

target_order = [];
targetCnt = [];
cue_text = '';

action_text = {'T', 'I', 'M', 'R', 'P', 'Po', 'Pi', 'Tr', 'A','A', 'F','F' };
trial_actions = Data.TrialActions;


for i = 1:length(trial_actions)
    num = randi(3);
    action = trial_actions(i);
    temp = ones(1,num)*action;
    target_order = [target_order, temp,0];
    tmp_txt = '';
    for j = 1:num
        tmp_txt(end+1:end+length(action_text{action})) = action_text{action};
    end
    
    CUE_txt{i} = tmp_txt;
    targetCnt = [targetCnt, ones(1,num)*i,0];
end

st = 0.;
en = -0.8;

handkin =  zeros(8,1);
ind = 1;  
for i = 1:length(target_order)
    
    numbins = Params.ActionBins(randi(length(Params.ActionBins)));
    
    go1 = linspace(st,en,numbins);
    go2 = linspace(en,st,numbins);
    wait = go2(end)*ones(1,Params.PauseBins);
    go      = [go1, go2(2:end), wait];
    
    numPauseBins =  Params.TransitionBins(randi(length(Params.TransitionBins)));
    pause = zeros(1,numPauseBins);
    
    if target_order(i) > 0
        idx = ind:ind + length(go)-1;
        if target_order(i) < 6 
            handkin(target_order(i),idx) = go;
        elseif target_order(i) == 6
            for j = 1:5
                handkin(j,idx) = go;
            end
        elseif target_order(i) == 7
            for j = 1:2
               handkin(j,idx) = go;
            end
        elseif target_order(i) == 8
            for j = 1:3
                handkin(j,idx) = go;
            end
        elseif target_order(i) == 9 || target_order(i) == 10
            handkin(6,idx) = go;
        elseif target_order(i) == 11 || target_order(i) == 12
            handkin(7,idx) = go;
        end
    else
        idx = ind:ind+length(pause)-1;
        handkin(:,idx) = zeros(size(handkin,1), length(idx));
    end    
    
    ind = idx(end) + 1;
    target(idx) = target_order(i)*ones(1,length(idx));
    target_cnt(idx) = targetCnt(i);
end

    
%% Instructed Delay
if ~Data.ErrorID && Params.InstructedDelayTime>0
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Instructed Delay';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
    
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;

    [xa,xb,xc] = doubleToUDP(Cursor.State(1)*80);
    [ya,yb,yc] = doubleToUDP(Cursor.State(2)*80); 
    [za,zb,zc] = doubleToUDP(Cursor.State(3)*80) ;
    if Params.handVis
        fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, Data.TargetID]);
        fwrite(Params.udp, [0,2,0])
    
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
            
            Data.TaskState(1,end+1) = 1;
            Data.HandState(:,end+1) = zeros(8,1);
            Data.Action(1,end+1) = 0; 
            
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

    Screen('TextSize', Params.WPTR, 80);
    
    l = length(cell2mat(CUE_txt));
    d = 0.9*Params.ScreenRectangle(3)/(l + 5*length(trial_actions));
     
    s = 100;
    for i = 1:length(trial_actions)
        text_st(i) = s;
        s = s + (length(CUE_txt{i})+ 5)*d;
        DrawFormattedText(Params.WPTR, CUE_txt{i}, text_st(i),'center', [255,255,255]);
    end
%     DrawFormattedText(Params.WPTR, cue_text, 'center','center', [255,255,255]);
    Screen('DrawingFinished', Params.WPTR);
    Screen('Flip', Params.WPTR);

    
    while ~done
        % Update Time & Position
        tim = GetSecs;
        
        % for pausing and quitting expt
        if CheckPause, [Neuro,Data,Params] = ExperimentPause(Params,Neuro,Data); end
        
        % Update Screen
        if (tim-Cursor.LastPredictTime) > 1/Params.ScreenRefreshRate
            % timeCUE
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
            
            Data.TaskState(1,end+1) = 2;
            Data.HandState(:,end+1) = zeros(8,1);
            Data.Action(1,end+1) = 0; 

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
%     fwrite(Params.udp, [5, Data.TargetID, 1]) % DISPLAY GO CUE
    k = 0;
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
                k = k+1;
                dT = tim-Cursor.LastUpdateTime;
                dT_vec(end+1) = dT;
                Cursor.LastUpdateTime = tim;
                
                Data.NeuralTime(1,end+1) = tim;
                [Neuro,Data] = NeuroPipeline(Neuro,Data,Params);              

            end

            % Draw Text Cue
            for i = 1:length(trial_actions)
                DrawFormattedText(Params.WPTR, CUE_txt{i}, text_st(i),'center', [255,255,255]);
            end

            i = target_cnt(k);
            if i > 0
            DrawFormattedText(Params.WPTR, CUE_txt{i}, text_st(i),'center', [255,0,0]);
            end
            
            Screen('DrawingFinished', Params.WPTR);
            Screen('Flip', Params.WPTR);

        % Send Cursor State to hand visualization        
            [ta,tb,tc] = doubleToUDP(handkin(1,k)*80);
            [va,vb,vc] = doubleToUDP(handkin(2,k)*80);
            [wa,wb,wc] = doubleToUDP(handkin(3,k)*80); 
            [xa,xb,xc] = doubleToUDP(handkin(4,k)*80);
            [ya,yb,yc] = doubleToUDP(handkin(5,k)*80); 
            [za,zb,zc] = doubleToUDP(handkin(6,k)*80);
            [wfa,wfb,wfc] = doubleToUDP(handkin(7,k)*80);
            [wra,wrb,wrc] = doubleToUDP(handkin(8,k)*80);

            fwrite(Params.udp, [14, ta,tb,tc,va,vb,vc, wa,wb,wc, Data.TargetID]);
            fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, Data.TargetID]);
            fwrite(Params.udp, [8, wfa,wfb,wfc,wra,wrb,wrc,0,0,0 Data.TargetID]);
            fwrite(Params.udp, [15, target(k),0,0,0,0,0,0,0,0,0,0]);
       
            Data.TaskState(1,end+1) = 3;
            Data.HandState(:,end+1) = handkin(:,k);
            Data.Action(1,end+1) = target(k);
            
    if k == length(handkin)
        done = 1;
    end
    
    end % Reach Target Loop
   
end % only complete if no errors
end
%% Inter trial interval
% blank screen at end of trial but continue collecting data
Screen('Flip', Params.WPTR);
if Params.InterTrialInterval>0
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Inter Trial Interval';
    if Params.ArduinoSync, PulseArduino(Params.ArduinoPtr,Params.ArduinoPin,length(Data.Events)); end
    
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;

    fwrite(Params.udp, [0,1,0])
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

            Data.TaskState(1,end+1) = 4;
            Data.HandState(:,end+1) = zeros(8,1);
            Data.Action(1,end+1) = 0;
            
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
    Cursor.IntendedState = [0,0,0,0,1]';
    
    fprintf('ERROR: %s\n', Data.ErrorStr)
    
    if Params.FeedbackSound
        sound(Params.ErrorSound,Params.ErrorSoundFs)
    end
    WaitSecs(Params.ErrorWaitTime);
end

end % RunTrial