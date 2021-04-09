dates = {'20210402'};
pause(2.0)
for dateInd = 1:numel(dates)
    date = dates{dateInd};
    taskDir = ['/media/sarah/VICTOR-2/bravo1/', date, '/GangulyServer/', date, '/Robot3DArrow/'];
    tmp = dir(taskDir);
    blockDirs = {tmp.name};
    blockDirs = blockDirs(3:end);
% for blockInd = 1:numel(blockDirs)
for blockInd = 2
    blkDir = [taskDir, blockDirs{blockInd}, '/BCI_Fixed'];
    tmp = dir(blkDir);
    files = {tmp.name};
    files = files(3:end);
    numTrials = 0;
for fileInd = 1:numel(files)
    fn = [blkDir, '/', files{fileInd}];
% 
% fn = '/media/sarah/VICTOR-2/bravo1/20210319/GangulyServer/20210319/Robot3DArrow/110258/BCI_Fixed/Data0001.mat';



load(fn)
Params = TrialData.Params;

Params.udp = udp("127.0.0.1", 5006);
fopen(Params.udp)
fwrite(Params.udp, [0,1,0])                  % reset robot
fwrite(Params.udp, [0,2,Params.UpdateRate])  % set update rate
fwrite(Params.udp, [0,3,Params.RobotMode])   % set robot mode
fwrite(Params.udp, [0,4,Params.RobotDirectionLines])   % set debug lines
%%

% fn = '/media/sarah/VICTOR-2/bravo1/20210219/GangulyServer/20210219/RobotRR/150252/BCI_Fixed/Data0001.mat';
load(fn)

TD = TrialData;

%% Instructed Delay
done = 0;
InTargetTotalTime = 0;
tstart = GetSecs;
while ~done,
    % Update Time & Position
    tim = GetSecs;

    % end if in start target for hold time
    if tim - tstart > Params.InstructedDelayTime,
        done = 1;
    end
end % Instructed Delay Loop

%% Cue Time
done = 0;
InTargetTotalTime = 0;
tstart = GetSecs;

ReachTargetPos = TD.TargetPosition;

% Send target position
[xa,xb,xc] = doubleToUDP(ReachTargetPos(1));
[ya,yb,yc] = doubleToUDP(ReachTargetPos(2)); 
[za,zb,zc] = doubleToUDP(ReachTargetPos(3)) ;

fwrite(Params.udp, [1, xa,xb,xc,ya,yb,yc,za,zb,zc, 0]);  
while ~done
    % Update Time & Position
    tim = GetSecs;

    % end if in start target for hold time
    if tim - tstart > Params.InstructedDelayTime
        done = 1;
    end
end % Instructed Delay Loop

%% Go to reach target

tstart  = GetSecs;
Cursor.Counter = 0;

    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;
            tim = GetSecs;
    
    Cursor.LastPredictTime = tim;
    Cursor.LastUpdateTime = tim;
    
    i = 0;

    while i < length(TD.ClickerState)
        % Update Time & Position
        tim = GetSecs;
        
                % Update Screen
        if (tim-Cursor.LastPredictTime) > 1/Params.ScreenRefreshRate
            i = i + 1;
            % time
            dt = tim - Cursor.LastPredictTime;
            TotalTime = TotalTime + dt;
            Cursor.LastPredictTime = tim;

             % grab and process neural data
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate)
                dT = tim-Cursor.LastUpdateTime;
                Cursor.LastUpdateTime = tim;
                
            end
            
            Params.TargetID =  TD.TargetID;
                  
            % counter only if correct target is hit, training mode for now
            
            Click_Decision = TD.ClickerState(i);
            
            if Click_Decision == TD.TargetID
                Cursor.Counter = Cursor.Counter+1;
            else
                Cursor.Counter = 0;
            end
            
            % decision for clikcing and finishng trial
            if Cursor.Counter == Params.ClickCounter
                done=1;                        
                %Cursor.State(1:2) = Params.ReachTargetPositions(Click_Decision,:);
                Data.SelectedTargetID = Click_Decision;
    
               fwrite(Params.udp, [0, 5, 0])
            end
            
            Cursor.TaskState = 3;     
            
            % draw the arrow
           fwrite(Params.udp, [2, Click_Decision, 0])
        end

    end
    
    %%
    fwrite(Params.udp, [0,1,0])
    pause(1.0)

end
end
end