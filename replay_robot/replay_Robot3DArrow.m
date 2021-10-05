% Reanimates data from Robot3DArrow task
% Run python script runRobotClient.py first to start sim env 


clear, close all
addpath('../task_helpers')
%% Specify data directory which contains bravo1 directory
dataDir     = '/media/sarah/Unused/';

% Specify dates, blocks, trials 
dates       = {'20210910'};
blocks      = 1;
trials      = 1:21;

%% Run all trials
  
for dateInd = 1:numel(dates)    % iterate dates
    date        = dates{dateInd};
    taskDir     = [dataDir,'/bravo1/', dates{dateInd},'/GangulyServer/', dates{dateInd}, '/Robot3DArrow/'];
    tmp         = dir(taskDir);
    blockDirs   = {tmp.name};
    blockDirs   = blockDirs(3:end);

for blockInd = 1:numel(blockDirs)   % iterate blocks

    blkDir      = [taskDir, blockDirs{blockInd}, '/BCI_Fixed'];
    tmp         = dir(blkDir);
    files       = {tmp.name};
    files       = files(3:end);
    numTrials   = 0;
    
    %  load task params from first file
    fn          = [blkDir, '/', files{1}];
    load(fn)
    Params      = TrialData.Params;

    Params.udp  = udp("127.0.0.1", 5006);
    fopen(Params.udp)
    fwrite(Params.udp, [0,1,0])                  % reset robot
    fwrite(Params.udp, [0,2,Params.UpdateRate])  % set update rate
    fwrite(Params.udp, [0,3,Params.RobotMode])   % set robot mode
    fwrite(Params.udp, [0,4,Params.RobotDirectionLines])   % set debug lines

for fileInd = 1:numel(files)    % iterate trials
    
    % load data
    fn = [blkDir, '/', files{fileInd}]

    load(fn)
    Params = TrialData.Params;

    % Set up Trial
    Params.udp = udp("127.0.0.1", 5006);
    fopen(Params.udp)

    TD          = TrialData;
    inTargetOld = 0;   
    InTargetTotalTime = 0;

%% Instructed Delay
    done = 0;
    tstart = GetSecs;
    while ~done
        % Update Time
        tim = GetSecs;

        % end if in start target for hold time
        if tim - tstart > Params.InstructedDelayTime
            done = 1;
        end
    end % Instructed Delay Loop

%% Cue Time
    done                = 0;
    tstart              = GetSecs;

    ReachTargetPos = TD.TargetPosition;

    % Send target position to sim
    [xa,xb,xc] = doubleToUDP(ReachTargetPos(1));
    [ya,yb,yc] = doubleToUDP(ReachTargetPos(2)); 
    [za,zb,zc] = doubleToUDP(ReachTargetPos(3)) ;
    fwrite(Params.udp, [1, xa,xb,xc,ya,yb,yc,za,zb,zc, 0]);  

while ~done
    % Update Time
    tim = GetSecs;

    % end if in start target for hold time
    if tim - tstart > Params.CueTime
        done = 1;
    end
end % Cue Delay Loop

%% Go to reach target

    % Set Gripper pos
    [xa,xb,xc] = doubleToUDP(0);
    [ya,yb,yc] = doubleToUDP(0); 
    [za,zb,zc] = doubleToUDP(-50) ;

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
            % time;
            dt = tim - Cursor.LastPredictTime;
            TotalTime = TotalTime + dt;
            Cursor.LastPredictTime = tim;           

            % Send discrete selection
            ClickToSend = TD.FilteredClickerState(i);
            fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, ClickToSend]);
            
            % Count if in target
            if ClickToSend == TD.TargetID
                Counter = Counter+1;
            else
               Counter = 0;
            end
            
            % Send correct selection feedback
            if Counter == Params.ClickCounter
                fwrite(Params.udp, [0, 5, 0])
            end
        end

    end

%% Inter trial interval
    
    pause(0.125)
    fwrite(Params.udp, [0,1,0])
    while ~done
    % Update Time
    tim = GetSecs;

    % end if in start target for hold time
%     if tim - tstart > Params.InterTrialInterval
    if tim - tstart > 1
        done = 1;
    end
end % Inter trial interval loop


end
end
end
