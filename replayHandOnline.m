% Reanimates data from Robot3DArrow task
% Run python script runRobotClient.py first to start sim env 


clear, close all
addpath('../task_helpers')
%% Specify data directory which contains bravo1 directory
dataDir     = '/media/sarah/Unused/';

% Specify dates, blocks, trials 
dates       = {'20220225'};
blocks      = 3;
trials      = 1:8;

%% Run all trials
  
correct  = 0;


for dateInd = 1:numel(dates)    % iterate dates
    date        = dates{dateInd};
    taskDir     = [dataDir,'/bravo1/', dates{dateInd},'/GangulyServer/', dates{dateInd}, '/HandOnline/'];
    tmp         = dir(taskDir);
    blockDirs   = {tmp.name};
    blockDirs   = blockDirs(3:end);

% for blockInd = 1:numel(blockDirs)   % iterate blocks
for blockInd = blocks
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

for fileInd = 1:numel(files)    % iterate trials
    
    % load data
    fn = [blkDir, '/', files{fileInd}]

    load(fn)
    Params = TrialData.Params;

    % Set up Trial
    Params.udp = udp("127.0.0.1", 5006);
    fopen(Params.udp)

    TD          = TrialData;
    
    
    fwrite(Params.udp, [0,16,TD.TargetID]);

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

    [xa,xb,xc] = doubleToUDP(0*80);
    [ya,yb,yc] = doubleToUDP(0*80); 
    [za,zb,zc] = doubleToUDP(0*80) ;
    if Params.handVis
        fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, TD.TargetID]);
        fwrite(Params.udp, [0,2,0])
    end
fwrite(Params.udp, [5, TD.TargetID, 0])
while ~done
    % Update Time
    tim = GetSecs;

    % end if in start target for hold time
    if tim - tstart > Params.CueTime
        done = 1;
    end
end % Cue Delay Loop

%% Go to reach target
 fwrite(Params.udp, [5, TD.TargetID, 1]) 
    % Set Gripper pos
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
            
            offset = find(TD.TaskState == 3,1);
            Cursor.State = TD.CursorState(:,i+offset);
            
            % Send discrete selection
            ClickToSend = TD.FilteredClickerState(i);
            [ta,tb,tc] = doubleToUDP(Cursor.State(1)*80);
            [va,vb,vc] = doubleToUDP(Cursor.State(2)*80);
            [wa,wb,wc] = doubleToUDP(Cursor.State(3)*80); 
            [xa,xb,xc] = doubleToUDP(Cursor.State(4)*80);
            [ya,yb,yc] = doubleToUDP(Cursor.State(5)*80); 
            [za,zb,zc] = doubleToUDP(Cursor.State(6)*80);

            fwrite(Params.udp, [14, ta,tb,tc,va,vb,vc, wa,wb,wc, TD.TargetID]);
            fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, TD.TargetID]);

            fwrite(Params.udp, [15, ClickToSend,0,0,0,0,0,0,0,0,0,0]);
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


correct = correct + 1*(TD.TargetID == TD.SelectedTargetID);
end
end
end
