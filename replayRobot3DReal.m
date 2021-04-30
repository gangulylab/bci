dates = {'20210326'};
pause(5.0)
for dateInd = 1:numel(dates)
    date = dates{dateInd};
    taskDir = ['/media/sarah/Excellente/bravo1/', date, '/GangulyServer/', date, '/Robot/'];
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
% for fileInd = 1
   for fileInd = 1:numel(files)
    fn = [blkDir, '/', files{fileInd}];
    fileInd
% 
% fn = '/media/sarah/VICTOR-2/bravo1/20210319/GangulyServer/20210319/Robot3DArrow/110258/BCI_Fixed/Data0001.mat';



load(fn)
Params = TrialData.Params;

Params.udp = udp("127.0.0.1", 5006);
fopen(Params.udp)
fwrite(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0,0,0])                  % reset robot
% fwrite(Params.udp, [0,2,Params.UpdateRate])  % set update rate
% fwrite(Params.udp, [0,3,Params.RobotMode])   % set robot mode
% fwrite(Params.udp, [0,4,Params.RobotDirectionLines])   % set debug lines
%%

% fn = '/media/sarah/VICTOR-2/bravo1/20210219/GangulyServer/20210219/RobotRR/150252/BCI_Fixed/Data0001.mat';
load(fn)

TD = TrialData;
TD.Params.AssistAlpha
inTargetOld = 0;

if TD.TargetID <5

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

stInd = find(TD.TaskState==3,1);
tstart  = GetSecs;
Cursor.Counter = 0;
Cursor.Center = [0;0;0]';

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
            i = i + 1
            % time
            dt = tim - Cursor.LastPredictTime;
            TotalTime = TotalTime + dt;
            Cursor.LastPredictTime = tim;

             % grab and process neural data
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate)
                dT = tim-Cursor.LastUpdateTime;
                Cursor.LastUpdateTime = tim;
                
            end
            
            
            Cursor.State = TD.CursorState(:,i+stInd-1);
           TargetID = InTargetRobot3D(Cursor,Params.ReachTargetPositions,Params.RobotTargetRadius, Params.RobotTargetDim, TD.TargetID);

%             decision for clikcing and finishng trial
%             if Cursor.Counter == Params.ClickCounter
%                 done=1;                        
%                 Cursor.State(1:2) = Params.ReachTargetPositions(Click_Decision,:);
%                 Data.SelectedTargetID = Click_Decision;
%     
%                fwrite(Params.udp, [0, 5, 0])
%             end
%           
            Cursor.State

            ClickToSend = TD.FilteredClickerState(i);
            [xa,xb,xc] = doubleToUDP(Cursor.State(1));
            [ya,yb,yc] = doubleToUDP(Cursor.State(2)); 
            [za,zb,zc] = doubleToUDP(0) ;
            
            fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, ClickToSend]);
            

            [xa,xb,xc] = doubleToUDP(Cursor.State(4));
            [ya,yb,yc] = doubleToUDP(Cursor.State(5)); 
            [za,zb,zc] = doubleToUDP(0) ;
            
            fwrite(Params.udp, [5, xa,xb,xc,ya,yb,yc, za,zb,zc, ClickToSend]);
            
            
            RunningMode_ClickDec = TD.FilteredClickerState(i);
            
           

        end

    end
    
    %%
    fwrite(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
    pause(1.0)

end
end
end
end