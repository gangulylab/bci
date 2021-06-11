dates = {'20210312'};

a{1} = [1];
% a{3} = [2,3];
a{4} = [3];

for dateInd = 1:numel(dates)
    date = dates{dateInd};
    taskDir = ['/media/sarah/VICTOR-2/bravo1/', date, '/GangulyServer/', date, '/Robot/'];
    taskDir = ['/home/sarah/Downloads/RobotR2GAuto/'];
    tmp = dir(taskDir);
    blockDirs = {tmp.name};
    blockDirs = blockDirs(3:end);
% for blockInd = 1:numel(blockDirs)
for blockInd = [1]
    blkDir = [taskDir, blockDirs{blockInd}, '/BCI_Fixed'];
    tmp = dir(blkDir);
    files = {tmp.name};
    files = files(3:end);
    numTrials = 0;
for fileInd = a{blockInd}
%    for fileInd = 1:numel(files)
    fn = [blkDir, '/', files{fileInd}];
    fileInd
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
TD.Params.AssistAlpha
inTargetOld = 0;

inTargetOld = 0;
fing = 0.0;
AutoComplete = 0;

%% Instructed Delay
done = 0;
InTargetTotalTime = 0;
tstart = GetSecs;
while ~done
    % Update Time & Position
    tim = GetSecs;

    % end if in start target for hold time
    if tim - tstart > Params.InstructedDelayTime
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

fwrite(Params.udp, [11, xa,xb,xc,ya,yb,yc,za,zb,zc, 0]);  
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
            
            
            Cursor.State = TD.CursorState(:,i+stInd-1);
           TargetID = InTargetRobot3D(Cursor,Params.ReachTargetPositions,Params.RobotTargetRadius, Params.RobotTargetDim, TD.TargetID);
            
            if TargetID == TD.TargetID
                fwrite(Params.udp, [0, 15, 1])
                task = "HERE"
            end

%             decision for clikcing and finishng trial
%             if Cursor.Counter == Params.ClickCounter
%                 done=1;                        
%                 Cursor.State(1:2) = Params.ReachTargetPositions(Click_Decision,:);
%                 Data.SelectedTargetID = Click_Decision;
%     
%                fwrite(Params.udp, [0, 5, 0])
%             end
%             
            RunningMode_ClickDec = TD.FilteredClickerState(i);
            ClickToSend = RunningMode_ClickDec;
            
            
                            if RunningMode_ClickDec == 7
                    Cursor.State(4:6) = [0;0;0];          
                else 

                A = Params.dA;
                B = Params.dB;
                
                U = zeros(3,1);
                U(1) = int8(RunningMode_ClickDec == 1) - int8(RunningMode_ClickDec == 3);
                U(2) = int8(RunningMode_ClickDec == 2) - int8(RunningMode_ClickDec == 4);
                U(3) = int8(RunningMode_ClickDec == 5) - int8(RunningMode_ClickDec== 6);
                
                
                Data.TargetPosition = ReachTargetPos;
                vTarget = (Data.TargetPosition'- Cursor.State(1:3));
                norm_vTarget = vTarget/norm(vTarget);
                
                AssistVel = Params.AssistAlpha*B*norm_vTarget;

                 
                if norm(vTarget) < 15
                    if fing < 0.9
                        fing = fing + .1;
                     
                        [xa,xb,xc] = doubleToUDP(fing*80);

                        fwrite(Params.udp, [6, xa,xb,xc,0,0,0,0,0,0,0]);
                        
                        countRemoveBlock = 0;
                        ClickToSend = 0;
    
                    else
                         countRemoveBlock = countRemoveBlock + 1;
                         ClickToSend = 0;
                         if countRemoveBlock > 16
                             ClickToSend = 0;
                            done = 1;
                         end
                     end

                elseif norm(vTarget) < 100
                    AutoComplete = 1;
                    Cursor.State(1:3) = Cursor.State(1:3) + norm_vTarget*norm(vTarget)/10;
                    ClickToSend = 0;
                else
                
                Cursor.State = A*Cursor.State + (1-Params.AssistAlpha)*B*U + AssistVel;
                Cursor.IntendedState = [0 0 0 0 0]';              
                
                % Stop robot at boundaries
                
                end
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
                end
                
            %%%%% UPDATE CURSOR STATE OR POSITION BASED ON DECODED
            %%%%% DIRECTION

            [xa,xb,xc] = doubleToUDP(Cursor.State(1));
            [ya,yb,yc] = doubleToUDP(Cursor.State(2)); 
            [za,zb,zc] = doubleToUDP(Cursor.State(3)) ;
            
            fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, ClickToSend]);

            
            if TargetID==TD.TargetID
                inTargetOld = 1;
                InTargetTotalTime = InTargetTotalTime + dt;
                if Params.RobotClicker
                    if RunningMode_ClickDec == 7
                        done = 1;
                        Data.SelectedTargetID = TargetID;
                        Data.SelectedTargetPosition = Params.ReachTargetPositions(TargetID,:); 
                        fwrite(Params.udp, [0, 6, 0])
                    end
                end
            else
                InTargetTotalTime = 0;
                if inTargetOld
                    fwrite(Params.udp, [0, 7, 0])
                end
                inTargetOld = 0;
            end

        end

    end
    
    %%
    fwrite(Params.udp, [0,1,0])
    pause(1.0)

end
end
end