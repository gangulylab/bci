fn = '/media/sarah/VICTOR-2/bravo1/20210219/GangulyServer/20210219/RobotRR/150252/BCI_Fixed/Data0001.mat';

load(fn)
Params = TrialData.Params;

Params.udp = udp("127.0.0.1", 5006);
fopen(Params.udp)
fwrite(Params.udp, [0,1,0])                  % reset robot
fwrite(Params.udp, [0,2,Params.UpdateRate])  % set update rate
fwrite(Params.udp, [0,3,Params.RobotMode])   % set robot mode
fwrite(Params.udp, [0,4,Params.RobotDirectionLines])   % set debug lines
%%

fn = '/media/sarah/VICTOR-2/bravo1/20210219/GangulyServer/20210219/RobotRR/150252/BCI_Fixed/Data0001.mat';
load(fn)

TD = TrialData;

ReachTargetPos = TD.Params.ReachTargetPositions(TD.Params.TargetOrder(1),:);

[xa,xb,xc] = doubleToUDP(ReachTargetPos(1));
[ya,yb,yc] = doubleToUDP(ReachTargetPos(2)); 
[za,zb,zc] = doubleToUDP(ReachTargetPos(3)) ;

fwrite(Params.udp, [11, xa,xb,xc,ya,yb,yc,za,zb,zc, 0]);

ind = find(TD.TaskState == 3);

CS = TD.CursorState(:,ind);

Data.TargetID  = TD.Params.TargetOrder(1);
cnt = 0;
target = 1;
intarget2 = 0;
for i = 1:length(CS)
    
    cs = CS(:,i);
    ClickToSend = TD.FilteredClickerState(i);
    [xa,xb,xc] = doubleToUDP(cs(1));
    [ya,yb,yc] = doubleToUDP(cs(2)); 
    [za,zb,zc] = doubleToUDP(cs(3)) ;

    fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, ClickToSend]);
    norm(cs(1:3) - ReachTargetPos)
    
    
    if norm(cs(1:3) - ReachTargetPos') < 100
        if target == 1
            fwrite(Params.udp, [0, 15, 1])
            cnt = cnt + 1;
        else
            if intarget2 == 0
             fwrite(Params.udp, [0, 14, 1])
             intarget2  = 1;
            end
        end
        
        if cnt > TD.Params.TargetHoldTime*8*2
            fwrite(Params.udp, [0, 15, 2])
            
            if target == 1
            ReachTargetPos = Params.ReachTargetPositions(1,:);

                [xa,xb,xc] = doubleToUDP(ReachTargetPos(1));
                [ya,yb,yc] = doubleToUDP(ReachTargetPos(2)); 
                [za,zb,zc] = doubleToUDP(ReachTargetPos(3)) ;

                fwrite(Params.udp, [12, xa,xb,xc,ya,yb,yc,za,zb,zc, 0]);
                cnt  = 0;
                target = 2;
            end
        end
    end
    
    pause(1/8 - .01)
    
end


%fwrite(Params.udp, [0,1,0])    