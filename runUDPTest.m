Params.udp = udpport("LocalPort", 43210);

write(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", 5006)                  % reset robot

% 
% [xa,xb,xc] = doubleToUDP(Params.StartPos(1));
% [ya,yb,yc] = doubleToUDP(Params.StartPos(2)); 
% [za,zb,zc] = doubleToUDP(Params.StartPos(3) - 256) ;
% 
% fwrite(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, 0]); % send pos
% fwrite(Params.udp, [0,2,Params.RobotMode,0,0,0,0,0,0,0,0,0])
% fwrite(Params.udp, [0,3,Params.RobotTargetRadius,0,0,0,0,0,0,0,0,0])
% fwrite(Params.udp, [0,4,Params.TargetHoldTime,0,0,0,0,0,0,0,0,0])
% 
% 
% Params.udpIn = udpport("LocalPort", 43210);