clear
% move robot to position

Params.udp = udpport("LocalPort", 43210);
Params.pythonPort = 5006;
write(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);                  % reset robot

%%

pos = [240, -55, 410];

pos = [240, -90, 420];
% pos = [240, -220, 410];
% pos = [200, 80, 410];
% pos = [190, 180, 470];
% pos = [200, 30, 470];

[xa,xb,xc] = doubleToUDP(pos(1));
[ya,yb,yc] = doubleToUDP(pos(2)); 
[za,zb,zc] = doubleToUDP(pos(3) - 256) ;

write(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, 0], "127.0.0.1", Params.pythonPort) ; % send pos
write(Params.udp, [0,2,5,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);       
write(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);                  % reset robot
