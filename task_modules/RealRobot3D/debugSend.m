clear
% move robot to position

Params.udp = udpport("LocalPort", 43210);
Params.pythonPort = 5006;
% write(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);                  % reset robot
% write(Params.udp, [0,5,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); % open file    
%%
pos = [0, 80, 410];

[xa,xb,xc] = doubleToUDP(pos(1));
[ya,yb,yc] = doubleToUDP(pos(2)); 
[za,zb,zc] = doubleToUDP(pos(3) - 256) ;

write(Params.udp, [4, xa,xb,xc,ya,yb,yc, za,zb,zc, 0], "127.0.0.1", Params.pythonPort) ; % send pos
write(Params.udp, [0,2,5,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);    


Params.k_v      = 0.7;
Params.k_i      = 40;

Params.r_v      = 0.8;
Params.r_i      = 100;

write(Params.udp, [0,26,Params.k_v*10,Params.k_i,Params.r_v*10,Params.r_i,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 
write(Params.udp, [0,16,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort); 


write(Params.udp, [0,1,0,0,0,0,0,0,0,0,0,0], "127.0.0.1", Params.pythonPort);                  % reset robot

pause(1.0)



pause()
clickOrder = [ones(8,1); 4*ones(8,1); 3*ones(8,1); 2*ones(8,1)];
% clickOrder = [ones(10,1)]
for i = 1:length(clickOrder)
    ClickToSend = clickOrder(i);
    write(Params.udp, [5, 0,0,0,0,0,0,0,0,0, ClickToSend], "127.0.0.1", Params.pythonPort);
    pause(0.2)
end