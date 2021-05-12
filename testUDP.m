clear all
% echoudp('on',5007)
% Params.udp = udp("127.0.0.1", 5007);

Params.udp = udp("0.0.0.0", 'LocalPort', 5007);
fopen(Params.udp)
% fwrite(Params.udp, [0,1,0])                  % reset robot
u = Params.udp;
a = 0;
while 1
get (u, 'Status');
rec = get (u, 'ValuesReceived')

if rec
    data = fread(Params.udp,1)
end
    pause(0.125)
    a = a+1;
end

