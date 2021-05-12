% u = udp('127.0.0.1','LocalHost', '','LocalPort',8824);
% fopen(u);
% 
% while(1)
%     A = fread(u,10);
% end
% 
% fclose(u)
% echoudp("on",3030);
% u = udpport('LocalPort', 8000, "EnablePortSharing", true)
% while(1)
%     u.NumDatagramsAvailable
%     A = read(u,10);
% end

% echoudp("on",3030);
% u = udpport("datagram","OutputDatagramSize",5);
% 
% write(u,1:20,"uint8","127.0.0.1",3030);


echoudp("on",3040);
u = udpport;

v =udpport("LocalPort", 43210);

