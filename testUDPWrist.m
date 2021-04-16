ang = [0:.01:pi/2];
ang = [0:.01:pi/2];
for i = 1:length(ang)
Cursor.State(1) = ang(i);
Cursor.State(2) = ang(i);
Cursor.State(3) = ang(i);

%%%%% UPDATE CURSOR STATE

[xa,xb,xc] = doubleToUDP(Cursor.State(1)*80);
[ya,yb,yc] = doubleToUDP(Cursor.State(2)*80);
[za,zb,zc] = doubleToUDP(Cursor.State(3)*80) ;

fwrite(Params.udp, [6, xa,xb,xc,ya,yb,yc, za,zb,zc, Data.TargetID]);
pause(.01)
end