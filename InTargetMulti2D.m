function n = InTargetMulti2D(Cursor,t)

x = Cursor.State(1:2)';

d = vecnorm([x - t]');

n = d < 50;


end