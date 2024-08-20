function target = InTargetMulti2D(Cursor, targets, radius)

dist = vecnorm([[Cursor.State(1),Cursor.State(2)] - targets]');
target = dist < radius;
a = 1;
end