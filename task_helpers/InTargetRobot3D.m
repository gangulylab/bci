function TargetID = InTargetRobot3D(Cursor,TargetPos, Radius, dim, ID)
% function inFlag = InTargetRadial(Cursor,TargetsVerts,InnerCircleRadius)
% function to tell if cursor is inside any of the targets
% 
% Inputs:
%   Cursor - includes .State (posx and posy)
%   TargetsVerts - Cell array of vertices defining each target
%     cell is length: # of targets
%     each cell is a matrix: rows are vertices, w/ x and y as cols
%   InnerCircleRadius - excludes selections within radius
% 
% Outputs:
%   TargetID - index of target that cursor is in, (0 if not in any)

TargetID = 0;


temp = Cursor.State(1:3) - Cursor.Center(1:3)';

for i=1:length(TargetPos)
dist = temp - [TargetPos(i,1); TargetPos(i,2); TargetPos(i,3)];

in = norm(dist) < Radius;
if in, TargetID = i; end
end

a = 1;


