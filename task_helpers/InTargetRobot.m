function TargetID = InTargetRobot(Cursor,TargetPos, Radius, dim, ID)
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

if dim == 2
    temp = [Cursor.State(1) - Cursor.Center(1) ;Cursor.State(2) - Cursor.Center(2)];

  for i=1:length(TargetPos)
    dist = temp - [TargetPos(i,1); TargetPos(i,2)];
    
    in = norm(dist) < Radius;
    if in, TargetID = i; end
    end

elseif dim == 1
    target = TargetPos(ID,:)';
    index = find(target);
    
    temp = [Cursor.State(1) - Cursor.Center(1) ;Cursor.State(2) - Cursor.Center(2)];
    dist = temp(index) - target(index);
    
    in = norm(dist) < Radius;
    if in, TargetID = ID; end
end

end % InTarget

