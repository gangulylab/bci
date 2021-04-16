colors = [rgb('crimson'); rgb('orange'); rgb('mediumspringgreen');...
    rgb('dodgerblue');  rgb('MediumBlue');rgb('Deeppink'); rgb('orangered'); rgb('lime'); rgb('mediumturquoise'); rgb('blueviolet') ];




TD = TrialData;

if TD.TargetID >6
    c = colors(TD.TargetID,:);
else
    c = colors(TD.TargetID,:);
end
pos = TD.CursorState;



subplotOrder = [1, 5, 4, 2, 3, 6];