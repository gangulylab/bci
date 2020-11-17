function modeVal = runningMode(v)

n(1) = sum(v == 1);
n(2) = sum(v == 2);
n(3) = sum(v == 3);
n(4) = sum(v == 4);


maxVal = max(n);
idx = find(n == maxVal);

if (length(idx) >  1)
    modeVal = 0;
elseif (maxVal <= 0.5* length(v))
    modeVal = 0;
else
    modeVal = idx;
    
end