function modeVal = RunningMode(v)

n(1) = sum(v == 1);
n(2) = sum(v == 2);
n(3) = sum(v == 3);
n(4) = sum(v == 4);
n(5) = sum(v == 5);
n(6) = sum(v == 6);
n(7) = sum(v == 0);


maxVal = max(n);
idx = find(n == maxVal);

if (length(idx) >  1) % tie
    modeVal = 0;
elseif (maxVal <= 0.5* length(v)) % 
    modeVal = 0;
elseif idx == 7
    modeVal = 0;
else
    modeVal = idx;
    
end