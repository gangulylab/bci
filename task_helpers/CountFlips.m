function cnt = CountFlips(x)

last = 0;
cnt  = 0;

for i = 1:length(x)
    val = x(i);
    
    if val ~= last && val ~=0
        cnt = cnt + 1;
        last = val;
    end
end

end