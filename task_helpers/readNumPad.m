function key = readNumPad()

map = [86, 89, 84, 81, 80, 88, 85];

[down, ~, code] = KbCheck();

if down % Key pressed
    index   = find(code,1);  
    if any(index == map) % valid key
        key = find(index == map);
    else % invalid key pressed
        key = 0;
    end
else % No key pressed
    key = 0;
end

end