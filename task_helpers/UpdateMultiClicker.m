function UpdateMultiClicker(Params, Neuro, Clicker)
% Applies discrete decoder to update click state (within cursor)
% Clicker Code:
% Click_Decision = {0-no click, 1-click_class1, 2-click_class2, ...}
% Cursor.ClickState = [ 1 x num_click_classes ]

global Cursor

if (Params.GenNeuralFeaturesFlag), 
   [~,~,B] = GetMouse(); % multi-click
   B = B([1,3,2]); % swap 'scroll click' and 'right click'
   Click_Decision = find(B);
   if isempty(Click_Decision),
       Click_Decision = 0;
   end
else,
   [ Click_Decision,~] = Clicker.Func(Neuro.NeuralFeatures);
end

% must click for X bins in a row
if Click_Decision, % clicking
    Cursor.ClickState(setdiff(1:Params.NumClickerClasses, Click_Decision)) = 0;
    Cursor.ClickState(Click_Decision) = Cursor.ClickState(Click_Decision) + 1;
    return;
else, % not clicking
    Cursor.ClickState = zeros(1,Params.NumClickerClasses);
end

end % UpdateMultiClicker