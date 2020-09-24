function [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params, Neuro, Clicker)
% Multistate Clicker update
% Checks the current decoded target
% If it is the correct target, it updates the Cursor.ClickState

global Cursor

if (Params.GenNeuralFeaturesFlag), 
   [~,~,B] = GetMouse(); % multi-click
   B = B([1,3,2]); % swap 'scroll click' and 'right click'
   Click_Decision = find(B);
   if isempty(Click_Decision),
       Click_Decision = 0;
   end
else,
   [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.NeuralFeatures);
end


% 
% % must click for X bins in a row
% if Click_Decision, % clicking
%     Cursor.ClickState(setdiff(1:Params.NumClickerClasses, Click_Decision)) = 0;
%     Cursor.ClickState(Click_Decision) = Cursor.ClickState(Click_Decision) + 1;
%     return;
% else, % not clicking
%     Cursor.ClickState = zeros(1,Params.NumClickerClasses);
% end

end % UpdateMultiClicker