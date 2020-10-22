function [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params, Neuro, Clicker)
% Multistate Clicker update
% Checks the current decoded target
% If it is the correct target, it updates the Cursor.ClickState

global Cursor

if Params.ControlMode == 2 %mouse    
    Click_Decision = randperm(5,1)-1;
    Click_Distance = 0;
else,
    if Params.SmoothDataFlag ==1
        [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.FilteredFeatures);
    else        
        [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.NeuralFeatures);
    end
end

end % UpdateMultiClicker