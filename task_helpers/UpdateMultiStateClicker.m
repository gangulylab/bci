function [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params, Neuro, Clicker)
% Multistate Clicker update
% Checks the current decoded target
% If it is the correct target, it updates the Cursor.ClickState

global Cursor

if Params.ControlMode == 2 %mouse
    Click_Decision = randperm(5,1)-1;
    Click_Distance = 0;
else,
    if Params.NeuralNetFlag == 1
        %if Params.SmoothDataFlag==1
        X = Neuro.FilteredFeatures;
        X = X(:);
        X = X(769:end);
        Decision_Prob = multistate_discrete_6Target(X);
        [aa bb]=max(Decision_Prob);
        if aa >= Params.NeuralNetSoftMaxThresh
            Click_Decision = bb;
            Click_Distance = aa;
        else
            Click_Decision = 0;
            Click_Distance = 0;
        end            
        %else
        %    [ Click_Decision,Click_Distance] = multilayer_perceptron(Neuro.NeuralFeatures);        
        %end
    else        
        if Params.SmoothDataFlag ==1
            [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.FilteredFeatures);
        else
            [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.NeuralFeatures);
        end        
    end
end

end % UpdateMultiClicker