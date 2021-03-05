function [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params, Neuro, Clicker)
% Multistate Clicker update
% Checks the current decoded target
% If it is the correct target, it updates the Cursor.ClickState

global Cursor

if Params.ControlMode == 2 %mouse
    val = 0;
    
    val = double(get(gcf,'CurrentCharacter'));
    pause(0.1)
    if val == 49
        Click_Decision = 6;
    elseif val == 50
        Click_Decision = 2;
    elseif val == 52
        Click_Decision = 3;
    elseif val == 54
        Click_Decision = 1;
    elseif val == 55
        Click_Decision = 5;
    elseif val == 56
        Click_Decision = 4;
    else
        Click_Decision = 0;
    end
%     if strcmp(Params.Task, 'RobotDistance')
%         p = randi(2);
%         if p == 1
%             Click_Decision = randi(5)-1;
%         else
%                 Click_Decision = Params.TargetID;
%         end
%     else
% %     Click_Decision = randperm(6,1)-1;
%         p = randi(4);
%         if p < 2
% %             Click_Decision = 7;
%             Click_Decision = randperm(7,1)-1;
%         else
%             if (Params.TargetID < 7)
%                 Click_Decision = Params.TargetID
%             else
%                 Click_Decision = randperm(8,1)-1;
%             end
%         end
%     end
% %     Click_Decision = 7;
    Click_Distance = 0;
else,
    if Params.NeuralNetFlag == 1
        %if Params.SmoothDataFlag==1
        X = Neuro.FilteredFeatures;
        X = X(:);
        if Params.Use3Features
            idx = [1:128 385:512 641:768]+128; 
            X = X(idx);
        else
            X = X(769:end);
        end
        %Decision_Prob = multistate_discrete_6Target(X);
        Decision_Prob = feval(Params.NeuralNetFunction,X);
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