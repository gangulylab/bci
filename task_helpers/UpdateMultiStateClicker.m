function [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params, Neuro, Clicker)
% Multistate Clicker update
% Checks the current decoded target
% If it is the correct target, it updates the Cursor.ClickState

global Cursor

if Params.ControlMode == 2 %mouse
    Click_Decision = randperm(5,1)-1;
    Click_Distance = 0;
else
    if Params.NeuralNetFlag == 1
        %if Params.SmoothDataFlag==1
        X = Neuro.FilteredFeatures;
        X = X(:);

         % get delta, beta and hG removing bad channels 
        X = X([257:512 1025:1280 1537:1792]);        
        bad_ch = [108 113 118];
        good_ch = ones(length(X),1);
        for ii=1:length(bad_ch)
            bad_ch_tmp = bad_ch(ii)*[1 2 3];
            good_ch(bad_ch_tmp)=0;
        end
        X = X(logical(good_ch));

        %Decision_Prob = multilayer_perceptron_Day1to7(X);
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
        
    elseif Params.ConvNeuralNetFlag == 1
        chtemp=[];
        chmap=Params.ChMapB2;
        X = Neuro.FilteredFeatures;
        X = X(:);
        feat_idx = [129:256 257:384 385:512 513:640 641:768 769:896];
        X = X(feat_idx);
        f1 = (X(1:128));
        f2 = (X(129:256));
        f3 = (X(257:384));
        f4 = (X(385:512));
        f5 = (X(513:640));
        f6 = (X(641:768));
        chtemp(:,:,1) = f1(chmap);
        chtemp(:,:,2) = f2(chmap);
        chtemp(:,:,3) = f3(chmap);
        chtemp(:,:,4) = f4(chmap);
        chtemp(:,:,5) = f5(chmap);
        chtemp(:,:,6) = f6(chmap);
        %act = squeeze(activations(Params.ConvNeuralNet.net,chtemp,12),'ExecutionEnvironment','CPU')
        act = predict(Params.ConvNeuralNet.net,chtemp,'ExecutionEnvironment','cpu');
        [aa bb]=max(act);
        if aa<  Params.ConvNeuralNetSoftMaxThresh
            Click_Decision = 0;
            Click_Distance = aa;
        else
            Click_Decision = bb;
            Click_Distance = aa;
        end
    %else
%         if Params.SmoothDataFlag ==1
%             [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.FilteredFeatures);
%         else
%             [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.NeuralFeatures);
%         end
    end
end

end % UpdateMultiClicker