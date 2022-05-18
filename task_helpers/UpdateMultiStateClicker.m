function [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params, Neuro, Clicker,bl_mean,bl_std)
% Multistate Clicker update
% Checks the current decoded target
% If it is the correct target, it updates the Cursor.ClickState

global Cursor

if Params.ControlMode == 2 %mouse
    %     Click_Decision  = Params.clickOrder(Params.index);
    Click_Decision = Params.TargetID;
    Click_Distance = 0;
    
else,
    if Params.NeuralNetFlag == 1
        if Params.ChPooling == 0
            %if Params.SmoothDataFlag==1
            X = Neuro.FilteredFeatures;
            X = X(:);
            X = X(129:end);% all features except delta phase
            %X = X(769:end);% only hG
            idx=[1:128 385:512 641:768];
            X=X(idx);
            
            % 2-norm
            if Params.Norm2
                X = X./norm(X);
            end
            
            if Params.AdaptiveBaseline
                m = bl_mean(idx+128);
                s = bl_std(idx+128);
                X = (X-m)./s;
            end
            
            %Decision_Prob = multilayer_perceptron_Day1to7(X);
            Decision_Prob = feval(Params.NeuralNetFunction,X);
            %disp(Decision_Prob);
            [aa bb]=max(Decision_Prob);
            if aa >= Params.NeuralNetSoftMaxThresh
                Click_Decision = bb;
                Click_Distance = aa;
            else
                Click_Decision = 0;
                Click_Distance = 0;
            end
            
        else
            chmap=Params.ChMap;
            [xx yy] = size(chmap);
            X = Neuro.FilteredFeatures;
            X = X(:);
            feat_idx = [129:256 513:640 769:896];
            X = X(feat_idx);
            
            % pooling
            f1 = (X(1:128));
            f2 = (X(129:256));
            f3 = (X(257:384));
            f1 = f1(chmap);
            f2 = f2(chmap);
            f3 = f3(chmap);
            pooled_data=[];
            for i=1:2:xx
                for j=1:2:yy
                    delta = (f1(i:i+1,j:j+1));delta=mean(delta(:));
                    beta = (f2(i:i+1,j:j+1));beta=mean(beta(:));
                    hg = (f3(i:i+1,j:j+1));hg=mean(hg(:));
                    pooled_data = [pooled_data; delta; beta ;hg];
                end
            end
            X = pooled_data;
            
            % 2-norm
            if Params.Norm2
                X = X./norm(X);
            end
            
            % eval the classifier
            Decision_Prob = feval(Params.NeuralNetFunction,X);
            
            if Params.NeuralBias
                Decision_Prob(Params.NeuralNetBiasDirection) = ...
                    Decision_Prob(Params.NeuralNetBiasDirection) - Params.NeuralNetBiasCorrection;
            end
            
            [aa bb]=max(Decision_Prob);
            if aa >= Params.NeuralNetSoftMaxThresh
                Click_Decision = bb;
                Click_Distance = aa;
            else
                Click_Decision = 0;
                Click_Distance = 0;
            end
            
        end
    elseif Params.NeuralNet2Flag==1
        chmap=Params.ChMap;
        [xx yy] = size(chmap);
        X = Neuro.FilteredFeatures;
        X = X(:);
        feat_idx = [129:256 513:640 769:896];
        X = X(feat_idx);
        % pooling
        f1 = (X(1:128));
        f2 = (X(129:256));
        f3 = (X(257:384));
        f1 = f1(chmap);
        f2 = f2(chmap);
        f3 = f3(chmap);
        pooled_data=[];
        for i=1:2:xx
            for j=1:2:yy
                delta = (f1(i:i+1,j:j+1));delta=mean(delta(:));
                beta = (f2(i:i+1,j:j+1));beta=mean(beta(:));
                hg = (f3(i:i+1,j:j+1));hg=mean(hg(:));
                pooled_data = [pooled_data; delta; beta ;hg];
            end
        end
        X = pooled_data;
        
        % 2-norm
        if Params.Norm2
            X = X./norm(X);
        end
        
        act = predict(Params.NeuralNet2.net, X' ,'ExecutionEnvironment','cpu')';
        
        if Params.NeuralBias
            act(Params.NeuralNetBiasDirection) = ...
                act(Params.NeuralNetBiasDirection) - Params.NeuralNetBiasCorrection;
        end
        
        [aa bb]=max(act);
        if aa >=  Params.NeuralNet2SoftMaxThresh
            Click_Decision = bb;
            Click_Distance = aa;
        else
            Click_Decision = 0;
            Click_Distance = 0;
        end
        
    elseif Params.ConvNeuralNetFlag == 1
        chtemp=[];
        chmap=Params.ChMap;
        X = Neuro.FilteredFeatures;
        X = X(:);
        feat_idx = [129:256 513:640 769:896];
        X = X(feat_idx);
        
        % 2-norm
        if Params.Norm2
            X = X./norm(X);
        end
        
        f1 = (X(1:128));
        f2 = (X(129:256));
        f3 = (X(257:384));
        chtemp(:,:,1) = f1(chmap);
        chtemp(:,:,2) = f2(chmap);
        chtemp(:,:,3) = f3(chmap);
        %act = squeeze(activations(Params.ConvNeuralNet.net,chtemp,20));
        act = predict(Params.ConvNeuralNet.net,chtemp,'ExecutionEnvironment','cpu');
        [aa bb]=max(act);
        if aa<  Params.ConvNeuralNetSoftMaxThresh
            Click_Decision = 0;
            Click_Distance = aa;
        else
            Click_Decision = bb;
            Click_Distance = aa;
        end
        %
        %     elseif Params.UseSVM ==1
        %         [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.FilteredFeatures);
        %disp(Click_Decision)
    elseif Params.NeuralNetEnsemble == 1
        chmap=Params.ChMap;
        [xx yy] = size(chmap);
        X = Neuro.FilteredFeatures;
        X = X(:);
        feat_idx = [129:256 513:640 769:896];
        X = X(feat_idx);
        
        % pooling
        f1 = (X(1:128));
        f2 = (X(129:256));
        f3 = (X(257:384));
        f1 = f1(chmap);
        f2 = f2(chmap);
        f3 = f3(chmap);
        pooled_data=[];
        for i=1:2:xx
            for j=1:2:yy
                delta = (f1(i:i+1,j:j+1));delta=mean(delta(:));
                beta = (f2(i:i+1,j:j+1));beta=mean(beta(:));
                hg = (f3(i:i+1,j:j+1));hg=mean(hg(:));
                pooled_data = [pooled_data; delta; beta ;hg];
            end
        end
        X = pooled_data;
        
        % 2-norm
        if Params.Norm2
            X = X./norm(X);
        end
        
        % eval the classifier
        nets = Params.NeuralNetFunction(1).(Params.NeuralNetName);    
        Decision_Prob=[];
        for iter=1:1%length(nets)
            Decision_Prob = [Decision_Prob nets{iter}(X)];
        end
        %Decision_Prob = mean(Decision_Prob,2);    
        
        
        if Params.NeuralBias
            Decision_Prob(Params.NeuralNetBiasDirection) = ...
                Decision_Prob(Params.NeuralNetBiasDirection) - Params.NeuralNetBiasCorrection;
        end
        
        [aa bb]=max(Decision_Prob);
        if aa >= Params.NeuralNetSoftMaxThresh
            Click_Decision = bb;
            Click_Distance = aa;
        else
            Click_Decision = 0;
            Click_Distance = 0;
        end
    end
    
    %     else
    %         if Params.SmoothDataFlag ==1
    %             [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.FilteredFeatures);
    %         else
    %             [ Click_Decision,Click_Distance] = Clicker.Func(Neuro.NeuralFeatures);
    %         end
end


end % UpdateMultiClicker