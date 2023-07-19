function [Click_Decision,Click_Distance] = UpdateMultiStateClicker(Params, Neuro, Clicker,bl_mean,bl_std)
% Multistate Clicker update
% Checks the current decoded target
% If it is the correct target, it updates the Cursor.ClickState

global Cursor

if Params.ControlMode == 2 %mouse
%     Click_Decision  = Params.clickOrder(Params.index);
%         Click_Decision = Params.TargetID;
%      Click_Decision = 1;

    Click_Distance = 0;
    Click_Decision = readNumPad();

else,
    if Params.biLSTMFlag == 1
        if Params.LSTM_Output_Method
            pred = activations(Params.LSTM,Neuro.LSTMFeatures,'fc_2','ExecutionEnvironment','cpu');
            if size(pred,1) > size(pred,2)
                pred=pred';
            end
            X = [pred;Params.lstm_output_pattern];
            R = corrcoef(X');
            pred = R(1,2:end);
            [aa bb]=max(pred);
            if aa >= Params.LSTM_Output_Method_Thresh
                Click_Decision = bb;
                Click_Distance = aa;
            else
                Click_Decision = 0;
                Click_Distance = 0;
            end

        else
            pred =  predict(Params.LSTM,Neuro.LSTMFeatures,'ExecutionEnvironment','cpu');

            %         pred = pred +  [0.1, 0, -0.3, 0, 0, 0, 0];
            [aa bb]=max(pred);
            if aa >=  Params.biLSTMSoftMaxThresh
                Click_Decision = bb;
                Click_Distance = aa;
            else
                Click_Decision = 0;
                Click_Distance = 0;
            end
        end


        %disp(['LSTM o/p '  num2str(Click_Decision)])
    end


    if Params.NeuralNetFlag == 1
        if Params.ChPooling == 0
            %if Params.SmoothDataFlag==1
            X = Neuro.FilteredFeatures;
            X = X(:);
            X = X(129:end);% all features except delta phase
            %X = X(769:end);% only hG
            if Params.Use3Features
                idx=[1:128 385:512 641:768];
            elseif Params.Use4Features
                idx=[1:128 385:512 513:640 641:768];
            end
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
            if Params.Use3Features
                feat_idx = [129:256 513:640 769:896];
                X = X(feat_idx);

                if Params.ChPooling
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
                end

            elseif Params.Use4Features
                feat_idx = [129:256 513:640 641:768 769:896];
                X = X(feat_idx);

                if Params.ChPooling
                    % pooling
                    f1 = (X(1:128));
                    f2 = (X(129:256));
                    f3 = (X(257:384));
                    f4 = (X(385:512));

                    f1 = f1(chmap);
                    f2 = f2(chmap);
                    f3 = f3(chmap);
                    f4 = f4(chmap);
                    pooled_data=[];
                    for i=1:2:xx
                        for j=1:2:yy
                            delta = (f1(i:i+1,j:j+1));delta=mean(delta(:));
                            beta = (f2(i:i+1,j:j+1));beta=mean(beta(:));
                            lg = (f3(i:i+1,j:j+1));lg=mean(lg(:));
                            hg = (f4(i:i+1,j:j+1));hg=mean(hg(:));
                            pooled_data = [pooled_data; delta; beta ;lg;hg];
                        end
                    end
                    X = pooled_data;
                end
            end

            % 2-norm
            if Params.Norm2
                X = X./norm(X);
            end

            % eval the classifier
            Decision_Prob = feval(Params.NeuralNetFunction,X);
            %Decision_Prob = Params.NeuralNet(X);

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