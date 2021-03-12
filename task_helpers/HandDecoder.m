function [ClickDecision] = HandDecoder(Data,Params)
%function [ClickDecision] = HandDecoder(Data,Params)

kinax = Data.TaskState;
idx = find(kinax==3);

if length(idx)<17
    l = 17 - length(idx);
    idx = [idx(1)-l:idx idx];    
end

features = Data.SmoothedNeuralFeatures(idx);
X = cell2mat(features);
X = X(129:end,:);
decision_distance = Params.HandOpenCloseDecoder*X(:);

if decision_distance<=0
    ClickDecision = 7;
else
    ClickDecision = 8;
end


end