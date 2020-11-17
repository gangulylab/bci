function [decision, distance_from_boundary] = multistate_discrete(X,model,dec_bound)
%
% INPUT:% 
% X - Feature vector (vector of length 896)
% model - SVM tensor contrasting each target vs. others
% dec_bound - whether classification should be based on max-vote or
%             <dec_bound> distance from decision boundary
%
% OUTPUT:
% DECISION - for 1 of the 4 possible targets
% DISTANCE_FROM_BOUNDARY - Distance from the decision boudary (More
% negative values : greater confidence in decision)

% make it a row vector
if size(X,1) ~= 1
    X=X';
end

% ignore delta phase
X = X(129:end);

% only hG 
X = X(641:end);


% evaluate neural features at this data point with all pairwise parallel
% classifiers
d1 = X*squeeze(model(1,:,:))';
d2 = X*squeeze(model(2,:,:))';
d3 = X*squeeze(model(3,:,:))';
d4 = X*squeeze(model(4,:,:))';

% store results
Dec = [d1;d2;d3;d4];
Dec_values=[sum(d1);sum(d2);sum(d3);sum(d4)];
Dec_thresh=sum(Dec'<0);


if dec_bound == 0    
    % make decision on max-vote strategy
    [aa bb]=max(Dec_thresh);
    if length(find(Dec_thresh==aa)) == 1
        decision = bb;
        distance_from_boundary = Dec_values(bb);
    else
        [aa1 bb1]=min(Dec_values);
        decision = bb1;
        distance_from_boundary = Dec_values(bb1);
    end    
else    
    % make decision  based on distance from decision boundary
    Dec_values(Dec_values>dec_bound)=0;
    [aa1 bb1]=min(Dec_values);
    if aa1~=0
        decision = bb1;
        distance_from_boundary = Dec_values(bb1);
    else
        decision = 0;
        distance_from_boundary = 0;
    end    
end

end

