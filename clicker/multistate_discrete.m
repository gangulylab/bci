function [decision, distance_from_boundary] = multistate_discrete(X,model,dec_bound)
%
% INPUT: 
% model - model weights
% X - Feature vector (vector of length 896)
%
% OUTPUT: 
% DECISION - for 1 of the 4 possible targets
% DISTANCE_FROM_BOUNDARY - Mean distance away from the classifier and
% towards the decoded class. 
% More negative values : greater confidence in decision

% make it a row vector
if size(X,1) ~= 1
    X=X';
end

% ignore delta phase
X = X(129:end);

% evaluate neural features at this data point with all pairwise parallel
% classifiers
d1 = X*squeeze(model(1,:,:))';
d2 = X*squeeze(model(2,:,:))';
d3 = X*squeeze(model(3,:,:))';
d4 = X*squeeze(model(4,:,:))';
%d5 = X*squeeze(model(5,:,:))';
%d6 = X*squeeze(model(6,:,:))';
%d7 = X*squeeze(model(7,:,:))';
%d8 = X*squeeze(model(8,:,:))';

% store results
Dec = [d1;d2;d3;d4];
Dec_values=[sum(d1);sum(d2);sum(d3);sum(d4)];
Dec_thresh=sum(Dec'<dec_bound);


% decision based on distance 
% 
% [aa1 bb1]=min(Dec_values);
% decision = bb1;
% distance_from_boundary = Dec_values(bb1);




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


end

