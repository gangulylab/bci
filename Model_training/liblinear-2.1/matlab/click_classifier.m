function [decision distance_from_boundary] = click_classifier(X)
% 
%
% function [decision distance_from_boundary] = click_classifier(X)
%
% INPUT: 
% X - Feature vector (vector of length 896)
%
% OUTPUT: 
% DECISION - -1 for click and +1 for no click
% DISTANCE_FROM_BOUNDARY - How far given data point is from the classifier.
% More negative values : greater confidence in click. More positive
% values : greater confidence in no click.

%load clicker_svm_mdl
load clicker_svm_mdl_819_OKOnly

if size(X,1) ~= 1
    X=X';
end

distance_from_boundary = X(129:end)*model.w'; % ignoring delta phase info
if distance_from_boundary <= -0.5
    decision = -1;
else
    decision =1;
end

% if distance_from_boundary>0
%     decision = 1;
% else
%     decision = -1;
% end






