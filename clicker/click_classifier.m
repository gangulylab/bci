function [decision, distance_from_boundary] = click_classifier(X,model)
%% function [decision, distance_from_boundary] = click_classifier(X);
%
% INPUT: 
% model - model parameters, includes .w (model weights)
% X - Feature vector (vector of length 896)
%
% OUTPUT: 
% DECISION - +1 for click and -1 for no click
% DISTANCE_FROM_BOUNDARY - How far given data point is from the classifier.
% More negative values : greater confidence in no click. More positive
% values : greater confidence in click.

if size(X,1) ~= 1
    X=X';
end

distance_from_boundary = X(129:end)*model.w'; % ignoring delta phase info
if distance_from_boundary>-0.2, % not clicking
    decision = 1;
else, % clicking
    decision = -1;
end

% distance_from_boundary




