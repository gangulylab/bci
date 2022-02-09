function beta_scalar = betaband_output(Params,Neuro)


% get the beta band features
if Params.ControlMode == 2
    X = randn(896,1);
else
    X = Neuro.FilteredFeatures;
end

X = X(385:512);


% get PC wts
beta_wts = Params.BetaWts.betawts;
beta_mean = Params.BetaMean.betamean;

% get projected value
if Params.ControlMode == 2
    beta_scalar =  10 + 2*randn(1);
else    
    X = X - beta_mean;
    beta_scalar = X'*beta_wts;
end

% check four boundaries
if beta_scalar < 0
    beta_scalar = 0;
elseif beta_scalar>Params.BetaBarValue
    beta_scalar = Params.BetaBarValue;
end


end