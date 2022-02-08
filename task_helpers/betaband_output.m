function beta_scalar = betaband_output(Params,Neuro)


% get the beta band features
X = Neuro.FilteredFeatures;
X = X(385:512);

% get PC wts
beta_wts = Params.BetaWts.betawts;
beta_mean = Params.BetaMean.betamean;

% get projected value
X = X - beta_mean;
beta_scalar = X'*beta_wts;

% check four boundaries
if beta_scalar < 0
    beta_scalar = 0;
elseif beta_scalar>Params.BetaBarValue
    beta_scalar = Params.BetaBarValue;
end


end