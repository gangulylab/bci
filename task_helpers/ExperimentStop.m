function ExperimentStop(Params, fromPause)
if ~exist('fromPause', 'var'), fromPause = 0; end

% Close Screen
Screen('CloseAll');

% Close Serial Port and audio
fclose('all');

if size(Params.KF.A,1)==5,
    dim = '';
else,
    dim = '_1D';
end

% Update parameters of full kalman filter in persistence folder
f1 = load(fullfile(Params.Persistencedir, sprintf('full_kf_params%s.mat',dim)));
f2 = load(fullfile(Params.Persistencedir, sprintf('kf_params%s.mat',dim)));
% KF = f1.KF; % initialize to full kf, then update params
% KF.P                                            = f2.KF.P;
% KF.R                                            = f2.KF.R;
% KF.ESS                                          = f2.KF.ESS;
% KF.S(Params.FeatureMask,:)                      = f2.KF.S;
% KF.T(Params.FeatureMask,Params.FeatureMask)     = f2.KF.T;
% KF.C(Params.FeatureMask,:)                      = f2.KF.C;
% KF.Q(Params.FeatureMask,Params.FeatureMask)     = f2.KF.Q;
% KF.Tinv(Params.FeatureMask,Params.FeatureMask)  = f2.KF.Tinv;
% KF.Qinv(Params.FeatureMask,Params.FeatureMask)  = f2.KF.Qinv;
KF=[];
save(fullfile(Params.Persistencedir, sprintf('full_kf_params%s.mat',dim)), ...
    'KF', '-v7.3', '-nocompression')

% Save persistence folder with Data
mkdir(fullfile(Params.Datadir,'persistence'))
copyfile(fullfile(Params.Persistencedir, sprintf('kf_params%s.mat',dim)),...
    fullfile(Params.Datadir, 'persistence', sprintf('kf_params%s.mat',dim)));


% quit
fprintf('Ending Experiment\n')
if fromPause, keyboard; end

end % ExperimentStop
