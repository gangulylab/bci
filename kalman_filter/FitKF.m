function KF = FitKF(Params,datadir,fitFlag,KF,TrialBatch,dimRedFunc)
% function KF = FitKF(Params,datadir,fitFlag,KF,TrialBatch)
% Uses all trials in given data directory to initialize matrices for kalman
% filter. Returns KF structure containing matrices: A,W,P,C,Q
% 
% datadir - directory containing trials to fit data on
% fitFlag - 0-fit on actual state,
%           1-fit on intended kinematics (refit algorithm)
%           2-fit on intended kinematics (smoothbatch algorithm)
% KF - kalman filter structure containing matrices: A,W,P,C,Q
% TrialBatch - cell array of filenames w/ trials to use in smooth batch
% dimRedFunc - function handle for dimensionality red. redX = dimRedFunc(X)

% Initialization of KF
if ~exist('KF','var'),
    KF = Params.KF;
    KF.Lambda = Params.CLDA.Lambda;
end

% If Initialization Mode = 3, manually choose datadir & fit KF
if KF.InitializationMode==3 && fitFlag==0,
    datadir = uigetdir(datadir,'Select Directory to Initialize Kalman Filter');
end

% If Initialization Mode = 4, load kf params from persistence folder
if KF.InitializationMode==4 && fitFlag==0,
    f=load(fullfile(Params.Persistencedir, 'kf_params.mat'));

    % TRY UNCOMMENTING AND RUNNING THIS BLOCK 

    %%% block start
%     KF.Lambda = Params.CLDA.Lambda;
%     KF.P = [];
%     KF.R  = [];
%     KF.ESS = [];
%     KF.S = [];
%     KF.T = [];
%     KF.C = [];
%     KF.Q = [];
%     KF.Tinv = [];
%     KF.Qinv = [];
%%% block END

    % ELSE UNCOMMENT AND RUN THIS BLOCK
    %%% block start

    %if all(Params.FeatureMask == f.FeatureMask), % load reduced KF
        KF.Lambda = Params.CLDA.Lambda;
        KF.P = f.KF.P;
        KF.R = f.KF.R;
        KF.ESS = f.KF.ESS;
        KF.S = f.KF.S;
        KF.T = f.KF.T;
        KF.C = f.KF.C;
        KF.Q = f.KF.Q;
        KF.Tinv = f.KF.Tinv;
        KF.Qinv = f.KF.Qinv;
    %%% block END
      
    %else, % load full KF
%         f=load(fullfile(Params.Persistencedir, 'full_kf_params.mat'));
%         KF.Lambda = Params.CLDA.Lambda;
%         KF.P = f.KF.P;
%         KF.R = f.KF.R;
%         KF.ESS = f.KF.ESS;
%         KF.S = f.KF.S(Params.FeatureMask,:);
%         KF.T = f.KF.T(Params.FeatureMask,Params.FeatureMask);
%         KF.C = f.KF.C(Params.FeatureMask,:);
%         KF.Q = f.KF.Q(Params.FeatureMask,Params.FeatureMask);
%         KF.Tinv = f.KF.Tinv(Params.FeatureMask,Params.FeatureMask);
%         KF.Qinv = f.KF.Qinv(Params.FeatureMask,Params.FeatureMask);
    %end
    fprintf('\n\nLoading Previous Kalman Filter:\n')
    return
end

% ouput to screen
fprintf('\n\nFitting Kalman Filter:\n')
switch fitFlag,
    case 0,
        fprintf('  Initial Fit\n')
        switch KF.InitializationMode,
            case {1,2,3},
                fprintf('  Data in %s\n', datadir)
            case 4,
                fprintf('  Using most recent KF.\n')
        end
    case 1,
        fprintf('  ReFit\n')
        fprintf('  Data in %s\n', datadir)
    case 2,
        fprintf('  Smooth Batch\n')
        fprintf('  Data in %s\n', datadir)
        fprintf('  Trials: {%s-%s}\n', TrialBatch{1},TrialBatch{end})
end

% grab data trial data
datafiles = dir(fullfile(datadir,'Data*.mat'));
if fitFlag==2, % if smooth batch, only use files TrialBatch
    names = {datafiles.name};
    idx = zeros(1,length(names))==1;
    for i=1:length(TrialBatch),
        idx = idx | strcmp(names,TrialBatch{i});
    end
    datafiles = datafiles(idx);
end

Tfull = [];
Xfull = [];
Y = [];
T = [];
for i=1:length(datafiles),
    % load data
    load(fullfile(datadir,datafiles(i).name)) %#ok<LOAD>
    % ignore inter-trial interval data
    if strcmp(TrialData.Events(1).Str, 'Inter Trial Interval'),
        tidx = TrialData.Time >= TrialData.Events(2).Time;
    else,
        tidx = TrialData.Time >= TrialData.Events(1).Time;
    end
    % grab cursor pos and time
    Tfull = cat(2,Tfull,TrialData.Time(tidx));
    if fitFlag==0, % fit on true kinematics
        Xfull = cat(2,Xfull,TrialData.CursorState(:,tidx));
    else, % refit on intended kinematics
        Xfull = cat(2,Xfull,TrialData.IntendedCursorState(:,tidx));
    end
    T = cat(2,T,TrialData.NeuralTime(tidx));
    Y = cat(2,Y,TrialData.NeuralFeatures{tidx});
end

% interpolate to get cursor pos and vel at neural times
if size(Xfull,2)>size(Y,2)
    X = interp1(Tfull',Xfull',T')';
else,
    X = Xfull;
end

% ignore "bad features"
Y = Y(Params.FeatureMask,:);

% if DimRed is on, reduce dimensionality of neural features
if exist('dimRedFunc','var'),
    Y = dimRedFunc(Y);
end

% full cursor state at neural times
D = size(X,2);

% if initialization mode returns shuffled weights
if fitFlag==0 && KF.InitializationMode==2, % return shuffled weights
    fprintf('  *Shuffled Weights\n')
    idx = randperm(size(Y,2));
    Y = Y(:,idx);
end

% fit kalman matrices
if KF.VelKF, % only use vel to fit C, set pos terms to 0
    C = (Y*X(3:end,:)') / (X(3:end,:)*X(3:end,:)');
    C = [zeros(size(C,1),2),C];
else,
    C = (Y*X') / (X*X');
end
Q = (1/D) * ((Y-C*X) * (Y-C*X)');

% update kalman matrices
switch fitFlag,
    case {0,1},
        % fit sufficient stats
        if KF.VelKF, % only use vel to fit C, set pos terms to 0
            X = X(3:end,:);
        end
        KF.R = X*X';
        KF.S = Y*X';
        KF.T = Y*Y';
        KF.ESS = D;
        KF.C = C;
        KF.Q = Q;
        KF.Tinv = inv(KF.T);
        KF.Qinv = inv(Q);
    case 2, % smooth batch
        try
        alpha = Params.CLDA.Alpha;
        KF.C = alpha*KF.C + (1-alpha)*C;
        KF.Q = alpha*KF.Q + (1-alpha)*Q;
        KF.Qinv = inv(KF.Q);
        catch
        end
end

end % FitKF
