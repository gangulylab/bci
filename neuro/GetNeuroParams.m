function Params = GetNeuroParams(Params)
% Params = GetNeuroParams(Params)
% specifically sets params related to the neural processing pipeline

Params.ZBufSize = 120; % secs
Params.ZscoreRawFlag = true;
Params.UpdateChStatsFlag = false;
Params.ZscoreFeaturesFlag = true;
Params.UpdateFeatureStatsFlag = false;
Params.SaveRaw = true;
Params.SaveProcessed = false;

Params.DimRed.Flag = false;
Params.DimRed.InitMode = 2; % 1-use imagined mvmts, 2-choose dir
Params.DimRed.InitAdapt = true;
Params.DimRed.InitFixed = ~Params.DimRed.InitAdapt;
Params.DimRed.Method = 1; % 1-pca, 2-fa
Params.DimRed.AvgTrialsFlag = false; % 0-cat imagined mvmts, 1-avg imagined mvmts
Params.DimRed.NumDims = 500;

Params.Fs = 1000;
Params.NumChannels = 128;
Params.NumFeatureBins = 1;
Params.BufferTime = 2; % secs longer for better phase estimation of low frqs
Params.BufferSamps = Params.BufferTime * Params.Fs;
RefModeStr = {'none','common_mean','common_median'};
Params.ReferenceMode = 2; % 0-no ref, 1-common mean, 2-common median
Params.ReferenceModeStr = RefModeStr{Params.ReferenceMode+1};

Params.ChMapFile = 'ECOG_Grid_8596-002135';

% filter bank - each element is a filter bank
% fpass - bandpass cutoff freqs
% feature - # of feature (can have multiple filters for a single feature
% eg., high gamma is composed of multiple freqs)
Params.FilterBank = [];
Params.FilterBank(end+1).fpass = [.5,4];    % delta
Params.FilterBank(end).buffer_flag = true;
Params.FilterBank(end).hilbert_flag = true;
Params.FilterBank(end).phase_flag = true;
Params.FilterBank(end).feature = 2;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [4,8];     % theta
Params.FilterBank(end).buffer_flag = true;
Params.FilterBank(end).hilbert_flag = true;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 3;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [8,13];    % alpha
Params.FilterBank(end).buffer_flag = true;
Params.FilterBank(end).hilbert_flag = true;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 4;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [13,19];   % beta1
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 5;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [19,30];   % beta2
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 5;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [30,36];   % low gamma1
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 6;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [36,42];   % low gamma2
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 6;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [42,50];   % low gamma3
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 6;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [70,77];   % high gamma1
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 7;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [77,85];   % high gamma2
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 7;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [85,93];   % high gamma3
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 7;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [93,102];  % high gamma4
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 7;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [102,113]; % high gamma5
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 7;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [113,124]; % high gamma6
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 7;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [124,136]; % high gamma7
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 7;
Params.FilterBank(end).spatial_filt_sz = 3;

Params.FilterBank(end+1).fpass = [136,150]; % high gamma8
Params.FilterBank(end).buffer_flag = false;
Params.FilterBank(end).hilbert_flag = false;
Params.FilterBank(end).phase_flag = false;
Params.FilterBank(end).feature = 7;
Params.FilterBank(end).spatial_filt_sz = 3;

% compute filter coefficients
for i=1:length(Params.FilterBank),
    [b,a] = butter(3,Params.FilterBank(i).fpass/(Params.Fs/2));
    Params.FilterBank(i).b = b;
    Params.FilterBank(i).a = a;
end


% get lpfilt
lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
    'PassbandFrequency',25,'PassbandRipple',0.2, ...
    'SampleRate',1e3);
Params.lpFilt = lpFilt;

% unique pwr feature + all phase features
Params.NumBuffer = sum([Params.FilterBank.buffer_flag]);
Params.NumHilbert = sum([Params.FilterBank.hilbert_flag]);
Params.NumPhase = sum([Params.FilterBank.phase_flag]);
Params.NumPower = length(unique([Params.FilterBank.feature]));
Params.NumFeatures = Params.NumPower + Params.NumPhase;

% make sure spatial filter sizes are consistent within features
features = [Params.FilterBank.feature];
for feature=Params.NumPhase+1:Params.NumFeatures,
    idx = features == feature;
    spatial_filt_szs = [Params.FilterBank(idx).spatial_filt_sz];
    if length(spatial_filt_szs)>1,
        assert(all(spatial_filt_szs==spatial_filt_szs(1)),...
            'Spatial Filter Sizes are Inconsistent')
    end
end

end % GetNeuroParams