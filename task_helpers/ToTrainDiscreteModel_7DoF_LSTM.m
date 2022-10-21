%% FINE TUNING LSTM MODEL FOR A BATCH UPDATE ON ARROW DATA

clc;clear
addpath '/home/ucsf/Projects/bci/task_helpers/'
addpath '/home/ucsf/Projects/bci/clicker/'

root_path = '/home/ucsf/Data/bravo1/';
foldernames = {'20221021'};

% filter bank hg
Params=[];
Params.Fs = 1000;
Params.FilterBank(1).fpass = [70,77];   % high gamma1
Params.FilterBank(end+1).fpass = [77,85];   % high gamma2
Params.FilterBank(end+1).fpass = [85,93];   % high gamma3
Params.FilterBank(end+1).fpass = [93,102];  % high gamma4
Params.FilterBank(end+1).fpass = [102,113]; % high gamma5
Params.FilterBank(end+1).fpass = [113,124]; % high gamma6
Params.FilterBank(end+1).fpass = [124,136]; % high gamma7
Params.FilterBank(end+1).fpass = [136,150]; % high gamma8
Params.FilterBank(end+1).fpass = [30,36]; % lg1
Params.FilterBank(end+1).fpass = [36,42]; % lg2
Params.FilterBank(end+1).fpass = [42,50]; % lg3
% compute filter coefficients
for i=1:length(Params.FilterBank),
    [b,a] = butter(3,Params.FilterBank(i).fpass/(Params.Fs/2));
    Params.FilterBank(i).b = b;
    Params.FilterBank(i).a = a;
end

% low pass filters
lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
    'PassbandFrequency',25,'PassbandRipple',0.2, ...
    'SampleRate',1e3);

% get all the folders
filepath = fullfile(root_path,foldernames{1},'Robot3DArrow');
folders = dir(filepath);
folders=folders(3:end-3);
%folders=folders(3:8);

% load the decoder.. use the decoder name run in the experiment 
load net_bilstm_20220824
net_bilstm = net_bilstm_20220824;


%%%%% get the training data
folders_train = folders;

% get the files
files_train=[];
for j=1:length(folders_train)
    subfolder = fullfile(folders_train(j).folder,folders_train(j).name,'BCI_Fixed');
    tmp = findfiles('mat',subfolder,1)';
    files_train =[files_train;tmp];
end

% get inidividual samples
[XTrain,XTest,YTrain,YTest] = get_lstm_features(files_train,Params,lpFilt);

% update the decoder
batch_size=128;
val_freq = floor(length(XTrain)/batch_size);
learning_rate = 1e-4;
options = trainingOptions('adam', ...
    'MaxEpochs',50, ...
    'MiniBatchSize',batch_size, ...
    'GradientThreshold',10, ...
    'Verbose',true, ...
    'ValidationFrequency',val_freq,...
    'Shuffle','every-epoch', ...
    'ValidationData',{XTest,YTest},...
    'ValidationPatience',5,...
    'Plots','training-progress',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.1,...    
    'LearnRateDropPeriod',30,...
    'InitialLearnRate',learning_rate,...
    'CheckpointPath','/home/ucsf/Projects/bci/lstm_models');

% train the model
% IMPORTANT -> VALIDATION SHOULD NOT BE WILDLY DIFFERENT THAT TRAINING
% LOSS. IF SO, REDUCE LEARNING RATE (2E-4 TO 1E-4 ETC.) AND RE RUN. IF IT
% CONVERGES TOO QUICKLY (<10 EPOCHS) RE RUN AGAIN AS IT IS STUCK IN LOCAL
% MINIMA
clear net
layers = net_bilstm.Layers;
net = trainNetwork(XTrain,YTrain,layers,options);


net_bilstm_20220824_update = net;
save net_bilstm_20220824_update net_bilstm_20220824_update

tmp=load('net_checkpoint__252__2022_10_21__14_58_38.mat')
net_bilstm_20220824_update = tmp.net;
save net_bilstm_20220824_update net_bilstm_20220824_update

%% FINE TUNING LSTM MODEL FOR A BATCH UPDATE ON ROBOT CENTER OUT DATA


clc;clear
addpath '/home/ucsf/Projects/bci/task_helpers/'
addpath '/home/ucsf/Projects/bci/clicker/'

root_path = '/home/ucsf/Data/bravo1/';
foldernames = {'20220902'};

% filter bank hg
Params=[];
Params.Fs = 1000;
Params.FilterBank(1).fpass = [70,77];   % high gamma1
Params.FilterBank(end+1).fpass = [77,85];   % high gamma2
Params.FilterBank(end+1).fpass = [85,93];   % high gamma3
Params.FilterBank(end+1).fpass = [93,102];  % high gamma4
Params.FilterBank(end+1).fpass = [102,113]; % high gamma5
Params.FilterBank(end+1).fpass = [113,124]; % high gamma6
Params.FilterBank(end+1).fpass = [124,136]; % high gamma7
Params.FilterBank(end+1).fpass = [136,150]; % high gamma8
Params.FilterBank(end+1).fpass = [30,36]; % lg1
Params.FilterBank(end+1).fpass = [36,42]; % lg2
Params.FilterBank(end+1).fpass = [42,50]; % lg3
% compute filter coefficients
for i=1:length(Params.FilterBank),
    [b,a] = butter(3,Params.FilterBank(i).fpass/(Params.Fs/2));
    Params.FilterBank(i).b = b;
    Params.FilterBank(i).a = a;
end

% low pass filters
lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
    'PassbandFrequency',25,'PassbandRipple',0.2, ...
    'SampleRate',1e3);

% get all the folders
filepath = fullfile(root_path,foldernames{1},'RealRobotBatch');
folders = dir(filepath);
folders=folders(3:end);
%folders=folders(3:8);

% load the decoder.. use the decoder name run in the experiment 
load net_bilstm_robot_20220824
net_bilstm = net_bilstm_robot_20220824;


%%%%% get the training data
folders_train = folders;

% get the files
files_train=[];
for j=1:length(folders_train)
    subfolder = fullfile(folders_train(j).folder,folders_train(j).name,'BCI_Fixed');
    tmp = findfiles('mat',subfolder,1)';
    files_train =[files_train;tmp];
end

% get inidividual samples
[XTrain,XTest,YTrain,YTest] = get_lstm_features_robotBatch(files_train,Params,lpFilt);

% update the decoder
batch_size=128;
val_freq = floor(length(XTrain)/batch_size);
learning_rate = 2e-4; % important parameter to tune
options = trainingOptions('adam', ...
    'MaxEpochs',50, ...
    'MiniBatchSize',batch_size, ...
    'GradientThreshold',10, ...
    'Verbose',true, ...
    'ValidationFrequency',val_freq,...
    'Shuffle','every-epoch', ...
    'ValidationData',{XTest,YTest},...
    'ValidationPatience',6,...
    'Plots','training-progress',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.1,...
    'OutputNetwork','best-validation-loss',...
    'LearnRateDropPeriod',30,...
    'InitialLearnRate',learning_rate);

% train the model
% IMPORTANT -> VALIDATION SHOULD NOT BE WILDLY DIFFERENT THAT TRAINING
% LOSS. IF SO, REDUCE LEARNING RATE (2E-4 TO 1E-4 ETC.) AND RE RUN. IF IT
% CONVERGES TOO QUICKLY (<10 EPOCHS) RE RUN AGAIN AS IT IS STUCK IN LOCAL
% MINIMA
clear net
layers = net_bilstm.Layers;
net = trainNetwork(XTrain,YTrain,layers,options);
net_bilstm_robot_20220824_update = net;
save net_bilstm_robot_20220824_update net_bilstm_robot_20220824_update



