%% FINE TUNING LSTM MODEL FOR A BATCH UPDATE ON ARROW DATA

clc;clear
addpath '/home/ucsf/Projects/bci/task_helpers/'
addpath '/home/ucsf/Projects/bci/clicker/'
addpath '/home/ucsf/Projects/bci/lstm_models/'

root_path = '/home/ucsf/Data/bravo1/';
foldernames = {'20231020'};
lstm_folder_path = '/home/ucsf/Projects/bci/lstm_models/';
clicker_path = '/home/ucsf/Projects/bci/clicker/';

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
folders=folders(3:end);
%folders=folders(3:8);

% load the decoder.. use the decoder name run in the experiment 
load net_bilstm_9DoF
net_bilstm = net_bilstm_9DoF;


%%%%% get the training data
folders_train = folders;

% get the files
files_train=[];
for j=[1:2, 8:11]%length(folders_train)
    subfolder = fullfile(folders_train(j).folder,folders_train(j).name,'BCI_Fixed');
    tmp = findfiles('mat',subfolder,1)';
    files_train =[files_train;tmp];
end

for j=[3:7]%length(folders_train)
    subfolder = fullfile(folders_train(j).folder,folders_train(j).name,'Imagined');
    tmp = findfiles('mat',subfolder,1)';
    files_train =[files_train;tmp];
end



% get inidividual samples
[XTrain,XTest,YTrain,YTest] = get_lstm_features_9DoF(files_train,Params,lpFilt);

% decoder fine tune params 
batch_size=32;
val_freq = floor(length(XTrain)/batch_size);
learning_rate = 1e-4;
patience = 6;
tmp_time = clock;
saved_folder=[];
for i=4:6
    saved_folder =[saved_folder num2str(round(tmp_time(i)))];
end
todays_date = [num2str(tmp_time(1)) num2str(tmp_time(2)) num2str(tmp_time(3))];
check_pt_foldername = fullfile(lstm_folder_path,todays_date,saved_folder);
if ~exist(check_pt_foldername)
    mkdir(check_pt_foldername)
end

options = trainingOptions('adam', ...
    'MaxEpochs',60, ...
    'MiniBatchSize',batch_size, ...
    'GradientThreshold',10, ...
    'Verbose',true, ...
    'ValidationFrequency',val_freq,...
    'Shuffle','every-epoch', ...
    'ValidationData',{XTest,YTest},...
    'ValidationPatience',patience,...
    'Plots','training-progress',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.1,...    
    'LearnRateDropPeriod',40,...
    'InitialLearnRate',learning_rate,...
    'CheckpointPath',check_pt_foldername);


% get rid of batch normalization layer as it throws an error due to matlab
% version issues
clear net
layers = net_bilstm.Layers;
idx=[];
for i = 1:length(layers)
    if ~strcmp(layers(i).Name,'batchnorm')
        idx=[idx i];        
    end
end
layers = layers(idx);

% train the model
% IMPORTANT -> VALIDATION SHOULD NOT BE WILDLY DIFFERENT THAT TRAINING
% LOSS. IF SO, REDUCE LEARNING RATE (2E-4 TO 1E-4 ETC.) AND RE RUN. IF IT
% CONVERGES TOO QUICKLY (<10 EPOCHS) RE RUN AGAIN AS IT IS STUCK IN LOCAL
% MINIMA

cd(check_pt_foldername) % if you want to  see the models being saved at each epoch
net = trainNetwork(XTrain,YTrain,layers,options);

% loading the best performing (goat) model 
saved_nets = dir(check_pt_foldername);
saved_nets = saved_nets(3:end);
[~,idx] = sort([saved_nets.datenum]);
saved_nets = saved_nets(idx);
goat_model = saved_nets(end-patience).name;
goat_model = load(fullfile(check_pt_foldername,goat_model));

% saving to clicker folder
cd(clicker_path)
net_bilstm_9DoF_update_20231020_3 = goat_model.net;
save net_bilstm_9DoF_update_20231020_3 net_bilstm_9DoF_update_20231020_3


% %% FINE TUNING LSTM MODEL FOR A BATCH UPDATE ON ROBOT CENTER OUT DATA
% 
% 
% clc;clear
% addpath '/home/ucsf/Projects/bci/task_helpers/'
% addpath '/home/ucsf/Projects/bci/clicker/'
% addpath '/home/ucsf/Projects/bci/lstm_models/'
% 
% root_path = '/home/ucsf/Data/bravo1/';
% foldernames = {'20230526'};
% lstm_folder_path = '/home/ucsf/Projects/bci/lstm_models/';
% clicker_path = '/home/ucsf/Projects/bci/clicker/';
% 
% % filter bank hg
% Params=[];
% Params.Fs = 1000;
% Params.FilterBank(1).fpass = [70,77];   % high gamma1
% Params.FilterBank(end+1).fpass = [77,85];   % high gamma2
% Params.FilterBank(end+1).fpass = [85,93];   % high gamma3
% Params.FilterBank(end+1).fpass = [93,102];  % high gamma4
% Params.FilterBank(end+1).fpass = [102,113]; % high gamma5
% Params.FilterBank(end+1).fpass = [113,124]; % high gamma6
% Params.FilterBank(end+1).fpass = [124,136]; % high gamma7
% Params.FilterBank(end+1).fpass = [136,150]; % high gamma8
% Params.FilterBank(end+1).fpass = [30,36]; % lg1
% Params.FilterBank(end+1).fpass = [36,42]; % lg2
% Params.FilterBank(end+1).fpass = [42,50]; % lg3
% % compute filter coefficients
% for i=1:length(Params.FilterBank),
%     [b,a] = butter(3,Params.FilterBank(i).fpass/(Params.Fs/2));
%     Params.FilterBank(i).b = b;
%     Params.FilterBank(i).a = a;
% end
% 
% % low pass filters
% lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
%     'PassbandFrequency',25,'PassbandRipple',0.2, ...
%     'SampleRate',1e3);
% 
% % get all the folders
% filepath = fullfile(root_path,foldernames{1},'RealRobotBatch');
% folders = dir(filepath);
% folders=folders(3:end);
% 
% 
% % load the decoder.. use the decoder name run in the experiment 
% load net_bilstm_robot_20220824
% net_bilstm = net_bilstm_robot_20220824;
% 
% 
% %%%%% get the training data
% folders_train = folders;
% 
% % get the files
% files_train=[];
% for j=1:length(folders_train)
%     subfolder = fullfile(folders_train(j).folder,folders_train(j).name,'BCI_Fixed');
%     tmp = findfiles('mat',subfolder,1)';
%     files_train =[files_train;tmp];
% end
% 
% % get inidividual samples
%  [XTrain,XTest,YTrain,YTest] = get_lstm_features_robotBatch(files_train,Params,lpFilt);
% %[XTrain,XTest,YTrain,YTest] = get_lstm_features(files_train,Params,lpFilt);
% 
% % decoder fine tune params 
% batch_size=32;
% val_freq = floor(length(XTrain)/batch_size);
% learning_rate = 1e-4;
% patience = 5;
% tmp_time = clock;
% saved_folder=[];
% for i=4:6
%     saved_folder =[saved_folder num2str(round(tmp_time(i)))];
% end
% todays_date = [num2str(tmp_time(1)) num2str(tmp_time(2)) num2str(tmp_time(3))];
% check_pt_foldername = fullfile(lstm_folder_path,todays_date,saved_folder);
% if ~exist(check_pt_foldername)
%     mkdir(check_pt_foldername)
% end
% 
% options = trainingOptions('adam', ...
%     'MaxEpochs',60, ...
%     'MiniBatchSize',batch_size, ...
%     'GradientThreshold',10, ...
%     'Verbose',true, ...
%     'ValidationFrequency',val_freq,...
%     'Shuffle','every-epoch', ...
%     'ValidationData',{XTest,YTest},...
%     'ValidationPatience',patience,...
%     'Plots','training-progress',...
%     'LearnRateSchedule','piecewise',...
%     'LearnRateDropFactor',0.1,...    
%     'LearnRateDropPeriod',40,...
%     'InitialLearnRate',learning_rate,...
%     'CheckpointPath',check_pt_foldername)
% %     'ExecutionEnvironment', 'cpu');
% 
% % get rid of batch normalization layer as it throws an error due to matlab
% % version issues
% clear net
% layers = net_bilstm.Layers;
% idx=[];
% for i = 1:length(layers)
%     if ~strcmp(layers(i).Name,'batchnorm')
%         idx=[idx i];        
%     end
% end
% layers = layers(idx);
% 
% % train the model
% % IMPORTANT -> VALIDATION SHOULD NOT BE WILDLY DIFFERENT THAT TRAINING
% % LOSS. IF SO, REDUCE LEARNING RATE (2E-4 TO 1E-4 ETC.) AND RE RUN. IF IT
% % CONVERGES TOO QUICKLY (<10 EPOCHS) RE RUN AGAIN AS IT IS STUCK IN LOCAL
% % MINIMA
% 
% cd(check_pt_foldername) % if you want to  see the models being saved at each epoch
% net = trainNetwork(XTrain,YTrain,layers,options);
% 
% % loading the best performing (goat) model 
% saved_nets = dir(check_pt_foldername);
% saved_nets = saved_nets(3:end);
% [~,idx] = sort([saved_nets.datenum]);
% saved_nets = saved_nets(idx);
% goat_model = saved_nets(end-patience).name;
% goat_model = load(fullfile(check_pt_foldername,goat_model));
% 
% % saving to clicker folder
% cd(clicker_path)
% net_bilstm_robot_20220824_update_20230526_p = goat_model.net;
% save net_bilstm_robot_20220824_update_20230526_p net_bilstm_robot_20220824_update_20230526_p
% 







