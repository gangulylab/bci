%% UPDATING THE PNP MLP WITH ROBOT DATA

clc;clear
close all

% PATH
addpath('/home/ucsf/Projects/bci/task_helpers')

% ROBOT DATA
%root_path= 'F:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate B3\20231127\Robot3D'
root_path = '/home/ucsf/Data/Bravo3/20231229/RealRobotBatch';
folders = {'140019', '140316', '140555','141354'};
files=[];
for ii=1:length(folders)
    folderpath = fullfile(root_path, folders{ii},'BCI_Fixed');     
    files = [files;findfiles('',folderpath)'];
end

%load all the data 
condn_data_overall = [load_data_for_MLP_TrialLevel_B3_bci(files,1)];

% split the data into validation and training sets and get training options
num_classes=7; % 7-actions
chk = 0;iter=0;
while chk==0
    test_idx = randperm(length(condn_data_overall),round(0.3*length(condn_data_overall)));
    test_idx=test_idx(:);
    I = ones(length(condn_data_overall),1);
    I(test_idx)=0;
    train_idx = find(I~=0);train_idx=train_idx(:);
    [options,XTrain,YTrain] = ...
        get_options(condn_data_overall,test_idx,train_idx,1e-4);
    ytest=options.ValidationData{2};
    options.Plots='training-progress';
    if length(unique(double(YTrain))) == num_classes && ...
            length(unique(double(ytest))) == num_classes
        chk=1;
    end
    iter=iter+1;
end
disp(['Data split done in ' num2str(iter) ' iterations'])

%%%%%%%% CODE TO BATCH UPDATE THE PNP DECODER %%%%%%%
cd('/home/ucsf/Projects/bci/clicker')
load net_PnP
layers=net_PnP.Layers;
net_PnP_RobotBatchUpdate = trainNetwork(XTrain,YTrain,layers,options);
save net_PnP_RobotBatchUpdate net_PnP_RobotBatchUpdate %CARRY THIS OVER TO GET PARAMS



%% RUN THE TASK AGAIN AFTER UPDATING THE DECODER NAME 

cd('/home/ucsf/Projects/bci')
clear;clc
% ExperimentStart('Robot3DArrow','Bravo3',4,1,0)
ExperimentStart('Robot3D','Bravo3',4,1,0)


