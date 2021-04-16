
%% TRAIN CNN


clc;clear
root_path = '/home/ucsf/Data/bravo1/20210324/CenterOut';
foldernames = {'140739'};
cd(root_path)

% load the data for each target
D1=[];
D2=[];
D3=[];
D4=[];
for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'Imagined');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = length(features)-30:length(features);
        temp = cell2mat(features(kinax));
        temp = temp(129:end,:);
        if TrialData.TargetID == 1
            D1 = [D1 temp];
        elseif TrialData.TargetID == 2
            D2 = [D2 temp];
        elseif TrialData.TargetID == 3
            D3 = [D3 temp];
        elseif TrialData.TargetID == 4
            D4 = [D4 temp];
        end
    end
end


% DATA FROM ONLINE TASKS 24th MAR
root_path = '/home/ucsf/Data/Bravo2/20210324/DiscreteArrow/';
cd(root_path)
clear foldernames
D = dir(root_path);
k=1;
for i=3:length(D)
    foldernames{k} = D(i).name;
    k=k+1;
end
    
cd(root_path)
foldernames = foldernames(1:3);

% online successful data
for k=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{k},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        %if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = [find(TrialData.TaskState==2) find(TrialData.TaskState==3)];
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax) ;
            %idx = length(kinax)+-2 - TrialData.Params.ClickCounter:length(kinax) ;
            %idx=idx(idx>0);
            %kinax = kinax(idx);
            temp = cell2mat(features(kinax));
            temp = temp(129:end,:);
            
            if TrialData.TargetID == 1
                D1 = [D1 temp];
            elseif TrialData.TargetID == 2
                D2 = [D2 temp];
            elseif TrialData.TargetID == 3
                D3 = [D3 temp];
            elseif TrialData.TargetID == 4
                D4 = [D4 temp];
            end
        %end
    end
end



% DATA FROM ONLINE TASKS 31 MAR
root_path = '/home/ucsf/Data/Bravo2/20210331/DiscreteArrow/';
cd(root_path)
clear foldernames
D = dir(root_path);
k=1;
for i=3:length(D)
    foldernames{k} = D(i).name;
    k=k+1;
end
    
cd(root_path)


% online successful data
for k=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{k},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        %if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = [find(TrialData.TaskState==2) find(TrialData.TaskState==3)];
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax) ;
            %idx = length(kinax)+-2 - TrialData.Params.ClickCounter:length(kinax) ;
            %idx=idx(idx>0);
            %kinax = kinax(idx);
            temp = cell2mat(features(kinax));
            temp = temp(129:end,:);
            
            if TrialData.TargetID == 1
                D1 = [D1 temp];
            elseif TrialData.TargetID == 2
                D2 = [D2 temp];
            elseif TrialData.TargetID == 3
                D3 = [D3 temp];
            elseif TrialData.TargetID == 4
                D4 = [D4 temp];
            end
        %end
    end
end

% PIXELATED CURSOR 31 mar
root_path = '/home/ucsf/Data/Bravo2/20210331/PixelatedCursor/';
cd(root_path)
clear foldernames
D = dir(root_path);
k=1;
for i=3:length(D)
    foldernames{k} = D(i).name;
    k=k+1;
end
    
cd(root_path)
foldernames = foldernames(1:2);

% online successful data
for k=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{k},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        %if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = [find(TrialData.TaskState==2) find(TrialData.TaskState==3)];
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax) ;
            %idx = length(kinax)+-2 - TrialData.Params.ClickCounter:length(kinax) ;
            %idx=idx(idx>0);
            %kinax = kinax(idx);
            temp = cell2mat(features(kinax));
            temp = temp(129:end,:);
            
            if TrialData.TargetID == 1
                D1 = [D1 temp];
            elseif TrialData.TargetID == 2
                D2 = [D2 temp];
            elseif TrialData.TargetID == 3
                D3 = [D3 temp];
            elseif TrialData.TargetID == 4
                D4 = [D4 temp];
            end
        %end
    end
end



% DATA FROM ONLINE TASKS 7 APR
root_path = '/home/ucsf/Data/Bravo2/20210407/DiscreteArrow/';
cd(root_path)
clear foldernames
D = dir(root_path);
k=1;
for i=3:length(D)
    foldernames{k} = D(i).name;
    k=k+1;
end
    
cd(root_path)


% online successful data
for k=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{k},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        %if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = [find(TrialData.TaskState==2) find(TrialData.TaskState==3)];
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax) ;
            %idx = length(kinax)+-2 - TrialData.Params.ClickCounter:length(kinax) ;
            %idx=idx(idx>0);
            %kinax = kinax(idx);
            temp = cell2mat(features(kinax));
            temp = temp(129:end,:);
            
            if TrialData.TargetID == 1
                D1 = [D1 temp];
            elseif TrialData.TargetID == 2
                D2 = [D2 temp];
            elseif TrialData.TargetID == 3
                D3 = [D3 temp];
            elseif TrialData.TargetID == 4
                D4 = [D4 temp];
            end
        %end
    end
end



clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210407/CenterOut';
foldernames = {'153243'};
cd(root_path)

% load the data for each target
D1=[];
D2=[];
D3=[];
D4=[];
for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'Imagined');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = length(features)-30:length(features);
        temp = cell2mat(features(kinax));
        temp = temp(129:end,:);
        if TrialData.TargetID == 1
            D1 = [D1 temp];
        elseif TrialData.TargetID == 2
            D2 = [D2 temp];
        elseif TrialData.TargetID == 3
            D3 = [D3 temp];
        elseif TrialData.TargetID == 4
            D4 = [D4 temp];
        end
    end
end




clear condn_data
% combing both onlien plus offline
idx=1;
condn_data{1}=[D1(idx:end,:) ]'; % right hand
condn_data{2}= [D2(idx:end,:)]'; % both feet
condn_data{3}=[D3(idx:end,:)]'; % left hand
condn_data{4}=[D4(idx:end,:)]';


cd('/home/ucsf/Projects/bci')
load('/home/ucsf/Projects/bci/clicker/ECOG_Grid_8596-002131.mat')
chmap = ecog_grid;


condn_data_train={};
condn_data_test={};
for i=1:length(condn_data)
    temp = condn_data{i};
    % resize into images
    tmp_resized=[];
    for j=1:128:size(temp,2)
        chdata = temp(:,j:j+127);
        chdata = zscore(chdata')';
        % reshape as a grid
        chtemp=[];
        for k=1:size(chdata,1)
            t = chdata(k,:);
            chtemp(:,:,k) = t(chmap);            
        end
        tmp_resized = cat(4,tmp_resized,chtemp);        
    end
    temp = permute(tmp_resized,[1 2 4 3]);
    % splitting the data
    l = round(0.97*size(temp,4));
    idx = randperm(size(temp,4),l);
    % setting aside training trials
    condn_data_train{i} = temp(:,:,:,idx);
    % setting aside testing trials
    I = ones(size(temp,4),1);
    I(idx)=0;
    I=logical(I);
    condn_data_test{i} = temp(:,:,:,I);
end

% getting the data ready
XTrain = [];
YTrain = [];
for i=1:length(condn_data_train)
    tmp = condn_data_train{i};
    XTrain = cat(4,XTrain,tmp);
    YTrain = [YTrain;i*ones(size(tmp,4),1)];
end
YTrain = categorical(YTrain);
XTest = [];
YTest = [];
for i=1:length(condn_data_test)
    tmp = condn_data_test{i};
    XTest = cat(4,XTest,tmp);
    YTest = [YTest;i*ones(size(tmp,4),1)];
end
YTest = categorical(YTest);


%%%%%% CNN construction %%%%%
layers = [
    imageInputLayer([8 16 6])
    
    convolution2dLayer(2,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2,'Stride',1)
    
    convolution2dLayer(2,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2,'Stride',1)
    
    %     convolution2dLayer(2,32,'Padding','same')
    %     batchNormalizationLayer
    % %    dropoutLayer(0.2)
    %     reluLayer
    %     maxPooling2dLayer(2,'Stride',1)
    %
    %     convolution2dLayer(3,32,'Padding','same')
    %     batchNormalizationLayer
    %     reluLayer
    %     maxPooling2dLayer(3,'Stride',1)
    
    fullyConnectedLayer(32)
    batchNormalizationLayer
    dropoutLayer(0.2)
    reluLayer
    
    fullyConnectedLayer(16)
    batchNormalizationLayer
    dropoutLayer(0.2)
    reluLayer
    
    fullyConnectedLayer(4)
    softmaxLayer
    classificationLayer];

options = trainingOptions('adam', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',20, ...
    'Shuffle','every-epoch', ...
    'Verbose',true, ...
    'Plots','training-progress',...    
    'MiniBatchSize',80,...
    'ValidationData',{XTest,YTest}, ...
    'ValidationFrequency',30,...
    'L2Regularization',1e-4);

%%%%%% CNN construction %%%%%
% build the classifier
net = trainNetwork(XTrain,YTrain,layers,options);

cd('/home/ucsf/Projects/bci/clicker')
save CNN_classifier_B2_new_C net


%% need imagined mvmt data and online data
clc;clear
close all
cd('E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker')
load 6DOF_Online_Data_3Feat
load('E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20201030\DiscreteArrow\114043\BCI_Fixed\Data0007.mat')
chmap=TrialData.Params.ChMap;
% split the data into training and testing data 
% set the ecog activity as images of the grid, 3 channels
% format is width, height, channels and number of samples
condn_data_train={};
condn_data_test={};
for i=1:length(condn_data)
    temp = condn_data{i};
    % resize into images
    tmp_resized=[];
    for j=1:128:size(temp,2)
        chdata = temp(:,j:j+127);
        chdata = zscore(chdata')';
        % reshape as a grid
        chtemp=[];
        for k=1:size(chdata,1)
            t = chdata(k,:);
            chtemp(:,:,k) = t(chmap);            
        end
        tmp_resized = cat(4,tmp_resized,chtemp);        
    end
    temp = permute(tmp_resized,[1 2 4 3]);
    % splitting the data
    l = round(0.95*size(temp,4));
    idx = randperm(size(temp,4),l);
    % setting aside training trials
    condn_data_train{i} = temp(:,:,:,idx);
    % setting aside testing trials
    I = ones(size(temp,4),1);
    I(idx)=0;
    I=logical(I);
    condn_data_test{i} = temp(:,:,:,I);
end
% getting the data ready
XTrain = [];
YTrain = [];
for i=1:length(condn_data_train)
    tmp = condn_data_train{i};
    XTrain = cat(4,XTrain,tmp);
    YTrain = [YTrain;i*ones(size(tmp,4),1)];
end
YTrain = categorical(YTrain);
XTest = [];
YTest = [];
for i=1:length(condn_data_test)
    tmp = condn_data_test{i};
    XTest = cat(4,XTest,tmp);
    YTest = [YTest;i*ones(size(tmp,4),1)];
end
YTest = categorical(YTest);
%%%%%% CNN construction %%%%%
layers = [
    imageInputLayer([8 16 3])
    convolution2dLayer(2,4,'Padding','same')
    batchNormalizationLayer
    reluLayer    
    maxPooling2dLayer(2,'Stride',1)
    convolution2dLayer(2,8,'Padding','same')
    batchNormalizationLayer
    reluLayer    
    maxPooling2dLayer(2,'Stride',1)
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(3,'Stride',1)    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(3,'Stride',1)
    fullyConnectedLayer(6)
    softmaxLayer
    classificationLayer];
options = trainingOptions('adam', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',10, ...
    'Shuffle','every-epoch', ...
    'Verbose',true, ...
    'Plots','training-progress',...    
    'MiniBatchSize',32,...
    'ValidationFrequency',30,...
    'L2Regularization',1e-4);
%%%%%% CNN construction %%%%%
% build the classifier
net = trainNetwork(XTrain,YTrain,layers,options);
%analyzeNetwork(net)
%classify(net,zeros(28,28,1))
% 
% XTrain=randn(16,8,3,500);
% XTrain(:,:,:,251:end) = XTrain(:,:,:,251:end)+1;
% YTrain = categorical(([ones(250,1);2*ones(250,1)]));
%%%%%% TESTING THE CNN ABOVE ON DATA FROM 19TH MAR 2021
% test this NN classifier on data from other time-periods when B1 was
% clicking
filepath='E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20210326\Robot3DArrow';
foldernames = dir(filepath);
foldernames = foldernames(3:end);
files=[];
for i=3:length(foldernames)
   folderpath = fullfile(filepath, foldernames(i).name, 'BCI_Fixed');
   files = [files;findfiles('',folderpath)'];  
end
% get the predictions
pred_acc = zeros(7);
pred_acc_cnn = zeros(7);
for i=1:length(files)
    disp(i)
    load(files{i})
    chmap=TrialData.Params.ChMap;
    tid = TrialData.TargetID;
    feat = TrialData.SmoothedNeuralFeatures;
    idx = find(TrialData.TaskState==3);
    feat = cell2mat(feat(idx));
    % augment the accuracy for the neural net during testing
    pred = TrialData.ClickerState;
    for j=1:length(pred)
        if pred(j)==0
            pred(j)=7;
        end
        pred_acc(tid,pred(j)) = pred_acc(tid,pred(j))+1;
    end
    % find predictions from the CNN
    feat_idx = [129:256 513:640 769:896];
    for j=1:size(feat,2)
        chtemp = [];
        tmp = feat(feat_idx,j);
        f1 = zscore(tmp(1:128));
        f2 = zscore(tmp(129:256));
        f3 = zscore(tmp(257:384));
        chtemp(:,:,1) = f1(chmap);
        chtemp(:,:,2) = f2(chmap);
        chtemp(:,:,3) = f3(chmap);
        %out  = classify(net,chtemp);
        act = squeeze(activations(net,chtemp,20));
        [aa out]=max(act);
        if aa< TrialData.Params.NeuralNetSoftMaxThresh
            out=7;
        end
        pred_acc_cnn(tid,out) = pred_acc_cnn(tid,out)+1;
    end    
end