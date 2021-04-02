%% CODE FOR SONOMA IN RETRAINING MODEL

% run from the BCI folder 

clc;clear
root_path = 'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker';
foldernames = {'20201218','20210108','20210115','20210128','20210201','20210212','20210219','20210226',...
    '20210305','20210312','20210319'};
cd(root_path)
addpath(fullfile('task_helpers'))

files=[];
for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'Robot3DArrow');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name,'BCI_Fixed');
        files = [files;findfiles('',filepath)'];
    end
end


% load the data for each target
D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
for i=1:length(files)
    disp(i)
    load(files{i});
    features  = TrialData.SmoothedNeuralFeatures;
    kinax = TrialData.TaskState;
    if kinax(1) == 0
        kinax=kinax(2:end);
    end
    idx = [find(kinax==2) find(kinax==3)];
    l = idx(end)-7:idx(end);
    idx=l(l>0);    
    kinax =  (idx);
    
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
    elseif TrialData.TargetID == 5
        D5 = [D5 temp];
    elseif TrialData.TargetID == 6
        D6 = [D6 temp];
    end
end

clear condn_data
%idx=641;
idx = [1:128 385:512 641:768];
condn_data{1}=[D1(idx,:) ]'; % right hand
condn_data{2}= [D2(idx,:)]'; % both feet
condn_data{3}=[D3(idx,:)]'; % left hand
condn_data{4}=[D4(idx,:)]'; % head
condn_data{5}=[D5(idx,:)]'; % mime up
condn_data{6}=[D6(idx,:)]'; % tongue in


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
        chdata = (chdata')';
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
    maxPooling2dLayer(2,'Stride',1)    
    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2,'Stride',1)
        
    fullyConnectedLayer(6)
    softmaxLayer
    classificationLayer];

%'ValidationData',{XTest,YTest},...
options = trainingOptions('adam', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',20, ...
    'Shuffle','every-epoch', ...
    'Verbose',true, ...
    'Plots','training-progress',...    
    'MiniBatchSize',64,...
    'ValidationFrequency',30,...
    'L2Regularization',1e-4,...
    'ValidationData',{XTest,YTest},...
    'ExecutionEnvironment','auto');

%%%%%% CNN construction %%%%%


% build the classifier
net = trainNetwork(XTrain,YTrain,layers,options);
% save this in the clicker folder
save CNN_classifier net

