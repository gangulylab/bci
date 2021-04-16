%% CODE FOR SONOMA IN RETRAINING MODEL

% run from the BCI folder 

clc;clear
root_path = '/home/ucsf/Data/bravo1';
foldernames = {'20210115','20210128','20210201','20210212','20210219','20210226',...
    '20210305','20210312','20210319','20210326','20210402','20210409','20210416'};
addpath(fullfile('task_helpers'))

files=[];
for i=length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'Robot3DArrow');
    D=dir(folderpath);
    for j=8:length(D)
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
    idx = [find(kinax==2)];
    %l = idx(end)-7:idx(end);
    %idx=l(l>0);    
    % kinax =  (idx);
     kinax = find(kinax==3);
     %l=length(kinax)+1;
     %kinax = kinax(l-TrialData.Params.ClickCounter:end);
     kinax = [idx kinax];
    
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
    %temp = temp(:,:,1,:);
    
    
    % splitting the data
    l = round(0.98*size(temp,4));
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
    
    convolution2dLayer(2,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2,'Stride',1)    
%     
    convolution2dLayer(3,32,'Padding','same')    
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(3,'Stride',1)   
    
    fullyConnectedLayer(128)
    batchNormalizationLayer
    reluLayer
    dropoutLayer(.5)
    
%     fullyConnectedLayer(64)
%     batchNormalizationLayer
%     reluLayer
%     dropoutLayer(.5)
        
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
%analyzeNetwork(net)
% save this in the clicker folder
cd('/home/ucsf/Projects/bci/clicker')
save CNN_classifier_Online_Apr16_2021_C net


%% TO TRAIN MULTI LAYER PERCEPTRON

clear N
idx = [1:128 385:514 641:768];
A = condn_data{1};B = condn_data{2};C = condn_data{3};D = condn_data{4};
E = condn_data{5};F = condn_data{6};

clear N
N = [A' B' C' D' E' F'];
%N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    5*ones(size(E,1),1);6*ones(size(F,1),1)];
%T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1)];


T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;
[aa bb]=find(T1==5);[aa(1) aa(end)]
T(aa(1):aa(end),5)=1;
[aa bb]=find(T1==6);[aa(1) aa(end)]
T(aa(1):aa(end),6)=1;

% 
% % All 12D of control
% clear N
% idx = [1:128 385:514 641:768];
% A = condn_data{1};B = condn_data{2};C = condn_data{3};D = condn_data{4};
% E = condn_data{5};F = condn_data{6};G = condn_data{7};H = condn_data{8};
% I = condn_data{9};J = condn_data{10};K = condn_data{11};L = condn_data{12};
% N = [A(:,idx)' B(:,idx)' C(:,idx)' D(:,idx)' E(:,idx)' F(:,idx)'...
%     G(:,idx)' H(:,idx)' I(:,idx)' J(:,idx)' K(:,idx)' L(:,idx)'];
% %N = [A' B' C' D'];
% T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
%     5*ones(size(E,1),1);6*ones(size(F,1),1);7*ones(size(G,1),1);8*ones(size(H,1),1);...
%     9*ones(size(I,1),1);10*ones(size(J,1),1);11*ones(size(K,1),1);12*ones(size(L,1),1)];
% %T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1)];
% 
% 
% T = zeros(size(T1,1),12);
% [aa bb]=find(T1==1);[aa(1) aa(end)]
% T(aa(1):aa(end),1)=1;
% [aa bb]=find(T1==2);[aa(1) aa(end)]
% T(aa(1):aa(end),2)=1;
% [aa bb]=find(T1==3);[aa(1) aa(end)]
% T(aa(1):aa(end),3)=1;
% [aa bb]=find(T1==4);[aa(1) aa(end)]
% T(aa(1):aa(end),4)=1;
% [aa bb]=find(T1==5);[aa(1) aa(end)]
% T(aa(1):aa(end),5)=1;
% [aa bb]=find(T1==6);[aa(1) aa(end)]
% T(aa(1):aa(end),6)=1;
% [aa bb]=find(T1==7);[aa(1) aa(end)]
% T(aa(1):aa(end),7)=1;
% [aa bb]=find(T1==8);[aa(1) aa(end)]
% T(aa(1):aa(end),8)=1;
% [aa bb]=find(T1==9);[aa(1) aa(end)]
% T(aa(1):aa(end),9)=1;
% [aa bb]=find(T1==10);[aa(1) aa(end)]
% T(aa(1):aa(end),10)=1;
% [aa bb]=find(T1==11);[aa(1) aa(end)]
% T(aa(1):aa(end),11)=1;
% [aa bb]=find(T1==12);[aa(1) aa(end)]
% T(aa(1):aa(end),12)=1;

% code to train a neural network
net = patternnet([128 128 128 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T','UseParallel','yes');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'multilayer_perceptron_6DoF_Online_Apr16_2021')
delete(gcp)
%genFunction(net,'multilayer_perceptron_4Dir_MimeUpTongueIn_OnlineData')

