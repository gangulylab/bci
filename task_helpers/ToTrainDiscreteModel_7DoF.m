%% Creating MLP classifier
clc;clear
% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20210709/Robot3DArrow';
% enter the folder names for the Task. These can be increased as more data
% is collected. For exaple: 

% foldernames = {'111533', '112105', '112436', '112747', '113524',
% '113817', '114114',  '115402'}; %morning

% foldernames = {'134509', '134821', '135145', '135525', '140741'}; 
% 
% foldernames = {'111602', '112055', '112427', '112811', '113337', '113552', '113924', '114205'...
%     '114800', '115132', '115413', '132949', '133554', '133844'};


% foldernames = {'132949', '133554', '133844', '134859', '135041', '135057', '135324', '135652', '140726', '144542', '144754', '145042'};

% foldernames = {'110604', '111123', '111649', '112201', '113524',
% '113754', '113909', '114318', '114537'}; % 

% foldernames = {'101338', '101857', '102342', '102825'};

% foldernames = {'133210', '133905', '134420'} 

%

cd(root_path)

%FOR IMAGINED MOVEMENT DATA, 
D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
D7=[];
for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'Imagined');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [find(TrialData.TaskState==2) find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        temp = temp(129:end,:);
        
        
%         %adaptive baseline training
%         idx_bl = find(TrialData.TaskState==1);
%         temp_bl = cell2mat(features(idx_bl));
%         temp_bl = temp_bl(129:end,:);
%         m = mean(temp_bl,2);
%         s = std(temp_bl')';
%         temp = (temp-m)./s;
        
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
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end

% FOR ONLINE DATA
for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [find(TrialData.TaskState==2) find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        temp = temp(129:end,:);
        
        
%         %adaptive baseline training
%         idx_bl = find(TrialData.TaskState==1);
%         temp_bl = cell2mat(features(idx_bl));
%         temp_bl = temp_bl(129:end,:);
%         m = mean(temp_bl,2);
%         s = std(temp_bl')';
%         temp = (temp-m)./s;
        
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
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end



clear condn_data
% combing delta beta and high gamma
idx=[1:128 385:512 641:768];
condn_data{1}=[D1(idx,:) ]'; 
condn_data{2}= [D2(idx,:)]'; 
condn_data{3}=[D3(idx,:)]'; 
condn_data{4}=[D4(idx,:)]'; 
condn_data{5}=[D5(idx,:)]'; 
condn_data{6}=[D6(idx,:)]'; 
condn_data{7}=[D7(idx,:)]'; 

A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};
E = condn_data{5};
F = condn_data{6};
G = condn_data{7};

clear N
N = [A' B' C' D' E' F' G'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    5*ones(size(E,1),1);6*ones(size(F,1),1);7*ones(size(G,1),1)];

T = zeros(size(T1,1),7);
[aa bb]=find(T1==1);
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);
T(aa(1):aa(end),4)=1;
[aa bb]=find(T1==5);
T(aa(1):aa(end),5)=1;
[aa bb]=find(T1==6);
T(aa(1):aa(end),6)=1;
[aa bb]=find(T1==7);
T(aa(1):aa(end),7)=1;


%%%%%%% CODE SNIPPET FOR TRAINING A MODEL FROM SCRATCH %%%%%
% training a simple MLP
% IMPORTANT, CLICK THE CONFUSION MATRIX BUTTON IN GUI TO VERIFY THAT THE
% TEST VALIDATION DOESN'T HAVE NaNs AND THAT PERFORMANCE IS REASONABLE
net = patternnet([128 128 128 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');
%%%%%%% END CODE SNIPPET %%%%%%%%%%%%%%%

%%% SAVING %%%%%%%
cd('/home/ucsf/Projects/bci/clicker')
% classifier name
classifier_name = 'MLP_Imag_Actions_0615_E';
genFunction(net,classifier_name); % make sure to update GetParams


% to restart exp run following lines
clear
clc
cd('/home/ucsf/Projects/bci')
 ExperimentStart('Robot3DArrow','bravo1',4,1,0)

 %ExperimentStart('Robot3D','bravo1',4,1,0)
% 
% ExperimentStart('RobotReachStop','bravo1',4,1,0)

%% TO CREATE CLASSIFIER BASED ON POOLING DATA



clc;clear
% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20220727/Robot3DArrow';
% enter the folder names for the Task. These can be increased as more data
% is collected. For exaple: 

% root_path = 'F:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20220323\Robot3DArrow\'



cd(root_path)

%FOR IMAGINED MOVEMENT DATA, 

%foldernames = {'111214', '111833', '112240', '112646','142116','141531','135832','135152'};
foldernames = {};


D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
D7=[];
for ii=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{ii},'Imagined');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        
        % get the pooled data
        new_temp=[];
        [xx yy] = size(TrialData.Params.ChMap);
        for k=1:size(temp,2)
            tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
            tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
            tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
            pooled_data=[];
            for i=1:2:xx
                for j=1:2:yy
                    delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
                    beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
                    hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
                    pooled_data = [pooled_data; delta; beta ;hg];
                end
            end
            new_temp= [new_temp pooled_data];
        end
        temp=new_temp;

        
        
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
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end

% FIXED ARROW
foldernames = {'105959', '110625', '111023', '111324', '111624'};

cd(root_path)
 
for ii=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{ii},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        
        
        % get the pooled data
        new_temp=[];
        [xx yy] = size(TrialData.Params.ChMap);
        for k=1:size(temp,2)
            tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
            tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
            tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
            pooled_data=[];
            for i=1:2:xx
                for j=1:2:yy
                    delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
                    beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
                    hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
                    pooled_data = [pooled_data; delta; beta ;hg];
                end
            end
            new_temp= [new_temp pooled_data];
        end
        temp=new_temp;
        
        
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
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end

size(D7)

% ROBOT BATCH
root_path = '/home/ucsf/Data/bravo1/20220720/RealRobotBatch';
foldernames = {};

for ii=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{ii},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        
        
        % get the pooled data
        new_temp=[];
        [xx yy] = size(TrialData.Params.ChMap);
        for k=1:size(temp,2)
            tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
            tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
            tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
            pooled_data=[];
            for i=1:2:xx
                for j=1:2:yy
                    delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
                    beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
                    hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
                    pooled_data = [pooled_data; delta; beta ;hg];
                end
            end
            new_temp= [new_temp pooled_data];
        end
        temp=new_temp;
        
        
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
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end
size(D7)

clear condn_data
% combing delta beta and high gamma
idx=[1:size(D1,1)];
condn_data{1}=[D1(idx,:) ]'; 
condn_data{2}= [D2(idx,:)]'; 
condn_data{3}=[D3(idx,:)]'; 
condn_data{4}=[D4(idx,:)]'; 
condn_data{5}=[D5(idx,:)]'; 
condn_data{6}=[D6(idx,:)]'; 
condn_data{7}=[D7(idx,:)]'; 


% 2norm
for i=1:length(condn_data)
   tmp = condn_data{i}; 
   for j=1:size(tmp,1)
       tmp(j,:) = tmp(j,:)./norm(tmp(j,:));
   end
   condn_data{i}=tmp;
end


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};
E = condn_data{5};
F = condn_data{6};
G = condn_data{7};

clear N
N = [A' B' C' D' E' F' G'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    5*ones(size(E,1),1);6*ones(size(F,1),1);7*ones(size(G,1),1)];

T = zeros(size(T1,1),7);
[aa bb]=find(T1==1);
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);
T(aa(1):aa(end),4)=1;
[aa bb]=find(T1==5);
T(aa(1):aa(end),5)=1;
[aa bb]=find(T1==6);
T(aa(1):aa(end),6)=1;
[aa bb]=find(T1==7);
T(aa(1):aa(end),7)=1;


%%%%%% NEW CODE SNIPPET FOR CREATING A DECODER FROM SCRATCH 05/25/2005 %%%%%

%if new from scratch
% clear net 
% net = patternnet([48 48 48]) ;
% net.performParam.regularization=0.2;
% 
% % if loading from pretrained
% %load net_new_7DoF_ZWrist_05252022 
% %net = net_new_7DoF_ZWrist_05252022.net 
% 
% 
% net = train(net,N,T');
% 
% % give it a name, this goes to  Params.NeuralNetFunctionName in line 86 and
% % 87 of GetParams
% net_new_7DoF_ZWrist_07012022D = net;
% 
% % save the weights
% cd('/home/ucsf/Projects/bci/clicker')
% save net_new_7DoF_ZWrist_07012022D net_new_7DoF_ZWrist_07012022D
% classifier_name = 'MLP_7DoF_ZWrist_07012022D'; % enter the name
% genFunction(net_new_7DoF_ZWrist_07012022D,classifier_name);



%%%%% CODE SNIPPET FOR UPDATING A PRETRAINED DECODER %%%%%
% USE 2 BLOCKS OF ONLINE DAA, EACH BLOCK WITH 21 TRIALS %%%
% cd('/home/ucsf/Projects/bci/clicker')
% load net_7DoF_PnP_2022Mar_2norm 
% % load pretrain_net_mlp % NEW PNP DECODER FOR BATCH UPDATE
% 
% %net_7DoF_PnP_2022Mar_2norm.divideParam.trainRatio=0.8;
% %net_7DoF_PnP_2022Mar_2norm.divideParam.valRatio=0.1;
% %net_7DoF_PnP_2022Mar_2norm.divideParam.testRatio=0.1;
% 
% 
% net_7DoF_PnP_2022Mar_2norm.trainParam.epochs=30;
% 
% net_7DoF_PnP_2022Mar_2norm = train(net_7DoF_PnP_2022Mar_2norm,N,T');
% 
% classifier_name = 'MLP_7DoF_PnP_2022Mar_2norm_0518_pm1'; % enter the name
% genFunction(net_7DoF_PnP_2022Mar_2norm,classifier_name); % make sure to update Params.NeuralNetFunction in GetParams with the new name of the classifier
% 
% %%%% limit the nunber of epochs the NN trains for
% net_7DoF_PnP_2022Mar_2norm.trainParam.epochs=25
% 
% %%%% lower the learning rate of the NN weight change
% net_7DoF_PnP_2022Mar_2norm.trainParam.sigma=5e-7

% 
% updating the ensemble decoder
% load net_7DoF_PnP4_ensemble_batch_0520A
% net_7DoF_PnP4_ensemble_batch_0520B=net_7DoF_PnP4_ensemble_batch_0520A;
% for i=2:length(net_7DoF_PnP4_ensemble_batch_0520B)
%     net = net_7DoF_PnP4_ensemble_batch_0520B{i};
%     net = train(net,N,T','useParallel','no');
%     net_7DoF_PnP4_ensemble_batch_0520B{i} = net;
% end
% cd('/home/ucsf/Projects/bci/clicker')
% save net_7DoF_PnP4_ensemble_batch_0520B net_7DoF_PnP4_ensemble_batch_0520B



% 
% %%%% IS USING THE ADAM OPTIMIZER %%%%
% 
% % first load the saved mlp network
% load('net_mlp_7DoF_Feb2022.mat')
% 
% % get the layers of the model 
% clear layers
% layers=net_mlp_7DoF_Feb2022.Layers;
% 
% % organize the data
% X = N;
% Y=categorical(T1);
% idx = randperm(length(Y),round(0.8*length(Y)));
% Xtrain = X(:,idx);
% Ytrain = Y(idx);
% I = ones(length(Y),1);
% I(idx)=0;
% idx1 = find(I~=0);
% Xtest = X(:,idx1);
% Ytest = Y(idx1);
% 
% 
% 
% %'ValidationData',{XTest,YTest},...
% options = trainingOptions('adam', ...
%     'InitialLearnRate',0.01, ...
%     'MaxEpochs',20, ...
%     'Verbose',true, ...
%     'Plots','training-progress',...
%     'MiniBatchSize',32,...
%     'ValidationFrequency',30,...
%     'ValidationPatience',6,...
%     'ExecutionEnvironment','cpu',...
%     'ValidationData',{Xtest',Ytest});

% % build the classifier
% net = trainNetwork(Xtrain',Ytrain,layers,options);
% net_mlp_7DoF_Feb2022 = net;
% 
% % save it
% save net_mlp_7DoF_Feb2022 net_mlp_7DoF_Feb2022

%%%%%%% CODE SNIPPET FOR TRAINING A MODEL FROM SCRATCH %%%%%
% training a simple MLP
% IMPORTANT, CLICK THE CONFUSION MATRIX BUTTON IN GUI TO VERIFY THAT THE
% TEST VALIDATION DOESN'T HAVE NaNs AND THAT PERFORMANCE IS REASONABLE
%  clear net
%  net = patternnet([64 64 64]) ;
%  net.performParam.regularization=0.2;
% 
% cd('/home/ucsf/Projects/bci/clicker')
% % % load net net
% % % 
%  net = train(net,N,T');
% % % 
% % 
% % 
% % %%%%%%% SAVING THE MODEL %%%%%%%%
% % % cd('/home/ucsf/Projects/bci/clicker')
% % % save net net
% % 
% % % classifier name
% net_new_7DoF_05252022 = net;
% classifier_name = 'net_new_7DoF_05252022';
% genFunction(net_new_7DoF_05252022,classifier_name); % make sure to update GetParams
% % % % 
% % % % 
% % % % % to restart exp run following lines
%  
% 
% % to restart exp run following lines
% % clear
% % clc
% % cd('/home/ucsf/Projects/bci')
% % ExperimentStart('Robot3DArrow','bravo1',4,1,0)
% %ExperimentStart('RealRobotBatch','bravo1',4,1,0)
% 
% 
%  ExperimentStart('Robot3DArrow','bravo1',4,1,0)
%  cd('/home/ucsf/Projects/bci')
%  ExperimentStart('RealRobotBatch','bravo1',4,1,0)
% ExperimentStart('Robot3D','bravo1',4,1,0)
%  %ExperimentStart('Robot3DArrow','bravo1',4,1,0)
% %  ExperimentStart('RobotR2GModeSwitch','bravo1',4,1,0)
%  
% %  

%% TO UPDATE CLASSIFIER USING POOLING ON 4 FEATURES with low gamma

clc;clear
% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20220727/Robot3DArrow';
% enter the folder names for the Task. These can be increased as more data
% is collected. For exaple: 

% root_path = 'F:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20220323\Robot3DArrow\'



cd(root_path)

%FOR IMAGINED MOVEMENT DATA, 

%foldernames = {'111214', '111833', '112240', '112646','142116','141531','135832','135152'};
foldernames = {};


D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
D7=[];
for ii=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{ii},'Imagined');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        
        %get the pooled data
        new_temp=[];
        [xx yy] = size(TrialData.Params.ChMap);
        for k=1:size(temp,2)
            tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
            tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
            tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
            tmp4 = temp(641:768,k);tmp4 = tmp4(TrialData.Params.ChMap);%low gamma
            pooled_data=[];
            for i=1:2:xx
                for j=1:2:yy
                    delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
                    beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
                    hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
                    lg = (tmp4(i:i+1,j:j+1));lg=mean(lg(:));
                    pooled_data = [pooled_data; delta; beta ;lg;hg];
                end
            end
            new_temp= [new_temp pooled_data];
        end
        temp_data=new_temp;

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
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end

% FIXED ARROW
foldernames = {'131235', '131806','132102', '132445', '132755'};

cd(root_path)
 
for ii=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{ii},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        
        
    %get the pooled data
        new_temp=[];
        [xx yy] = size(TrialData.Params.ChMap);
        for k=1:size(temp,2)
            tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
            tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
            tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
            tmp4 = temp(641:768,k);tmp4 = tmp4(TrialData.Params.ChMap);%low gamma
            pooled_data=[];
            for i=1:2:xx
                for j=1:2:yy
                    delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
                    beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
                    hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
                    lg = (tmp4(i:i+1,j:j+1));lg=mean(lg(:));
                    pooled_data = [pooled_data; delta; beta ;lg;hg];
                end
            end
            new_temp= [new_temp pooled_data];
        end
        temp_data=new_temp;
        temp = temp_data;
        
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
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end

size(D7)

% ROBOT BATCH
root_path = '/home/ucsf/Data/bravo1/20220722/RealRobotBatch';
foldernames = {};

for ii=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{ii},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        
        
%get the pooled data
        new_temp=[];
        [xx yy] = size(TrialData.Params.ChMap);
        for k=1:size(temp,2)
            tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
            tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
            tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
            tmp4 = temp(641:768,k);tmp4 = tmp4(TrialData.Params.ChMap);%low gamma
            pooled_data=[];
            for i=1:2:xx
                for j=1:2:yy
                    delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
                    beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
                    hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
                    lg = (tmp4(i:i+1,j:j+1));lg=mean(lg(:));
                    pooled_data = [pooled_data; delta; beta ;lg;hg];
                end
            end
            new_temp= [new_temp pooled_data];
        end
        temp_data=new_temp;
        temp = temp_data;
        
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
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end
size(D7)



% REAL ROBOT PATH
root_path = '/home/ucsf/Data/bravo1/20220720/RealRobotPath';
foldernames = {};

for ii=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{ii},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        
        
        % get the pooled data
        new_temp=[];
        [xx yy] = size(TrialData.Params.ChMap);
        for k=1:size(temp,2)
            tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
            tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
            tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
            tmp4 = temp(641:768,k);tmp4 = tmp4(TrialData.Params.ChMap);%low gamma
            pooled_data=[];
            for i=1:2:xx
                for j=1:2:yy
                    delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
                    beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
                    hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
                    lg = (tmp4(i:i+1,j:j+1));lg=mean(lg(:));
                    pooled_data = [pooled_data; delta; beta ;lg;hg];
                end
            end
            new_temp= [new_temp pooled_data];
        end
        temp_data=new_temp;
        temp = temp_data;
        
        % select bins in trial for each directions
        D1 = [D1 temp(:,TrialData.CorrectDecode == 1)];
        D2 = [D2 temp(:,TrialData.CorrectDecode == 2)];
        D3 = [D3 temp(:,TrialData.CorrectDecode == 3)];
        D4 = [D4 temp(:,TrialData.CorrectDecode == 4)];
        D5 = [D5 temp(:,TrialData.CorrectDecode == 5)];
        D6 = [D6 temp(:,TrialData.CorrectDecode == 6)];
        D7 = [D7 temp(:,TrialData.CorrectDecode == 7)];
    end
end
size(D7)

clear condn_data
% combing delta beta and high gamma
idx=[1:size(D1,1)];
condn_data{1}=[D1(idx,:) ]'; 
condn_data{2}= [D2(idx,:)]'; 
condn_data{3}=[D3(idx,:)]'; 
condn_data{4}=[D4(idx,:)]'; 
condn_data{5}=[D5(idx,:)]'; 
condn_data{6}=[D6(idx,:)]'; 
condn_data{7}=[D7(idx,:)]'; 


% 2norm
for i=1:length(condn_data)
   tmp = condn_data{i}; 
   for j=1:size(tmp,1)
       tmp(j,:) = tmp(j,:)./norm(tmp(j,:));
   end
   condn_data{i}=tmp;
end


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};
E = condn_data{5};
F = condn_data{6};
G = condn_data{7};

clear N
N = [A' B' C' D' E' F' G'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    5*ones(size(E,1),1);6*ones(size(F,1),1);7*ones(size(G,1),1)];

T = zeros(size(T1,1),7);
[aa bb]=find(T1==1);
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);
T(aa(1):aa(end),4)=1;
[aa bb]=find(T1==5);
T(aa(1):aa(end),5)=1;
[aa bb]=find(T1==6);
T(aa(1):aa(end),6)=1;
[aa bb]=find(T1==7);
T(aa(1):aa(end),7)=1;


%%%%%% UPDATING THE DECODER WEIGHTS %%%%%

%if new from scratch
%clear net 
%net = patternnet([64 64 64]) ;
%net.performParam.regularization=0.2;

% if loading from pretrained
load net_7DoF_PnP_lg 
net = net_7DoF_PnP_lg;


net = train(net,N,T');

% save the weights
cd('/home/ucsf/Projects/bci/clicker')
classifier_name = 'MLP_7DoF_PnP_2022July_lg_0727_PM'; % enter the name
genFunction(net,classifier_name);



% clc;clear
% cd('/home/ucsf/Projects/bci')
%  %ExperimentStart('Robot3DArrow','bravo1',4,1,0)
% %  cd('/home/ucsf/Projects/bci')
%   ExperimentStart('RealRobotBatch','bravo1',4,1,0)
% ExperimentStart('Robot3D','bravo1',4,1,0)





