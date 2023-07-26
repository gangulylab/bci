%% TO CREATE CLASSIFIER BASED ON POOLING DATA



clc;clear
% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20230519/HandImagined';
% enter the folder names for the Task. These can be increased as more data
% is collected. For exaple: 

foldernames = {'110913', '110535', '110010', '105633', '105053', '104503'};

cd(root_path)

%FOR IMAGINED MOVEMENT DATA, 
D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
D7=[];
D8=[];
D9=[];
D10=[];
D11=[];
D12=[];
D13=[];
D14=[];
for ii=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{ii},'Imagined');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        
%         % get the pooled data
%         new_temp=[];
%         [xx yy] = size(TrialData.Params.ChMap);
%         for k=1:size(temp,2)
%             tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
%             tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
%             tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
%             pooled_data=[];
%             for i=1:2:xx
%                 for j=1:2:yy
%                     delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
%                     beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
%                     hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
%                     pooled_data = [pooled_data; delta; beta ;hg];
%                 end
%             end
%             new_temp= [new_temp pooled_data];
%         end
%         temp=new_temp;

        
        
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
        elseif TrialData.TargetID == 8
            D8 = [D8 temp];
        elseif TrialData.TargetID == 9
            D9 = [D9 temp];
        elseif TrialData.TargetID == 10
            D10 = [D10 temp];
        elseif TrialData.TargetID == 11
            D11 = [D11 temp];
        elseif TrialData.TargetID == 12
            D12 = [D12 temp];      
        elseif TrialData.TargetID == 13
            D13 = [D13 temp];      
        elseif TrialData.TargetID == 14
            D14 = [D14 temp];      
        end
    end
end

% FIXED
root_path = '/home/ucsf/Data/bravo1/20230519/HandOnline';
foldernames = {'112006', '112551', '113226', '113919'};

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
        
%         % get the pooled data
%         new_temp=[];
%         [xx yy] = size(TrialData.Params.ChMap);
%         for k=1:size(temp,2)
%             tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
%             tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
%             tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
%             pooled_data=[];
%             for i=1:2:xx
%                 for j=1:2:yy
%                     delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
%                     beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
%                     hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
%                     pooled_data = [pooled_data; delta; beta ;hg];
%                 end
%             end
%             new_temp= [new_temp pooled_data];
%         end
%         temp=new_temp;

        
        
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
        elseif TrialData.TargetID == 8
            D8 = [D8 temp];
        elseif TrialData.TargetID == 9
            D9 = [D9 temp];
        elseif TrialData.TargetID == 10
            D10 = [D10 temp];
        elseif TrialData.TargetID == 11
            D11 = [D11 temp];
        elseif TrialData.TargetID == 12
            D12 = [D12 temp];      
        elseif TrialData.TargetID == 13
            D13 = [D13 temp];      
        elseif TrialData.TargetID == 14
            D14 = [D14 temp];      
        end
    end
end

clear condn_data
% combing delta beta and high gamma
%idx=[1:size(D1,1)];
idx=[129:256 513:640 769:896];
%idx1=1:516;
condn_data{1}=[D1(idx,:) ]'; 
condn_data{2}= [D2(idx,:)]'; 
condn_data{3}=[D3(idx,:)]'; 
condn_data{4}=[D4(idx,:)]'; 
condn_data{5}=[D5(idx,:)]'; 
condn_data{6}=[D6(idx,:)]'; 
condn_data{7}=[D7(idx,:)]'; 
condn_data{8}=[D8(idx,:)]'; 
condn_data{9}=[D9(idx,:)]'; 
condn_data{10}=[D10(idx,:)]'; 
condn_data{11}=[D11(idx,:)]'; 
condn_data{12}=[D12(idx,:)]'; 
% condn_data{13}=[D13(idx,:)]'; 
% condn_data{14}=[D14(idx,:)]'; 

A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};
E = condn_data{5};
F = condn_data{6};
G = condn_data{7};
H = condn_data{8};
I = condn_data{9};
J = condn_data{10};
K = condn_data{11};
L = condn_data{12};
% M = condn_data{13};
% N1 = condn_data{14};


clear N
%N = [A' B' C' D' E' F' G' H' I' J' K' L' M' N1'];
N = [A' B' C' D' E' F' G' H' I' J' K' L' ];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    5*ones(size(E,1),1);6*ones(size(F,1),1);7*ones(size(G,1),1);8*ones(size(H,1),1);...
    9*ones(size(I,1),1);10*ones(size(J,1),1);11*ones(size(K,1),1);12*ones(size(L,1),1);...
    ];

T = zeros(size(T1,1),length(condn_data));
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
[aa bb]=find(T1==8);
T(aa(1):aa(end),8)=1;
[aa bb]=find(T1==9);
T(aa(1):aa(end),9)=1;
[aa bb]=find(T1==10);
T(aa(1):aa(end),10)=1;
[aa bb]=find(T1==11);
T(aa(1):aa(end),11)=1;
[aa bb]=find(T1==12);
T(aa(1):aa(end),12)=1;
% [aa bb]=find(T1==13);
% T(aa(1):aa(end),13)=1;
% [aa bb]=find(T1==14);
% T(aa(1):aa(end),14)=1;


%%%%% CODE SNIPPET FOR UPDATING A PRETRAINED DECODER %%%%%
% USE 2 BLOCKS OF ONLINE DAA, EACH BLOCK WITH 21 TRIALS %%%
cd('/home/ucsf/Projects/bci/clicker')
% load pretrain_net_mlp % NEW PNP DECODER FOR BATCH UPDATE
net_mlp = patternnet([64 64 64]) ;
net_mlp.performParam.regularization=0.2;
net_mlp.divideParam.trainRatio=0.8;
net_mlp.divideParam.valRatio=0.1;
net_mlp.divideParam.testRatio=0.1;
net_mlp = train(net_mlp,N,T','UseGPU','yes');
classifier_name = 'MLP_Hand_12DoF_CL2_05192023_noPool'; % enter the name
genFunction(net_mlp,classifier_name); % make sure to update Params.NeuralNetFunction in GetParams with the new name of the classifier


% 
% %%%% IS USING THE ADAM OPTIMIZER %%%%
% 
% % first load the saved mlp network
% load('net_mlp_7DoF_Feb2022.mat')
% 
% % get the layers of the model 
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
%     'MaxEpochs',100, ...
%     'Verbose',true, ...
%     'Plots','training-progress',...
%     'MiniBatchSize',32,...
%     'ValidationFrequency',32,...
%     'ValidationPatience',5,...
%     'ExecutionEnvironment','GPU',...
%     'ValidationData',{Xtest',Ytest});
% 
% % build the classifier
% net = trainNetwork(Xtrain',Ytrain,layers,options);
% net_mlp_7DoF_Feb2022 = net;
% 
% % save it
% save net_mlp_7DoF_Feb2022 net_mlp_7DoF_Feb2022

%%%%%%% CODE SNIPPET FOR TRAINING A MODEL FROM SCRATCH %%%%%
%%
% training a simple MLP
% IMPORTANT, CLICK THE CONFUSION MATRIX BUTTON IN GUI TO VERIFY THAT THE
% TEST VALIDATION DOESN'T HAVE NaNs AND THAT PERFORMANCE IS REASONABLE
%  clear net
%  net = patternnet([64 64 64]) ;
%  net.performParam.regularization=0.2;
% 
% % cd('/home/ucsf/Projects/bci/clicker')
% % load net net
% % 
%  net = train(net,N,T');
% % 
% 
% 
% %%%%%%% SAVING THE MODEL %%%%%%%%
% % cd('/home/ucsf/Projects/bci/clicker')
% % save net net
% 
% % classifier name
%  classifier_name = 'MLP_20210922';
%  genFunction(net,classifier_name); % make sure to update GetParams
% % % 
% % % 
% % % % to restart exp run following lines
  clear
  clc
  cd('/home/ucsf/Projects/bci')
 ExperimentStart('HandOnline','bravo1',4,1,0)
% ExperimentStart('Robot3D','bravo1',4,1,0)
%  %ExperimentStart('Robot3DArrow','bravo1',4,1,0)
% %  ExperimentStart('RobotR2GModeSwitch','bravo1',4,1,0)
%  
% %  
