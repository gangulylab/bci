
%% CREATING A  CLASSIFIER (either new or batch update, select apporp folder)
% USES POOLED DATA

% if creating new classifier use only imagined mvmt data folders
% if doing batch updae, use both imagined and fixed folders

clc;clear

% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20230726/Robot3DArrow';

% enter the folder names for the Task. These can be increased as more data
% is collected. For example: 

foldernames = {'132810', '133457', '133936', '134402'};

cd(root_path)

%FOR IMAGINED MOVEMENT DATA 
D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
D7=[];
D8=[];
D9=[];
for ii=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{ii},'Imagined');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
        temp = cell2mat(features(kinax));
        
       % get data in the 4 features
       temp=temp([129:256 513:640 641:768 769:896],:);
        
        
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
        end
    end
end

% FIXED CONTROL
foldernames = {'135617', '140041', '140527', '140914'};
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
        
        
       % get data in the 4 features
       temp=temp([129:256 513:640 641:768 769:896],:);
        
        
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
        end
    end
end

size(D7)
%ROBOT BATCH
root_path = '/home/ucsf/Data/bravo1/20220216/RealRobotBatch';
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
        
       % get data in the 4 features
       temp=temp([129:256 513:640 641:768 769:896],:); 
               
        
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
        end
    end
end
size(D7)

clear condn_data
% combing delta beta and high gamma
idx=[1:size(D1,1)];
condn_data{1}=[D1(idx,:) ]'; 
condn_data{2}=[D2(idx,:)]'; 
condn_data{3}=[D3(idx,:)]'; 
condn_data{4}=[D4(idx,:)]'; 
condn_data{5}=[D5(idx,:)]'; 
condn_data{6}=[D6(idx,:)]'; 
condn_data{7}=[D7(idx,:)]'; 
condn_data{8}=[D8(idx,:)]'; 
condn_data{9}=[D9(idx,:)]'; 

A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};
E = condn_data{5};
F = condn_data{6};
G = condn_data{7};
H = condn_data{8};
I = condn_data{9};

clear N
N = [A' B' C' D' E' F' G' H' I'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    5*ones(size(E,1),1);6*ones(size(F,1),1);7*ones(size(G,1),1);8*ones(size(H,1),1);...
    9*ones(size(I,1),1)];

T = zeros(size(T1,1),9);
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





%%%%% CODE SNIPPET FOR UPDATING A PRETRAINED MODEL %%%%%
cd('/home/ucsf/Projects/bci/clicker')
load net_9DoF_July2023

net_9DoF_July2023.divideParam.trainRatio=0.85;
net_9DoF_July2023.divideParam.valRatio=0.15;
net_9DoF_July2023.divideParam.testRatio=0;
net_9DoF_July2023 = train(net_9DoF_July2023,N,T','useParallel','yes');
classifier_name = 'MLP_9DoF_07262023Update2'; % enter the name
genFunction(net_9DoF_July2023,classifier_name); % make sure to update Params.NeuralNetFunction in GetParams with the new name of the classifier
delete(gcp)


%%%%% CODE SNIPPET FOR TRAINING AND SAVING THE DECODER FROM SCRATCH %%%%%
% cd('/home/ucsf/Projects/bci/clicker')
% clear net
% net = patternnet([64 64 64]);
% net.performParam.regularization=0.2;
% net.divideParam.trainRatio=0.8;
% net.divideParam.valRatio=0.1;
% net.divideParam.testRatio=0.1;
% net = train(net,N,T','useParallel','yes');
% classifier_name = 'MLP_PreTrained_9DoF_02092022_PM2'; % enter the name
% genFunction(net,classifier_name); % make sure to update Params.NeuralNetFunction in GetParams with the new name of the classifier
% delete(gcp)

%%% STARTING THE EXPERIMENT %%%%
clear
clc
cd('/home/ucsf/Projects/bci')
ExperimentStart('Robot3DArrow','bravo1',4,1,0)
%  ExperimentStart('RobotLateralR2G','bravo1',4,1,0)
% ExperimentStart('Robot3D','bravo1',4,1,0)
%  %ExperimentStart('Robot3DArrow','bravo1',4,1,0)
% %  ExperimentStart('RobotR2GModeSwitch','bravo1',4,1,0)
%
% %  



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
