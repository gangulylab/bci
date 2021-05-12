%% Creating MLP classifier
clc;clear
% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20210512/Robot3DArrow';
% enter the folder names for the Task. These can be increased as more data
% is collected. For exaple: 

foldernames = {'110329', '110946', '111317', '111724', '112038', '112501', '113856', '114424', '114722', '134708', '135326', '140043', '140148', '140730'};
% foldernames = {'134708', '135326', '140043', '140148', '140730'};
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


% training a simple MLP
% IMPORTANT, CLICK THE CONFUSION MATRIX BUTTON IN GUI TO VERIFY THAT THE
% TEST VALIDATION DOESN'T HAVE NaNs AND THAT PERFORMANCE IS REASONABLE
net = patternnet([128 128 128 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');
cd('/home/ucsf/Projects/bci/clicker')
% classifier name
classifier_name = 'MLP_Lips_RtThumb_LtThumb_RtMiddle_Tongue_LfMiddle_SquezeHands_2';
% generates the MLP as function in the clikcer folder. Make sure to update
% the MLP classifier name in GetParams.m in the Neural Network section. 
genFunction(net,classifier_name);



