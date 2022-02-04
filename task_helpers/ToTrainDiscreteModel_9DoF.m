
%% CREATING A  CLASSIFIER (either new or batch update, select apporp folder)
% USES POOLED DATA

% if creating new classifier use only imagined mvmt data folders
% if doing batch updae, use both imagined and fixed folders

clc;clear

% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20220204/Robot3DArrow';

% enter the folder names for the Task. These can be increased as more data
% is collected. For example: 

foldernames = {'133312', '134225', '134655','135121'};

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
        elseif TrialData.TargetID == 8
            D8 = [D8 temp];
        elseif TrialData.TargetID == 9
            D9 = [D9 temp];
        end
    end
end

% FIXED CONTROL
 foldernames = {'140159','140607','141121','141718'};

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
        elseif TrialData.TargetID == 8
            D8 = [D8 temp];
        elseif TrialData.TargetID == 9
            D9 = [D9 temp];
        end
    end
end

size(D7)
root_path = '/home/ucsf/Data/bravo1/20210813/RealRobotBatch';
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


%%%%% CODE SNIPPET FOR TRAINING AND SAVING THE DECODER %%%%%
cd('/home/ucsf/Projects/bci/clicker')
clear net
net = patternnet([64 64 64]);
net.performParam.regularization=0.2;
net.divideParam.trainRatio=0.8;
net.divideParam.valRatio=0.1;
net.divideParam.testRatio=0.1;
net = train(net,N,T','useParallel','yes');
classifier_name = 'MLP_PreTrained_9DoF_02042022_PM2'; % enter the name
genFunction(net,classifier_name); % make sure to update Params.NeuralNetFunction in GetParams with the new name of the classifier
delete(gcp)

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
  ExperimentStart('Robot3DArrow','bravo1',4,1,0)
%  ExperimentStart('RobotLateralR2G','bravo1',4,1,0)
% ExperimentStart('Robot3D','bravo1',4,1,0)
%  %ExperimentStart('Robot3DArrow','bravo1',4,1,0)
% %  ExperimentStart('RobotR2GModeSwitch','bravo1',4,1,0)
%  
% %  
