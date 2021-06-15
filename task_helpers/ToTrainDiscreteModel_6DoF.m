%% Creating MLP classifier
clc;clear
% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20210611/Robot3DArrow';
% enter the folder names for the Task. These can be increased as more data
% is collected. For exaple: 




foldernames = {'133803','134610', '135244', '135725', '140307', '141405'};

% foldernames = {'133417', '134037', '134340'};


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
        kinax = [ find(TrialData.TaskState==3)];
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
        end
    end
end




% FOR ONLINE DATA
%  foldernames = {'135314','140629'};

for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = [ find(TrialData.TaskState==3)];
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

A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};
E = condn_data{5};
F = condn_data{6};


clear N
N = [A' B' C' D' E' F' ];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    5*ones(size(E,1),1);6*ones(size(F,1),1)];

T = zeros(size(T1,1),6);
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


%%%%%%% CODE SNIPPET FOR TRAINING A MODEL FROM SCRATCH %%%%%
% training a simple MLP
% IMPORTANT, CLICK THE CONFUSION MATRIX BUTTON IN GUI TO VERIFY THAT THE
% TEST VALIDATION DOESN'T HAVE NaNs AND THAT PERFORMANCE IS REASONABLE
net = patternnet([128 128 128 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');
%%%%%%% END CODE SNIPPET %%%%%%%%%%%%%%%


%%%%%%% UPDATING THE WEIGHTS OF A PRETRAINED MODEL %%%%%%%%
cd('/home/ucsf/Projects/bci/clicker')
load net % pretrained model
net = train(net,N,T'); % updates the weights


%%% SAVING %%%%%%%
cd('/home/ucsf/Projects/bci/clicker')
% classifier name
classifier_name = 'MLP_Imag_Actions_0611_D';
genFunction(net,classifier_name); % make sure to update GetParams





% to restart exp run following lines
% clear
% clc
% cd('/home/ucsf/Projects/bci')
%  ExperimentStart('Robot3DArrow','bravo1',4,1,0)
 
 %ExperimentStart('Robot3D','bravo1',4,1,0)
% 
% ExperimentStart('RobotReachStop','bravo1',4,1,0)



%% TO CREATE CLASSIFIER BASED ON POOLING DATA


clc;clear
% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20210611/Robot3DArrow';
% enter the folder names for the Task. These can be increased as more data
% is collected. For exaple: 

foldernames = {'110732', '111455','111935', '112357'};

% foldernames = {'133417', '134037', '134340'};

cd(root_path)

%FOR IMAGINED MOVEMENT DATA, 
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
        end
    end
end




% FOR ONLINE DATA
%  foldernames = {'135314','140629'};

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
        end
    end
end



clear condn_data
% combing delta beta and high gamma
idx=[1:size(D1,1)];
condn_data{1}=[D1(idx,:) ]'; 
condn_data{2}= [D2(idx,:)]'; 
condn_data{3}=[D3(idx,:)]'; 
condn_data{4}=[D4(idx,:)]'; 
condn_data{5}=[D5(idx,:)]'; 
condn_data{6}=[D6(idx,:)]'; 

A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};
E = condn_data{5};
F = condn_data{6};


clear N
N = [A' B' C' D' E' F' ];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    5*ones(size(E,1),1);6*ones(size(F,1),1)];

T = zeros(size(T1,1),6);
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


%%%%%%% CODE SNIPPET FOR TRAINING A MODEL FROM SCRATCH %%%%%
% training a simple MLP
% IMPORTANT, CLICK THE CONFUSION MATRIX BUTTON IN GUI TO VERIFY THAT THE
% TEST VALIDATION DOESN'T HAVE NaNs AND THAT PERFORMANCE IS REASONABLE
clear net
net = patternnet([96 96 96]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');



%%%%%%% SAVING THE MODEL %%%%%%%%
cd('/home/ucsf/Projects/bci/clicker')

% classifier name
classifier_name = 'MLP_Imag_Actions_0611_PA';
genFunction(net,classifier_name); % make sure to update GetParams


% % to restart exp run following lines
% clear
% clc
% cd('/home/ucsf/Projects/bci')
%  ExperimentStart('Robot3DArrow','bravo1',4,1,0)
 


