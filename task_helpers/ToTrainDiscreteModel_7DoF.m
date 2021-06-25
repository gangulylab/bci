%% Creating MLP classifier
clc;clear
% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20210616/Robot3DArrow';
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

foldernames = {'111251', '111821'};

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
%  ExperimentStart('Robot3DArrow','bravo1',4,1,0)
 
 %ExperimentStart('Robot3D','bravo1',4,1,0)
% 
% ExperimentStart('RobotReachStop','bravo1',4,1,0)

%% TO CREATE CLASSIFIER BASED ON POOLING DATA



clc;clear
% enter the root path from the Data folder
root_path = '/home/ucsf/Data/bravo1/20210625/Robot3DArrow';
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

% foldernames = {'133244', '133928','134357'};

foldernames = {'111134', '112108', '112805', '113645', '114239', '132902', '134133', '142139'};
cd(root_path)

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
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end




% FOR ONLINE DATA
%   foldernames = {'135435','135630','135830','140530','142530','142723'};
 
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
clear net
net = patternnet([96 96 96]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');



%%%%%%% SAVING THE MODEL %%%%%%%%
cd('/home/ucsf/Projects/bci/clicker')

% classifier name
classifier_name = 'MLP_Imag_Actions_0625_7DoF_PM2';
genFunction(net,classifier_name); % make sure to update GetParams
% 
% 
% % to restart exp run following lines
clear
clc
cd('/home/ucsf/Projects/bci')
% ExperimentStart('RobotStop','bravo1',4,1,0)
%ExperimentStart('Robot3D','bravo1',4,1,0)
 %ExperimentStart('Robot3DArrow','bravo1',4,1,0)
%  
