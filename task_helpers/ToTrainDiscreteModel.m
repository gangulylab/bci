
%% CREATING A MLP FOR B3, NO POOLING

clc;clear
close all


% get the imagined data from cursor center out 
clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210203/CenterOut';
foldernames = {'145654','150544'};
cd(root_path)


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
        kinax = length(features)-20:length(features);
        temp = cell2mat(features(kinax));

        % get delta, beta and hG removing bad channels 
        temp = temp([257:512 1025:1280 1537:1792],:);        
        bad_ch = [108 113 118];
        good_ch = ones(size(temp,1),1);
        for ii=1:length(bad_ch)
            bad_ch_tmp = bad_ch(ii)*[1 2 3];
            good_ch(bad_ch_tmp)=0;
        end
        temp = temp(logical(good_ch),:);

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


% online data
root_path = '/home/ucsf/Data/Bravo2/20210203/DiscreteArrow';
foldernames = {'154256'};
cd(root_path)

for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = find(TrialData.TaskState==3);
        temp = cell2mat(features(kinax));

        % get delta, beta and hG removing bad channels
        temp = temp([257:512 1025:1280 1537:1792],:);
        bad_ch = [108 113 118];
        good_ch = ones(size(temp,1),1);
        for ii=1:length(bad_ch)
            bad_ch_tmp = bad_ch(ii)*[1 2 3];
            good_ch(bad_ch_tmp)=0;
        end
        temp = temp(logical(good_ch),:);

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
condn_data{4}=[D4(idx:end,:)]'; % head


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};

clear N
N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    ];
T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;

% code to train a neural network
net = patternnet([256 256 256 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T','UseGPU','yes','UseParallel','yes');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_4Dir_Imagined_20210210_Days12_Online_AllFeat')



%%



clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210203/CenterOut';
foldernames = {'145654','150544'};
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
        kinax = length(features)-20:length(features);
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

%Day2
root_path = '/home/ucsf/Data/Bravo2/20210210/CenterOut';
foldernames = {'143225'};
cd(root_path)

for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'Imagined');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = length(features)-20:length(features);
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

% Day 1 online data
root_path = '/home/ucsf/Data/Bravo2/20210203/DiscreteArrow';
foldernames = {'154256'};
cd(root_path)

for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = find(TrialData.TaskState==3);
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



%Day2 online data
root_path = '/home/ucsf/Data/Bravo2/20210210/DiscreteArrow';
foldernames = {'151341','151744','154405'};
cd(root_path)

for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = find(TrialData.TaskState==3);
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
condn_data{4}=[D4(idx:end,:)]'; % head


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};

clear N
N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    ];
T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;

% code to train a neural network
net = patternnet([256 256 256 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T','UseGPU','yes','UseParallel','yes');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_4Dir_Imagined_20210210_Days12_Online_AllFeat')

%% DAY 3 REINIT

clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210217/CenterOut';
foldernames = {'142905'};
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
        kinax = length(features)-25:length(features);
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
condn_data{4}=[D4(idx:end,:)]'; % head


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};

clear N
N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    ];
T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;


% code to train a neural network
net = patternnet([256 256 256 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_4Dir_Imagined_20210217_Day3_AllFeat')


% DAY 3 ONLINE DATA, LOTS OF SUCCESSFUL TRIALS

clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210217/DiscreteArrow';
D = dir(root_path);
j=1;
for i=3:length(D)
    foldernames{j} = D(i).name;
    j=j+1;
end    
foldernames

cd(root_path)
D1=[];
D2=[];
D3=[];
D4=[];

for k=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{k},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = find(TrialData.TaskState==3);
            idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax);
            kinax = kinax(idx);
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
end


clear condn_data
% combing both onlien plus offline
idx=1;
condn_data{1}=[D1(idx:end,:) ]'; % right hand
condn_data{2}= [D2(idx:end,:)]'; % both feet
condn_data{3}=[D3(idx:end,:)]'; % left hand
condn_data{4}=[D4(idx:end,:)]'; % head


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};

clear N
N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    ];
T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;


% code to train a neural network
net = patternnet([256 256 256 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_4Dir_Online_20210217_Day3_AllFeat')



%% %% DAY 4 REINIT

clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210224/CenterOut';
foldernames = {'143622','161431'};
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
        kinax = length(features)-25:length(features);
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



%Day4 online data
root_path = '/home/ucsf/Data/Bravo2/20210224/DiscreteArrow/';
cd(root_path)
clear foldernames
D = dir(root_path);
k=1;
for i=3:length(D)
    foldernames{k} = D(i).name;
    k=k+1;
end
    
cd(root_path)

for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = find(TrialData.TaskState==3);
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
condn_data{4}=[D4(idx:end,:)]'; % head


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};

clear N
N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    ];
T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;


% code to train a neural network
net = patternnet([256 256 256 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_4Dir_Imagined_Online_20210224_Day4_AllFeat')

%% TRAIN USING ONLINE DATA ACROSS DAYS
%%%% RETRAINING USING ONLINE DATA

% DAY 3 ONLINE DATA, LOTS OF SUCCESSFUL TRIALS

clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210217/DiscreteArrow';
D = dir(root_path);
j=1;
for i=3:length(D)
    foldernames{j} = D(i).name;
    j=j+1;
end    
foldernames

cd(root_path)
D1=[];
D2=[];
D3=[];
D4=[];

for k=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{k},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = find(TrialData.TaskState==3);
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax);
            idx = (length(kinax)+1 - TrialData.Params.ClickCounter-5):length(kinax) ;
            idx=idx(idx>0);
            kinax = kinax(idx);
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
end

%DAY 4
root_path = '/home/ucsf/Data/Bravo2/20210224/DiscreteArrow/';
cd(root_path)
clear foldernames
D = dir(root_path);
k=1;
for i=3:length(D)
    foldernames{k} = D(i).name;
    k=k+1;
end
    
cd(root_path)

for k=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{k},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = find(TrialData.TaskState==3);
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax);
            idx = (length(kinax)+1 - TrialData.Params.ClickCounter-5):length(kinax) ;
            idx=idx(idx>0);
            kinax = kinax(idx);
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
end


% day 5

root_path = '/home/ucsf/Data/Bravo2/20210303/DiscreteArrow/';
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
        if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = find(TrialData.TaskState==3);
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax) ;
            idx = (length(kinax)+1 - TrialData.Params.ClickCounter-5):length(kinax) ;
            idx=idx(idx>0);
            kinax = kinax(idx);
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
end



% DAY 6

root_path = '/home/ucsf/Data/Bravo2/20210310/DiscreteArrow/';
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
        if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = find(TrialData.TaskState==3);
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax) ;
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax) ;
            idx = (length(kinax)+1 - TrialData.Params.ClickCounter-5):length(kinax) ;
            idx=idx(idx>0);
            kinax = kinax(idx);
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
end




% DAY 7

root_path = '/home/ucsf/Data/Bravo2/20210317/DiscreteArrow/';
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
        if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = find(TrialData.TaskState==3);
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax) ;
            idx = (length(kinax)+1 - TrialData.Params.ClickCounter-5):length(kinax) ;
            idx=idx(idx>0);
            kinax = kinax(idx);
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
end




% code to train a neural network

clear condn_data
% combing both onlien plus offline
idx=1;
condn_data{1}=[D1(idx:end,:) ]'; % right hand
condn_data{2}= [D2(idx:end,:)]'; % both feet
condn_data{3}=[D3(idx:end,:)]'; % left hand
condn_data{4}=[D4(idx:end,:)]'; % head


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};

clear N
N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    ];
T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;
net = patternnet([256 256 256 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');
cd('/home/ucsf/Projects/bci/clicker')
%genFunction(net,'MLP_4Dir_Online_20210303_Day3to7_AllFeat_Online')


%% DAY 5 REINIT


clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210303/CenterOut';
foldernames = {'142521'};
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
        kinax = length(features)-25:length(features);
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



root_path = '/home/ucsf/Data/Bravo2/20210303/DiscreteArrow/';
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
        if TrialData.SelectedTargetID == TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = find(TrialData.TaskState==3);
            %idx = length(kinax)+1 - TrialData.Params.ClickCounter:length(kinax) ;
            idx = length(kinax)+-2 - TrialData.Params.ClickCounter:length(kinax) ;
            idx=idx(idx>0);
            kinax = kinax(idx);
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
        
        if TrialData.TargetID==3 && TrialData.SelectedTargetID ~= TrialData.TargetID
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = find(TrialData.TaskState==3); 
            kinax = kinax(4:end);
            temp = cell2mat(features(kinax));
            temp = temp(129:end,:);
            D3 = [D3 temp];
        end
    end
end
% 
% %   all online data
% for i=1:length(foldernames)
%     folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
%     D=dir(folderpath);
%     for j=3:length(D)
%         filepath=fullfile(folderpath,D(j).name);
%         load(filepath)
%         features  = TrialData.SmoothedNeuralFeatures;
%         kinax = find(TrialData.TaskState==3);
%         temp = cell2mat(features(kinax));
%         temp = temp(129:end,:);
%         if TrialData.TargetID == 1
%             D1 = [D1 temp];
%         elseif TrialData.TargetID == 2
%             D2 = [D2 temp];
%         elseif TrialData.TargetID == 3
%             D3 = [D3 temp];
%         elseif TrialData.TargetID == 4
%             D4 = [D4 temp];
%         end
%     end
% end



clear condn_data
% combing both onlien plus offline
idx=1;
condn_data{1}=[D1(idx:end,:) ]'; % right hand
condn_data{2}= [D2(idx:end,:)]'; % both feet
condn_data{3}=[D3(idx:end,:)]'; % left hand
condn_data{4}=[D4(idx:end,:)]'; % head


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};

idx = randperm(size(C,1),220);
C = C(idx,:);

clear N
N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    ];
T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;
net = patternnet([256 256 256 ]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_4Dir_Online_20210303_Day5_Online')


%% DAY 6 10TH MARCH USING IMAGINED WORDS FOR 4 DIRECTIONS


clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210310/CenterOut';
foldernames = {'150625'};
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




%%% USING ONLINE DATA


root_path = '/home/ucsf/Data/Bravo2/20210310/DiscreteArrow/';
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
            kinax = find(TrialData.TaskState==3);
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


clear condn_data
% combing both onlien plus offline
idx=641;
condn_data{1}=[D1(idx:end,:) ]'; % right hand
condn_data{2}= [D2(idx:end,:)]'; % both feet
condn_data{3}=[D3(idx:end,:)]'; % left hand
condn_data{4}=[D4(idx:end,:)]'; % head


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};


clear N
N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    ];
T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;
net = patternnet([256 256 256]) ;
net.performParam.regularization=0.2;
net = train(net,N,T');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_4Dir_Words_Online_hG_20210303_Day6')

%% DAY 8 MAR 24, 4 DIR IMAGINED ACTIONS WITH INSTRUCTIONS TO SACCADE AND FIXATE



clc;clear
root_path = '/home/ucsf/Data/Bravo2/20210324/CenterOut';
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


% DATA FROM ONLINE TASKS
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


% PIXELATED CURSOR
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




clear condn_data
% combing both onlien plus offline
idx=1;
condn_data{1}=[D1(idx:end,:) ]'; % right hand
condn_data{2}= [D2(idx:end,:)]'; % both feet
condn_data{3}=[D3(idx:end,:)]'; % left hand
condn_data{4}=[D4(idx:end,:)]'; % head


A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};

parpool(4)

clear N
N = [A' B' C' D'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    ];
T = zeros(size(T1,1),4);
[aa bb]=find(T1==1);[aa(1) aa(end)]
T(aa(1):aa(end),1)=1;
[aa bb]=find(T1==2);[aa(1) aa(end)]
T(aa(1):aa(end),2)=1;
[aa bb]=find(T1==3);[aa(1) aa(end)]
T(aa(1):aa(end),3)=1;
[aa bb]=find(T1==4);[aa(1) aa(end)]
T(aa(1):aa(end),4)=1;
net = patternnet([128 128 128]) ;
net.performParam.regularization=0.2;
net = train(net,N,T','useParallel','yes');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_4Dir_Actions_AllOnline_20210331c')
delete(gcp)


