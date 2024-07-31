%% CREATING A MLP FOR B3, NO POOLING

clc;clear
close all



% IMAGINED 
clc;clear
root_path = '/home/ucsf/Data/bravo1/20240731/Robot3DArrow';
foldernames = {'140330', '140715', '141041', '141406', '141727', '142055'};
cd(root_path)


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
        %kinax = length(features)-20:length(features);
        kinax = find(TrialData.TaskState==3);
        temp = cell2mat(features(kinax));

        % get delta, beta and hG removing bad channels 
        temp = temp([257:512 1025:1280 1537:1792],:);        
        bad_ch = [108 113 118];
        good_ch = ones(size(temp,1),1);
        for ii=1:length(bad_ch)
            %bad_ch_tmp = bad_ch(ii)*[1 2 3];
            bad_ch_tmp = bad_ch(ii)+(256*[0 1 2]);
            good_ch(bad_ch_tmp)=0;
        end
        temp = temp(logical(good_ch),:);

        % 2-norm
        for ii=1:size(temp,2)
            temp(:,ii) = temp(:,ii)./norm(temp(:,ii));
        end

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


% ONLINE DATA AS WELL
root_path = '/home/ucsf/Data/bravo1/20240731/Robot3DArrow';
foldernames = {'142859', '143131', '143357'};
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
            %bad_ch_tmp = bad_ch(ii)*[1 2 3];
             bad_ch_tmp = bad_ch(ii)+(256*[0 1 2]);
            good_ch(bad_ch_tmp)=0;
        end
        temp = temp(logical(good_ch),:);

        % 2-norm
        for ii=1:size(temp,2)
            temp(:,ii) = temp(:,ii)./norm(temp(:,ii));
        end

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



% ROBOT DATA AS WELL 
root_path = '/home/ucsf/Data/bravo1/20240731/RealRobotBatch';
foldernames = {'144519'};
cd(root_path)

for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = find(TrialData.TaskState==3);
        features =  features(kinax);
        %l = min(length(features),8);
        %features=features(1:l);
        temp = cell2mat(features);

        % get delta, beta and hG removing bad channels
        temp = temp([257:512 1025:1280 1537:1792],:);
        bad_ch = [108 113 118];
        good_ch = ones(size(temp,1),1);
        for ii=1:length(bad_ch)
            bad_ch_tmp = bad_ch(ii)*[1 2 3];
            good_ch(bad_ch_tmp)=0;
        end
        temp = temp(logical(good_ch),:);

        % 2-norm
        for ii=1:size(temp,2)
            temp(:,ii) = temp(:,ii)./norm(temp(:,ii));
        end

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
% combing both onlien plus offline
idx=1;
condn_data{1}=[D1(idx:end,:) ]'; % right thumb
condn_data{2}= [D2(idx:end,:)]'; % both feet
condn_data{3}=[D3(idx:end,:)]'; % left pinch
condn_data{4}=[D4(idx:end,:)]'; % head
condn_data{5}=[D5(idx:end,:)]'; % lips
condn_data{6}=[D6(idx:end,:)]'; % tong
condn_data{7}=[D7(idx:end,:)]'; % both hands

A = condn_data{1};
B = condn_data{2};
C = condn_data{3};
D = condn_data{4};
E= condn_data{5};
F= condn_data{6};
G= condn_data{7};


clear N
N = [A' B' C' D' E' F' G'];
T1 = [ones(size(A,1),1);2*ones(size(B,1),1);3*ones(size(C,1),1);4*ones(size(D,1),1);...
    5*ones(size(E,1),1);6*ones(size(F,1),1);7*ones(size(G,1),1)]
T = zeros(size(T1,1),7);
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
[aa bb]=find(T1==7);[aa(1) aa(end)]
T(aa(1):aa(end),7)=1;

%%%%%%%% CODE TO UPDATE THE PNP DECODER %%%%%%%
% cd('/home/ucsf/Projects/bci/clicker')
% % load net_B3B1_PnP_trfLearn_patternet_0717 % this is the original PnP... if updating and saving, use the name below
% % net = net_B3B1_PnP_trfLearn_patternet;
% % net = train(net,N,T','UseParallel','no');
% % genFunction(net,'net_B3B1_PnP_trfLearn_patternet_Update1')
% 
% % some bugs in codes above, could not get the val and test
% % modified version below
% load net_B3B1_PnP_trfLearn_patternet_0717 % this is the original PnP net net = net_B3B1_PnP_trfLearn_patternet; 
% net = net_B3B1_PnP_trfLearn_patternet;
% trainRatio = 0.7; 
% valRatio = 0.15; 
% testRatio = 0.15; 
% [trainInd,valInd,testInd]=dividerand(length(T),trainRatio,valRatio,testRatio);
% 
% net.divideFcn='divideind';
% net.divideParam.trainInd=trainInd;
% net.divideParam.valInd=valInd;
% net.divideParam.testInd=testInd;
% net = train(net, N, T', 'UseParallel', 'no'); 
% genFunction(net, 'net_B3B1_PnP_trfLearn_patternet_20240731_Update1') 


% backup
% load net_B3B1_PnP_trfLearn_patternet_0717 % this is the original PnP net net = net_B3B1_PnP_trfLearn_patternet; 
% net = net_B3B1_PnP_trfLearn_patternet;
% 
% 
% % Split dataset 
% trainRatio = 0.7; 
% valRatio = 0.15; 
% testRatio = 0.15; 
% totalLength=length(T);
% 
% % shuffle idx
% randnIdx = randperm(totalLength);
% 
% trainSize=round(totalLength*trainRatio);
% valSize=round(totalLength*valRatio);
% testSize=totalLength - trainSize - valSize;
% 
% % set idx of different set manually
% trainInd=randnIdx(1:trainSize);
% valInd=randnIdx(trainSize+1:trainSize+valSize);
% testInd=randnIdx(trainSize+valSize+1:end);
% 
% 
% % Set params 
% net.divideFcn = 'divideind'; 
% net.divideParam.trainInd = trainInd; 
% net.divideParam.valInd = valInd; 
% net.divideParam.testInd = testInd; 
% net = train(net, N, T', 'UseParallel', 'no'); 
% genFunction(net, 'net_B3B1_PnP_trfLearn_patternet_20240731_Update1') 

% 
% 
% 
% 
% 
% 
% 



% 






% net_mlp_pnp_update1 = net;
% save net_mlp_pnp_update1 net_mlp_pnp_update1
%%%%%%%%%%%%% END SECTION %%%%%%%%%%%%

% %%%% CODE TO TRAIN A NEURAL NETWORK FROM SCRATCH
clear net
net = patternnet([120]);
net.performParam.regularization=0.2;
net.divideParam.trainRatio=0.80;
net.divideParam.valRatio=0.10;
net.divideParam.testRatio=0.1;
net = train(net,N,T','UseParallel','no');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_7Dir_B1_20240731_CL3_NoPooling')
save net net
% %%%%%%%%%%%%%%%%%%%%%%% END SECTION %%%%%
% 
% 
% 

cd('/home/ucsf/Projects/bci')
% clear;clc
% ExperimentStart('Robot3DArrow','bravo1',4,1,0)
%ExperimentStart('Robot3D','Bravo3',4,1,0)


