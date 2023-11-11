%% CREATING A MLP FOR B3, NO POOLING

clc;clear
close all



% IMAGINED 
clc;clear
root_path = '/home/ucsf/Data/Bravo3/20231110/Robot3DArrow';
foldernames = {'151210', '152025', '152352', '152755', '153205'};
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
root_path = '/home/ucsf/Data/Bravo3/20231110/Robot3DArrow';
foldernames = {'153829', '154338', '154618', '154906', '155654', '155935', '160204'};
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
        elseif TrialData.TargetID == 5
            D5 = [D5 temp];
        elseif TrialData.TargetID == 6
            D6 = [D6 temp];
        elseif TrialData.TargetID == 7
            D7 = [D7 temp];
        end
    end
end



% ROBOT DATA AS WELL --> FIRST 5 BINS OR 1.0S
% root_path = '/home/ucsf/Data/Bravo3/20230420/Robot3D';
% foldernames = {'113952', '114341', '114820', '115318', '120050', '120326'};
% cd(root_path)
% 
% for i=1:length(foldernames)
%     folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
%     D=dir(folderpath);
%     for j=3:length(D)
%         filepath=fullfile(folderpath,D(j).name);
%         load(filepath)
%         features  = TrialData.SmoothedNeuralFeatures;
%         kinax = find(TrialData.TaskState==3);
%         features =  features(kinax);
%         l = min(length(features),8);
%         features=features(1:l);
%         temp = cell2mat(features);
% 
%         % get delta, beta and hG removing bad channels
%         temp = temp([257:512 1025:1280 1537:1792],:);
%         bad_ch = [108 113 118];
%         good_ch = ones(size(temp,1),1);
%         for ii=1:length(bad_ch)
%             bad_ch_tmp = bad_ch(ii)*[1 2 3];
%             good_ch(bad_ch_tmp)=0;
%         end
%         temp = temp(logical(good_ch),:);
% 
%         if TrialData.TargetID == 1
%             D1 = [D1 temp];
%         elseif TrialData.TargetID == 2
%             D2 = [D2 temp];
%         elseif TrialData.TargetID == 3
%             D3 = [D3 temp];
%         elseif TrialData.TargetID == 4
%             D4 = [D4 temp];
%         elseif TrialData.TargetID == 5
%             D5 = [D5 temp];
%         elseif TrialData.TargetID == 6
%             D6 = [D6 temp];
%         elseif TrialData.TargetID == 7
%             D7 = [D7 temp];
%         end
%     end
% end




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
    5*ones(size(E,1),1);6*ones(size(F,1),1);7*ones(size(G,1),1);];
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
cd('/home/ucsf/Projects/bci/clicker')
% load net_mlp_pnp % this is the original PnP... if updating and saving, use the name below
% net = net_mlp_pnp;
% net = train(net,N,T','UseParallel','no');
% genFunction(net,'MLP_7Dir_B3_PnP_04042023_NoPooling_Update2')
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
genFunction(net,'MLP_7Dir_B3_20231110_CL3_NoPooling')
save net net
% %%%%%%%%%%%%%%%%%%%%%%% END SECTION %%%%%
% 
% 
% 

cd('/home/ucsf/Projects/bci')
% clear;clc
% ExperimentStart('Robot3DArrow','Bravo3',4,1,0)
%ExperimentStart('Robot3D','Bravo3',4,1,0)


