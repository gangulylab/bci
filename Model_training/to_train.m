clear;clc
close all

root='E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\';

folders = {'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200911\CenterOut\145330\Imagined\',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200911\CenterOut\145855\Imagined\',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200911\CenterOut\150113\Imagined',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200911\CenterOut\150229\Imagined',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200911\CenterOut\150556\Imagined',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200911\CenterOut\153125\Imagined',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200911\CenterOut\155038\Imagined',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200918\CenterOut\145458\Imagined',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200918\CenterOut\151515\Imagined',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200918\CenterOut\154012\Imagined',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200918\CenterOut\155929\Imagined',...
    'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20200925\CenterOut\145022\Imagined'};



% collecting all the training data
files=[];
for i=1:length(folders)
    cd(folders{i})
    temp = findfiles('',pwd);
    files=[files;temp'];
end




% load the data for each target
D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
D7=[];
D8=[];
for i=1:length(files)
    disp(i)
    load(files{i});
    features  = TrialData.NeuralFeatures;
    kinax = length(features)+[-20:0];
    temp = cell2mat(features(kinax));
    temp = temp(129:end,:);
    for j=1:size(temp,1)
        temp(j,:) = smooth(temp(j,:),3);
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
    elseif TrialData.TargetID == 8
        D8 = [D8 temp];
    end
end

condn_data{1}=D1'; % right arm
condn_data{2}=D2'; % right thigh
condn_data{3}=D3'; % both feet
condn_data{4}=D4'; % left thigh
condn_data{5}=D5'; % left arm
condn_data{6}=D6'; % left shoulder
condn_data{7}=D7'; % head
condn_data{8}=D8'; % right shoulder


% only the 4 cardinal directions
clear condn_data
condn_data{1}=D1'; % right hand
condn_data{2}=D3'; % both feet
condn_data{3}=D5'; % left hand
condn_data{4}=D7'; % head



% TACKING ON THE ONLINE DISCRETE DATA SESSION 1
filepath = 'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20201002\RadialTask\';
% discrete arrow
folders = {'145058','145910','150443','151528','152108'};
% pixelated
%folders1 = {'153238','153653','154847'};
folders1 = '';
% together
folders =[ folders folders1];
files=[];
for i=1:length(folders)
    fullpath = [filepath folders{i} '\BCI_Fixed\'];
    files = [files findfiles('', fullpath)];
end

% get the trials

% load the data for each target
D11=[];
D22=[];
D33=[];
D44=[];
for i=1:length(files)
    disp(i)
    load(files{i});
    features  = TrialData.NeuralFeatures;
    kinax = [4:length(features)];
    %kinax = length(features)+[-10:0];
    temp = cell2mat(features(kinax));
    temp = temp(129:end,:);
    for j=1:size(temp,1)
        temp(j,:) = smooth(temp(j,:),3);
    end
    if TrialData.TargetID == 1
        D11 = [D11 temp];
    elseif TrialData.TargetID == 2
        D22 = [D22 temp];
    elseif TrialData.TargetID == 3
        D33 = [D33 temp];
    elseif TrialData.TargetID == 4
        D44 = [D44 temp];
    end
end
% % only the 4 cardinal directions
% condn_data{1}=D11'; % right hand
% condn_data{2}=D22'; % both feet
% condn_data{3}=D33'; % left hand
% condn_data{4}=D44'; % head


% TACKING ON THE ONLINE DISCRETE DATA SESSION 2 ARROW TASK
filepath = 'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20201009\RadialTask\';
% discrete arrow
folders = {'134342','134836','135450','135947','142059'};
folders1 = '';
% together
folders =[ folders folders1];
files=[];
for i=1:length(folders)
    fullpath = [filepath folders{i} '\BCI_Fixed\'];
    files = [files findfiles('', fullpath)];
end

% get the trials

% load the data for each target
D2_1=[];
D2_2=[];
D2_3=[];
D2_4=[];
for i=1:length(files)
    disp(i)
    load(files{i});
    features  = TrialData.SmoothedNeuralFeatures;
    kinax = TrialData.TaskState;
    if kinax(1) == 0
        kinax=kinax(2:end);
    end
    kinax = find(kinax==3);
    temp = cell2mat(features(kinax));
    temp = temp(129:end,:);
    if TrialData.TargetID == 1
        D2_1 = [D2_1 temp];
    elseif TrialData.TargetID == 2
        D2_2 = [D2_2 temp];
    elseif TrialData.TargetID == 3
        D2_3 = [D2_3 temp];
    elseif TrialData.TargetID == 4
        D2_4 = [D2_4 temp];
    end
end

% TACKING ON THE ONLINE DISCRETE DATA SESSION 3 ARROW TASK
filepath = 'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20201016\DiscreteArrow\';
% discrete arrow
folders = {'142924','144138','144924'};
folders1 = '';
% together
folders =[ folders folders1];
files=[];
for i=1:length(folders)
    fullpath = [filepath folders{i} '\BCI_Fixed\'];
    files = [files findfiles('', fullpath)];
end

% get the trials

% load the data for each target
D3_1=[];
D3_2=[];
D3_3=[];
D3_4=[];
for i=1:length(files)
    disp(i)
    load(files{i});
    features  = TrialData.SmoothedNeuralFeatures;
    kinax = TrialData.TaskState;
    if kinax(1) == 0
        kinax=kinax(2:end);
    end
    kinax = find(kinax==3);
    temp = cell2mat(features(kinax));
    temp = temp(129:end,:);
    if TrialData.TargetID == 1
        D3_1 = [D3_1 temp];
    elseif TrialData.TargetID == 2
        D3_2 = [D3_2 temp];
    elseif TrialData.TargetID == 3
        D3_3 = [D3_3 temp];
    elseif TrialData.TargetID == 4
        D3_4 = [D3_4 temp];
    end
end



% TACKING ON THE ONLINE DISCRETE DATA SESSION 5 ARROW TASK
filepath = 'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20201030\DiscreteArrow\';
% discrete arrow
folders = {'133902','134541','134838','142426'};
folders1 = '';
% together
folders =[ folders folders1];
files=[];
for i=1:length(folders)
    fullpath = [filepath folders{i} '\BCI_Fixed\'];
    files = [files findfiles('', fullpath)];
end

% load the data for each target
D5_1=[];
D5_2=[];
D5_3=[];
D5_4=[];
for i=1:length(files)
    disp(i)
    load(files{i});
    features  = TrialData.SmoothedNeuralFeatures;
    kinax = TrialData.TaskState;
    if kinax(1) == 0
        kinax=kinax(2:end);
    end
    kinax = find(kinax==3);
    temp = cell2mat(features(kinax));
    temp = temp(129:end,:);
    if TrialData.TargetID == 1
        D5_1 = [D5_1 temp];
    elseif TrialData.TargetID == 2
        D5_2 = [D5_2 temp];
    elseif TrialData.TargetID == 3
        D5_3 = [D5_3 temp];
    elseif TrialData.TargetID == 4
        D5_4 = [D5_4 temp];
    end
end

clear condn_data
% combing both onlien plus offline
idx=641;
condn_data{1}=[D1(idx:end,:) D11(idx:end,:) D2_1(idx:end,:) D3_1(idx:end,:) D5_1(idx:end,:) ]'; % right hand
condn_data{2}= [D3(idx:end,:) D22(idx:end,:) D2_2(idx:end,:) D3_2(idx:end,:) D5_2(idx:end,:)]'; % both feet and right foot here
condn_data{3}=[D5(idx:end,:) D33(idx:end,:) D2_3(idx:end,:) D3_3(idx:end,:) D5_3(idx:end,:)]'; % left hand
condn_data{4}=[D7(idx:end,:) D44(idx:end,:) D2_4(idx:end,:) D3_4(idx:end,:) D5_4(idx:end,:)]'; % head
% 
% 
% % looking at day specific differences per target
% condn_data{1}=[D1(641:768,:) ]'; % day 1
% condn_data{2}= [D11(641:768,:) ]'; % day 2
% condn_data{3}=[D2_1(641:768,:) ]'; % day 3
% condn_data{4}=[D3_1(641:768,:) ]'; % day 4



% MUTLIC-CLASS CLASSIFIERS AND ASSIGNING MAX SCORES
% partition the data
Conf_Matrix_Overall=[];
for iter=1:1
    condn_data_train={};
    condn_data_test={};
    for i=1:length(condn_data)
        temp = condn_data{i};
        l = round(1*size(temp,1));
        idx = randperm(size(temp,1),l);
        % setting aside training trials
        condn_data_train{i} = temp(idx,:);
        % setting aside testing trials
        I = ones(size(temp,1),1);
        I(idx)=0;
        I=logical(I);
        condn_data_test{i} = temp(I,:);
    end
    
    % build the classifiers
    svm_model={};D=[];
    for i=1:length(condn_data_train)
        disp(i)
        A = [condn_data_train{i}];
        A=[A ];
        for j=i:length(condn_data_train)
            if i==j
                svm_model{i,j}=0;
            else
                B = [condn_data_train{j}];
                B=[B ];
                [res_acc, model,pval] = svm_linear(A',B',1,1);
                %[model] = svm_nonlinear(A',B',1);                
                %no pruning 
                svm_model{i,j} = mean(model,1);
                svm_model{j,i} = -mean(model,1);                
                % pruning weights
%                 temp_wts = mean(model,1);
%                 [aa bb]=sort(abs(temp_wts));
%                 len = round(0.5*length(temp_wts));
%                 temp_wts(bb(1:len))=0;
%                 svm_model{i,j} = temp_wts;
%                 svm_model{j,i} = -temp_wts;
            end
        end
    end
    
    
%     % non-linear
%     clear model
%     model = svm_model;
%     cd('E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker')
%     save clicker_svm_mdl_Days123_5Partial_OnlyArrow_RBF model -v7.3
    
    % collate the condition specific classifiers
    models=[];
    for i=1:size(svm_model,1)
        temp=[];
        for j=1:size(svm_model,2)
            if length(svm_model{i,j}) > 1
                temp=[temp;svm_model{i,j}];
            end
        end
        models(i,:,:)=temp;
    end
    model = models;
    cd('E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker')
    save clicker_svm_mdl_Day0_hG_4Dir model -v7.3
    
    % now test them based on max vote strategy
    % i.e. test each data point on whether it came from one of the 8 classes
    Conf_Matrix=[];
    num_trials=[];
    for i=1:length(condn_data_test)
        temp = condn_data_test{i};
        % smooth data moving average
        %         for k=1:size(temp,2)
        %             temp(:,k) = smooth(temp(:,k),5);
        %         end
        predicted_classes=zeros(1,length(condn_data_test));
        for j=1:size(temp,1)
            d1 = temp(j,:)*squeeze(models(1,:,:))';
            d2 = temp(j,:)*squeeze(models(2,:,:))';
            d3 = temp(j,:)*squeeze(models(3,:,:))';
            d4 = temp(j,:)*squeeze(models(4,:,:))';
            %d5 = temp(j,:)*squeeze(models(5,:,:))';
            %d6 = temp(j,:)*squeeze(models(6,:,:))';
            %d7 = temp(j,:)*squeeze(models(7,:,:))';
            %d8 = temp(j,:)*squeeze(models(8,:,:))';
            %Dec = [d1;d2;d3;d4;d5;d6;d7;d8];
            Dec = [d1;d2;d3;d4;];
            Dec_values=[sum(d1);sum(d2);sum(d3);sum(d4)];
            %Dec_values=[sum(d1);sum(d2);sum(d3);sum(d4);sum(d5);sum(d6);sum(d7);sum(d8)];
            Dec_thresh=sum(Dec'<0);
            
            %%on max vote strategy
            [aa bb]=max(Dec_thresh);
            if length(find(Dec_thresh==aa)) == 1
                predicted_classes(bb)=predicted_classes(bb)+1;
            else
                [aa1 bb1]=min(Dec_values);
                predicted_classes(bb1)=predicted_classes(bb1)+1;
            end
            
            %%on max distance from decision boundary
%             Dec_values(Dec_values>-2)=0;
%             [aa1 bb1]=min(Dec_values);
%             if aa1~=0
%                 decision = bb1;
%                 predicted_classes(bb1)=predicted_classes(bb1)+1;
%             end
        end
        Conf_Matrix(i,:) = predicted_classes./sum(predicted_classes);
        num_trials = [num_trials sum(predicted_classes)/j];
    end
    figure;imagesc(Conf_Matrix);colormap bone;caxis([0 .8]);
    Conf_Matrix;
    diag(Conf_Matrix);
    Conf_Matrix_Overall(iter,:,:)=Conf_Matrix;
end
Conf_Matrix = squeeze(mean(Conf_Matrix_Overall,1));
diag(Conf_Matrix)
figure;stem(ans,'LineWidth',1)

xticklabels({'Right hand','Right thigh','Both feet','Left thigh',...
    'Left arm','Left shoulder','Head','Right shoulder'})
yticklabels({'Right hand','Right thigh','Both feet','Left thigh',...
    'Left arm','Left shoulder','Head','Right shoulder'})
set(gcf,'Color','w')
set(gca,'FontSize',14)
hh=hline(1/4);
set(hh,'LineWidth',2)
set(gcf,'Color','w')
set(gca,'FontSize',14)
figure;imagesc(Conf_Matrix);colormap bone;caxis([0.05 .8]);
xticks([1 2 3 4])
yticks([1 2 3 4])
xticklabels({'Right arm','Left arm','Right leg','Left leg'})
yticklabels({'Right arm','Left arm','Right leg','Left leg'})
set(gcf,'Color','w')
set(gca,'FontSize',14)
title('Confusion Matrix')
colorbar



body_part={'Right hand','Right thigh','Both feet','Left thigh',...
    'Left hand','Left shoulder','Head','Right shoulder'};
figure;ha=tight_subplot(8,1);
for i=1:8
    axes(ha(i))
    stem(Conf_Matrix(i,:),'LineWidth',2);
    xlim([0.5 8.5])
    ylim([0 1])
    yticks ''
    if i~=8
        xticks ''
    else
        h = (gca);
        h.XTickLabel = body_part;
    end
    ylabel(body_part{i});
    set(gca,'FontSize',10)
    hline(1/8)
end
set(gcf,'Color','w')


% plotting models on top of each other to show channels
idx=2;
figure;
hold on
cols = parula(6);k=1;
temp_wts=[];
for i=1:128:768
    stem(squeeze(models(idx,2,i:i+128-1)),'Color',[cols(k,:) .5],'LineWidth',1);
    k=k+1;
    temp_wts = [temp_wts squeeze(models(idx,1,i:i+128-1))];
end

figure;
for i=1:128
    subplot(16,8,i)
    stem(abs(temp_wts(i,:)),'LineWidth',1);
    axis tight
    if i~=128
        xticks ''       
    end
     yticks ''
end
xticklabels({'delta','theta','alpha','beta','lg','hg'})


% test a model online 
filepath='E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20201030\DiscreteArrow\';
cd(filepath)
folders =  {'141625'};
files=[];
for i=1:length(folders)
    fullpath = [filepath folders{i} '\BCI_Fixed\'];
    files = [files findfiles('', fullpath)];
end

cd('E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker')
%load('clicker_svm_mdl_Day3Smoothed.mat')
%load 'clicker_svm_mdl_Days0123_5Partial_OnlyArrow.mat'
load clicker_svm_mdl_Days123_5Partial_OnlyArrow_RBF.mat

ax=[];
for i=1:length(files)
    clear TrialData
    load(files{i})    
    %TrialData.TargetID
    idx = find(TrialData.TaskState==3);
    feat = TrialData.SmoothedNeuralFeatures;
    feat=feat(idx);
    feat = cell2mat(feat);
    ClickState=[];
    for j=1:size(feat,2)
        [Click_Decision,distance_boundary] = multistate_discrete_local_RBF(feat(:,j),model,1);
        ClickState = [ClickState Click_Decision];
    end
    a1=sum(TrialData.ClickerState == TrialData.TargetID);
    a2=sum(ClickState == TrialData.TargetID) ;   
    ax= [ax -a1+a2];
    disp([TrialData.TargetID -a1+a2])    
end
sum(ax>0)

