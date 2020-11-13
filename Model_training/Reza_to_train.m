clear;clc
close all

%% experiment info: TargetID: 1: right target;2:down target; 3: left target; 4: top target; 
expts = [];
expts(end+1).yymmdd = '20201030';
expts(end).Imagined_hhmmss = {'150244','150407','151019'};

%% collecting all the training data

% load the data for each target
D1=[];
D2=[];
D3=[];
D4=[];

for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.Imagined_hhmmss)
        display (['Session:',num2str(expt.Imagined_hhmmss{1,j})])
        % go through datafiles in fixed blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,'GangulyServer',...
            yymmdd,'CenterOut',expt.Imagined_hhmmss{1,j},'Imagined');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        for k=1:length(files)
            display (['Trial:',num2str(k)])
            load(fullfile(datadir,files(k).name));
            % Nikhilesh lines
            features  = TrialData.SmoothedNeuralFeatures;
            kinax = length(features)+[-25:0];
            temp = cell2mat(features(kinax));
            temp = temp(129:end,:);
            % only hG
            % X = X(769:end);
            
            
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

condn_data{1}=D1'; % right hand
condn_data{2}=D2'; % left leg
condn_data{3}=D3'; % left hand
condn_data{4}=D4'; % head

%% MUTLIC-CLASS CLASSIFIERS AND ASSIGNING MAX SCORES

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
    save clicker_svm_mdl_Day5 model -v7.3

end

delete(gcp)


