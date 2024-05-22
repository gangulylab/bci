%% CREATING A MLP FOR B3, NO POOLING

clc;clear
close all



% IMAGINED
clc;clear
root_path = '/home/ucsf/Data/Bravo1/20240522/HandImagined';
foldernames = {'103410', '104032', '104409', '105436', '105817', '110148', '110624'};
cd(root_path)


%SET BAD CHANNELS
Params.SetBadChannels = [];


D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
D7=[];
D8=[];
D9=[];
D10=[];
D11=[];
D12=[];
for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'Imagined');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name);
        load(filepath)
        features  = TrialData.SmoothedNeuralFeatures;
        kinax = TrialData.TaskState;

        %%%%% if screen update rate is different from decoder update rate

        if length(TrialData.Time) ~= length(TrialData.NeuralTime)
            time_exp = TrialData.Time;
            time_neural = TrialData.NeuralTime;
            time_exp = time_exp-time_exp(1);
            time_neural = time_neural-time_neural(1);
            time_exp = time_exp(find(kinax==3));
            kinax_neural=[];diff_times=[];
            for k=1:length(time_exp) % find the index of closed neural time
                [aa bb]= min(abs(time_neural - time_exp(k)));
                kinax_neural= [kinax_neural bb];
            end
            kinax_neural=unique(kinax_neural);
            temp = cell2mat(features(kinax_neural));
        else
            kinax=find(kinax==3);
            temp = cell2mat(features(kinax));
        end
        %%%%%%%%%%%%%%%%%%

        % get delta, beta and hG removing bad channels
        temp = temp([257:512 1025:1280 1537:1792],:);
        bad_ch = [228, 13, 14 Params.SetBadChannels];
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
        elseif TrialData.TargetID == 8
            D8 = [D8 temp];
        elseif TrialData.TargetID == 9
            D9 = [D9 temp];
        elseif TrialData.TargetID == 10
            D10 = [D10 temp];
        elseif TrialData.TargetID == 11
            D11 = [D11 temp];
        elseif TrialData.TargetID == 12
            D12 = [D12 temp];
        end
    end
end


% ONLINE DATA AS WELL
% <<<<<<< Updated upstream
% root_path = '/home/ucsf/Data/Bravo3/20231011/HandOnline';
% foldernames = {'150717', '151014', '151304', '151821', '152025', '152237','152647','152902','153129'};
% cd(root_path)
% 
% for i=1:length(foldernames)
%     folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
%     D=dir(folderpath);
%     for j=3:length(D)
%         filepath=fullfile(folderpath,D(j).name);
%         load(filepath)
%         features  = TrialData.SmoothedNeuralFeatures;
%         kinax = TrialData.TaskState;
% 
%         %%%%% if screen update rate is different from decoder update rate
% 
%         if length(TrialData.Time) ~= length(TrialData.NeuralTime)
%             time_exp = TrialData.Time;
%             time_neural = TrialData.NeuralTime;
%             time_exp = time_exp-time_exp(1);
%             time_neural = time_neural-time_neural(1);
%             time_exp = time_exp(find(kinax==3));
%             kinax_neural=[];diff_times=[];
%             for k=1:length(time_exp) % find the index of closed neural time
%                 [aa bb]= min(abs(time_neural - time_exp(k)));
%                 kinax_neural= [kinax_neural bb];
%             end
%             kinax_neural=unique(kinax_neural);
%             temp = cell2mat(features(kinax_neural));
%         else
%             kinax=find(kinax==3);
%             temp = cell2mat(features(kinax));
%         end
%         %%%%%%%%%%%%%%%%%%
% 
% 
%         % get delta, beta and hG removing bad channels
%         temp = temp([257:512 1025:1280 1537:1792],:);
%         bad_ch = [108 113 118 Params.SetBadChannels];
%         good_ch = ones(size(temp,1),1);
%         for ii=1:length(bad_ch)
%             %bad_ch_tmp = bad_ch(ii)*[1 2 3];
%             bad_ch_tmp = bad_ch(ii)+(256*[0 1 2]);
%             good_ch(bad_ch_tmp)=0;
%         end
%         temp = temp(logical(good_ch),:);
% 
%         % 2-norm
%         for ii=1:size(temp,2)
%             temp(:,ii) = temp(:,ii)./norm(temp(:,ii));
%         end
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
%         elseif TrialData.TargetID == 8
%             D8 = [D8 temp];
%         elseif TrialData.TargetID == 9
%             D9 = [D9 temp];
%         elseif TrialData.TargetID == 10
%             D10 = [D10 temp];
%         elseif TrialData.TargetID == 11
%             D11 = [D11 temp];
%         elseif TrialData.TargetID == 12
%             D12 = [D12 temp];
%         end
%     end
% end
% 
% =======
% root_path = '/home/ucsf/Data/Bravo1/20240522/HandOnline';
% foldernames = {};
% cd(root_path)
% 
% for i=1:length(foldernames)
%     folderpath = fullfile(root_path, foldernames{i},'BCI_Fixed');
%     D=dir(folderpath);
%     for j=3:length(D)
%         filepath=fullfile(folderpath,D(j).name);
%         load(filepath)
%         features  = TrialData.SmoothedNeuralFeatures;
%         kinax = TrialData.TaskState;
% 
%         %%%%% if screen update rate is different from decoder update rate
% 
%         if length(TrialData.Time) ~= length(TrialData.NeuralTime)
%             time_exp = TrialData.Time;
%             time_neural = TrialData.NeuralTime;
%             time_exp = time_exp-time_exp(1);
%             time_neural = time_neural-time_neural(1);
%             time_exp = time_exp(find(kinax==3));
%             kinax_neural=[];diff_times=[];
%             for k=1:length(time_exp) % find the index of closed neural time
%                 [aa bb]= min(abs(time_neural - time_exp(k)));
%                 kinax_neural= [kinax_neural bb];
%             end
%             kinax_neural=unique(kinax_neural);
%             temp = cell2mat(features(kinax_neural));
%         else
%             kinax=find(kinax==3);
%             temp = cell2mat(features(kinax));
%         end
%         %%%%%%%%%%%%%%%%%%
% 
% 
%         % get delta, beta and hG removing bad channels
%         temp = temp([257:512 1025:1280 1537:1792],:);
%         bad_ch = [228, 13, 14 Params.SetBadChannels];
%         good_ch = ones(size(temp,1),1);
%         for ii=1:length(bad_ch)
%             %bad_ch_tmp = bad_ch(ii)*[1 2 3];
%             bad_ch_tmp = bad_ch(ii)+(256*[0 1 2]);
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
%         elseif TrialData.TargetID == 8
%             D8 = [D8 temp];
%         elseif TrialData.TargetID == 9
%             D9 = [D9 temp];
%         elseif TrialData.TargetID == 10
%             D10 = [D10 temp];
%         elseif TrialData.TargetID == 11
%             D11 = [D11 temp];
%         elseif TrialData.TargetID == 12
%             D12 = [D12 temp];
%         end
%     end
% end
% 

clear condn_data
idx=1;
condn_data{1}=[D1(idx:end,:) ]';
condn_data{2}= [D2(idx:end,:)]';
condn_data{3}=[D3(idx:end,:)]';
condn_data{4}=[D4(idx:end,:)]';
condn_data{5}=[D5(idx:end,:)]';
condn_data{6}=[D6(idx:end,:)]';
condn_data{7}=[D7(idx:end,:)]';
condn_data{8}=[D8(idx:end,:)]';
condn_data{9}=[D9(idx:end,:)]';
condn_data{10}=[D10(idx:end,:)]';
condn_data{11}=[D11(idx:end,:)]';
condn_data{12}=[D12(idx:end,:)]';


N=[];
T1=[];
for i=1:length(condn_data)
    tmp=condn_data{i};
    N = [N tmp'];
    T1 = [T1;i*ones(size(tmp,1),1)];
end


T = zeros(size(T1,1),length(condn_data));
for i=1:length(condn_data)
    [aa bb]=find(T1==i);
    T(aa(1):aa(end),i)=1;
end


%%%% CODE TO TRAIN A NEURAL NETWORK FROM SCRATCH
cd('/home/ucsf/Projects/bci/clicker')
load net_mlp_hand
clear net
net = net_mlp_hand;
%net = patternnet([128 128 ]) ;
%net.performParam.regularization=0.3;
net = train(net,N,T','UseParallel','no');
cd('/home/ucsf/Projects/bci/clicker')
genFunction(net,'MLP_Hand_20240522_CL1_NoPooling_MimeHand')


%%%%%%%%%%%%%%%%%%%%%%% END SECTION %%%%%


% cd('/home/ucsf/Projects/bci')
% clear
% % %ExperimentStart('Robot3DArrow','Bravo3',4,1,0)
% ExperimentStart('HandOnline','Bravo3',4,1,0)


