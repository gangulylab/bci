%% ERPs
% load a particular day'a data.
% Normalize the length to be constant across trials
% plot ERPs by target


clc;clear
root_path = 'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker';
foldernames = {'20201218','20210108','20210115','20210128','20210201','20210212','20210219','20210226',...
    '20210305','20210312','20210319','20210402','20210326','20210409'};
cd(root_path)

files=[];
for i=length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'Robot3DArrow');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name,'BCI_Fixed');
        files = [files;findfiles('',filepath)'];
    end
end


% load the data for each target
D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
D7=[];
for i=1:length(files)
    disp(i)
    load(files{i});
    features  = TrialData.SmoothedNeuralFeatures;
    features = cell2mat(features);
    features = features(129:256,:);
    fs = TrialData.Params.UpdateRate;
    kinax = TrialData.TaskState;
    state1 = find(kinax==1);
    state2 = find(kinax==2);
    state3 = find(kinax==3);
    state4 = find(kinax==4);
    tmp_data = features(:,state3);
    idx1= ones(length(state1),1);
    idx2= 2*ones(length(state2),1);
    idx3= 3*ones(length(state3),1);
    idx4= 4*ones(length(state4),1);
    
    % interpolate
    tb = (1/fs)*[1:size(tmp_data,2)];
    t=(1/fs)*[1:10];
    tb = tb*t(end)/tb(end);
    tmp_data1 = interp1(tb,tmp_data',t,'spline')';
    idx3 = interp1(tb,idx3,t,'spline');
    
    % now stick all the data together
    data = [features(:,[state1 state2]) tmp_data1 features(:,[state4])];

    % now get the ERPs
    if TrialData.TargetID == TrialData.SelectedTargetID
        if TrialData.TargetID == 1
            D1 = cat(3,D1,data);            
        elseif TrialData.TargetID == 2
            D2 = cat(3,D2,data);
        elseif TrialData.TargetID == 3
            D3 = cat(3,D3,data);
        elseif TrialData.TargetID == 4
            D4 = cat(3,D4,data);
        elseif TrialData.TargetID == 5
            D5 = cat(3,D5,data);
        elseif TrialData.TargetID == 6
            D6 = cat(3,D6,data);
        end        
    end    
end

% plot the ERPs now
chMap=TrialData.Params.ChMap;
figure
ha=tight_subplot(8,16);
d = 1;
set(gcf,'Color','w')
tim = cumsum([length(idx1) length(idx2) length(idx3) length(idx4)]);
for i = 1:size(D2,1)
    [x y] = find(chMap==d);
    if x == 1
        axes(ha(y));
        %subplot(8, 16, y)
    else
        s = 16*(x-1) + y;
        axes(ha(s));
        %subplot(8, 16, s)
    end
    hold on
    erps =  squeeze(D2(i,:,:));
    plot(erps, 'color', [0.4 0.4 0.4 ]);
    ylabel (num2str(i))
    axis tight
    ylim([-2 2])
    plot(mean(erps,2),'r','LineWidth',1.5)
    %set(gca,'LineWidth',1)
    %vline([time(2:4)])   
    h=vline(tim);
    %set(h,'LineWidth',1)
    set(h,'Color','b')
    h=hline(0);
    set(h,'LineWidth',1.5)    
    if i~=102
        yticklabels ''
        xticklabels ''
    else
        %xticks([tim])
        %xticklabels({'S1','S2','S3','S4'})
    end
    d = d+1;
end

%% ERPS WITH channel STATS
% find trials where it takes less than 2s for the correct action to be
% decoded


clc;clear
root_path = '/home/ucsf/Data/bravo1/20210430/';
foldernames = {'113936','114912'};
cd(root_path)

files=[];
for i=1:length(foldernames)
    folderpath = fullfile(root_path, 'DiscreteArrow',foldernames{i});
    D=dir(folderpath);
    for j=3
        filepath=fullfile(folderpath,D(j).name);
        files1 = dir(filepath);
        for k=3:length(files1)
            files=[files; fullfile(filepath, files1(k).name)];
        end        
    end
end


% load the data for each target
D1=[];
D2=[];
D3=[];
D4=[];
D5=[];
D6=[];
time_to_target=zeros(2,6);
for i=1:length(files)
    disp(i)
    load(files{i});
    features  = TrialData.SmoothedNeuralFeatures;
    features = cell2mat(features);
    features = features(769:end,:);
    %features = features(513:640,:);
    fs = TrialData.Params.UpdateRate;
    kinax = TrialData.TaskState;
    state1 = find(kinax==1);
    state2 = find(kinax==2);
    state3 = find(kinax==3);
    state4 = find(kinax==4);
    tmp_data = features(:,state3);
    idx1= ones(length(state1),1);
    idx2= 2*ones(length(state2),1);
    idx3= 3*ones(length(state3),1);
    idx4= 4*ones(length(state4),1);
    
    % interpolate
    tb = (1/fs)*[1:size(tmp_data,2)];
    t=(1/fs)*[1:10];
    tb = tb*t(end)/tb(end);
    tmp_data1 = interp1(tb,tmp_data',t,'spline')';
    idx3 = interp1(tb,idx3,t,'spline');
    
    % now stick all the data together
    trial_dur = (length(state3)-TrialData.Params.ClickCounter)*(1/fs);
    data = [features(:,[state1 state2]) tmp_data1 features(:,[state4])];
    
    % correction
    if length(state1)<8
        data  =[data(:,1) data];
    end
    
    % store the time to target data
    time_to_target(2,TrialData.TargetID) = time_to_target(2,TrialData.TargetID)+1;
    if trial_dur<=2
        time_to_target(1,TrialData.TargetID) = time_to_target(1,TrialData.TargetID)+1;
    end

    % now get the ERPs
    if TrialData.TargetID == TrialData.SelectedTargetID && trial_dur<=2
        if TrialData.TargetID == 1
            D1 = cat(3,D1,data);            
        elseif TrialData.TargetID == 2
            D2 = cat(3,D2,data);
        elseif TrialData.TargetID == 3
            D3 = cat(3,D3,data);
        elseif TrialData.TargetID == 4
            D4 = cat(3,D4,data);
        elseif TrialData.TargetID == 5
            D5 = cat(3,D5,data);
        elseif TrialData.TargetID == 6
            D6 = cat(3,D6,data);
        end        
    end    
end


time_to_target(1,:)./time_to_target(2,:)

% plot the ERPs now
chMap=TrialData.Params.ChMap;
figure
ha=tight_subplot(8,16);
d = 1;
set(gcf,'Color','w')
tim = cumsum([length(idx1) length(idx2) length(idx3) length(idx4)]);
for i = 1:size(D1,1)
    [x y] = find(chMap==d);
    if x == 1
        axes(ha(y));
        %subplot(8, 16, y)
    else
        s = 16*(x-1) + y;
        axes(ha(s));
        %subplot(8, 16, s)
    end
    hold on
    erps =  squeeze(D1(i,:,:));
    plot(erps, 'color', [0.4 0.4 0.4 ]);
    ylabel (num2str(i))
    axis tight
    ylim([-2 2])
    plot(mean(erps,2),'r','LineWidth',1.5)
    %set(gca,'LineWidth',1)
    %vline([time(2:4)])   
    h=vline(tim);
    %set(h,'LineWidth',1)
    set(h,'Color','b')
    h=hline(0);
    set(h,'LineWidth',1.5)    
    if i~=102
        yticklabels ''
        xticklabels ''
    else
        %xticks([tim])
        %xticklabels({'S1','S2','S3','S4'})
    end
    d = d+1;
end
% 
% 
% bootstrapped confidnce intervals for a particular channel
% chdata = ch100;
% 
% zscore the data to the first 8 time-bins
% tmp_data=chdata(1:8,:);
% m = mean(tmp_data(:));
% s = std(tmp_data(:));
% chdata = (chdata -m)./s;
%  
% 
% m = mean(chdata,2);
% mb = sort(bootstrp(1000,@mean,chdata'));
% figure;
% hold on
% plot(m,'b')
% plot(mb(25,:),'--b')
% plot(mb(975,:),'--b')
% hline(0)
% 
% % shuffle the data and see the results
% tmp_mean=[];
% for i=1:1000
%     %tmp = circshift(chdata,randperm(size(chdata,1),1));
%     tmp = chdata;
%     tmp(randperm(numel(chdata))) = tmp;
%     tmp_data=tmp(1:8,:);
%     m = mean(tmp_data(:));
%     s = std(tmp_data(:));
%     tmp = (tmp -m)./s;
%     tmp_mean(i,:) = mean(tmp,2);
% end
% 
% tmp_mean = sort(tmp_mean);
% plot(tmp_mean(25,:),'--r')
% plot(tmp_mean(975,:),'--r')



% plot the ERPs with bootstrapped C.I. shading
chMap=TrialData.Params.ChMap;
figure
ha=tight_subplot(8,16);
d = 1;
set(gcf,'Color','w')
tim = cumsum([length(idx1) length(idx2) length(idx3) length(idx4)]);
for i = 1:size(D5,1)
    [x y] = find(chMap==i);
    if x == 1
        axes(ha(y));
        %subplot(8, 16, y)
    else
        s = 16*(x-1) + y;
        axes(ha(s));
        %subplot(8, 16, s)
    end
    hold on
    erps =  squeeze(D5(i,:,:));
    
    chdata = erps;
    % zscore the data to the first 8 time-bins
    tmp_data=chdata(1:8,:);
    m = mean(tmp_data(:));
    s = std(tmp_data(:));
    chdata = (chdata -m)./s;
    
    % get the confidence intervals
    m = mean(chdata,2);
    mb = sort(bootstrp(1000,@mean,chdata'));
    tt=1:32;
    [fillhandle,msg]=jbfill(tt,mb(25,:),mb(975,:)...
        ,[0.3 0.3 0.7],[0.3 0.3 0.7],1,.2);
    hold on
    plot(m,'b')
    %plot(mb(25,:),'--b')
    %plot(mb(975,:),'--b')
    %hline(0)
    
    % shuffle the data for null confidence intervals
    tmp_mean=[];
    for j=1:1000
        %tmp = circshift(chdata,randperm(size(chdata,1),1));
        tmp = chdata;
        tmp(randperm(numel(chdata))) = tmp;
        tmp_data=tmp(1:8,:);
        m = mean(tmp_data(:));
        s = std(tmp_data(:));
        tmp = (tmp -m)./s;
        tmp_mean(j,:) = mean(tmp,2);
    end
    
    tmp_mean = sort(tmp_mean);
    %plot(tmp_mean(25,:),'--r')
    %plot(tmp_mean(975,:),'--r')
    [fillhandle,msg]=jbfill(tt,tmp_mean(25,:),tmp_mean(975,:)...
        ,[0.7 0.3 0.3],[0.7 0.3 0.3],1,.2);
    
    
    % statistical test
    % if the mean is outside confidence intervals in state 3
    m = mean(chdata,2);
    idx=13:27;
    mstat = m((idx));
    pval=[];
    for j=1:length(idx)
        pval(j) = (sum(abs(mstat(j)) >= abs(tmp_mean(:,idx(j)))))./(size(tmp_mean,1));
    end
    
    res=sum((1-pval)<=0.05);
    if res>=7
        suc=1;
    else
        suc=0;
    end
    
    % beautify
    ylabel (num2str(i))
    axis tight
    ylim([-2 2])    
    %set(gca,'LineWidth',1)
    %vline([time(2:4)])   
    h=vline(tim);
    %set(h,'LineWidth',1)
    set(h,'Color','k')
    h=hline(0);
    set(h,'LineWidth',1.5)    
    if i~=102
        yticklabels ''
        xticklabels ''
    else
        %xticks([tim])
        %xticklabels({'S1','S2','S3','S4'})
    end
    
    if suc==1
        box on
        set(gca,'LineWidth',2)
        set(gca,'XColor','g')
        set(gca,'YColor','g')
    end
    d = d+1;
end

%% putting it all together ERPs with stats across days

clc;clear
root_path = 'E:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker';
foldernames = {'20201218','20210115','20210128','20210201','20210212','20210219','20210226',...
    '20210305','20210312','20210319','20210402','20210326','20210409','20210416'};
cd(root_path)

T1=[];
T2=[];
T3=[];
T4=[];
T5=[];
T6=[];
time_to_target_overall=[];
for i=length(foldernames)
    disp(['processing folder ' num2str(i) ' of ' num2str(length(foldernames))])
    
    % get all the files
    files=[];
    folderpath = fullfile(root_path, foldernames{i},'Robot3DArrow');
    D=dir(folderpath);
    for j=3:length(D)
        filepath=fullfile(folderpath,D(j).name,'BCI_Fixed');
        files = [files;findfiles('',filepath)'];
    end
    
    % now load all the data
    [D1,D2,D3,D4,D5,D6,time_to_target] = load_erp_data(files);
    time_to_target = time_to_target(1,:)./time_to_target(2,:);
    time_to_target_overall(i,:) = time_to_target;
    % chmap file
    load(files{1})
    chMap=TrialData.Params.ChMap;
    
    % get the mask of sig. channels for each movement type  and save ERPs
    clear res1 res2 res3 res res5 res6
    res1 = erp_stats(D1,chMap,foldernames{i},'T1',0);
    res2 = erp_stats(D2,chMap,foldernames{i},'T2',0);
    res3 = erp_stats(D3,chMap,foldernames{i},'T3',0);
    res4 = erp_stats(D4,chMap,foldernames{i},'T4',0);
    res5 = erp_stats(D5,chMap,foldernames{i},'T5',0);
    res6 = erp_stats(D6,chMap,foldernames{i},'T6',0);
    
    % collate across days
    T1(i,:,:) = res1;
    T2(i,:,:) = res2;
    T3(i,:,:) = res3;
    T4(i,:,:) = res4;
    T5(i,:,:) = res5;
    T6(i,:,:) = res6;
    
end

% save the results
save ERP_data

% plotting the sig channels over days using all the days
figure
ha=tight_subplot(7,2);
for i=1:14
    axes(ha(i))
    if i==4 || i==6 || i==8 || i==14
        box off
        axis off
    else
        imagesc(squeeze(T1(i,:,:)));
        colormap parula
        caxis([0 1])
        xticks ''
        yticks ''
        box on
    end
end
set(gcf,'Color','w')

% plotting the sig channels over days using only useful days
idx = [1:3 5 7 9:13];
data = T6(idx,:,:);
figure
ha=tight_subplot(5,2);
for i=1:10
    axes(ha(i))
    imagesc(squeeze(data(i,:,:)));
    colormap parula
    caxis([0 1])
    xticks ''
    yticks ''
    box on
    
end
set(gcf,'Color','w')

figure
ha=tight_subplot(3,2);
axes(ha(1))
imagesc(squeeze(sum(T1(idx,:,:),1)));colormap bone; caxis([0 10])
xticks ''
yticks ''
axes(ha(2))
imagesc(squeeze(sum(T2(idx,:,:),1)));colormap bone; caxis([0 10])
xticks ''
yticks ''
axes(ha(3))
imagesc(squeeze(sum(T3(idx,:,:),1)));colormap bone; caxis([0 10])
xticks ''
yticks ''
axes(ha(4))
imagesc(squeeze(sum(T4(idx,:,:),1)));colormap bone; caxis([0 10])
xticks ''
yticks ''
axes(ha(5))
imagesc(squeeze(sum(T5(idx,:,:),1)));colormap bone; caxis([0 10])
xticks ''
yticks ''
axes(ha(6))
imagesc(squeeze(sum(T6(idx,:,:),1)));colormap bone; caxis([0 10])
xticks ''
yticks ''
set(gcf,'Color','w')



figure;imagesc(squeeze(sum(T1,1)));colormap bone
caxis([0 13])
figure;imagesc(squeeze(sum(T2,1)));colormap bone
caxis([0 13])
figure;imagesc(squeeze(sum(T3,1)));colormap bone
caxis([0 13])
figure;imagesc(squeeze(sum(T4,1)));colormap bone
caxis([0 13])
figure;imagesc(squeeze(sum(T5,1)));colormap bone
caxis([0 13])
figure;imagesc(squeeze(sum(T6,1)));colormap bone
caxis([0 13])





