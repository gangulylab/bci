%% PERFORMANCE MEASURE USING MAX VOTE STRATEGY FOR DECODES 

clc;clear

root_path='F:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker';

foldernames = {'20220603'};
cd(root_path)

files=[];
for i=1:length(foldernames)
    folderpath = fullfile(root_path, foldernames{i},'Robot3DArrow');
    D=dir(folderpath);
    for j=1:length(D)
        filepath=fullfile(folderpath,D((j)).name,'BCI_Fixed');
        if exist(filepath)
            disp(filepath)
            files = [files;findfiles('',filepath)'];
        end
    end
end


% look at the decodes per direction to get a max vote
T=zeros(7);
tim_to_target=[];
num_suc=[];
num_fail=[];
for i=1:length(files)
    disp(i)
    indicator=1;
    try
        load(files{i});
    catch ME
        warning('Not able to load file, skipping to next')
        indicator = 0;
    end
    if indicator
        kinax = TrialData.TaskState;
        clicker_state = TrialData.FilteredClickerState;
        idx = TrialData.TargetID;
        t(1) = sum(clicker_state ==1);
        t(2) = sum(clicker_state ==2);
        t(3) = sum(clicker_state ==3);
        t(4) = sum(clicker_state ==4);
        t(5) = sum(clicker_state ==5);
        t(6) = sum(clicker_state ==6);
        t(7) = sum(clicker_state ==7);
        [aa bb]=max(t);
        T(idx,bb) = T(idx,bb)+1;
        if TrialData.TargetID == TrialData.SelectedTargetID
            tim_to_target = [tim_to_target length(clicker_state)-TrialData.Params.ClickCounter];
            num_suc = [num_suc 1];
        else%if TrialData.SelectedTargetID ==0%~= TrialData.SelectedTargetID %&& TrialData.SelectedTargetID~=0
            tim_to_target = [tim_to_target length(clicker_state)];
            num_fail = [num_fail 1];
        end
    end
end

% get acc.
for i=1:size(T)
    T(i,:) = T(i,:)./sum(T(i,:));
end
figure;imagesc(T)
colormap bone
caxis([0 1])
xticks([1:7])
yticks([1:7])
xticklabels({'Rt thumb','Both Feet','Lt. thumb','Grasp Ext.', 'Rt wrist','Tong','Both middle'})
yticklabels({'Rt thumb','Both Feet','Lt. thumb','Grasp Ext.', 'Rt wrist','Tong','Both middle'})
set(gcf,'Color','w')
set(gca,'FontSize',12)
colorbar

% bit rate calculations
tim_to_target = tim_to_target.*(1/TrialData.Params.UpdateRate);
B = log2(7-1) * (sum(num_suc)-sum(num_fail)) / sum(tim_to_target)
title(['Accuracy of ' num2str(100*mean(diag(T))) '%' ' and bitrate of ' num2str(B)])
