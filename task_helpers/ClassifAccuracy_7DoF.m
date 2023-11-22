%% GETTING ACCURACY AT BIN LEVEL AND TRIAL LEVEL FROM 3D ARROW TASK


clc;clear
root_path = '/home/ucsf/Data/Bravo3/20231115/Robot3DArrow';
folders = {'144845', '145418', '145745', '150314', '150646'};
addpath('/home/ucsf/Projects/bci/task_helpers')

files=[];
for ii=1:length(folders)
    folderpath = fullfile(root_path, folders{ii},'BCI_Fixed');
    %cd(folderpath)
    files = [files;findfiles('',folderpath)'];
end


acc=zeros(7); % trial level decoding accuracy
acc1=zeros(7); % bin level decoding accuracy 
for i=1:length(files)
    file_loaded=1;
    try
        load(files{i});
    catch
        file_loaded=0;
    end

    if file_loaded
        out = TrialData.ClickerState;
        out1 = TrialData.FilteredClickerState;
        tid = TrialData.TargetID;
        decodes=[];
        for ii=1:7
            decodes(ii) = sum(out==ii);
        end
        [aa bb]=max(decodes);
        acc(tid,bb) = acc(tid,bb)+1; % trial level
        for j=1:length(out)
            if out(j)>0
                acc1(tid,out(j)) = acc1(tid,out(j))+1; % bin level
            end
        end
    end
end

for i=1:7
    acc(i,:) = acc(i,:)/sum(acc(i,:));
    acc1(i,:) = acc1(i,:)/sum(acc1(i,:));
end

clc
figure;
imagesc(acc1)
colormap bone
caxis([0 1])
title(['Bin level decoding acc. ' num2str(100*mean(diag(acc1)))])
disp('Bin level decoding accuracy')
disp(acc1)






