function [trial_data,coeff,score,latent] = load_data_for_MLP_TrialLevel_B3_bci(files,trial_type)
%function [condn_data] = load_data_for_MLP(files)


trial_data=[];
tmp_data=[];
for ii=1:length(files)
    disp(ii)
    file_loaded=1;
    try
        load(files{ii});
    catch
        file_loaded=0;
    end
    if file_loaded
        features  = TrialData.SmoothedNeuralFeatures;
        temp = cell2mat(features);
        kinax = find(TrialData.TaskState==3);
        temp = cell2mat(features(kinax));

        kinax1 = find(TrialData.TaskState==1);
        temp_state1 = cell2mat(features(kinax1));

        % relative to state 1:
        %temp = temp - mean(temp_state1,2);


        % get delta, beta and hG removing bad channels
        temp = temp([257:512 1025:1280 1537:1792],:); %delta, beta, hG
        %temp = temp([1537:1792],:);% only hG
        bad_ch = [108 113 118];
        good_ch = ones(size(temp,1),1);
        for iii=1:length(bad_ch)
            %bad_ch_tmp = bad_ch(iii)*[1 2 3];
            bad_ch_tmp = bad_ch(iii)+(256*[0 1 2]);
            good_ch(bad_ch_tmp)=0;
        end
        temp = temp(logical(good_ch),:);

        %
        %         % get the pooled data
        %         new_temp=[];
        %         TrialData.Params.ChMap = ecog_grid;
        %         [xx yy] = size(TrialData.Params.ChMap);
        %         for k=1:size(temp,2)
        %             tmp1 = temp(129:256,k);tmp1 = tmp1(TrialData.Params.ChMap);
        %             tmp2 = temp(513:640,k);tmp2 = tmp2(TrialData.Params.ChMap);
        %             tmp3 = temp(769:896,k);tmp3 = tmp3(TrialData.Params.ChMap);
        %             pooled_data=[];
        %             for i=1:2:xx
        %                 for j=1:2:yy
        %                     delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
        %                     beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
        %                     hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
        %                     pooled_data = [pooled_data; delta; beta ;hg];
        %                 end
        %             end
        %             new_temp= [new_temp pooled_data];
        %         end
        %         temp=new_temp;
        %

        % 2-norm
        for i=1:size(temp,2)
            temp(:,i) = temp(:,i)./norm(temp(:,i));
        end

        trial_data(ii).neural = temp;
        trial_data(ii).targetID = TrialData.TargetID;
        trial_data(ii).trial_type = trial_type;
        % store for PCA
        tmp_data = [tmp_data temp];
    end
end

% return PCs
%[coeff,score,latent] = pca(tmp_data','Centered','off');


end




