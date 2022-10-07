function [XTrain,XTest,YTrain,YTest] = get_lstm_features_robotBatch(files,Params,lpFilt,chmap)


% get the raw data
D1={};
D2={};
D3={};
D4={};
D5={};
D6={};
D7={};
for i=1:length(files)    

    try
        load(files{i})
        file_loaded = true;
        warning('off')
    catch
        file_loaded=false;
        disp(['Could not load ' files{i}]);
        files_not_loaded=[files_not_loaded;files(i)];
    end
    if file_loaded

        idx00 = find(TrialData.TaskState==1) ;
        idx0 = find(TrialData.TaskState==2) ;
        idx = find(TrialData.TaskState==3) ;
        idx=[idx00 idx0 idx];
        raw_data = cell2mat(TrialData.BroadbandData(idx)');
        idx1 = find(TrialData.TaskState==4) ;
        raw_data4 = cell2mat(TrialData.BroadbandData(idx1)');
        id = TrialData.TargetID;
        s = size(raw_data,1);
        %         if s>7800
        %             s=7800;
        %         end
        data_seg={};
        if s<800 % for really quick decisions just pad data from state 4
            len = 800-s;
            tmp = raw_data4(1:len,:);
            raw_data = [raw_data;tmp];
            data_seg = raw_data;
        elseif s>800 && s<1000 % if not so quick, prune to data to 600ms
            raw_data = raw_data(1:800,:);
            data_seg = raw_data;
        elseif s>1000% for all other data length, have to parse the data in overlapping chuncks of 600ms, 50% overlap
            %bins =1:400:s; % originally only for state 3
            bins = 1:500:s;
            jitter = round(100*rand(size(bins)));
            bins=bins+jitter;
            raw_data = [raw_data;raw_data4];
            for k=1:length(bins)-1
                tmp = raw_data(bins(k)+[0:999],:);
                data_seg = cat(2,data_seg,tmp);
            end
        end

        feat_stats = TrialData.FeatureStats;
        feat_stats.Mean = feat_stats.Mean(769:end);
        feat_stats.Var = feat_stats.Var(769:end);
        clear feat_stats1
        feat_stats1(1:length(data_seg)) = feat_stats;

        if id==1
            D1 = cat(2,D1,data_seg);
            %D1f = cat(2,D1f,feat_stats1);
        elseif id==2
            D2 = cat(2,D2,data_seg);
            %D2f = cat(2,D2f,feat_stats1);
        elseif id==3
            D3 = cat(2,D3,data_seg);
            %D3f = cat(2,D3f,feat_stats1);
        elseif id==4
            D4 = cat(2,D4,data_seg);
            %D4f = cat(2,D4f,feat_stats1);
        elseif id==5
            D5 = cat(2,D5,data_seg);
            %D5f = cat(2,D5f,feat_stats1);
        elseif id==6
            D6 = cat(2,D6,data_seg);
            %D6f = cat(2,D6f,feat_stats1);
        elseif id==7
            D7 = cat(2,D7,data_seg);
            %D7f = cat(2,D7f,feat_stats1);
        end      
    end
end


% run it through the lstm processing pipeline
condn_data_new = [];
Y=[];
jj=1;
len=1000;
condn_data1 = zeros(len,128,length(D1));
k=1;
for i=1:length(D1)
    %disp(k)
    tmp = D1{i};
    %tmp1(:,1,:)=tmp;
    tmp1=tmp;
    condn_data1(:,:,k) = tmp1;
    k=k+1;
    Y = [Y ;1];
end

disp('Processing action 1')
for ii=1:size(condn_data1,3)
     %disp(ii)

    tmp = squeeze(condn_data1(:,:,ii));

    tmp = extract_lstm_features(tmp,Params,lpFilt);

    % store
    condn_data_new(:,:,jj) = tmp;
    jj=jj+1;
end

condn_data2 = zeros(len,128,length(D2));
k=1;
for i=1:length(D2)
    %disp(k)
    tmp = D2{i};
    %tmp1(:,1,:)=tmp;
    tmp1=tmp;
    condn_data2(:,:,k) = tmp1;
    k=k+1;
    Y = [Y ;2];
end
disp('Processing action 2')
for ii=1:size(condn_data2,3)
    %disp(ii)

    tmp = squeeze(condn_data2(:,:,ii));

    tmp = extract_lstm_features(tmp,Params,lpFilt);

    % store
    condn_data_new(:,:,jj) = tmp;
    jj=jj+1;
end

condn_data3 = zeros(len,128,length(D3));
k=1;
for i=1:length(D3)
    %disp(k)
    tmp = D3{i};
    %tmp1(:,1,:)=tmp;
    tmp1=tmp;
    condn_data3(:,:,k) = tmp1;
    k=k+1;
    Y = [Y ;3];
end
disp('Processing action 3')
for ii=1:size(condn_data3,3)
    %disp(ii)

    tmp = squeeze(condn_data3(:,:,ii));

    tmp = extract_lstm_features(tmp,Params,lpFilt);

    % store
    condn_data_new(:,:,jj) = tmp;
    jj=jj+1;
end

condn_data4 = zeros(len,128,length(D4));
k=1;
for i=1:length(D4)
    %disp(k)
    tmp = D4{i};
    %tmp1(:,1,:)=tmp;
    tmp1=tmp;
    condn_data4(:,:,k) = tmp1;
    k=k+1;
    Y = [Y ;4];
end
disp('Processing action 4')
for ii=1:size(condn_data4,3)
    %disp(ii)

    tmp = squeeze(condn_data4(:,:,ii));

    tmp = extract_lstm_features(tmp,Params,lpFilt);


    % store
    condn_data_new(:,:,jj) = tmp;
    jj=jj+1;
end

condn_data5 = zeros(len,128,length(D5));
k=1;
for i=1:length(D5)
    %disp(k)
    tmp = D5{i};
    %tmp1(:,1,:)=tmp;
    tmp1=tmp;
    condn_data5(:,:,k) = tmp1;
    k=k+1;
    Y = [Y ;5];
end
disp('Processing action 5')
for ii=1:size(condn_data5,3)
    %disp(ii)

    tmp = squeeze(condn_data5(:,:,ii));


    tmp = extract_lstm_features(tmp,Params,lpFilt);

    % store
    condn_data_new(:,:,jj) = tmp;
    jj=jj+1;
end

condn_data6 = zeros(len,128,length(D6));
k=1;
for i=1:length(D6)
    %disp(k)
    tmp = D6{i};
    %tmp1(:,1,:)=tmp;
    tmp1=tmp;
    condn_data6(:,:,k) = tmp1;
    k=k+1;
    Y = [Y ;6];
end
disp('Processing action 6')
for ii=1:size(condn_data6,3)
    %disp(ii)

    tmp = squeeze(condn_data6(:,:,ii));

    tmp = extract_lstm_features(tmp,Params,lpFilt);

    % store
    condn_data_new(:,:,jj) = tmp;
    jj=jj+1;
end
condn_data7 = zeros(len,128,length(D7));
k=1;
for i=1:length(D7)
    %disp(k)
    tmp = D7{i};
    %tmp1(:,1,:)=tmp;
    tmp1=tmp;
    condn_data7(:,:,k) = tmp1;
    k=k+1;
    Y = [Y ;7];
end
disp('Processing action 7')
for ii=1:size(condn_data7,3)
    %disp(ii)

    tmp = squeeze(condn_data7(:,:,ii));

    tmp = extract_lstm_features(tmp,Params,lpFilt);


    % store
    condn_data_new(:,:,jj) = tmp;
    jj=jj+1;
end

% get rid of artifacts, any channel with activity >15SD, set it to near zero
for i=1:size(condn_data_new,3)
    xx=squeeze(condn_data_new(:,1:128,i));
    I = abs(xx)>15;
    I = sum(I);
    [aa bb]=find(I>0);
    xx(:,bb) = 1e-5*randn(size(xx(:,bb)));
    condn_data_new(:,1:128,i)=xx;

    xx=squeeze(condn_data_new(:,129:256,i));
    I = abs(xx)>15;
    I = sum(I);
    [aa bb]=find(I>0);
    xx(:,bb) = 1e-5*randn(size(xx(:,bb)));
    condn_data_new(:,129:256,i)=xx;

    %     xx=squeeze(condn_data_new(:,257:384,i));
    %     I = abs(xx)>15;
    %     I = sum(I);
    %     [aa bb]=find(I>0);
    %     xx(:,bb) = 1e-5*randn(size(xx(:,bb)));
    %     condn_data_new(:,257:384,i)=xx;
end

% normalize the data to be between 0 and 1
for i=1:size(condn_data_new,3)
    tmp=squeeze(condn_data_new(:,:,i));
    tmp1=tmp(:,1:128);
    tmp1 = (tmp1 - min(tmp1(:)))/(max(tmp1(:))-min(tmp1(:)));

    tmp2=tmp(:,129:256);
    tmp2 = (tmp2 - min(tmp2(:)))/(max(tmp2(:))-min(tmp2(:)));

    %     tmp3=tmp(:,257:384);
    %     tmp3 = (tmp3 - min(tmp3(:)))/(max(tmp3(:))-min(tmp3(:)));

    %tmp = [tmp1 tmp2 tmp3];
    tmp = [tmp1 tmp2 ];
    condn_data_new(:,:,i)=tmp;
end




% split into training and validation datasets
idx = randperm(size(condn_data_new,3),round(0.8*size(condn_data_new,3)));
I = zeros(size(condn_data_new,3),1);
I(idx)=1;

XTrain={};
XTest={};
YTrain=[];
YTest=[];
for i=1:size(condn_data_new,3)
    tmp = squeeze(condn_data_new(:,:,i));
    if I(i)==1
        XTrain = cat(1,XTrain,tmp');
        YTrain = [YTrain Y(i)];
    else
        XTest = cat(1,XTest,tmp');
        YTest = [YTest Y(i)];
    end
end

% shuffle
idx  = randperm(length(YTrain));
XTrain = XTrain(idx);
YTrain = YTrain(idx);

YTrain = categorical(YTrain');
YTest = categorical(YTest');

% data augmentation: introduce random noise plus some mean shift to each
% channel for about 50k samples
aug_idx = randperm(length(XTrain));
for i=1:length(aug_idx)
    %disp(i)
    tmp = XTrain{aug_idx(i)}';
    t_id=categorical(YTrain(aug_idx(i)));

    % hG
    tmp1 = tmp(:,1:128);
    % add variable noise
    %var_noise=randsample(400:1200,size(tmp1,2))/1e3;
    var_noise=0.7;
    add_noise=randn(size(tmp1)).*std(tmp1).*var_noise;
    tmp1n = tmp1 + add_noise;
    % add variable mean offset between 5 and 25%
    m=mean(tmp1);
    add_mean =  m*.25;
    %add_mean=randsample(0:500,size(tmp1,2))/1e3;
    flip_sign = rand(size(add_mean));
    flip_sign(flip_sign>0.5)=1;
    flip_sign(flip_sign<=0.5)=-1;
    add_mean=add_mean.*flip_sign+m;
    tmp1m = tmp1n + add_mean;
    %tmp1m = (tmp1m-min(tmp1m(:)))/(max(tmp1m(:))-min(tmp1m(:)));
    %  figure;plot(tmp1(:,3));hold on;plot(tmp1m(:,3))

    % lmp
    tmp2 = tmp(:,129:256);
    % add variable noise
    var_noise=0.7;
    %var_noise=randsample(400:1200,size(tmp2,2))/1e3;
    add_noise=randn(size(tmp2)).*std(tmp2).*var_noise;
    tmp2n = tmp2 + add_noise;
    % add variable mean offset between 5 and 25%
    m=mean(tmp2);
    add_mean =  m*.35;
    %add_mean=randsample(0:500,size(tmp2,2))/1e3;
    flip_sign = rand(size(add_mean));
    flip_sign(flip_sign>0.5)=1;
    flip_sign(flip_sign<=0.5)=-1;
    add_mean=add_mean.*flip_sign+m;
    tmp2m = tmp2n + add_mean;
    % tmp2m = (tmp2m-min(tmp2m(:)))/(max(tmp2m(:))-min(tmp2m(:)));
    
    %tmp=[tmp1m tmp2m tmp3m]';
    tmp=[tmp1m tmp2m]';

    XTrain=cat(1,XTrain,tmp);
    YTrain = cat(1,YTrain,t_id);
end







