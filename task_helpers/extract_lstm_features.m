function out = extract_lstm_features(tmp,Params,lpFilt,chmap)


% do an optional laplacian ref
%tmp = laplacian_ref(tmp',chmap)';

%get hG through filter bank approach
filtered_data=zeros(size(tmp,1),size(tmp,2),8);
for i=1:8 % only hg
    filtered_data(:,:,i) =  ((filter(...
        Params.FilterBank(i).b, ...
        Params.FilterBank(i).a, ...
        tmp)));
end
tmp_hg = squeeze(mean(filtered_data.^2,3)); % hgEnvelope


% LFO low pass filtering
tmp_lp = filter(lpFilt,tmp);

% % get lg thru filter bank approach
% filtered_data=zeros(size(tmp,1),size(tmp,2),3);
% for i=9:11 % only lg
%     filtered_data(:,:,i) =  ((filter(...
%         Params.FilterBank(i).b, ...
%         Params.FilterBank(i).a, ...
%         tmp)));
% end
% tmp_lg = squeeze(mean(filtered_data.^2,3));

% downsample the data
%     tmp_lp = resample(tmp_lp,200,800);
%     tmp_hg = resample(tmp_hg,200,800)*5e2;

% decimate the data, USE AN OPTIONAL SMOOTHING INFO HERE
%     tmp_hg1=[];
%     tmp_lp1=[];
%     for i=1:size(tmp_hg,2)
%         tmp_hg1(:,i) = decimate(tmp_hg(:,i),20)*5e2;
%         tmp_lp1(:,i) = decimate(tmp_lp(:,i),20);
%     end

% resample
%tmp_lg1=resample(tmp_lg,80,size(tmp_lg,1))*5e2;
tmp_hg1=resample(tmp_hg,size(tmp_hg,1)/10,size(tmp_hg,1))*5e2;
tmp_lp1=resample(tmp_lp,size(tmp_hg,1)/10,size(tmp_lp,1));

% laplacian referencing 
%tmp_hg1 = laplacian_ref(tmp_hg1',chmap)';
%tmp_lp1 = laplacian_ref(tmp_lp1',chmap)';

% make new data structure
tmp = [tmp_hg1 tmp_lp1 ];
out=tmp;



