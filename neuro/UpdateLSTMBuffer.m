function Neuro = UpdateLSTMBuffer(Neuro)
% function Neuro = UpdateLSTMBuffer(Neuro)

% update the buffer 
samps = size(Neuro.BroadbandData,1);  
if samps>Neuro.LSTMBufferSize
    samps = Neuro.LSTMBufferSize;
end
Neuro.LSTMBuffer = circshift(Neuro.LSTMBuffer,-samps,2);
Neuro.LSTMBuffer(:,(end-samps+1):end) = Neuro.BroadbandData(end-samps+1:end,:)';

% do the hg filtering
filtered_data=zeros(size(Neuro.LSTMBuffer',1),size(Neuro.LSTMBuffer',2),8);
for i=9:16%hg features
    filtered_data(:,:,i) =  ((filter(...
        Neuro.FilterBank(i).b, ...
        Neuro.FilterBank(i).a, ...
        Neuro.LSTMBuffer')));
end
tmp_hg = squeeze(mean(filtered_data.^2,3));

% low pass filtering
tmp_lp = filter(Neuro.lpFilt,Neuro.LSTMBuffer');

% down sampling
tmp_hg = resample(tmp_hg,80,800)*5e2;
tmp_lp = resample(tmp_lp,80,800);

% removing errors in the data
I = abs(tmp_hg>15);
I = sum(I);
[aa bb]=find(I>0);
tmp_hg(:,bb) = 1e-5*randn(size(tmp_hg(:,bb)));

I = abs(tmp_lp>15);
I = sum(I);
[aa bb]=find(I>0);
tmp_lp(:,bb) = 1e-5*randn(size(tmp_lp(:,bb)));

% normalizing between 0 and 1
tmp_hg = (tmp_hg - min(tmp_hg(:)))/(max(tmp_hg(:))-min(tmp_hg(:)));
tmp_lp = (tmp_lp - min(tmp_lp(:)))/(max(tmp_lp(:))-min(tmp_lp(:)));

% concatenating into LSTM features
Neuro.LSTMFeatures = [tmp_hg tmp_lp]';

end
