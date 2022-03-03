function Neuro = UpdateLSTMBuffer(Neuro)
% function Neuro = UpdateLSTMBuffer(Neuro)

% update the buffer 
samps = size(Neuro.BroadbandData{end},1);
Neuro.LSTMBuffer = circshift(Neuro.LSTMBuffer,-samps,2);
Neuro.LSTMBuffer(:,(end-samps+1):end) = Neuro.BroadbandData{end}';

% do the hg filtering
filtered_data=zeros(size(Neuro.LSTMBuffer',1),size(Neuro.LSTMBuffer',2),8);
for i=1:length(Neuro.FilterBank)
    filtered_data(:,:,i) =  ((filter(...
        Neuro.FilterBank(i).b, ...
        Neuro.FilterBank(i).a, ...
        Neuro.LSTMBuffer')));
end
tmp_hg = squeeze(mean(filtered_data.^2,3));

% low pass filtering
tmp_lp = filter(Neuro.lpFilt,Neuro.DataBuf');

% down sampling
tmp_lp = resample(tmp_lp,200,800);
tmp_hg = resample(tmp_hg,200,800)*5e2;

% concatenating into LSTM features


end