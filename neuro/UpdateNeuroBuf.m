function Neuro = UpdateNeuroBuf(Neuro)
% Neuro = UpdateNeuroBuf(Neuro)
% efficiently replaces old data in circular buffer with new filtered
% signals
samps = size(Neuro.BroadbandData{end},1);

Neuro.DataBuf = circshift(Neuro.DataBuf,-samps,2);
Neuro.DataBuf(:,(end-samps+1):end) = Neuro.BroadbandData{end}';

end % UpdateNeuroBuf
