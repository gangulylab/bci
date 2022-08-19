function Neuro = UpdateNeuroBuf(Neuro)
% Neuro = UpdateNeuroBuf(Neuro)
% efficiently replaces old data in circular buffer with new filtered
% signals


samps = Neuro.NumSamps;
Neuro.FilterDataBuf = circshift(Neuro.FilterDataBuf,-samps);
Neuro.FilterDataBuf((end-samps+1):end,:,:) = Neuro.FilteredData(:,:,1:Neuro.NumBuffer);

end % UpdateNeuroBuf
