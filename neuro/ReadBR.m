function Neuro = ReadBR(Neuro)
% function Neuro = ReadBR(Neuro)
% reads in all data in the data buffer. intended to be run at each
% timestep. if data exceeds limit in Params, data is truncated.
%
% Neuro -
%   .TimeStamp - returns timestamp of data collection
%   .BroadbandData - returns matrix of neural data [ samples x channels ]
%   .NumSamps - number of samples of neural data per channel

% The data looks like { [] [] [ samples x channels ] }
% read data from blackrock
[timestamp, X] = cbmex('trialdata',1); % read buffer
% try
neural_data = double(horzcat(X{1:Neuro.NumChannels,3}));
% bad_ch = sort([108 113 118]);
% good_ch = ones(size(neural_data,2),1);
% good_ch(bad_ch) = 0;
% good_ch = logical(good_ch);
% neural_data = neural_data(:,good_ch);

% catch
%     X{125,3} = downsample(X{125,3},30);
%     neural_data = double(horzcat(X{1:Neuro.NumChannels,3}));
% end

% limit to buffer size
num_samps = size(neural_data,1);
if num_samps > Neuro.BufferSamps,
    neural_data = neural_data(1:Neuro.BufferSamps,:);
    num_samps = Neuro.BufferSamps;
end

% put in Neuro structure
Neuro.TimeStamp = timestamp;
Neuro.BroadbandData = neural_data;
Neuro.NumSamps = num_samps;

end % ReadBR

