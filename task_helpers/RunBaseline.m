function Neuro = RunBaseline(Params,Neuro)
% Neuro = RunBaseline(Params,Neuro)
% Baseline period - grab neural features to get baseline for z-scoring
% 
% Neuro - contains NeuralStats structure w/ mean and var of each chan

global Cursor
if isempty(Cursor),
    Cursor.LastUpdateTime = GetSecs;
end

msg = '';
tstart  = GetSecs;
tlast = tstart;
done = 0;
while ~done,
    % Update Time & Position
    tim = GetSecs;
    
    % for pausing and quitting expt
    if CheckPause, ExperimentPause(Params); end
    
    % Grab data every Xsecs
    if (tim-tlast) > 1/Params.ScreenRefreshRate,
        % time
        tlast = tim;
        
        % grab and process neural data
        if Params.BLACKROCK && ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate),
            Cursor.LastUpdateTime = tim;
            Cursor.LastPredictTime = tim;
            Neuro = NeuroPipeline(Neuro, [], Params);
        end
        
        % update command line with progress
        fprintf(repmat('\b',1,length(msg)-1))
        msg = sprintf('Computing Baseline: %i%s\n', round(100*(tim-tstart)/Params.BaselineTime),'%%');
        fprintf(1,msg)
        
        % update screen with progress
        tex = sprintf('Computing Baseline: %.1f%% \n', 100*(tim-tstart)/Params.BaselineTime);
      %  DrawFormattedText(Params.WPTR, tex,'center','center',255);
       % Screen('Flip', Params.WPTR);
        
    end
    
    % end if takes too long
    if (tim - tstart) > Params.BaselineTime,
        done = 1;
    end
end

fprintf('Done\n\n')

end % RunBaseline