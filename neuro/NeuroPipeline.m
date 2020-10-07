function varargout = NeuroPipeline(Neuro,Data,Params),
% function Neuro = NeuroPipeline(Neuro)
% function [Neuro,Data] = NeuroPipeline(Neuro,Data)
% Neuro processing pipeline. To change processing, edit this function.

% process neural data
if Neuro.Blackrock,
    Neuro = ReadBR(Neuro);
    Neuro = RefNeuralData(Neuro);
    if Neuro.UpdateChStatsFlag,
        Neuro = UpdateChStats(Neuro);
    end
    if Neuro.ZscoreRawFlag,
        Neuro = ZscoreChannels(Neuro);
    end
    Neuro = ApplyFilterBank(Neuro);
    Neuro = UpdateNeuroBuf(Neuro);
    Neuro = CompNeuralFeatures(Neuro);
    if Neuro.UpdateFeatureStatsFlag,
        Neuro = UpdateFeatureStats(Neuro);
    end
    if Neuro.ZscoreFeaturesFlag,
        Neuro = ZscoreFeatures(Neuro);
    end
    
    % smoothing option
    if Neuro.SmoothDataFlag
        Neuro = SmoothNeuro(Neuro);
    end
    
end

% override neural data if generating neural features
if Params.GenNeuralFeaturesFlag,
    Neuro.NeuralFeatures = VelToNeuralFeatures(Params);
end

% dimensionality reduction on neural features
Neuro.MaskedFeatures = Neuro.NeuralFeatures(Neuro.FeatureMask);
if Neuro.DimRed.Flag,
    %Neuro.NeuralFactors = Neuro.DimRed.F(Neuro.NeuralFeatures);
    Neuro.NeuralFactors = Neuro.DimRed.F(Neuro.MaskedFeatures);
end


varargout{1} = Neuro;

% if Data exists and is not empty, fill structure
if exist('Data','var') && ~isempty(Data),
    if Neuro.Blackrock,
        Data.NeuralTimeBR(1,end+1) = Neuro.TimeStamp;
        Data.NeuralSamps(1,end+1) = Neuro.NumSamps;
        if Neuro.SaveRaw,
            Data.BroadbandData{end+1} = Neuro.BroadbandData;
            Data.Reference{end+1} = Neuro.Reference;
        end
        if Neuro.SaveProcessed,
            Data.ProcessedData{end+1} = Neuro.FilteredData;
        end
    end
    
    Data.NeuralFeatures{end+1} = Neuro.NeuralFeatures;
    Data.SmoothedNeuralFeatures{end+1} = Neuro.FilteredFeatures;
    if Neuro.DimRed.Flag,
        Data.NeuralFactors{end+1} = Neuro.NeuralFactors;
    end
    
    varargout{2} = Data;
end

end % NeuroPipeline