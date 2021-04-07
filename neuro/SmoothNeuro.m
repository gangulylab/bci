function Neuro = SmoothNeuro(Neuro)
% function Neuro = SmoothNeuro(Neuro)

Neuro.FeatureDataBuffer = circshift(Neuro.FeatureDataBuffer,1);
Neuro.FeatureDataBuffer(1,:) = Neuro.NeuralFeatures';
Neuro.FilteredFeatures = (mean(Neuro.FeatureDataBuffer,1))';

end