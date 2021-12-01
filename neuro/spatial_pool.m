function data_temporal_spatial = spatial_pool(inputData,chmap)
    gridMapData = zeros(8,16,200);
    for j = 1:size(inputData,2)
       oneChannel = inputData(:,j); 
       [row, column, ~] = find(chmap == j);
       gridMapData(row, column,:)=oneChannel;
    end
    dlX = dlarray(gridMapData,'SSCB'); % Transform to DL format
    dlY = avgpool(dlX, [2 2], 'Stride',2); % Perform avgpooling
    poolingAfter = extractdata(dlY); % extract from DL format to regular
    counter = 1;
    data_temporal_spatial = zeros(200,32);
    for j = 1:size(poolingAfter,1)
       for k = 1:size(poolingAfter,2)
           temp_data = squeeze(poolingAfter(j,k,:));
           data_temporal_spatial(:,counter)=temp_data;
           counter = counter + 1;
       end
    end
    
end

