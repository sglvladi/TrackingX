MeasurementScans = MeasurementListX.empty();
measTagGen = UuidTagGeneratorX();
timestamp = datetime;
for i = 1:size(DataList,2)
    measurements = MeasurementX.empty();
    measurementVectors = DataList{i};
%     measurementVectors = measurementVectors(1:2,:);
    for j = 1:size(measurementVectors,2)
        tag = measTagGen.generate();
        measurements(end+1) = MeasurementX(measurementVectors(1:2,j), timestamp, tag);
        measurements(end).Metadata.LatLon = measurementVectors(3:4,j);
    end
    %mtags = measTagGen.generate(size(measurementVectors,2));
    if isempty(measurements)
        MeasurementScans(i) = MeasurementListX(measurementVectors,timestamp);
    else
        MeasurementScans(i) = MeasurementListX(measurements);%Vectors, timestamp, mtags);
    end
    timestamp = timestamp + seconds(2);
end
MeasurementScans_hard = MeasurementScans;