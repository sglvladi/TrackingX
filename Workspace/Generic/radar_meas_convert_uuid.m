MeasurementScans_new = MeasurementListX.empty();
for scan = MeasurementScans_hard
    measurements = MeasurementX.empty();
    for measurement = scan.Measurements
        latlon = measurement.Metadata.LatLon;
        if~(latlon(2,:)>lat_lim(2) || latlon(1,:)<lat_lim(1) || latlon(2,:)<lon_lim(1) || latlon(2,:)>lon_lim(2))
            measurements(end+1) = measurement;
        end
    end
    MeasurementScans_new (end+1) = MeasurementListX(measurements);
end
MeasurementScans = MeasurementScans_new;
    