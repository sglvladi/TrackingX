clear all;
%% Reader config
%config.Filename = "D:\Dropbox\University of Liverpool\PhD\Workspace\Track Analytics\NelsonRadarAisTracking\data\ais_2017-01-23_2017-01-28.csv";
config.Filename = "C:\Users\sglvladi\Documents\GitHub\Stone-Soup-TA\examples\track_analytics\data\exact_earth\exactEarth_historical_data_2017-08-10.csv";
config.StateFields = {'Longitude', 'Latitude'};
config.TimeField = 'Time';
config.TimeFormat = 'yyyyMMdd_HHmmss';
%config.TimeFormat = 'yyyy-MM-dd''T''HH:mm:ssX';
config.TimeZone = 'UTC';

%% Initialise reader
reader = CsvMeasurementReaderX(config);

%% Read from reader
measurements = reader.read();