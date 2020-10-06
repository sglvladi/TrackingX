clear all;
%% Reader config
config.Filename = "D:\Dropbox\University of Liverpool\PhD\Workspace\Track Analytics\StoneSoup\examples\track_analytics\data\exact_earth\id\exactEarth_historical_data_2017-08-10_id.csv";
config.StateFields = {'Longitude', 'Latitude'};
config.TimeField = 'Time';
config.TimeFormat = 'yyyyMMdd_HHmmSS';
% config.TimeFormat = 'yyyy-MM-dd''T''HH:mm:ssX';
config.TimeZone = 'UTC';

%% Initialise reader
reader = CsvMeasurementReaderX(config);

%% Read from reader
measurements = reader.read();