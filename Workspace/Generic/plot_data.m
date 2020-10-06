RadarName = 'staddon_detection';
RadarCoords.lat = 50.346069;
RadarCoords.lon = -4.113670;

RadarCoords_lat = 50.339026;
RadarCoords_lon = -4.126307;

% Load dataset
load(strcat(RadarName,'.mat'));
% load('jpda_mon_no_detection.mat');

% lon_lim = [-4.2637   -4.1024];
% lat_lim = [50.2699   50.3853];

lon_lim_obs = [-4.16737200000000,-4.14897200000000];
lat_lim_obs = [50.3537860000000,50.3617750000000];
% radar_meas_convert;


% Map plot
figure('units','normalized','outerposition',[0 0 .5 0.8])
ax(1) = gca;
plot_google_map('Axis',ax(1),'APIKey','AIzaSyBXKujdtXRZiqya1soVS9pxBzYR4g7aGvM','Resize',3,'Scale',2,'MapType','satellite');
axis(ax(1),[lon_lim(1) lon_lim(2) lat_lim(1) lat_lim(2)])
plots = [];
hold on;

N = size(MeasurementScans,2); % Simulation length

Tracks = [TrackList, DeletedTracks];
t = [368, 912, 640, 673, 146, 594, 342, 182, 912, 505, 672, 195, 594, 342, ...
    798, 267, 54, 794, 826, 125, 767, 524, 63, 200, 754, 244, 74, 777, 969, ...
    77, 51, 74, 828, 307, 861, 97, 807, 931, 665, 404, 177, 938, 225, 392, ...
    985, 585, 608, 319, 641, 915, 661, 538, 517, 586, 27, 164];
for k=2:N
    
    MeasurementList = MeasurementScans(k);
    if(MeasurementList.NumMeasurements>0)
        % Convert measurements to LLA and plot them
    %     [lat,lon,~] = ned2geodetic(MeasurementList.Vectors(2,:),...
    %                                MeasurementList.Vectors(1,:),...
    %                                0,...
    %                                RadarCoords.lat,...
    %                                RadarCoords.lon,...
    %                                0,...
    %                                referenceEllipsoid('wgs84'));
        latlon = [MeasurementList.Metadata.LatLon];
        lat = latlon(1,:);
        lon = latlon(2,:);
        plots(end+1) = plot(ax(1), lon,lat,'y.','MarkerSize', 4);
    end
end

% Plot all existing tracks
for j=1:numel(Tracks)
    track = Tracks{j};
%     if find(t==track.Tag.ID,1)
%         continue
%     end
    % Convert track trajectory to LLA and plot it
    means = [track.Trajectory.Mean];
    [lat,lon,~] = ned2geodetic(means(3,:),...
                               means(1,:),...
                               0,...
                               RadarCoords.lat,...
                               RadarCoords.lon,...
                               0,...
                               referenceEllipsoid('wgs84'));
    traj_length = size(lon,2);
    if(traj_length>NumPersistFrames)
        start = traj_length-NumPersistFrames;
    else
        start = 1;
    end
    plots(end+1) = plot(ax(1), lon(:,start),lat(:,start),'go','MarkerSize',15,'LineWidth',2);
    plots(end+1) = plot(ax(1), lon(:,start:end),lat(:,start:end),'-.w','LineWidth',2);
    plots(end+1) = plot(ax(1), lon(:,end),lat(:,end),'rx','MarkerSize',15,'LineWidth',2);

    % Convert track velocity to LLA and plot it
%     [lat_vel,lon_vel,~] = ned2geodetic(track.State.Mean(4,end),...
%                                        track.State.Mean(2,end),...
%                                        0,...
%                                        lat(:,end),...
%                                        lon(:,end),...
%                                        0,...
%                                        referenceEllipsoid('wgs84'));
%     lat_vel = lat_vel-lat(:,end);
%     lon_vel = lon_vel-lon(:,end);
%     plots(end+1) = quiver(ax(1), lon(:,end),lat(:,end),20*lon_vel,20*lat_vel,'r','LineWidth',1.5);
%     if ShowTrackInfo
%     text(ax(1), lon(:,end), lat(:,end), num2str(track.Tag.ID),'FontSize',18,'Color','g') 
%     end
    
end

% Plot region of low PD
x = [-3422.85261310741,-3191.58441229017,-2993.55608122595,-3228.08558459284, -3422.85261310741];
y = [1.472966334316646e+03,1.123217741935625e+03,1.243952901078573e+03,1.579667558371468e+03, 1.472966334316646e+03];
[y,x,~] = ned2geodetic(y,...
                       x,...
                       0,...
                       RadarCoords.lat,...
                       RadarCoords.lon,...
                       0,...
                       referenceEllipsoid('wgs84'));
plots(end+1) = plot(ax(1), x, y, 'm-', 'LineWidth', 2);

% Plot region of high clutter
x = [-469.710602736631,-2572.92295458513,-2516.31234127506,-1320.74163107899,-835.307363807233, -469.710602736631];
y = [-4915.41679013513,-4065.50122858976,-2858.60918927224,-1380.19032924249,-1261.08679763585, -4915.41679013513];
%         x = [-2858.62556879160,-2856.97056906902,-108.925126527319,-108.988225201227,-2858.62556879160];
%         y = [-4044.74840882411,-974.655039469575,-975.424226884065,-4045.51804181679,-4044.74840882411];
[y,x,~] = ned2geodetic(y,...
                       x,...
                       0,...
                       RadarCoords.lat,...
                       RadarCoords.lon,...
                       0,...
                       referenceEllipsoid('wgs84'));
plots(end+1) = plot(ax(1), x, y, 'r-', 'LineWidth', 2);

plot(ax(1), RadarCoords_lon,RadarCoords_lat,...
                               '-s','MarkerSize',10,...
                               'MarkerEdgeColor','red',...
                               'MarkerFaceColor',[1 .6 .6]);
                           
h1 = plot(NaN, NaN, 'o','MarkerSize', 5, 'MarkerFaceColor', 'y', 'MarkerEdgeColor','k');
h2 = plot(NaN, NaN, 's','MarkerSize',10,...
                               'MarkerEdgeColor','red',...
                               'MarkerFaceColor',[1 .6 .6]);
h3 = plot(NaN, NaN, 'm-', 'LineWidth', 2);
h4 = plot(NaN, NaN, 'r-', 'LineWidth', 2);
h5 = plot(NaN, NaN, 'go','MarkerSize',15,'LineWidth',2);
h6 = plot(NaN, NaN, 'rx','MarkerSize',15,'LineWidth',2);
legend([h1, h2, h5, h6], 'Detections', 'Obscured Region', 'Track Birth', 'Track Death');
% legend([h1, h2, h3, h4], 'Detections', 'Radar Position', 'Obscured Region', 'High-clutter region');
% xlabel('Longitude');
% ylabel('Latitude');
set(ax(1),'YTickLabel',[]);
set(ax(1),'XTickLabel',[]);
