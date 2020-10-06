RadarName = 'staddon_clutter';
RadarCoords.lat = 50.346069;
RadarCoords.lon = -4.113670;

RadarCoords_lat = 50.339026;
RadarCoords_lon = -4.126307;

% Load dataset
load(strcat(RadarName,'.mat'));
load('gnn_llr_clutter_a_fakephd.mat');

lon_lim = [-4.176846509809199,-4.103030078951917];
lat_lim = [50.297532151330770,50.337022151330770];

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


%% PLot detections
for k=2:N
    
%     MeasurementList = MeasurementScans(k);
%     if(MeasurementList.NumMeasurements>0)
%         % Convert measurements to LLA and plot them
%     %     [lat,lon,~] = ned2geodetic(MeasurementList.Vectors(2,:),...
%     %                                MeasurementList.Vectors(1,:),...
%     %                                0,...
%     %                                RadarCoords.lat,...
%     %                                RadarCoords.lon,...
%     %                                0,...
%     %                                referenceEllipsoid('wgs84'));
%         latlon = [MeasurementList.Metadata.LatLon];
%         lat = latlon(1,:);
%         lon = latlon(2,:);
%         plots(end+1) = plot(ax(1), lon,lat,'y.','MarkerSize', 4);
%     end
end

%% Plot all existing tracks

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

    plots(end+1) = plot(ax(1), lon(:,1),lat(:,1),'go','MarkerSize',15,'LineWidth',2);
    t = [51806, 6423, 8579, 82699, 37841];
    td = [28571, 57088, 46505, 50816];
    if (find(t == track.Tag.ID,1))
        plots(end+1) = plot(ax(1), lon, lat,'-','LineWidth',4,'Color', [0.9856, 0.7372, 0.2537]);
        plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',1);
    elseif(find(td==track.Tag.ID,1))
        % Delayed initiation
        plots(end+1) = plot(ax(1), lon, lat,'-','LineWidth',4,'Color', [0.0265, 0.6137, 0.8135]);
        plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',1);
    else
        plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',2);
    end
    plots(end+1) = plot(ax(1), lon(:,end),lat(:,end),'rx','MarkerSize',15,'LineWidth',2);
%     text(ax(1), lon(:,1), lat(:,1), num2str(track.Tag.ID),'FontSize',18,'Color','g');

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
%     end
    
end

%% Plot region of high clutter
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
%% Draw fragmented tracks
pos = [-4.129656720082580,50.304822261284470,0.004444050429163,0.002011064814816]; 
rectangle(ax(1), 'Position', pos, 'LineWidth', 2, 'EdgeColor', 'c');
pos = [-4.167167518620260,50.307747446469660,0.003464852876975,0.002513831018518]; 
rectangle(ax(1), 'Position', pos, 'LineWidth', 2, 'EdgeColor', 'c');
pos = [-4.169427205279155,50.312455166377070,0.003163561322454,0.002056770833327]; 
rectangle(ax(1), 'Position', pos, 'LineWidth', 2, 'EdgeColor', 'c');
pos = [-4.129656720082580,50.302994020543736,0.002485655324787,0.001508298611107]; 
rectangle(ax(1), 'Position', pos, 'LineWidth', 2, 'EdgeColor', 'c');
% x = [-4.12945334365527,-4.12917944224207,-4.12739908305628,-4.12753603376287,-4.12945334365527];
% y = [50.3032914000914,50.3022941778692,50.3023772797210,50.3034576037951,50.3032914000914];
% plot(ax(1), x,y, 'c-','LineWidth', 2);       
%% PLot false tracks
pos = [-4.130334626080249,50.326669738136324,0.004293404651903,0.004753425925927];
rectangle(ax(1), 'Position', pos, 'LineWidth', 2, 'EdgeColor', 'm');

%% Legend
h1 = plot(NaN, NaN, 'o','MarkerSize', 5, 'MarkerFaceColor', 'y', 'MarkerEdgeColor','k');
h2 = plot(NaN, NaN, 's','MarkerSize',10,...
                               'MarkerEdgeColor','red',...
                               'MarkerFaceColor',[1 .6 .6]);
h3 = plot(NaN, NaN, 'ms', 'MarkerSize',15,'LineWidth', 2);
h4 = plot(NaN, NaN, 'r-', 'LineWidth', 2);
h5 = plot(NaN, NaN, 'go','MarkerSize',15,'LineWidth',2);
h6 = plot(NaN, NaN, 'rx','MarkerSize',15,'LineWidth',2);
h7 = plot(NaN, NaN, 'cs','MarkerSize',15,'LineWidth',2);
h8 = plot(NaN, NaN, 'b-','MarkerSize',15,'LineWidth',4);
h9 = plot(NaN, NaN, '-','LineWidth',4,'Color', [0.9856, 0.7372, 0.2537]);
h10 = plot(NaN, NaN, '-', 'LineWidth',4, 'Color', [0.0265, 0.6137, 0.8135]);
legend([h4, h5, h6, h3, h7, h9, h10], 'High-clutter Region', 'Track Birth', 'Track Death','False Tracks', 'Fragmentation Examples', 'Bonus Tracks', 'Delayed Initiation', 'Location','northeast');
% legend([h4, h5, h6], 'High-clutter Region', 'Track Birth', 'Track Death','Location','northeast');
% legend([h1, h3, h5, h6, h9], 'Detections', 'Obscured Region', 'Track Birth', 'Track Death', 'Delayed Initiation', 'Location','southeast');
% legend([h1, h2, h3, h4], 'Detections', 'Radar Position', 'Obscured Region', 'High-clutter region');
% legend([h1, h3, h5, h6, h7, h8], 'Detections', 'Obscured Region', 'Track Birth', 'Track Death', 'Fragmentation Examples', 'Jump-track','Location','southeast');
% legend([h1, h2, h3, h4], 'Detections', 'Radar Position', 'Obscured Region', 'High-clutter region');
% xlabel('Longitude');
% ylabel('Latitude');
set(ax(1),'YTickLabel',[]);
set(ax(1),'XTickLabel',[]);
