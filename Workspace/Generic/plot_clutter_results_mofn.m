RadarName = 'staddon_clutter';
RadarCoords.lat = 50.346069;
RadarCoords.lon = -4.113670;

RadarCoords_lat = 50.339026;
RadarCoords_lon = -4.126307;

% Load dataset
load(strcat(RadarName,'.mat'));
filename = 'gnn_mon_clutter_16_20.mat';
load(filename);

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
t = [368, 912, 640, 673, 146, 594, 342, 182, 912, 505, 672, 195, 594, 342, ...
    798, 267, 54, 794, 826, 125, 767, 524, 63, 200, 754, 244, 74, 777, 969, ...
    77, 51, 74, 828, 307, 861, 97, 807, 931, 665, 404, 177, 938, 225, 392, ...
    985, 585, 608, 319, 641, 915, 661, 538, 517, 586, 27, 164];

%% PLot detections
% for k=2:N
%     
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
% end
% x = [-469.710602736631,-2572.92295458513,-2516.31234127506,-1320.74163107899,-835.307363807233, -469.710602736631];
% y = [-4915.41679013513,-4065.50122858976,-2858.60918927224,-1380.19032924249,-1261.08679763585, -4915.41679013513];
% %         x = [-2858.62556879160,-2856.97056906902,-108.925126527319,-108.988225201227,-2858.62556879160];
% %         y = [-4044.74840882411,-974.655039469575,-975.424226884065,-4045.51804181679,-4044.74840882411];
% [y,x,~] = ned2geodetic(y,...
%                        x,...
%                        0,...
%                        RadarCoords.lat,...
%                        RadarCoords.lon,...
%                        0,...
%                        referenceEllipsoid('wgs84'));
% plots(end+1) = plot(ax(1), x, y, 'r-', 'LineWidth', 2);
% h1 = plot(NaN, NaN, 'o','MarkerSize', 5, 'MarkerFaceColor', 'y', 'MarkerEdgeColor','k');
% h4 = plot(NaN, NaN, 'r-', 'LineWidth', 2);
% legend([h1, h4], 'Detections', 'High-clutter Region','Location','northeast');
% 
% set(ax(1),'YTickLabel',[]);
% set(ax(1),'XTickLabel',[]);
%% Plot all existing tracks

for j=1:numel(Tracks)
    track = Tracks{j};
%     if find(t==track.Tag.ID,1)
%         continue
%     end

    means = [track.Trajectory.Mean];
    if strcmp(filename, 'gnn_mon_clutter_24_30.mat')
        t = [656944, 236783, 967201];
        if(find(t==track.Tag.ID,1))
            continue
        end
        
        if(track.Tag.ID == 663897)
            means = means(:, 100:end);
        end
        if(track.Tag.ID == 260895)
            means = means(:, 100:end);
        end
        if(track.Tag.ID == 993690)
            means = means(:, 50:end);
        end
    elseif strcmp(filename, 'gnn_mon_clutter_16_20.mat')
        t = [324379];
        if(find(t==track.Tag.ID,1))
            continue
        end
        
        if(track.Tag.ID == 472513)
            means = means(:, 50:end);
        end
        if(track.Tag.ID == 378804)
            means = means(:, 50:end);
        end
    end
    % Convert track trajectory to LLA and plot it
    [lat,lon,~] = ned2geodetic(means(3,:),...
                               means(1,:),...
                               0,...
                               RadarCoords.lat,...
                               RadarCoords.lon,...
                               0,...
                               referenceEllipsoid('wgs84'));
    traj_length = size(lon,2);

    plots(end+1) = plot(ax(1), lon(:,1),lat(:,1),'go','MarkerSize',15,'LineWidth',2);
    
    t = [472513, 378804, 492969];
%     tf = [
    if(find(t==track.Tag.ID,1))
        % Delayed initiation
        plots(end+1) = plot(ax(1), lon, lat,'-','LineWidth',4,'Color', [0.0265, 0.6137, 0.8135]);
        plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',1);
    elseif(track.Tag.ID==946079)
        plots(end+1) = plot(ax(1), lon, lat,'b-','LineWidth',4);
        plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',1);
    else
        plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',2);
    end 
%     plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',2);

%     text(ax(1), lon(:,1), lat(:,1), num2str(track.Tag.ID),'FontSize',18,'Color','g')
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
%     end
    
end

%% Plot jump track
% Convert track trajectory to LLA and plot it
% means = [jump_track.Trajectory.Mean];
% [lat,lon,~] = ned2geodetic(means(3,:),...
%                            means(1,:),...
%                            0,...
%                            RadarCoords.lat,...
%                            RadarCoords.lon,...
%                            0,...
%                            referenceEllipsoid('wgs84'));
% traj_length = size(lon,2);
% 
% plots(end+1) = plot(ax(1), lon(:,1),lat(:,1),'go','MarkerSize',15,'LineWidth',3);
% plots(end+1) = plot(ax(1), lon, lat,'-b','LineWidth',4);
% plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',1);
% 
% %     text(ax(1), lon(:,end), lat(:,end), num2str(track.Tag.ID),'FontSize',18,'Color','g')
% plots(end+1) = plot(ax(1), lon(:,end),lat(:,end),'rx','MarkerSize',15,'LineWidth',2);

%% Plot region of low PD
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
pos = [-4.165811706624922,50.304959379340030,0.002335009547527,0.001416886574070]; 
rectangle(ax(1), 'Position', pos, 'LineWidth', 2, 'EdgeColor', 'c');
x = [-4.16950252816779,-4.16701687284300,-4.16468186329547,-4.16724284150889,-4.16950252816779];
y = [50.3142377010993,50.3099870413771,50.3105355135993,50.3147861733215,50.3142377010993];
plot(ax(1), x,y, 'c-','LineWidth', 2);
x = [-4.16658235867053,-4.16414724790782,-4.16222292442139,-4.16457242597750,-4.16658235867053];
y = [50.3096604413594,50.3057670085599,50.3062977172654,50.3101529840629,50.3096604413594];
plot(ax(1), x,y, 'c-','LineWidth', 2);
x = [-4.16400395729781,-4.16317540552288,-4.16121701041850,-4.16197023930480,-4.16400395729781];
y = [50.3055535575808,50.3047308492474,50.3052336154511,50.3062848538771,50.3055535575808];
plot(ax(1), x,y, 'c-','LineWidth', 2);
x = [-4.16264814530247,-4.15940926109138,-4.15737554309838,-4.16046378153220,-4.16264814530247];
y = [50.3040909649882,50.3002516594326,50.3008915436919,50.3047765552660,50.3040909649882];
plot(ax(1), x,y, 'c-','LineWidth', 2);
x = [-4.12945334365527,-4.12917944224207,-4.12739908305628,-4.12753603376287,-4.12945334365527];
y = [50.3032914000914,50.3022941778692,50.3023772797210,50.3034576037951,50.3032914000914];
plot(ax(1), x,y, 'c-','LineWidth', 2);
%% PLot false tracks
pos = [-4.133422864514073,50.325389969617810,0.009038746635585,0.006535960648144];
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
% h9 = plot(NaN, NaN, '-','LineWidth',4,'Color', [0.9856, 0.7372, 0.2537]);
h9 = plot(NaN, NaN, '-', 'LineWidth',4, 'Color', [0.0265, 0.6137, 0.8135]);
legend([h4, h5, h6, h3, h7, h8, h9], 'High-clutter Region', 'Track Birth', 'Track Death', 'False Tracks', 'Fragmentation Examples', 'Jump Track', 'Delayed Initiation', 'Location','northeast');
% legend([h1, h2, h3, h4], 'Detections', 'Radar Position', 'Obscured Region', 'High-clutter region');
% legend([h4, h5, h6, h7, h8, h9], 'High-clutter Region', 'Track Birth', 'Track Death', 'Fragmentation Examples', 'Jump Track', 'Bonus Tracks','Location','northeast');
% legend([h1, h2, h3, h4], 'Detections', 'Radar Position', 'Obscured Region', 'High-clutter region');
% xlabel('Longitude');
% ylabel('Latitude');
set(ax(1),'YTickLabel',[]);
set(ax(1),'XTickLabel',[]);
