RadarName = 'staddon';
RadarCoords.lat = 50.346069;
RadarCoords.lon = -4.113670;

RadarCoords_lat = 50.339026;
RadarCoords_lon = -4.126307;

% Load dataset
load(strcat(RadarName,'.mat'));

lon_lim = [-4.2637   -4.1024];
lat_lim = [50.3516   50.3640];

lon_lim_obs = [-4.1674   -4.1490];
lat_lim_obs = [50.3516   50.3640];
% radar_meas_convert;
figure('units','normalized','outerposition',[0 0 1 1])
[ha, pos] = tight_subplot(1,2,[.0 .0],[.0 .0],[.0 .0]);
for i = 1:2 
    ax = ha(i);
    
    if i==1
        plot_google_map('Axis',ax,'APIKey','AIzaSyBXKujdtXRZiqya1soVS9pxBzYR4g7aGvM','Resize',3,'Scale',2,'MapType','satellite');
        axis(ax,[lon_lim(1) lon_lim(2) lat_lim(1) lat_lim(2)])
        plots = [];
        hold on;

        N = size(MeasurementScans,2); % Simulation length

        for k=2:N

            MeasurementList = MeasurementScans(k);
            % Convert measurements to LLA and plot them
            latlon = [MeasurementList.Metadata.LatLon];
            lat = latlon(1,:);
            lon = latlon(2,:);
            plot(ax, lon,lat,'y.','MarkerSize', 2);
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
        plot(ax, x, y, 'm-', 'LineWidth', 2);

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
        plot(ax, x, y, 'r-', 'LineWidth', 2);

        plot(ax, RadarCoords_lon,RadarCoords_lat,...
                                       '-s','MarkerSize',10,...
                                       'MarkerEdgeColor','red',...
                                       'MarkerFaceColor',[1 .6 .6]);
        x = [lon_lim_obs(1), lon_lim_obs(2), lon_lim_obs(2), lon_lim_obs(1), lon_lim_obs(1)];
        y = [lat_lim_obs(2), lat_lim_obs(2), lat_lim_obs(1), lat_lim_obs(1), lat_lim_obs(2)];
        plot(ax, x, y, 'w-', 'LineWidth', 2);
%         rectangle(ax, 'Position', [lon_lim_obs(1), lat_lim_obs(1), ...
%                                    lon_lim_obs(2)-lon_lim_obs(1), ...
%                                    lat_lim_obs(2)-lat_lim_obs(1)], ...
%                       'EdgeColor', 'w', 'LineWidth', 2)
%         plot(ax, [-4.167515391073365,-4.097638729563496], [50.361870612297780,50.379278070444585],...
%              'w--', 'LineWidth', 2);
%         plot(ax, [-4.16738343137302, -4.09740837419749], [50.3537539650738, 50.2830389777679],...
%              'w--', 'LineWidth', 2);
%         plot(ax, [-4.14892471073231, -4.09733683909138], [50.3618050077183, 50.3701207750788],...
%              'w--', 'LineWidth', 2);
%         plot(ax, [-4.14882330020060, -4.09677422728014], [50.3536549416168, 50.3336381730257],...
%              'w--', 'LineWidth', 2);
        h1 = plot(ax, NaN, NaN, 'o','MarkerSize', 5, 'MarkerFaceColor', 'y', 'MarkerEdgeColor','k');
        h2 = plot(ax, NaN, NaN, 's','MarkerSize',10,...
                                       'MarkerEdgeColor','red',...
                                       'MarkerFaceColor',[1 .6 .6]);
        h3 = plot(ax, NaN, NaN, 'm-', 'LineWidth', 2);
        h4 = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2);
        legend(ax, [h1, h2, h3, h4], 'Detections', 'Radar Position', 'Obscured Region', 'High-clutter region');
        xlabel('Longitude');
        ylabel('Latitude');
%         axis(ax,[lon_lim(1) lon_lim(2) lat_lim(1) lat_lim(2)])
    else
        plot_google_map('Axis',ax,'APIKey','AIzaSyBXKujdtXRZiqya1soVS9pxBzYR4g7aGvM','Resize',3,'Scale',2,'MapType','satellite');
        axis(ax,[lon_lim_obs(1) lon_lim_obs(2) lat_lim_obs(1) lat_lim_obs(2)])
        hold on;

        N = size(MeasurementScans,2); % Simulation length

        for k=2:N

            MeasurementList = MeasurementScans(k);
            % Convert measurements to LLA and plot them
            latlon = [MeasurementList.Metadata.LatLon];
            lat = latlon(1,:);
            lon = latlon(2,:);
            plot(ax, lon,lat,'y.','MarkerSize', 2);
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
        plot(ax, x, y, 'm-', 'LineWidth', 2);

        h1 = plot(ax, NaN, NaN, 'o','MarkerSize', 5, 'MarkerFaceColor', 'y', 'MarkerEdgeColor','k');
        h2 = plot(ax, NaN, NaN, 'm-', 'LineWidth', 2);
        legend(ax, [h1, h2, h2], 'Detections', 'Obscured Region');
%         xlabel('Longitude');
%         ylabel('Latitude');
    end
%     plot(ax, randn(10,i)); 
%     set(ax,'YTickLabel',[]);
end


