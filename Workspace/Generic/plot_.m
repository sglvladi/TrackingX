% TrackList_b = TrackList;
% TrackList = [TrackList_b, DeletedTracks];
TrackList = TrackList_b;
% Map plot
figure('units','normalized','outerposition',[0 0 .5 1])
ax(1) = gca;
if PlotLatLon
    plot_google_map('Axis',ax(1),'APIKey','AIzaSyBXKujdtXRZiqya1soVS9pxBzYR4g7aGvM','Resize',3,'Scale',2,'MapType','satellite');
    axis(ax(1),[lon_lim(1) lon_lim(2) lat_lim(1) lat_lim(2)])
end

% Delete all plots (other than map)
% for i = 1:numel(plots)
%     delete(plots(i))
% end
plots = [];
hold on;

if(MeasurementList.NumMeasurements>0)
    % Convert measurements to LLA and plot them
    [lat,lon,~] = ned2geodetic(MeasurementList.Vectors(2,:),...
                               MeasurementList.Vectors(1,:),...
                               0,...
                               RadarCoords.lat,...
                               RadarCoords.lon,...
                               0,...
                               referenceEllipsoid('wgs84'));
    plots(end+1) = plot(ax(1), lon,lat,'y*','MarkerSize', 10);
end

% Plot region of low PD
x = [-3560.53119790131, -3560.53119790131,  -3080.76025084906,  -3080.76025084906,  -3367.20269177091,  -3367.20269177091,  -3560.53119790131];
y = [1349.54778996881,  1590.19881386064,   1590.19881386064,   1108.03444256592,   1108.03444256592,   1349.54778996881,   1349.54778996881];
[y,x,~] = ned2geodetic(y,...
                       x,...
                       0,...
                       RadarCoords.lat,...
                       RadarCoords.lon,...
                       0,...
                       referenceEllipsoid('wgs84'));
plots(end+1) = plot(ax(1), x, y, 'm-');

% Plot region of high clutter
x = [-2858.62556879160,-2856.97056906902,-108.925126527319,-108.988225201227,-2858.62556879160];
y = [-4044.74840882411,-974.655039469575,-975.424226884065,-4045.51804181679,-4044.74840882411];
[y,x,~] = ned2geodetic(y,...
                       x,...
                       0,...
                       RadarCoords.lat,...
                       RadarCoords.lon,...
                       0,...
                       referenceEllipsoid('wgs84'));
plots(end+1) = plot(ax(1), x, y, 'r-');

% Plot all existing tracks
for j=1:numel(TrackList)
    track = TrackList{j};
    % Convert track trajectory to LLA and plot it
    means = [track.Trajectory.Mean];
    if PlotLatLon
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
        plots(end+1) = plot(ax(1), lon(:,start:end),lat(:,start:end),'-.w','LineWidth',2);
        plots(end+1) = plot(ax(1), lon(:,end),lat(:,end),'ws','MarkerSize',15);

        % Convert track velocity to LLA and plot it
        [lat_vel,lon_vel,~] = ned2geodetic(track.State.Mean(4,end),...
                                           track.State.Mean(2,end),...
                                           0,...
                                           lat(:,end),...
                                           lon(:,end),...
                                           0,...
                                           referenceEllipsoid('wgs84'));
        lat_vel = lat_vel-lat(:,end);
        lon_vel = lon_vel-lon(:,end);
        plots(end+1) = quiver(ax(1), lon(:,end),lat(:,end),20*lon_vel,20*lat_vel,'r','LineWidth',1.5);
        if ShowTrackInfo
            speed_kmph = sqrt(track.State.Mean(4,end)^2+track.State.Mean(2,end)^2)*3.6;
            speed_knot = speed_kmph/1.852;
            plots(end+1) = text(ax(1), lon(:,end)+0.001,lat(:,end)+0.00027,strcat("Sog:",num2str(speed_knot,2)," kt"),'FontSize',8,'Color','w');
            plots(end+1) = text(ax(1), lon(:,end)+0.001,lat(:,end)-0.00027,strcat("PoE:",num2str(track.ExistenceProbability*100,3)," %"),'FontSize',8,'Color','w');
        end
    else
        traj_length = size(lon,2);
        if(traj_length>NumPersistFrames)
            start = traj_length-NumPersistFrames;
        else
            start = 1;
        end
        plots(end+1) = plot(ax(1), means(1,start:end),means(3,start:end),'-.k','LineWidth',2);
        plots(end+1) = plot(ax(1), means(1,end),means(3,end),'ks','MarkerSize',15);

        x_vel = track.State.Mean(3,end)*cos(track.State.Mean(4,end));
        y_vel = track.State.Mean(3,end)*sin(track.State.Mean(4,end));
        plots(end+1) = quiver(ax(1), means(1,end),means(2,end),20*x_vel,20*y_vel,'r','LineWidth',1.5);

        if ShowTrackInfo
            speed_kmph = track.State.Mean(3,end)*3.6;
            speed_knot = speed_kmph/1.852;
            plots(end+1) = text(ax(1), means(1,end)+60,means(2,end)+50,strcat("Sog:",num2str(speed_knot,2)," kt"),'FontSize',8,'Color','k');
            plots(end+1) = text(ax(1), means(1,end)+60,means(2,end)-50,strcat("PoE:",num2str(track.ExistenceProbability*100,3)," %"),'FontSize',8,'Color','k');
        end
    end
end

% Add axis labels
%         xlabel('Longitude')
%         ylabel('Latitude')
title(datestr(timestamp_k));
%         set(gca,'visible','off')
pause(0.0001)