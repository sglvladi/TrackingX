RadarName = 'staddon_clutter';
RadarCoords.lat = 50.346069;
RadarCoords.lon = -4.113670;

RadarCoords_lat = 50.339026;
RadarCoords_lon = -4.126307;

% Load dataset
load(strcat(RadarName,'.mat'));
filename = 'smc_phd_clutter_fakellr.mat';
load(filename);
% load('smc_phd_no_clutter.mat');
% load('smc_phd_clutter_fakellr.mat');

% lon_lim = [-4.2637   -4.1024];
% lat_lim = [50.2699   50.3853];

lon_lim = [-4.176846509809199,-4.103030078951917];
lat_lim = [50.297532151330770,50.337022151330770];
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

%% Plot all existing tracks

for j=1:numel(Tracks)
    
    track = Tracks{j};
    if(strcmp(filename, 'smc_phd_clutter_fakellr_0.mat'))
        t = [183531,328072];
        if(find(t==track.Tag.ID,1))
            continue
        end
    else
    end
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
    t = [766616, 128226, 592651, 237840];
    if (find(t == track.Tag.ID,1))
        plots(end+1) = plot(ax(1), lon, lat,'-','LineWidth',4,'Color', [0.9856, 0.7372, 0.2537]);
        plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',1);
    else
        plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',2);
    end
%     plots(end+1) = plot(ax(1), lon, lat,'-.w','LineWidth',2);
    plots(end+1) = plot(ax(1), lon(:,end),lat(:,end),'rx','MarkerSize',15,'LineWidth',2);
%     text(ax(1), lon(:,end), lat(:,end), num2str(track.Tag.ID),'FontSize',18,'Color','g');
end

% Plot extra track
means = [-4.12069846817424,50.3000236150299;-4.12101955629663,50.2999210786118;-4.12152794582374,50.2999381680148;-4.12192930597672,50.2999381680148;-4.12227715144264,50.2999210786118;-4.12259823956503,50.2999039892088;-4.12310662909214,50.2999210786118;-4.12348123190159,50.2999381680148;-4.12484585642174,50.3000065256269;-4.12588939281949,50.3000919726419;-4.12698644390432,50.3000748832389;-4.12741456140083,50.3000748832389;-4.12768213483616,50.3000236150299;-4.12800322295854,50.3000065256269;-4.12811025233267,50.3000065256269;-4.12843134045506,50.3000236150299;-4.12859188451625,50.3000236150299;-4.12979596497520,50.2999894362238;-4.13033111184585,50.2999552574178;-4.13078598668589,50.2999381680148;-4.13094653074709,50.2999723468208;-4.13142816293067,50.3000065256269;-4.13150843496126,50.2999894362238;-4.13169573636599,50.2999552574178;-4.13207033917544,50.2999723468208;-4.13241818464136,50.2999894362238;-4.13273927276375,50.2999723468208;-4.13289981682494,50.2999552574178;-4.13314063291673,50.2999039892088;-4.13343496369559,50.2999381680148;-4.13383632384857,50.2999894362238;-4.13410389728389,50.2999894362238;-4.13455877212394,50.3000065256269;-4.13466580149807,50.2999723468208;-4.13501364696399,50.2999723468208;-4.13525446305578,50.2999381680148;-4.13576285258289,50.2998868998058;-4.13616421273587,50.2999723468208;-4.13664584491946,50.3000236150299;-4.13702044772891,50.3000065256269;-4.13768938131721,50.3000065256269;-4.13833155756199,50.3000748832389;-4.13883994708910,50.3000748832389;-4.13913427786795,50.3000748832389;-4.13926806458562,50.3001090620449;-4.13942860864681,50.3000919726419;-4.13996375551745,50.3000919726419;-4.14028484363984,50.3000919726419;-4.14060593176223,50.3001432408509;-4.14098053457168,50.3000748832389;-4.14138189472466,50.3000577938359;-4.14178325487764,50.3000577938359;-4.14253246049655,50.3000236150299;-4.14376329829903,50.3000065256269;-4.14411114376495,50.2998868998058;-4.14429844516968,50.2998698104028;-4.14453926126147,50.2998014527907;-4.14496737875798,50.2998185421937;-4.14520819484977,50.2998868998058;-4.14582361375101,50.2999381680148;-4.14595740046867,50.3000065256269;-4.14622497390400,50.2999552574178;-4.14646578999579,50.2999381680148;-4.14689390749230,50.2999381680148;-4.14726851030175,50.2999210786118;-4.14745581170648,50.2999210786118;-4.14780365717240,50.2999381680148;-4.14804447326419,50.2999723468208;-4.14823177466892,50.3000065256269;-4.14857962013483,50.2999723468208;-4.14879367888309,50.2999723468208;-4.14906125231841,50.2999894362238;-4.14940909778433,50.3000407044329;-4.14964991387612,50.3000407044329;-4.14991748731145,50.2999723468208;-4.15055966355622,50.3001603302539;-4.15104129573980,50.3000748832389;-4.15144265589278,50.3002628666720;-4.15200456010696,50.3004166712991;-4.15283403775646,50.3003654030901;-4.15334242728357,50.3003654030901;-4.15371703009302,50.3003824924931;-4.15414514758954,50.3004337607021;-4.15438596368133,50.3004337607021]';
plots(end+1) = plot(ax(1), means(1,1),means(2,1),'go','MarkerSize',15,'LineWidth',2);
plots(end+1) = plot(ax(1), means(1,:), means(2,:),'-','LineWidth',4,'Color', [0.9856, 0.7372, 0.2537]);
plots(end+1) = plot(ax(1), means(1,:), means(2,:),'-.w','LineWidth',1);
plots(end+1) = plot(ax(1), means(1,end),means(2,end),'rx','MarkerSize',15,'LineWidth',2);
means = [-4.15651854566730,50.3007257705160;-4.15696081580833,50.3007386100374;-4.15728246681998,50.3007771286017;-4.15760411783164,50.3007642890803;-4.15776494333747,50.3008156471659;-4.15846855492547,50.3008926842945;-4.15877010274890,50.3009183633373;-4.15921237288993,50.3009312028587;-4.15963453984273,50.3009697214230;-4.15999639723084,50.3009825609444;-4.16029794505427,50.3010467585515;-4.16045877056010,50.3010724375943;-4.16074021519530,50.3010595980729;-4.16132320765393,50.3010467585515;-4.16146392997153,50.3010724375943;-4.16232836706536,50.3011366352014;-4.16279074039461,50.3011879932871;-4.16323301053564,50.3012393513728;-4.16345414560616,50.3012136723300;-4.16453971777050,50.3013677465870;-4.16530363892319,50.3013677465870;-4.16578611544067,50.3015346603655;-4.16606756007587,50.3015860184511]';
plots(end+1) = plot(ax(1), means(1,1),means(2,1),'go','MarkerSize',15,'LineWidth',2);
plots(end+1) = plot(ax(1), means(1,:), means(2,:),'-','LineWidth',4,'Color', [0.9856, 0.7372, 0.2537]);
plots(end+1) = plot(ax(1), means(1,:), means(2,:),'-.w','LineWidth',1);
plots(end+1) = plot(ax(1), means(1,end),means(2,end),'rx','MarkerSize',15,'LineWidth',2);

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
pos = [-4.166866227065740,50.308021682580770,0.003389529988345,0.002102476851853]; 
rectangle(ax(1), 'Position', pos, 'LineWidth', 2, 'EdgeColor', 'c');
pos = [-4.162648145302467,50.302582666377070,0.002937592656566,0.001919652777772]; 
rectangle(ax(1), 'Position', pos, 'LineWidth', 2, 'EdgeColor', 'c');
pos = [-4.157676834652895,50.299840305265950,0.004293404651904,0.001508298611114]; 
rectangle(ax(1), 'Position', pos, 'LineWidth', 2, 'EdgeColor', 'c');
%% PLot false tracks
pos = [-4.129054136973541,50.326578326099290,0.002937592656566,0.001873946759261];
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
legend([h4, h5, h6, h3, h7, h9], 'High-clutter Region', 'Track Birth', 'Track Death','False Tracks', 'Fragmentation Examples', 'Bonus Tracks','Location','northeast');
% legend([h1, h3, h5, h6], 'Detections', 'Obscured Region', 'Track Birth', 'Track Death','Location','southeast');
% legend([h1, h2, h3, h4], 'Detections', 'Radar Position', 'Obscured Region', 'High-clutter region');
% legend([h1, h3, h5, h6, h7, h9], 'Detections', 'Obscured Region', 'Track Birth', 'Track Death', 'Fragmentation Examples', 'Bonus Tracks','Location','southeast');
% legend([h1, h2, h3, h4], 'Detections', 'Radar Position', 'Obscured Region', 'High-clutter region');
% xlabel('Longitude');
% ylabel('Latitude');
set(ax(1),'YTickLabel',[]);
set(ax(1),'XTickLabel',[]);
