%% Plot settings
ShowPlots = 1;
ShowPredict = 0;
ShowData = 1;
SkipFrames = 0;

%% Recording Settings
Record = 0;
clear F
clear M

%% Clutter settings
%lambdas = [0.0001,10,50,100,150,200,250,300,350,400,450,500];
%lambda_count = lambda_count + 1;
lambdaV = 100; %lambdas(lambda_count); % mean number of clutter points 
fprintf('STARTING NEW CYCLE for lambda = %d',lambdaV);
V_bounds = [0 10 0 10];
%V_bounds = [-1800 2400 -1500 1700];
%V_bounds = [-2000 2500 -1000 1700]; %[-700 -400 -700 400]; %[-2500 200 -3000 3000]; %[-2.5 .2 -3 3]; [-2 -.800 2 3] [-.700 -.400 -.700 .400]; % [x_min x_max y_min y_max]
V = (abs(V_bounds(2)-V_bounds(1))*abs(V_bounds(4)-V_bounds(3)));

%% Constant Velocity model
Params_cv.xDim = 4;
Params_cv.q = .0001;
CVmodel = ConstantVelocityModelX(Params_cv);

%% Constant Heading Model
Params_ch.q_vel = 2;
Params_ch.q_head = 0.3;
CHmodel = ConstantHeadingModelX(Params_ch);

%% 2D linear-Gaussian Observation Model
Params_po.xDim = 4;
Params_po.yDim = 2;
Params_po.r = .1;
POmodel = PositionalObsModelX(Params_po);

%% Polar2Cart Observation Model
% Params_meas.xDim = 4;
% Params_meas.R = [0.000304617419786709,0;0,25];%[0.00121846967914683, 0; 0, 100]; %[pi/3600 0 ; 0 0.05];
% POmodel = Polar2CartGaussModelX(Params_meas);

%% Assign PF parameter values
Params_pf.k               = 1;                   % initial iteration number
Params_pf.Np              = 5000;                 % number of particles
Params_pf.resampling_strategy = 'systematic_resampling';
Params_pf.DynModel = CVmodel;
Params_pf.ObsModel = POmodel;
Params_pf.particles_init = zeros(Params_pf.DynModel.Params.xDim, Params_pf.Np);
Params_pf.gen_x0 = @(x,Np) [mvnrnd(repmat(x',Np,1),CVmodel.Params.Q)]';
% Params_pf.gen_x0 = @(x,Np) [mvnrnd(repmat(x',Np,1),diag([CHmodel.Params.q_vel^2, CHmodel.Params.q_vel^2, 0, 0])), 10*rand(Np,1), 2*pi*rand(Np,1)]';
initFilter = ParticleFilterX(Params_pf);

%% Set TrackNum
TrackNum = 3;
TrueTracks = 3;
TrackIds = [];

%% Generate DataList                       (TrackNum, x_true, y_true, R, R_clutter, lambdaV, Iter) 
%[DataList,x1,y1] = gen_obs_cluttered_multi2(TrueTracks, x_true, y_true, POmodel.Params.r, V_bounds, lambdaV, 1);
[DataList, newGroundTruth] = dataGen(POmodel, GroundTruth, TrueTracks, lambdaV, V_bounds, [], 1);
% 
% %% Get GroundTruth
% for i=1:TrueTracks
%     GroundTruth{i} = [x_true(:,i), y_true(:,i)]; % ith target's GroundTruth
% end

%% Initiate PDAF parameters
Params_jpdaf.DataList = DataList{1}(:,:);  
Params_jpdaf.TrackList = [];
Params_jpdaf.pDetect = 0.9;
Params_jpdaf.pGate = 0.998; %0.999;
Params_jpdaf.gateLevel = chi2inv(Params_jpdaf.pGate,2);%14;

%% Instantiate JPDAF
jpdaf = JPDAFilterX(Params_jpdaf);
jpdaf.Params.TrackList = [];
jpdaf.Params.DataList = DataList{1}(:,:); 

%% Assign PHD parameter values
Params.k               = 1;                   % initial iteration number
Params.Np              = 50000;              % number of particles
Params.resampling_strategy = 'systematic_resampling'; % resampling strategy
Params.DynModel = CVmodel;
Params.ObsModel = POmodel;
Params.gen_x0 = @(Np)[abs(V_bounds(2)-V_bounds(1))*rand(Np,1)+V_bounds(1),abs(V_bounds(4)-V_bounds(3))*rand(Np,1)+V_bounds(3), zeros(Np,1), zeros(Np,1)] + mvnrnd(zeros(Np,4),CVmodel.Params.Q); % Uniform position and heading, Gaussian speedCVmodel.Params.Q)
%Params.gen_x0 = @(Np)[abs(V_bounds(2)-V_bounds(1))*rand(Np,1)+V_bounds(1),abs(V_bounds(4)-V_bounds(3))*rand(Np,1)+V_bounds(3), mvnrnd(zeros(Np,1),10), 2*pi*rand(Np,1)]; % Uniform position and heading, Gaussian speed
%Params.gen_x0 = @(Np) [(V_bounds(2)-V_bounds(1))*rand(Np,1),(V_bounds(4)-V_bounds(3))*rand(Np,1), mvnrnd(zeros(Np,1), CVmodel.Params.q^2), 2*pi*rand(Np,1)]; % Uniform position and heading, Gaussian speed
Params.particles_init = Params.gen_x0(Params.Np)'; % Generate inital particles as per gen_x0
Params.w_init = repmat(1/Params.Np, Params.Np, 1)'; % Uniform weights
Params.pBirth = 0.05;
Params.pDeath = 0.005;
Params.Jk = 500;
Params.pDetect = 0.9;
Params.lambda = lambdaV/V;
Params.pConf = 0.95;
Params.NpConf = Params_pf.Np;
Params.type = 'search';
Params.birth_strategy = 'mixture';

%% Instantiate PHD filter
myphd = SMC_PHDFilterX(Params);

%% Instantiate Log to store output
N=size(DataList,2);
Logs = [];
Log.xV = NaN(4,N);          %estmate        % allocate memory
Log.PV = NaN(4,4,N);
Log.ExistProbV = NaN(1,N);
Log.sV = NaN(2,N);          %actual
Log.zV = NaN(2,N);
Log.eV = NaN(2,N);

if(ShowPlots)
    img = imread('maze.png');

    % set the range of the axes
    % The image will be stretched to this.
    min_x = 0;
    max_x = 10;
    min_y = 0;
    max_y = 10;

    % make data to plot - just a line.
    x = min_x:max_x;
    y = (6/8)*x;

    figure('units','normalized','outerposition',[0 0 .5 1])
    ax(1) = gca;
    figure('units','normalized','outerposition',[.5 0 .5 1])
    ax(2) = gca;

    TrackIds = [];
end

if(exist('exec_time_phd','var'))
    exec_time_phd(:,end+1) = zeros(2,1);
else
    exec_time_phd = zeros(2,1);
end
for i = 1:N
    tic;
    fprintf('\nIteraration: %d/%d\n', i, N);
    
    %% Remove null measurements   
    DataList_k = DataList{i}(:,:);
    DataList_k( :, ~any(DataList_k,1) ) = [];
    
    %% Change JPDA and PHD filter parameters
    jpdaf.Params.DataList = DataList_k; % New observations
    myphd.Params.lambda = lambdaV/(pi*2500);%(pi/2*14);%V;
    myphd.Params.y = DataList_k; % New observations
    myphd.Params.k = 1; % Time index
    jpdaf.Params.k = 1; % Time index
%     for t = 1:length(jpdaf.Params.TrackList)
%         if(i==1)
%             myphd.Params.k = i; % Time index
%             jpdaf.Params.k = i; % Time index
%         else
%             myphd.Params.k = Dt(i-1); % Time index
%             jpdaf.Params.k = Dt(i-1); % Time index
%         end
%     end
    
    t1 = tic;
    %% 1) Predict the confirmed tracks
    jpdaf.Predict();
    %% 2) Update the confirmed tracks
    jpdaf.Update();
    %% 3) Perform track management
    t2 = tic;
    [jpdaf, myphd, TrackIds] = ExistProbPHDSearchX(jpdaf, myphd, initFilter, TrackIds);
    exec_time_phd(1,end) = exec_time_phd(1,end) + toc(t2);
    exec_time_phd(2,end) = exec_time_phd(2,end) + toc(t1);
    
    %% Store Logs
    % Expand Logs to fit new track
    TrackNums(i) = jpdaf.Params.nTracks;
    while(length(Logs)<length(TrackIds))
        Logs{end+1} = Log;
    end
    for j=1:jpdaf.Params.nTracks
        t = find(TrackIds==jpdaf.Params.TrackList{j}.trackId);
        Logs{t}.xV(:,i) = jpdaf.Params.TrackList{j}.TrackObj.Params.x;
        Logs{t}.PV(:,:,i) = weightedcov(jpdaf.Params.TrackList{j}.TrackObj.Params.particles',jpdaf.Params.TrackList{j}.TrackObj.Params.w);
        Logs{t}.ExistProbV(i) = jpdaf.Params.TrackList{j}.ExistProb;
    end
    %TrackNum_log(i) = TrackNum;
    if (ShowPlots)
        if(i==1 || rem(i,SkipFrames+1)==0)
            % Plot data
            cla(ax(1));
             % Flip the image upside down before showing it
            %imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            if(ShowData)
                Data = POmodel.obs_inv(0,DataList{i});
                data = plot(ax(1), Data(1,:),Data(2,:),'g*','MarkerSize', 5);
            end
            for j=1:TrueTracks
                h2 = plot(ax(1), x_true(1:i,j),y_true(1:i,j),'b.-','LineWidth',1);
                if j==2
                    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                h2 = plot(ax(1), x_true(i,j),y_true(i,j),'bo','MarkerSize', 10);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            end
            for j=1:length(Logs)
                colour = 'r';
%                 if(j==2)
%                    colour = 'c';
%                 elseif (j==3)
%                    colour = 'm';
%                 end
                h4 = plot(ax(1),Logs{j}.xV(1,:),Logs{j}.xV(2,:),'-r');
                if(~isnan(Logs{j}.xV(1,i)))
                    %h4 = plot(ax(1), Logs{j}.xV(1,:),Logs{j}.xV(2,:),strcat(colour,'.-'),'LineWidth',1);
                    h2 = plot_gaussian_ellipsoid(Logs{j}.xV(1:2,i), Logs{j}.PV(1:2,1:2,i), 1, [], ax(1));
                    %quiver(ax(1), Logs{j}.xV(1,i),Logs{j}.xV(2,i),30*Logs{j}.xV(3,i)*cos(Logs{j}.xV(4,i)),30*Logs{j}.xV(3,i)*sin(Logs{j}.xV(4,i)), 'Color', 'b', 'Linewidth', 1, 'ShowArrowHead','off', 'MaxHeadSize', 1)
                    set(h2,'color',colour);
                    set(h2,'LineWidth',1);
                    %plot(ax(1),jpdaf.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
                    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    text(ax(1),Logs{j}.xV(1,i)+0.25,Logs{j}.xV(2,i)-0.25,int2str(j));
                    text(ax(1),Logs{j}.xV(1,i)+0.25,Logs{j}.xV(2,i)-.35,num2str(Logs{j}.ExistProbV(i), 2));
                end
            end
                % set the y-axis back to normal.
            %set(ax(1),'ydir','normal');
            str = sprintf('Visualisation of tracking process');
            title(ax(1),str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
%            h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
%            set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
            axis(ax(1),V_bounds)
            %pause(0.01)
            
            % Plot PHD
            cla(ax(2), 'reset');
            [bandwidth,density,X,Y]=kde2d(myphd.Params.particles(1:2,:)');
            %contour3(X,Y,density,50);
            surf(ax(2),X,Y,density);
            shading interp
            colormap(ax(2), jet(3000))
            hold on;
            plot(ax(2), myphd.Params.particles(1,:), myphd.Params.particles(2,:), '.')
            hold on;
            Data = POmodel.obs_inv(0,DataList{i});
            plot(ax(2), Data(1,:), Data(2,:), 'y*');
            axis(ax(2),V_bounds);
            str = sprintf('Visualisation of PHD search track density');
            xlabel(ax(2),'X position (m)')
            ylabel(ax(2),'Y position (m)')
            zlabel(ax(2),'Intensity')
            title(ax(2),str)
            pause(0.01)
            
            if(Record)
                F(i) = getframe(ax(1));
            end
        end
    end
end

if(exist('ospa_vals_phd','var'))
    ospa_vals_phd(:,:,end+1)= zeros(N,3);
else
    ospa_vals_phd(:,:)= zeros(N,3);
end

ospa_c= 100;
ospa_p= 1;
for k=1:N
    trueX = [];
    for i=1:TrueTracks
        if(~isnan(newGroundTruth{k}(1,i))||~isnan(newGroundTruth{k}(2,i)))
            trueX = [trueX, [newGroundTruth{k}(1,i);newGroundTruth{k}(2,i)]];
        end
    end
    estX = [];
    for i=1:numel(Logs)
        if(~isnan(Logs{i}.xV(1:2,k)))
            estX = [estX, Logs{i}.xV(1:2,k)];
        end
    end
    [ospa_vals_phd(k,1,end), ospa_vals_phd(k,2,end), ospa_vals_phd(k,3,end)]= ospa_dist(trueX,estX,ospa_c,ospa_p);
end

figure; ospa= gcf; hold on;
subplot(3,1,1); plot(1:N,ospa_vals_phd(:,1,end),'k'); grid on; set(gca, 'XLim',[1 N]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist');
subplot(3,1,2); plot(1:N,ospa_vals_phd(:,2,end),'k'); grid on; set(gca, 'XLim',[1 N]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc');
subplot(3,1,3); plot(1:N,ospa_vals_phd(:,3,end),'k'); grid on; set(gca, 'XLim',[1 N]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card');
xlabel('Time');

% F = F(2:end);
% vidObj = VideoWriter(sprintf('phd_search.avi'));
% vidObj.Quality = 100;
% vidObj.FrameRate = 10;
% open(vidObj);
% writeVideo(vidObj, F);
% close(vidObj);


