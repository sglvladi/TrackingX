% phd_test.m
% ====================================================>
% This is a test script which demonstrates the usage of the "SMC_PHD" class.
% 
% SETUP:
%  * Before running the simulation, open "2_crossing_targets.mat" or "3_roaming_targets.mat" datasets, from the "datasets" folder
%  * The datasets have been extracted by simulating the motion of differential robots in a 2D-plane (x,y)
%  * The "gen_obs_cluttered_multi3" function takes as an input the ground truth data, including information about the measurement noise and clutter rate
%     and then produces 1xNk cell array of corrupted and cluttered measurements, Nk being the total number of timesteps

load('3_robots.mat');
tot_elapsed = 0;
% Plot settings
ShowPlots = 1;              % Set to 0 to hide plots
ShowPrediction = 0;         % Set to 0 to skip showing prediction
ShowUpdate = 1;             % Set to 0 to skip showing update
TrackNum = size(x_true,2);

tot_ellapsed = 0;

% Instantiate a Dynamic model
dyn = ConstantVelocityModelX_2D('VelocityErrVariance',0.0001);

% Instantiate an Observation model
obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',0.02,'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% n_x = 4;      % state dimensions
% q = 0.01;     % std of process noise 
% n_y = 2;      % measurement dimensions
% r = 0.1;      % std of measurement noise
lambdaV = 50; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]

% Generate observations (Poisson number with rate of lambdaV, positions are uniform over surveillance region)
numTrueTracks = 3;
[DataList,x1,y1] = gen_obs_cluttered_multi2(numTrueTracks, x_true, y_true, sqrt(obs.ObsErrVariance), V_bounds, lambdaV, 1); 
N=size(DataList,2); % timesteps 

% Instantiate PF filter
kf = KalmanFilterX(ssm);
obs_covar= obs.covariance();
kf.StateCovar = dyn.covariance() + blkdiag(obs_covar(1,1), 0, obs_covar(2,2),0);
mypf = ParticleFilterX(ssm);

% Initiate PDAF parameters
Params_jpdaf.Clusterer = NaiveClustererX();
Params_jpdaf.Gater = EllipsoidalGaterX(2,'ProbOfGating',0.998)';
Params_jpdaf.ProbOfDetect = 0.9;
mypdaf = JointProbabilisticDataAssocX(Params_jpdaf);
mypdaf.MeasurementList = DataList{1}(:,:); 

% Initiate Track Initiator
% Initiate PDAF parameters
Params_pdaf.Clusterer = NaiveClustererX();
Params_pdaf.Gater = EllipsoidalGaterX(2,'ProbOfGating',0.998)';
Params_pdaf.ProbOfDetect = 0.9;
pdaf = ProbabilisticDataAssocX(Params_pdaf);
pdaf.MeasurementList = DataList{1}(:,:); 
config_ti.InitFilter = kf;
config_ti.DataAssociator = pdaf;
config_ti.ConfirmThreshold = [8,10];
config_ti.DeleteThreshold = [5,10];
CovarThreshold = 4*kf.StateCovar;
config_ti.CustomDeleteConditionFcn = @(x,t) t.Filter.StateCovar(1,1)>CovarThreshold(1,1) || t.Filter.StateCovar(3,3)>CovarThreshold(3,3);
myti = MofN_TrackInitiatorX(config_ti);

TrackList = [];

% Create figure windows
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
end

% START OF SIMULATION
% ===================>
for k=1:N
    fprintf('Iteration = %d/%d\n================>\n',k,N);
    if(k==92)
    
    end
    % Extract DataList at time k
    tempDataList = DataList{k}(:,:);
    tempDataList( :, ~any(tempDataList,1) ) = [];
    
    % Process JPDAF
    mypdaf.MeasurementList = tempDataList;
    mypdaf.TrackList = TrackList;
    for j=1:numel(TrackList)
        mypdaf.TrackList{j}.Filter.predict();
    end
    mypdaf.associate();    
    mypdaf.updateTracks();
    
    if(exist('data_plot','var'))
        delete(data_plot);
    end
    data_plot = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
    
    % Perform Track initiation
%     myti.MeasurementList = tempDataList; % New observations
%     myti.TrackList = TrackList;
%     myti.AssocWeightsMatrix = mypdaf.AssocWeightsMatrix;
    [TrackList, TentativeTrackList] = myti.initiateTracks(mypdaf.TrackList, tempDataList, mypdaf.AssocWeightsMatrix);
    tic;
    
    % Plot prediction step results
    if(ShowPlots && ShowPrediction)
        % Plot data
        cla(ax(1));
         % Flip the image upside down before showing it
        imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

        % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
        hold on;
        h2 = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
        for j=1:TrackNum
            h2 = plot(ax(1), x_true(1:k,j),y_true(1:k,j),'b.-','LineWidth',1);
            if j==2
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            h2 = plot(ax(1), x_true(k,j),y_true(k,j),'bo','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end

        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Robot positions (Prediction)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis(ax(1),V_bounds)
            
        % Plot PHD
        cla(ax(2), 'reset');
        [bandwidth,density,X,Y]=kde2d(myphd.PredParticles([1,3],:)');
        %contour3(X,Y,density,50);
        h = surf(ax(2),X,Y,density);
        shading interp
        colormap(ax(2), jet(3000))
        set(h, 'edgecolor','none')
        hold on;
        plot(ax(2), myphd.PredParticles(1,:), myphd.PredParticles(3,:), '.')
        hold on;
        plot(ax(2), myphd.MeasurementList(1,:), myphd.MeasurementList(2,:), 'y*');
        axis(ax(2), [V_bounds 0 10]);
        str = sprintf('PHD intensity (Prediction)');
        xlabel(ax(2),'X position (m)')
        ylabel(ax(2),'Y position (m)')
        zlabel(ax(2),'Intensity')
        title(ax(2),str)
        pause(0.01)
    end
        
    % Update PHD filter
    %myphd.update();
    ellapsed = toc;
    tot_ellapsed = tot_ellapsed + ellapsed;
    fprintf("Ellapsed time: %f\n", ellapsed);
%     for j=1:numel(TrackList)
%         fprintf("Track %d: %f\n", TrackList{j}.TrackID, TrackList{j}.ProbOfExist);
%     end
    disp("");
    % Plot update step results
    if(ShowPlots && ShowUpdate)
        % Plot data
        cla(ax(1));
         % Flip the image upside down before showing it
        imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

        % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
        hold on;
        data_plot = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
        for j=1:numTrueTracks
            h2 = plot(ax(1), x_true(1:k,j),y_true(1:k,j),'b.-','LineWidth',1);
            if j==2
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            h2 = plot(ax(1), x_true(k,j),y_true(k,j),'bo','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
        for j=1:numel(TrackList)
            h2 = plot(ax(1), TrackList{j}.Filter.StateMean(1),TrackList{j}.Filter.StateMean(3),'.','LineWidth',1);
            if j==2
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            h2 = plot_gaussian_ellipsoid(TrackList{j}.Filter.StateMean([1 3]), TrackList{j}.Filter.StateCovar([1 3],[1 3]),'r',1,20,ax(1));
%             set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
%             text(ax(1),TrackList{j}.Filter.StateMean(1)+0.25,TrackList{j}.Filter.StateMean(3)-0.25,int2str(TrackList{j}.TrackID));
%             text(ax(1),TrackList{j}.Filter.StateMean(1)+0.25,TrackList{j}.Filter.StateMean(3)-.35,num2str(TrackList{j}.ProbOfExist,2));
        end
        for j=1:numel(TentativeTrackList)
            h2 = plot(ax(1), TentativeTrackList{j}.Filter.StateMean(1),TentativeTrackList{j}.Filter.StateMean(3),'.','LineWidth',1);
            if j==2
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            h2 = plot_gaussian_ellipsoid(TentativeTrackList{j}.Filter.StateMean([1 3]), TentativeTrackList{j}.Filter.StateCovar([1 3],[1 3]),'g',1,20,ax(1));
%             set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
%             text(ax(1),TrackList{j}.Filter.StateMean(1)+0.25,TrackList{j}.Filter.StateMean(3)-0.25,int2str(TrackList{j}.TrackID));
%             text(ax(1),TrackList{j}.Filter.StateMean(1)+0.25,TrackList{j}.Filter.StateMean(3)-.35,num2str(TrackList{j}.ProbOfExist,2));
        end
        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Robot positions (Update)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis(ax(1),V_bounds)
            
        % Plot PHD
%         cla(ax(2), 'reset');
%         [bandwidth,density,X,Y]=kde2d(myphd.Particles([1,3],:)');
%         %contour3(X,Y,density,50);
%         h = surf(ax(2),X,Y,density);        
%         shading interp
%         colormap(ax(2), jet(3000))
%         %set(h, 'edgecolor','none')
%         hold on;
%         plot(ax(2), myphd.Particles(1,:), myphd.Particles(3,:), '.')
%         hold on;
%         plot(ax(2), myphd.MeasurementList(1,:), myphd.MeasurementList(2,:), 'y*');
%         axis(ax(2), [V_bounds]);
%         str = sprintf('PHD intensity (Update)');
%         xlabel(ax(2),'X position (m)')
%         ylabel(ax(2),'Y position (m)')
%         zlabel(ax(2),'Intensity')
%         title(ax(2),str)
         pause(0.01)
    end
end