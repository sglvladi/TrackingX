% phd_test.m
% ====================================================>
% This is a test script which demonstrates the usage of the "SMC_PHD" class.

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
obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',0.04,'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

lambdaV = 50; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]

% Generate observations (Poisson number with rate of lambdaV, positions are uniform over surveillance region)
numTrueTracks = 3;
[DataList,x1,y1] = gen_obs_cluttered_multi2(numTrueTracks, x_true, y_true, 0.2, V_bounds, lambdaV, 1); 
N=size(DataList,2); % timesteps 

% Assign PHD parameter values
config.NumParticles = 20000;              % number of particles
config.Model = ssm;
q = dyn.covariance();
config.BirthIntFcn = @(Np) [(V_bounds(2)-V_bounds(1))*rand(Np,1), mvnrnd(zeros(Np,1), 0.05'),(V_bounds(4)-V_bounds(3))*rand(Np,1),mvnrnd(zeros(Np,1), 0.05')]'; % Uniform position and heading, Gaussian speed
config.PriorDistFcn = @ (Np) deal(config.BirthIntFcn(Np), repmat(1/Np, Np, 1)');
config.BirthScheme = {'Mixture', 0.01};
config.ProbOfDeath = 0.005;
config.ProbOfDetection = 0.9;
config.ClutterRate = lambdaV/V;
config.ResamplingScheme = 'Multinomial';

% Instantiate PHD filter
myphd = SMC_PHDFilterX(config);

% Instantiate PF filter
mypf = ParticleFilterX(ssm);

% Instantiate JPDAF (main tracker)
Params_pdaf.Clusterer = NaiveClustererX();
Params_pdaf.Gater = EllipsoidalGaterX(2,'ProbOfGating',0.998)';
Params_pdaf.ProbOfDetect = 0.9;
mypdaf = JointProbabilisticDataAssocX(Params_pdaf);
mypdaf.MeasurementList = DataList{1}(:,:); 

% Instantiate a Track Initiator
config_ti.Filter = mypf;
config_ti.PHDFilter = myphd;
config_ti.ProbOfGating = 0.998;
config_ti.ConfirmThreshold = 0.8;
myti = PhdTrackInitiatorX(config_ti);

% Instantiate Existance probability calculator
myepc = ExistProbCalculatorX(config.ProbOfDetection,config.ProbOfDeath,0.998);

% Instantiate a Track Deleter
mytd = ExistProbTrackDeleterX(0.1);

TrackList = {};

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
    
    % Compute new existence probabilities
    mypdaf.TrackList = myepc.calculate(mypdaf.TrackList,mypdaf.AssocWeightsMatrix);
    
    if(exist('data_plot','var'))
        delete(data_plot);
    end
    data_plot = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
    
    % Perform Track initiation
    myti.MeasurementList = tempDataList; % New observations
    % Compute rhi as given by Eq. (16) in [2] 
    numMeas = size(myti.MeasurementList,2);
    if(mypdaf.AssocWeightsMatrix>-1)  % Check if beta exists
        rhi = prod(1-mypdaf.AssocWeightsMatrix(:,2:end),1);
    else
        rhi = ones(1, numMeas);
    end
    myti.UnassocWeights = rhi;
    TrackList = [mypdaf.TrackList, myti.initiateTracks()];
    
    % Perform Track deletion
    [TrackList, DeletedTracks] = mytd.deleteTracks(TrackList);
    
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
        
    fprintf("Estimated number of targets: %f\n", myphd.NumTargets);
    %fprintf("Ellapsed time: %f\n", ellapsed);
    for j=1:numel(TrackList)
        fprintf("Track %d: %f\n", TrackList{j}.TrackID, TrackList{j}.ProbOfExist);
    end
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
            h2 = plot_gaussian_ellipsoid(TrackList{j}.Filter.StateMean([1 3]), TrackList{j}.Filter.StateCovar([1 3],[1 3]),1,20,ax(1));
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            text(ax(1),TrackList{j}.Filter.StateMean(1)+0.25,TrackList{j}.Filter.StateMean(3)-0.25,int2str(TrackList{j}.TrackID));
            text(ax(1),TrackList{j}.Filter.StateMean(1)+0.25,TrackList{j}.Filter.StateMean(3)-.35,num2str(TrackList{j}.ProbOfExist,2));
        end
        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Robot positions (Update)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis(ax(1),V_bounds)
            
        % Plot PHD
        cla(ax(2), 'reset');
        [bandwidth,density,X,Y]=kde2d(myphd.Particles([1,3],:)');
        %contour3(X,Y,density,50);
        h = surf(ax(2),X,Y,density);        
        shading interp
        colormap(ax(2), jet(3000))
        %set(h, 'edgecolor','none')
        hold on;
        plot(ax(2), myphd.Particles(1,:), myphd.Particles(3,:), '.')
        hold on;
        plot(ax(2), myphd.MeasurementList(1,:), myphd.MeasurementList(2,:), 'y*');
        axis(ax(2), [V_bounds]);
        str = sprintf('PHD intensity (Update)');
        xlabel(ax(2),'X position (m)')
        ylabel(ax(2),'Y position (m)')
        zlabel(ax(2),'Intensity')
        title(ax(2),str)
        pause(0.01)
    end
end