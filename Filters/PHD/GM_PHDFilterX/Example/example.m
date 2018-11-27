% phd_test.m
% ====================================================>
% This is a test script which demonstrates the usage of the "SMC_PHD" class.
% 
% SETUP:
%  * Before running the simulation, open "2_crossing_targets.mat" or "3_roaming_targets.mat" datasets, from the "datasets" folder
%  * The datasets have been extracted by simulating the motion of differential robots in a 2D-plane (x,y)
%  * The "gen_obs_cluttered_multi3" function takes as an input the ground truth data, including information about the measurement noise and clutter rate
%     and then produces 1xNk cell array of corrupted and cluttered measurements, Nk being the total number of timesteps

% Plot settings
ShowPlots = 1;              % Set to 0 to hide plots
ShowPrediction = 0;         % Set to 0 to skip showing prediction
ShowUpdate = 1;             % Set to 0 to skip showing update

lambdaV = 50; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]

% Instantiate a Transitionamic model
transition_model = ConstantVelocityX('NumDims',2,'VelocityErrVariance',0.0001);

% Instantiate a Measurement model
measurement_model = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',0.01,'Mapping',[1 3]);
%measurement_model = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[0.001,0.02],'Mapping',[1 3]);

% Instantiate a clutter model
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,'Limits',[V_bounds(1:2);V_bounds(3:4)]);

% Instantiate birth model
numBirthComponents = 10;
BirthComponents.Means = [(V_bounds(2)-V_bounds(1))*rand(numBirthComponents,1)';
                         mvnrnd(zeros(numBirthComponents,1), 0.5)';
                         (V_bounds(4)-V_bounds(3))*rand(numBirthComponents,1)';
                         mvnrnd(zeros(numBirthComponents,1), 0.5)'];
BirthComponents.Covars = repmat(diag([10 1 10 1]),1,1,numBirthComponents);
BirthComponents.Weights = 0.01*ones(1,numBirthComponents);
birth_distribution = GaussianMixtureX(BirthComponents.Means,BirthComponents.Covars, BirthComponents.Weights);
birth_model = DistributionBasedBirthModelX('Distribution', birth_distribution,...
                                           'BirthIntensity', 0.000001);
% Compile the State-Space model
ssm = StateSpaceModelX(transition_model,measurement_model,'Clutter',clutter_model, 'Birth', birth_model);

% Extract the ground truth data from the example workspace
load('example.mat');
NumIter = size(GroundTruth,2);

% Set BirthIntensity
NumTracks = 3;

% Generate DataList
meas_simulator = MeasurementSimulatorX('Model',ssm);
meas_simulator.DetectionProbability = 1;
[DataList, nGroundTruth] = meas_simulator.simulate(GroundTruth);

% Assign PHD parameter values
config.Model = ssm;
[priorParticles, priorWeights] = ssm.Birth.random(50000);
config.Filter = KalmanFilterX('Model',ssm);
config.StatePrior = GaussianMixtureStateX(birth_distribution);
config.StatePrior.Weights = config.StatePrior.Weights*100;
config.SurvivalProbability = 0.99;
config.DetectionProbability = 0.9; %meas_simulator.DetectionProbability;

% Instantiate PHD filter
myphd = GM_PHDFilterX(config);

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
for k=1:NumIter
    fprintf('Iteration = %d/%d\n================>\n',k,NumIter);
    
    % Extract DataList at time k
    tempDataList = DataList{k}(:,:);     
    
    % Change PHD filter parameters
    myphd.MeasurementList = tempDataList; % New observations
    
    % Predict PHD filter
    myphd.predict();
        
    % Update PHD filter
    myphd.update();
    
    fprintf("Estimated number of targets: %f\n", sum(myphd.StatePosterior.Weights));
    % Plot update step results
    if(ShowPlots && ShowUpdate)
        % Plot data
        cla(ax(1));
         % Flip the image upside down before showing it
        imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

        % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
        hold on;
        h2 = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
        for j=1:NumTracks
            h2 = plot(ax(1), x_true(1:k,j),y_true(1:k,j),'b.-','LineWidth',1);
            h2 = plot(ax(1), x_true(k,j),y_true(k,j),'bo','MarkerSize', 10);
        end
        
        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Robot positions (Update)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis(ax(1),V_bounds)
            
        % Plot PHD
        particles = zeros(2,0);
        numParticles = 10000;
        for i=1:myphd.StatePosterior.NumComponents
            numParts = round(myphd.StatePosterior.Weights(i)/sum(myphd.StatePosterior.Weights)*numParticles);
            particles = [particles, mvnrnd(myphd.StatePosterior.Means([1,3],i)',myphd.StatePosterior.Covars([1,3],[1,3],i)',numParts)'];
        end
        %gm = gmdistribution(myphd.Components.Means([1,3],:)',myphd.Components.Covars([1,3],[1,3],:),myphd.Components.Weights/sum(myphd.Components.Weights));
        cla(ax(2), 'reset');
        %fsurf(ax(2),@(x,y)pdf(gm,[x y]),[0 10 0 10])
        [bandwidth,density,X,Y]=kde2d(particles');
        %contour3(X,Y,density,50);
        h = surf(ax(2),X,Y,density);        
        %camlight
        shading interp
        colormap(ax(2), jet(3000))
        %set(h, 'edgecolor','none')
        hold on;
        plot(ax(2), particles(1,:), particles(2,:), '.')
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