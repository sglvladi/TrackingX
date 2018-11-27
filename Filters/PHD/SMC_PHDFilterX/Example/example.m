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

lambdaV = 10; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]

% Instantiate a Transitionamic model
transition_model = ConstantVelocityX('NumDims',2,'VelocityErrVariance',0.0001);

% Instantiate a Measurement model
measurement_model = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',0.02,'Mapping',[1 3]);
%measurement_model = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[0.001,0.02],'Mapping',[1 3]);

% Instantiate a clutter model
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,'Limits',[V_bounds(1:2);V_bounds(3:4)]);

% Instantiate birth model
birth_model = DistributionBasedBirthModelX('Distribution', UniformDistributionX([V_bounds(1:2); ...
                                                                                [-0.1 0.1 ];...
                                                                                V_bounds(3:4);...
                                                                                [-0.1 0.1 ]]),...
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
config.StatePrior = ParticleStateX(priorParticles,10*priorWeights);
config.BirthScheme = {'Expansion', 5000};
config.SurvivalProbability = 0.99;
config.DetectionProbability = 0.9; %meas_simulator.DetectionProbability;

% Instantiate PHD filter
myphd = SMC_PHDFilterX(config);

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
        p_i = myphd.IntensityPerHypothesis(2:end);
%         disp(p_i);
        ValidHypothesisInds = find(p_i>0.8)+1;
        for j = 1:length(ValidHypothesisInds)
            hypInd = ValidHypothesisInds(j);
            dist = ParticleDistributionX(myphd.StatePrediction.Particles,myphd.weightsPerHypothesis_(hypInd,:));
            dist.resample(1000);
            plot(ax(1), dist.Particles(1,:), dist.Particles(3,:), '.');
            text(ax(1), dist.Mean(1,:), dist.Mean(3,:), num2str(p_i(hypInd-1)),'FontSize',20,'Color','k');
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
        [bandwidth,density,X,Y]=kde2d(myphd.StatePosterior.Particles([1,3],:)');
        %contour3(X,Y,density,50);
        h = surf(ax(2),X,Y,density);        
        shading interp
        colormap(ax(2), jet(3000))
        %set(h, 'edgecolor','none')
        hold on;
%         plot(ax(2), myphd.Particles(1,:), myphd.Particles(3,:), '.')
%         hold on;
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