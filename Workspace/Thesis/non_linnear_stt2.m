% Plot settings
SHOW_PLOTS = 1;
SHOW_UPDATE = 1;
SHOW_ARENA = 0;
SHOW_PREDICT = 0;
NUM_SIMS = 1;

%% Model parameter shortcuts
lambdaV = 200; % Expected number of clutter measurements over entire surveillance region
V = 25*15;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 25 0 15]; % [x_min x_max y_min y_max]
P_D = 1;    % Probability of detection
timestep_duration = duration(0,0,1);

%% Extract the GroundTruth data from the example workspace
load('maneuvering_robot.mat');
NumIter = size(TrueTrack.Trajectory,2);

%% Models
% Instantiate a Transitionamic model
transition_model = ConstantVelocityX('NumDims',2,'VelocityErrVariance',0.0002);

% Instantiate an Observation model
% measurement_model = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,...
%                                     'MeasurementErrVariance',0.1^2,...
%                                     'Mapping',[1 3]);
measurement_model = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',...
                                            [(pi/90)^2,0.1^2],'Mapping',[1 3]);

% clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,...
%                                             'Limits',[V_bounds(1:2);...
%                                                       V_bounds(3:4)]); 
                                                  
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,...
                                            'Limits',[0, pi/2;...
                                                      0, 25]); 
detection_model = ConstantDetectionProbabilityX('DetectionProbability',P_D);
% Compile the State-Space model
model = StateSpaceModelX(transition_model,measurement_model,'Clutter',clutter_model, 'Detection', detection_model);

%% Simulation
% Data Simulator
dataSim = SingleTargetMeasurementSimulatorX(model);

% Instantiate a FilterList to store each filter
FilterList = [];
BaseFilterList = [];
% Instantiate all filters
%FilterList{end+1} = KalmanFilterX('Model',ssm);
BaseFilterList{end+1} = ExtendedKalmanFilterX('Model',model);
BaseFilterList{end+1} = UnscentedKalmanFilterX('Model', model);
BaseFilterList{end+1} = ParticleFilterX('Model',model);
% FilterList{end+1} = ExtendedParticleFilterX(ssm);
% FilterList{end+1} = UnscentedParticleFilterX(ssm);

% Store the number of the filters we want to instantiate
numFilters = 3;

% Log Containers
Logs = cell(1, numFilters); % 4 tracks
for i=1:numFilters
    Logs{i}.stateMean = zeros(4,NumIter,NUM_SIMS);       
    Logs{i}.stateCovar = zeros(4,4,NumIter,NUM_SIMS);
    Logs{i}.positionalError = zeros(1,NumIter,NUM_SIMS);
    Logs{i}.execTime = zeros(1,NUM_SIMS);
end

% Create figure windows
if(SHOW_PLOTS)
    if(SHOW_ARENA)
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
    end
    
    figure('units','normalized','outerposition',[0 0 .5 1])
    ax(1) = gca;
end

for simIter = 1:NUM_SIMS
    
    fprintf('Sim %d/%d\n', simIter, NUM_SIMS);
    
    % Simulate some measurements from ground-truth data
    MeasurementScans = dataSim.simulate(TrueTrack);
%     measurement_model = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',...
%                                             [(pi/180)^2,0.1^2],'Mapping',[1 2]);
    measurements = [MeasurementScans.Vectors];
    
    % Reset all filters to prior
    measurement = MeasurementScans(1).getMeasurements(1);
    timestamp = measurement.Timestamp;
    true_mean = TrueTrack.Trajectory(1).Vector;
    xPrior = true_mean;
    PPrior = 10*transition_model.covar();
    dist = GaussianDistributionX(xPrior,PPrior);
    StatePrior = ParticleStateX(dist,10000,timestamp);
    
    for filterInd=1:numFilters
        % Initiate a track using the generated prior
%         track = TrackX(StatePrior, TagX(1));
%         TrackList{1} = track;
%         TrackList{1}.addprop('Filter');
%         TrackList{1}.Filter = copy(base_filter);
%         TrackList{1}.Filter.initialise('Model',model,'StatePrior',StatePrior);
        % Data Associator
        config.ClutterModel = clutter_model;
        config.Clusterer = NaiveClustererX();
        config.Gater = EllipsoidalGaterX(2,'GateLevel',10)';
        config.DetectionModel = detection_model;
        filter = BaseFilterList{filterInd};
        filter.initialise('Model',model,'StatePrior',StatePrior);
        FilterList{filterInd} = ProbabilisticDataAssocX(config);
        FilterList{filterInd}.TrackList{1} = TrackX(copy(filter.StatePrior), TagX(1));
        FilterList{filterInd}.TrackList{1}.addprop('Filter');
        FilterList{filterInd}.TrackList{1}.Filter = copy(filter);
%         FilterList{filterInd}.TrackList{1}.Filter.initialise('Model',model,'StatePrior',StatePrior);
        
    end
    
    % START OF SIMULATION
    % ===================>
    for k = 1:NumIter

        % Update measurements
        MeasurementList = MeasurementScans(k);

        % Iterate and time all filters
        for i=1:numFilters
            tic;
            assocFilter = FilterList{i};
            assocFilter.MeasurementList = MeasurementList;
            assocFilter.predictTracks();
            assocFilter.associate();    
            assocFilter.updateTracks();
            Logs{i}.execTime(1,simIter) = Logs{i}.execTime(1,simIter) + toc;
        end

        % Store Logs
        for i=1:numFilters
            truth = TrueTrack.Trajectory(k).Vector;
            Logs{i}.stateMean(:,k,simIter) = FilterList{i}.TrackList{1}.State.Mean;       
            Logs{i}.stateCovar(:,:,k,simIter) = FilterList{i}.TrackList{1}.State.Covar;
            Logs{i}.positionalError(1,k,simIter) = ((truth(1,:) - FilterList{i}.TrackList{1}.State.Mean(1))^2 + (truth(3,:) - FilterList{i}.TrackList{1}.State.Mean(3))^2);
        end

      % Plot update step results
        if(SHOW_PLOTS && SHOW_UPDATE)
            % Plot data
            cla(ax(1));

            if(SHOW_ARENA)
                 % Flip the image upside down before showing it
                imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));
            end

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
            hold on;
%             meas = measurement_model.finv(measurements(:,1:k));
            meas = measurement_model.finv(MeasurementList.Vectors);
            true_means = [TrueTrack.Trajectory(1:k).Vector];
            h1 = plot(true_means(1,1:k), true_means(3,1:k),'.-k', meas(1,:), meas(3,:), 'rx');
            for i=1:size(h1,1)
                h2 = h1(i);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            end
            for i=1:numFilters
                h2 = plot(Logs{i}.stateMean(1,k,simIter), Logs{i}.stateMean(3,k,simIter), 'o', 'MarkerSize', 10);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
                %plot(pf.Params.particles(1,:), pf.Params.particles(2,:), 'b.', 'MarkerSize', 10);
                plot(Logs{i}.stateMean(1,1:k,simIter), Logs{i}.stateMean(3,1:k,simIter), '.-', 'MarkerSize', 10);
            end
            %legend('KF','UKF', 'PF', 'UPF')
            %legend('KF','EKF', 'UKF', 'PF');%, 'EPF', 'UPF')
            legend('EKF', 'UKF', 'PF');
            
            if(SHOW_ARENA)
                % set the y-axis back to normal.
                set(ax(1),'ydir','normal');
            end

            str = sprintf('Robot positions (Update)');
            title(ax(1),str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
            axis(ax(1),V_BOUNDS)
            pause(0.01);
        end
      %s = f(s) + q*randn(3,1);                % update process 
    end
end

figure
for i=1:numFilters
    hold on;
    plot(sqrt(sum(Logs{i}.positionalError,3)/NUM_SIMS), '.-');
end
%legend('KF','EKF', 'UKF', 'PF');%, 'EPF', 'UPF');
legend('EKF', 'UKF', 'PF');%, 'EPF', 'UPF');

figure
bars = zeros(1, numFilters);
% c = {'KF','EKF', 'UKF', 'PF'};%, 'EPF', 'UPF'};
% c = categorical(c, {'KF','EKF', 'UKF', 'PF'},'Ordinal',true); %, 'EPF', 'UPF'
c = {'EKF', 'UKF', 'PF'};%, 'EPF', 'UPF'};
c = categorical(c, {'EKF', 'UKF', 'PF'},'Ordinal',true); %, 'EPF', 'UPF'
for i=1:numFilters
    bars(i) =  mean(Logs{i}.execTime,2);
end
bar(c, bars);
% END OF SIMULATION
% ===================>

