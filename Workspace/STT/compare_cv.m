% Plot settings
SHOW_PLOTS = 1;
SHOW_UPDATE = 1;
SHOW_ARENA = 0;
SHOW_PREDICT = 0;
NUM_SIMS = 1;
V_BOUNDS = [0 25 0 15];

%% Extract the GroundTruth data from the example workspace
load('maneuvering_robot.mat');
NumIter = size(TrueTrack.Trajectory,2);

%% Models
% Instantiate a Transitionamic model
transition_model = ConstantVelocityX('NumDims',2,'VelocityErrVariance',0.0001);
% transition_model = ConstantHeadingX('VelocityErrVariance',(0.02)^2, 'HeadingErrVariance',0.2^2);

% Instantiate an Observation model
% measurement_model = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',1^2,'Mapping',[1 2]);
%measurement_model = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[0.001,0.02],'Mapping',[1 3]);
measurement_model = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',...
                                            [(pi/180)^2,0.1^2],'Mapping',[1 3]);
% Compile the State-Space model
model = StateSpaceModelX(transition_model,measurement_model);

%% Simulation
% Data Simulator
dataSim = SingleTargetMeasurementSimulatorX(model);

% Instantiate a FilterList to store each filter
FilterList = [];

% Instantiate all filters
%FilterList{end+1} = KalmanFilterX('Model',ssm);
FilterList{end+1} = ExtendedKalmanFilterX('Model',model);
FilterList{end+1} = UnscentedKalmanFilterX('Model', model);
FilterList{end+1} = ParticleFilterX('Model',model);
% FilterList{end+1} = ExtendedParticleFilterX(ssm);
% FilterList{end+1} = UnscentedParticleFilterX(ssm);

% Store the number of the filters we want to instantiate
numFilters = numel(FilterList);

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
%     xPrior = [true_mean(1);true_mean(3);0.2; pi/4];
    xPrior = [true_mean(1);0.2;true_mean(3); 0.2];
    PPrior = 30*transition_model.covar();
    dist = GaussianDistributionX(xPrior,PPrior);
    StatePrior = ParticleStateX(dist,5000,timestamp);
    for filterInd=1:numFilters
        % Setup prior
        FilterList{filterInd}.initialise('StatePrior',StatePrior);
    end
    
    % START OF SIMULATION
    % ===================>
    for k = 1:NumIter

        % Update measurements
        MeasurementList = MeasurementScans(k);

        % Iterate and time all filters
        for i=1:numFilters
            tic;
            filter = FilterList{i};
            filter.MeasurementList = MeasurementList;
%             prior = filter.StatePosterior;
            filter.predict();
            filter.update();
%             prediction = filter.predict(prior, MeasurementList.Timestamp);
%             posterior = filter.update(prediction, MeasurementList);
            Logs{i}.execTime(1,simIter) = Logs{i}.execTime(1,simIter) + toc;
        end

        % Store Logs
        for i=1:numFilters
            truth = TrueTrack.Trajectory(k).Vector;
            Logs{i}.stateMean(:,k,simIter) = FilterList{i}.StatePosterior.Mean;       
            Logs{i}.stateCovar(:,:,k,simIter) = FilterList{i}.StatePosterior.Covar;
            Logs{i}.positionalError(1,k,simIter) = ((truth(1,:) - FilterList{i}.StatePosterior.Mean(1))^2 + (truth(3,:) - FilterList{i}.StatePosterior.Mean(3))^2);
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
            meas = measurement_model.finv(measurements(:,1:k));
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

