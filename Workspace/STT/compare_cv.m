% Plot settings
SHOW_PLOTS = 0;
SHOW_UPDATE = 1;
SHOW_ARENA = 0;
SHOW_PREDICT = 0;
NUM_SIMS = 10;
V_BOUNDS = [0 25 0 15];

% Extract the ground truth data from the example workspace
% load('example.mat');
% x_true = truth(1,:);
% y_true = truth(2,:);
% Instantiate a Transitionamic model
dyn = ConstantVelocityX('NumDims',2,'VelocityErrVariance',0.001);

% Instantiate an Observation model
%obs = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',0.02,'Mapping',[1 3]);
obs = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[0.001,0.04],'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% Instantiate a FilterList to store each filter
FilterList = [];

% Instantiate all filters
%FilterList{end+1} = KalmanFilterX('Model',ssm);
FilterList{end+1} = ExtendedKalmanFilterX('Model',ssm);
FilterList{end+1} = UnscentedKalmanFilterX('Model',ssm);
FilterList{end+1} = ParticleFilterX('Model',ssm);
% FilterList{end+1} = ExtendedParticleFilterX(ssm);
% FilterList{end+1} = UnscentedParticleFilterX(ssm);

% Store the number of the filters we want to instantiate
numFilters = numel(FilterList);

% Log Containers
Logs = cell(1, numFilters); % 4 tracks
N = size(x_true,1)-2;
for i=1:numFilters
    Logs{i}.stateMean = zeros(4,N,NUM_SIMS);       
    Logs{i}.stateCovar = zeros(4,4,N,NUM_SIMS);
    Logs{i}.positionalError = zeros(1,N,NUM_SIMS);
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
    
    fprintf('Sim %d/%d', simIter, NUM_SIMS);
    
    measurementsCellArray = dataGen(obs,GroundTruth,1);
    % Generate ground truth and measurements
    for k = 1:N
        % Generate new measurement from ground truth
        truth(:,k) = [x_true(k+2,1); y_true(k+2,1)];   
        measurement(:,k) = measurementsCellArray{k}; 
    end
    
    priorStateMean = obs.finv(measurement(:,1));
    measErrCovar = ssm.Measurement.covar();
    stateErrCovar = ssm.Transition.covar();
    priorStateCovar = 10*stateErrCovar;%+ blkdiag(measErrCovar(1,1),0,measErrCovar(2,2),0);
    
    
    % Reset all filters to prior
    for filterInd=1:numFilters
        if(isa(FilterList{filterInd},'KalmanFilterX'))
            Prior.State = GaussianStateX(priorStateMean,priorStateCovar);
        elseif(isa(FilterList{i},'ParticleFilterX'))
            Prior.State = ParticleStateX(priorStateMean,priorStateCovar,5000);
        end
        FilterList{filterInd}.initialise('Prior',Prior);
    end
    
    % START OF SIMULATION
    % ===================>
    for k = 1:N

        % Update measurements
        for i=1:numFilters
            FilterList{i}.MeasurementList = measurement(:,k);
        end

        % Iterate and time all filters
        for i=1:numFilters
            tic;
            FilterList{i}.predict();
            FilterList{i}.update();
            Logs{i}.execTime(1,simIter) = Logs{i}.execTime(1,simIter) + toc;
        end

        % Store Logs
        for i=1:numFilters
            Logs{i}.stateMean(:,k,simIter) = FilterList{i}.Posterior.State.Mean;       
            Logs{i}.stateCovar(:,:,k,simIter) = FilterList{i}.Posterior.State.Covar;
            Logs{i}.positionalError(1,k,simIter) = ((truth(1,k) - FilterList{i}.Posterior.State.Mean(1))^2 + (truth(2,k) - FilterList{i}.Posterior.State.Mean(3))^2);
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
            meas = obs.finv(measurement(:,k));
            h2 = plot(ax(1),meas(1),meas(3),'k*','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            h2 = plot(ax(1),truth(1,1:k),truth(2,1:k),'b.-','LineWidth',1);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            h2 = plot(ax(1),truth(1,k),truth(2,k),'bo','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

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

