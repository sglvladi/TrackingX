% Script used to compare the estimation performance between EKF, UKF and PF
% (and EPF, UPF) using non-linear (both dynamic and observation) models
% =========================================================================

% Plot settings
SHOW_PLOTS = 0;
SHOW_UPDATE = 1;
SHOW_ARENA = 0;
SKIP_STEPS = 1;
NUM_SIMS = 20;
V_BOUNDS = [0 25 0 15];

% Load dataset
load('maneuvering_robot.mat');

% Instantiate a Dynamic model
dyn = ConstantHeadingModelX('VelocityErrVariance',(0.03)^2, 'HeadingErrVariance',0.07);

% Instantiate an Observation model
% obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',0.2,'Mapping',[1 2]);
obs = Polar2CartGaussModelX('NumStateDims',4,'RangeErrVariance',(0.1)^2,...
                          'ThetaErrVariance',(pi/45)^2,'Mapping',[1 2]);
                      
% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% Instantiate a FilterList to store each filter
FilterList = [];

% Instantiate all filters
FilterList{end+1} = ExtendedKalmanFilterX(ssm);
FilterList{end+1} = UnscentedKalmanFilterX(ssm);
FilterList{end+1} = ParticleFilterX(ssm,'NumParticles',10000);
FilterList{end+1} = ExtendedParticleFilterX(ssm,'NumParticles',10000);
FilterList{end+1} = UnscentedParticleFilterX(ssm,'NumParticles',10000);

% Store the number of the filters we want to instantiate
numFilters = numel(FilterList);

[measurementsCellArray,~,inds] = dataGen(obs,GroundTruth,1,0,[],[],SKIP_STEPS+1);

% Log Containers
Logs = cell(1, numFilters); % 4 tracks
N = numel(measurementsCellArray);
for i=1:numFilters
    Logs{i}.stateMean = zeros(4,N,NUM_SIMS);       
    Logs{i}.stateCovar = zeros(4,4,N,NUM_SIMS);
    Logs{i}.positionalError = zeros(1,N,NUM_SIMS);
    Logs{i}.RMSE = zeros(1,N,NUM_SIMS);
    Logs{i}.execTime = zeros(1,NUM_SIMS);
end
measurementError = zeros(1,N,NUM_SIMS);
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
    
    disp(strcat({'Running simulation '},num2str(simIter),{' of '},num2str(NUM_SIMS)));  
    
    [measurementsCellArray,~,inds] = dataGen(obs,GroundTruth,1,0,[],[],SKIP_STEPS+1);
    N = numel(measurementsCellArray);
    % Generate ground truth and measurements
    for k = 1:N
        % Generate new measurement from ground truth
        truth(:,k) = [x_true(inds(k),1); y_true(inds(k),1)];   
        measurement(:,k) = measurementsCellArray{k}; 
    end
    
    speed_xy = diff(truth')';
    speed_abs = sqrt(speed_xy(1,:).^2+speed_xy(2,1).^2);
    speed_init = speed_abs(1);
    headings = mod(atan2(speed_xy(1,:), speed_xy(2,:))*180/pi, 360);
    head_init = headings(1)*pi/180;
    priorStateMean = [truth(1,1);truth(2,1); speed_init; head_init];%obs.heval_inv([measurement(1,1); measurement(2,1)]);
    measErrCovar = ssm.Obs.covariance();
    stateErrCovar = ssm.Dyn.covariance();
    priorStateCovar = stateErrCovar + blkdiag(measErrCovar(1,1),measErrCovar(2,2),0,0);
    
    % Reset all filters to prior
    for filterInd=1:numFilters
        if(isa(FilterList{filterInd},'KalmanFilterX'))
            FilterList{filterInd}.initialise('PriorStateMean',priorStateMean,'PriorStateCov',priorStateCovar);
        elseif(isa(FilterList{i},'ParticleFilterX'))
            FilterList{filterInd}.initialise('PriorDistFcn',@(N)deal(mvnrnd(priorStateMean(:,ones(1,N))',priorStateCovar)',1/N*ones(1,N)));
        end
        Logs{filterInd}.stateMean(:,1,simIter) = FilterList{filterInd}.StateMean;       
        Logs{filterInd}.stateCovar(:,:,1,simIter) = FilterList{filterInd}.StateCovar;
        Logs{filterInd}.positionalError(1,1,simIter) = ((truth(1,1) - FilterList{filterInd}.StateMean(1))^2 + (truth(2,1) - FilterList{filterInd}.StateMean(2))^2);
    end
    meas2state = obs.heval_inv(measurement(:,1));
    measurementError(1,1,simIter) = ((truth(1,1) - meas2state(1))^2 + (truth(2,1) - meas2state(2))^2);
    
    % START OF SIMULATION
    % ===================>
    for k = 2:N

        % Update measurements
        for i=1:numFilters
            FilterList{i}.Measurement = measurement(:,k);
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
            Logs{i}.stateMean(:,k,simIter) = FilterList{i}.StateMean;       
            Logs{i}.stateCovar(:,:,k,simIter) = FilterList{i}.StateCovar;
            Logs{i}.positionalError(1,k,simIter) = ((truth(1,k) - FilterList{i}.StateMean(1))^2 + (truth(2,k) - FilterList{i}.StateMean(2))^2);
        end
        meas2state = obs.heval_inv(measurement(:,k));
        measurementError(1,k,simIter) = ((truth(1,k) - meas2state(1))^2 + (truth(2,k) - meas2state(2))^2);

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
            %h2 = plot(ax(1),measurement(1,k),measurement(2,k),'k*','MarkerSize', 10);
            meas = obs.heval_inv([measurement(1,k);measurement(2,k)]);
            h2 = plot(ax(1),meas(obs.Mapping(1)),meas(obs.Mapping(2)),'k*','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            h2 = plot(ax(1),truth(1,1:k),truth(2,1:k),'b.-','LineWidth',1);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            h2 = plot(ax(1),truth(1,k),truth(2,k),'bo','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

            for i=1:numFilters
                h2 = plot(Logs{i}.stateMean(1,k,simIter), Logs{i}.stateMean(2,k,simIter), 'o', 'MarkerSize', 10);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
                %plot(pf.Params.particles(1,:), pf.Params.particles(2,:), 'b.', 'MarkerSize', 10);
                plot(Logs{i}.stateMean(1,1:k,simIter), Logs{i}.stateMean(2,1:k,simIter), '.-', 'MarkerSize', 10);
            end
            %legend('KF','UKF', 'PF', 'UPF')
            legend('EKF', 'UKF', 'PF', 'EPF', 'UPF')

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
    
    % Store RMSE
    for i=1:numFilters
        Logs{i}.RMSE(1,:,simIter) = sqrt(mean((truth - FilterList{i}.StateMean([1,3],:)).^2));
    end
    measurementRMSE(1,:,simIter) = sqrt(mean((truth - measurement).^2));
end

figure
%disp(mean(sum(measurementRMSE,3)/NUM_SIMS)*100);
disp(mean(sqrt(mean((measurementError*100).^2,3))));
for i=1:numFilters
    hold on;
    plot(sqrt(mean(Logs{i}.positionalError,3)), '.-');
    disp(mean(sqrt(mean((Logs{i}.positionalError*100).^2,3))));
    %disp(mean(sum(Logs{i}.RMSE,3)/NUM_SIMS)*100);
end
legend('EKF','UKF', 'PF', 'EPF', 'UPF');

figure
bars = zeros(1, numFilters);
c = {'EKF', 'UKF', 'PF', 'EPF', 'UPF'};
c = categorical(c, {'EKF', 'UKF', 'PF', 'EPF', 'UPF'},'Ordinal',true); 
for i=1:numFilters
    bars(i) =  mean(Logs{i}.execTime,2);
end
bar(c, bars);
% END OF SIMULATION
% ===================>

