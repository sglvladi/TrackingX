%% Plot settings
ShowPlots = 1;
ShowArena = 1;
SkipFrames = 0;

% Instantiate a Dynamic model
dyn = ConstantVelocityModelX_2D('VelocityErrVariance',0.0001);

% Instantiate an Observation model
obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',0.04,'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% Instantiate a Kalman Filter object
kf = KalmanFilterX(ssm);

% Extract the ground truth data from the example workspace
%load('example.mat');
NumIter = size(GroundTruth,2);

% Set NumTracks
NumTracks = 2;

% Generate DataList
%[DataList,x1,y1] = gen_obs_cluttered_multi2(NumTracks, x_true, y_true, sqrt(obs.ObsErrVariance), [0 10 0 10], 50, 1);

% Initiate TrackList
for i=1:NumTracks
    Params_kf.priorStateMean = [GroundTruth{1}(1,i); 0; GroundTruth{1}(2,i); 0];
    Params_kf.priorStateCovar = dyn.covariance(); %blkdiag(POmodel.Params.R(1)/2, 2^2, 2*pi);%CVmodel.Params.Q(1);
    Params_kf.Model = ssm;
    TrackList{i} = TrackX();
    TrackList{i}.addprop('Filter');
    TrackList{i}.Filter = KalmanFilterX(Params_kf);
    
    %TrackList{i}.Filter = ExtendedParticleFilterX(ssm);
    %TrackList{i}.Filter.initialise('NumParticles', 5000, 'PriorDistFcn',@(N)deal(mvnrnd(Params_kf.priorStateMean(:,ones(1,N))',Params_kf.priorStateCovar)',1/N*ones(1,N)));
    %TrackList{i}.TrackObj = ParticleFilterX(Params_pf);%UKalmanFilterX(Params_kf, CHmodel, POmodel);% 
end

%% Initiate PDAF parameters
Params_pdaf.Clusterer = NaiveClustererX();
Params_pdaf.Gater = EllipsoidalGaterX(2,'GateLevel',10)';
Params_pdaf.ProbOfDetect = 0.9;

%% Instantiate PDAF
pdaf = JointProbabilisticDataAssocX(Params_pdaf);
pdaf.TrackList = TrackList;
pdaf.MeasurementList = DataList{1}(:,:); 

%% Instantiate Log to store output
N=size(DataList,2);
Logs = cell(1, NumTracks); % 4 tracks
N = size(x_true,1)-2;
for i=1:NumTracks
    Logs{i}.Estimates.StateMean = zeros(4,N);
    Logs{i}.Estimates.StateCovar = zeros(4,4,N);
    Logs{i}.Groundtruth.State = zeros(2,N);
end

% Create figure windows
if(ShowPlots)
    if(ShowArena)
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

exec_time = 0;
for i = 1:N
    % Remove null measurements   
    pdaf.MeasurementList = DataList{i}(:,:);
    pdaf.MeasurementList( :, ~any(pdaf.MeasurementList,1) ) = [];
    
    for j=1:NumTracks
        pdaf.TrackList{j}.Filter.predict();
    end
   
    pdaf.associate();    
    pdaf.updateTracks();
        
    %store Logs
    for j=1:NumTracks
        Logs{j}.Estimates.StateMean(:,i) = pdaf.TrackList{j}.Filter.StateMean;
        Logs{j}.Estimates.StateCovar(:,:,i) = pdaf.TrackList{j}.Filter.StateCovar;
        Logs{j}.Groundtruth.State(:,i) = GroundTruth{i}(:,j);
    end

    if (ShowPlots)
        if(i==1 || rem(i,SkipFrames+1)==0)
            % Plot data
            clf;
             % Flip the image upside down before showing it
            imagesc([min_x max_x], [min_y max_y], flipud(img));

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            for j=1:NumTracks
                h2 = plot(Logs{j}.Groundtruth.State(1,1:i),Logs{j}.Groundtruth.State(2,1:i),'b.-','LineWidth',1);
                if j==2
                    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                h2 = plot(Logs{j}.Groundtruth.State(1,i),Logs{j}.Groundtruth.State(2,i),'bo','MarkerSize', 10);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
                h3 = plot(Logs{j}.Estimates.StateMean(1,1:i),Logs{j}.Estimates.StateMean(3,1:i),'r.-','LineWidth',1);
                h4=plot_gaussian_ellipsoid(Logs{j}.Estimates.StateMean([1 3],i), Logs{j}.Estimates.StateCovar([1 3],[1 3],i));
            end
            h2 = plot(DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);          
            % set the y-axis back to normal.
            set(gca,'ydir','normal');
            str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
            title(str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
            axis([0 10 0 10])
            pause(0.01)
        end
    end
end

    
ospa_vals= zeros(N,3);
ospa_c= 1;
ospa_p= 1;
for k=1:N
    trueX = [x_true(k,:);y_true(k,:)];
    estX = zeros(2,NumTracks);
    for i=1:NumTracks
        estX(:,i) = Logs{i}.Estimates.StateMean([1 3],k);
    end
    [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(trueX,estX,ospa_c,ospa_p);
end

figure
subplot(2,2,[1 2]), plot(1:k,ospa_vals(1:k,1));
subplot(2,2,3), plot(1:k,ospa_vals(1:k,2));
subplot(2,2,4), plot(1:k,ospa_vals(1:k,3));