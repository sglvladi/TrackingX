function [jpdaf, jpdaf_init,TrackIds] = LogLikRatio_TM(jpdaf, jpdaf_init, thresholds, filter,TrackIds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track_Maintenance - Performs track initiation, confirmation and deletion 
% Input:
%   TrackList        - List of Tracks
%   ValidationMatrix - Matrix showing all valid data associations
% Output:
%   TrackList    - Updated list of Tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Confirmed tracks
    P_TM = thresholds(1); %0.1;                         % Track miss probability (from Blackman & Popoli)
    P_FC = thresholds(2); %10^-6;                     % False confirm probability (from Blackman & Popoli)
    THD = thresholds(3); %-5.9;                         % Deletion threshold diff (from Blackman & Popoli)
    %Vmax = 0.4;
    
    % High and low thresholds for confirmation and deletion
    %  of tentative tracks
    gamma_low = log(P_TM/(1-P_FC));
    gamma_high = log((1-P_TM)/P_FC);
    
    % 1) Compute LogLikelihoodRatios
    for t = 1:jpdaf.Params.nTracks
        
        % Get the index of the cluster which track belongs to
        % (Needed to extract lambda (clutter rate)
        for j=1:size(jpdaf.Params.ClusterList,2)
            Cluster = jpdaf.Params.ClusterList{j};
            if (ismember(t, Cluster.TrackIndList)~=0)
                clusterInd = j;
                break;
            end
        end
        
        % Update LPR
        jpdaf.Params.TrackList{t}.LPR = jpdaf.Params.TrackList{t}.LPR + jpdaf.Params.AssocWeightsMatrix(t,1)*log(1-jpdaf.Params.pDetect*jpdaf.Params.pGate);
        ValidDataInd = find(jpdaf.Params.ValidationMatrix(t,:));
        for j=1:numel(ValidDataInd)
            dataInd = ValidDataInd(j);
            jpdaf.Params.TrackList{t}.LPR = jpdaf.Params.TrackList{t}.LPR + jpdaf.Params.AssocWeightsMatrix(t,dataInd+1)*log(jpdaf.Params.pDetect*jpdaf.Params.LikelihoodMatrix(t,dataInd)/jpdaf.Params.ClusterList{clusterInd}.lambda);
        end
        % Update maximum LPR
        if (jpdaf.Params.TrackList{t}.LPR>jpdaf.Params.TrackList{t}.LPR_max)
            jpdaf.Params.TrackList{t}.LPR_max = jpdaf.Params.TrackList{t}.LPR;
        end
    end
    
     % 2) Delete tracks
    del_tracks = 0;
    del_flag = 0;
    for t = 1:jpdaf.Params.nTracks
        % Check against thresholds
        if (isnan(jpdaf.Params.TrackList{t}.LPR) || jpdaf.Params.TrackList{t}.LPR < jpdaf.Params.TrackList{t}.LPR_max+THD)
            jpdaf.Params.TrackList{t} = [];
            del_tracks = del_tracks + 1;
            del_flag = 1;
        end
    end
    if(del_flag)
        jpdaf.Params.TrackList = jpdaf.Params.TrackList(~cellfun('isempty',jpdaf.Params.TrackList));
        jpdaf.Params.nTracks = jpdaf.Params.nTracks - del_tracks;
    end
    
    
    %% Tentative tracks 
    % Get indices of all measurements where there exist no possible association
    UnassocMeasInd = find(~any(jpdaf.Params.ValidationMatrix,1));
    
    jpdaf_init.Params.DataList = jpdaf.Params.DataList(:,UnassocMeasInd); % Only use unassociated measurements
    
    % 1) Predict the confirmed tracks
    jpdaf_init.Predict();
    
    % 2) Update the confirmed track
    jpdaf_init.Update();
    
    % 1) Compute LogLikelihoodRatios
    for t = 1:jpdaf_init.Params.nTracks
        
        % Get the index of the cluster which track belongs to
        % (Needed to extract lambda (clutter rate)
        for j=1:size(jpdaf_init.Params.ClusterList,2)
            Cluster = jpdaf_init.Params.ClusterList{j};
            if (ismember(t, Cluster.TrackIndList)~=0)
                clusterInd = j;
                break;
            end
        end
        
        % Update LPR
        jpdaf_init.Params.TrackList{t}.LPR = jpdaf_init.Params.TrackList{t}.LPR + jpdaf_init.Params.AssocWeightsMatrix(t,1)*log(1-jpdaf_init.Params.pDetect*jpdaf_init.Params.pGate);
        ValidDataInd = find(jpdaf_init.Params.ValidationMatrix(t,:));
        for j=1:numel(ValidDataInd)
            dataInd = ValidDataInd(j);
            jpdaf_init.Params.TrackList{t}.LPR = jpdaf_init.Params.TrackList{t}.LPR + jpdaf_init.Params.AssocWeightsMatrix(t,dataInd+1)*log(jpdaf_init.Params.pDetect*jpdaf_init.Params.LikelihoodMatrix(t,dataInd)/jpdaf_init.Params.ClusterList{clusterInd}.lambda);
        end
        % Update maximum LPR
        if (jpdaf_init.Params.TrackList{t}.LPR>jpdaf_init.Params.TrackList{t}.LPR_max)
            jpdaf_init.Params.TrackList{t}.LPR_max = jpdaf_init.Params.TrackList{t}.LPR;
        end
    end
    
    % 2) Delete tracks
    del_tracks = 0;
    del_flag = 0;
    for t = 1:jpdaf_init.Params.nTracks
        if (isnan(jpdaf_init.Params.TrackList{t}.LPR) || jpdaf_init.Params.TrackList{t}.LPR < gamma_low)
            jpdaf_init.Params.TrackList{t} = [];
            del_tracks = del_tracks + 1;
            del_flag = 1;
        elseif (jpdaf_init.Params.TrackList{t}.LPR > gamma_high)
            disp("Confirmig new track");
            jpdaf.Params.TrackList{end+1}.TrackObj = copy(jpdaf_init.Params.TrackList{t}.TrackObj);
            jpdaf.Params.TrackList{end}.LPR = jpdaf_init.Params.TrackList{t}.LPR;
            jpdaf.Params.TrackList{end}.LPR_max = jpdaf_init.Params.TrackList{t}.LPR_max;
            newTrackId = 1;
            while(ismember(newTrackId,TrackIds))
                newTrackId = newTrackId + 1;
            end
            jpdaf.Params.TrackList{end}.trackId = newTrackId;
            TrackIds = [TrackIds, newTrackId];
            jpdaf.Params.nTracks = jpdaf.Params.nTracks + 1;
            jpdaf_init.Params.TrackList{t} = [];
            del_tracks = del_tracks + 1;
            del_flag = 1;
        end
    end
    if(del_flag)
        jpdaf_init.Params.TrackList = jpdaf_init.Params.TrackList(~cellfun('isempty',jpdaf_init.Params.TrackList));
        jpdaf_init.Params.nTracks = jpdaf_init.Params.nTracks - del_tracks;
    end
        
    UnassocMeasInd = find(~any(jpdaf_init.Params.ValidationMatrix,1));
    % Initiate new tracks
    for i=1:numel(UnassocMeasInd)
        %disp("Initating new track");
        if(isa(filter,'KalmanFilterX'))
            filter.Params.x = [jpdaf_init.Params.DataList(1,UnassocMeasInd(i)); jpdaf_init.Params.DataList(2,UnassocMeasInd(i));0;0];
            %filter.Params.P = weightedcov(phd.Params.NewTracks{j}.particles, phd.Params.NewTracks{j}.w);
        else
            filter.Params.particles = filter.Params.gen_x0([jpdaf_init.Params.DataList(1,UnassocMeasInd(i)); jpdaf_init.Params.DataList(2,UnassocMeasInd(i));0;0],filter.Params.Np);
            filter.Params.x = sum(filter.Params.w.*filter.Params.particles,2);
        end
        %pf.gen_x0 = @(Np) mvnrnd(repmat([DataList(1,MeasInd), DataList(2,MeasInd), 0, 0],Np,1),diag([q^2, q^2, 0.3^2, 4*pi^2]));
        %s.x = [DataList(1,MeasInd); DataList(2,MeasInd); 0; 0]; %initial state
        %s.P = diag([q^2, q^2, (Vmax^2/3), (Vmax^2/3)]);
        jpdaf_init.Params.TrackList{end+1}.TrackObj = copy(filter);
        jpdaf_init.Params.TrackList{end}.LPR = 1;
        jpdaf_init.Params.TrackList{end}.LPR_max = jpdaf_init.Params.TrackList{end}.LPR;
        jpdaf_init.Params.nTracks = jpdaf_init.Params.nTracks + 1;  
    end
end