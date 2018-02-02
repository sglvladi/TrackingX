function [jpdaf, jpdaf_init,TrackIds] = MoutOfN_TM(jpdaf, jpdaf_init, thresholds, filter, TrackIds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track_Maintenance - Performs track initiation, confirmation and deletion 
% Input:
%   TrackList        - List of Tracks
%   ValidationMatrix - Matrix showing all valid data associations
% Output:
%   TrackList    - Updated list of Tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % High and low thresholds for confirmation and deletion
    %  of tentative tracks
    confThresh  = thresholds(1);       
    delThresh   = thresholds(2);
    sigmaDelThresh = thresholds(3);

    %% Confirmed tracks
    if jpdaf.Params.nTracks>0

        % Update Miss counters for all tracks
        for t=1:jpdaf.Params.nTracks

            % Update counts
            if isempty(find(jpdaf.Params.ValidationMatrix(t,:), 1))
                jpdaf.Params.TrackList{t}.missedCounts = jpdaf.Params.TrackList{t}.missedCounts + 1;
            else
                if(jpdaf.Params.TrackList{t}.missedCounts>0)
                    jpdaf.Params.TrackList{t}.missedCounts = jpdaf.Params.TrackList{t}.missedCounts -1;
                end
            end 
        end
        
        % Delete Tracks
        del_tracks = 0;
        del_flag = 0;
        for t = 1:jpdaf.Params.nTracks
            if(isa(jpdaf.Params.TrackList{t}.TrackObj,'KalmanFilterX'))
                P = jpdaf.Params.TrackList{t}.TrackObj.Params.P; 
            else
                P  = weightedcov(jpdaf.Params.TrackList{t}.TrackObj.Params.particles',jpdaf.Params.TrackList{t}.TrackObj.Params.w);
            end
            if(jpdaf.Params.TrackList{t}.missedCounts > delThresh|| P(1,1)>sigmaDelThresh^2 || P(2,2)>sigmaDelThresh^2)
                jpdaf.Params.TrackList{t} = [];
                del_tracks = del_tracks + 1;
                del_flag = 1;
            end
        end
        if(del_flag)
            jpdaf.Params.TrackList = jpdaf.Params.TrackList(~cellfun('isempty',jpdaf.Params.TrackList));
            jpdaf.Params.nTracks = jpdaf.Params.nTracks - del_tracks;
        end
    end
    
    %% Tentative tracks 
    
    % Get indices of all measurements where there exist no possible association
    UnassocMeasInd = find(~any(jpdaf.Params.ValidationMatrix,1));
    
    jpdaf_init.Params.DataList = jpdaf.Params.DataList(:,UnassocMeasInd); % Only use unassociated measurements
    
    % 1) Predict the confirmed tracks
    jpdaf_init.Predict();
    
    % 2) Update the confirmed track
    jpdaf_init.Update();
    
    % Update Detection and Miss counters for all tracks
    for t=1:jpdaf_init.Params.nTracks
        
        % Update counts
        if isempty(find(jpdaf_init.Params.ValidationMatrix(t,:), 1))
            if(jpdaf_init.Params.TrackList{t}.detectedCounts>0)
                jpdaf_init.Params.TrackList{t}.detectedCounts = jpdaf_init.Params.TrackList{t}.detectedCounts - 1;
            end
            jpdaf_init.Params.TrackList{t}.missedCounts = jpdaf_init.Params.TrackList{t}.missedCounts + 1;
        else
            jpdaf_init.Params.TrackList{t}.detectedCounts = jpdaf_init.Params.TrackList{t}.detectedCounts + 1;
            if(jpdaf_init.Params.TrackList{t}.missedCounts>0)
                jpdaf_init.Params.TrackList{t}.missedCounts = jpdaf_init.Params.TrackList{t}.missedCounts - 1;
            end
        end
    end
    
    % Delete or Confirm Tentative Tracks
    del_tracks = 0;
    del_flag = 0;
    for t = 1:jpdaf_init.Params.nTracks
        if(isa(jpdaf_init.Params.TrackList{t}.TrackObj,'KalmanFilterX'))
            P = jpdaf_init.Params.TrackList{t}.TrackObj.Params.P(1:2,1:2); 
        else
            P  = weightedcov(jpdaf_init.Params.TrackList{t}.TrackObj.Params.particles',jpdaf_init.Params.TrackList{t}.TrackObj.Params.w);
        end
        
        if(jpdaf_init.Params.TrackList{t}.missedCounts > delThresh || P(1,1)>sigmaDelThresh^2 || P(2,2)>sigmaDelThresh^2) 
            jpdaf_init.Params.TrackList{t} = []; 
            del_tracks = del_tracks + 1;
            del_flag = 1;
        elseif(jpdaf_init.Params.TrackList{t}.detectedCounts > confThresh)
            disp("Confirmig new track");
            jpdaf.Params.TrackList{end+1}.TrackObj = copy(jpdaf_init.Params.TrackList{t}.TrackObj);
            jpdaf.Params.TrackList{end}.missedCounts = 0;
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
%         filter.Params.particles = filter.Params.gen_x0([jpdaf_init.Params.DataList(1,UnassocMeasInd(i)); jpdaf_init.Params.DataList(2,UnassocMeasInd(i));0;0],filter.Params.Np);
%         filter.Params.x = sum(filter.Params.w.*filter.Params.particles,2);
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
        jpdaf_init.Params.TrackList{end}.missedCounts = 0;
        jpdaf_init.Params.TrackList{end}.detectedCounts = 0;
        jpdaf_init.Params.nTracks = jpdaf_init.Params.nTracks + 1;  
    end
end