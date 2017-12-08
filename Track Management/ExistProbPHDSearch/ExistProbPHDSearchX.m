function [jpdaf, phd, TrackIds] = ExistProbPHDSearchX(jpdaf, phd, filter, TrackIds)
    %%ExistProbPHDSearchX Perform track management using existence
    %                     probabilities and a PHD filter search track.
 
    % 3) Compute Existence Probabilities
    for t = 1:jpdaf.Params.nTracks
        jpdaf.Params.TrackList{t}.ExistProb = (1-phd.Params.pDeath)*jpdaf.Params.TrackList{t}.ExistProb;
        c = sum(jpdaf.Params.AssocWeightsMatrix(t,:))*jpdaf.Params.TrackList{t}.ExistProb + jpdaf.Params.AssocWeightsMatrix(t,1)*(1-jpdaf.Params.TrackList{t}.ExistProb)/((1-jpdaf.Params.pDetect*jpdaf.Params.pGate));
        jpdaf.Params.TrackList{t}.ExistProb = (sum(jpdaf.Params.AssocWeightsMatrix(t,:))*jpdaf.Params.TrackList{t}.ExistProb)/c;% 1-(pf.betta(1)*(1-pf.ExistProb))/c;
    end
    % 4) Predict the PHD filter
    phd.Predict();
    % 5) Compute rhi as given by Eq. (16) in [2] 
    rhi = zeros(1, size(jpdaf.Params.DataList,2));
    for m = 1:size(jpdaf.Params.DataList,2)
        rhi_tmp = 1;
        if(jpdaf.Params.AssocWeightsMatrix>-1) % Check if beta exists
            for t = 1:jpdaf.Params.nTracks
                rhi_tmp = rhi_tmp*(1-jpdaf.Params.AssocWeightsMatrix(t,m+1));
            end
        end
        rhi(m) = rhi_tmp;
    end
    phd.Params.rhi = rhi;%==1;
    % 6) Update PHD filter
    phd.Update();
    % 7) Initiate tracks
    for j = 1:length(phd.Params.NewTracks)
        disp("Initating new track");
        if(isa(filter,'KalmanFilterX'))
            filter.Params.x = sum(phd.Params.NewTracks{j}.particles.*phd.Params.NewTracks{j}.w,2);
            filter.Params.P = weightedcov(phd.Params.NewTracks{j}.particles, phd.Params.NewTracks{j}.w);
        else
            filter.Params.particles = phd.Params.NewTracks{j}.particles;
            filter.Params.w = phd.Params.NewTracks{j}.w;
            filter.Params.x = sum(phd.Params.NewTracks{j}.particles.*phd.Params.NewTracks{j}.w,2);
        end
        jpdaf.Params.TrackList{end+1}.TrackObj = copy(filter);
        jpdaf.Params.TrackList{end}.ExistProb = phd.Params.NewTracks{j}.ExistProb;
        newTrackId = 1;
        while(ismember(newTrackId,TrackIds))
            newTrackId = newTrackId + 1;
        end
        jpdaf.Params.TrackList{end}.trackId = newTrackId;
        TrackIds = [TrackIds, newTrackId];
        jpdaf.Params.nTracks = jpdaf.Params.nTracks + 1;
    end
    % 8) Delete tracks
    del_tracks = 0;
    del_flag = 0;
    for t = 1:jpdaf.Params.nTracks
        if(jpdaf.Params.TrackList{t}.ExistProb<0.1)
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