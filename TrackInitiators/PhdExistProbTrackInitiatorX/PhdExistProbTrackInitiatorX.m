classdef PhdExistProbTrackInitiatorX < TrackInitiatorX 
% PhdExistProbTrackInitiatorX Abstract class
%
% Summary of PhdExistProbTrackInitiatorX:
% This is the base class for all TrackingX resamplers.
% Any custom defined Track Initiator should be derived from this PhdExistProbTrackInitiatorX base class. 
%
% TrackInitiatorX Properties:
%   None
%
% TrackInitiatorX Methods:
%   + TrackInitiatorX - Constructor method
%
% (+) denotes public properties/methods
%
% March 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
       TrackList
       MeasurementList
       AssocWeightsMatrix
       PHDFilter
       Filter
       ProbOfGating
       ProbOfConfirm
    end
        
    methods
        function this = PhdExistProbTrackInitiatorX(varargin)
        % PHDEXISTPROBTRACKINITIATORX Constructor method
        %   
        % DESCRIPTION: 
        % * PhdExistProbTrackInitiatorX('PHDFilter',PHDFilter,'ProbOfGating',ProbOfGating) 
        %   returns a PhdExistProbTrackInitiatorX object handle configured with the provided
        %   PHD filter instance PHDFilter and probability of gating ProbOfGating
            
            if(nargin==0)
                error('Not enough input arguments');
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.PHDFilter = config.PHDFilter;
                    this.Filter = config.Filter;
                    this.ProbOfGating = config.ProbOfGating;
                    this.ProbOfConfirm = config.ProbOfConfirm;
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = varargin{1};
            this.PHDFilter = config.PHDFilter;
            this.Filter = config.Filter;
            this.ProbOfGating = config.ProbOfGating;
            this.ProbOfConfirm = config.ProbOfConfirm;
        end
        
        function [TrackList] = initiateTracks(this,varargin)
        % INITIATETRACKS Perform track initiation
        %   
        % DESCRIPTION: 
        % * [TrackList] = initiateTracks(this,TrackList,MeasurementList) 
        %   utilises the provided lists of tracks (TrackList) and measurements
        %   (MeasurementList) 
        %
        % INPUTS:
        % * TrackList - (1-by-NumTracks) array of TrackX objects
        % * MeasurementList - (NumObsDims-by-NumMeas) matrix of
        %   measurements
        %
        % See also EllipsoidalGaterX/EllipsoidalGaterX
        
            if(nargin==1)
                TrackList = this.TrackList;
                MeasurementList = this.MeasurementList;
            end
            numTracks = numel(TrackList);
            numMeas   = size(MeasurementList,2);

            % 3) Compute Existence Probabilities
            for trackInd = 1:numTracks
                TrackList{trackInd}.ProbOfExist = (1 - this.PHDFilter.ProbOfDeath)*TrackList{trackInd}.ProbOfExist;
                %Lj_total = (1-this.PHDFilter.ProbOfDetection*this.ProbOfGating) ...
                %            + this.PHDFilter.ProbOfDetection*this.ProbOfGating/TrackList{trackInd}.ClutterDensity*sum(TrackList{trackInd}.Filter.MeasLikelihood,2)';
                %Lj = this.AssocWeightsMatrix(trackInd,1)*Lj_total...
                %     + Lj_total*sum(this.AssocWeightsMatrix(trackInd,find(TrackList{trackInd}.ValidationMatrix)));
                 denom = sum(this.AssocWeightsMatrix(trackInd,:))*TrackList{trackInd}.ProbOfExist + this.AssocWeightsMatrix(trackInd,1)*(1-TrackList{trackInd}.ProbOfExist)/((1-this.PHDFilter.ProbOfDetection*this.ProbOfGating));
                 TrackList{trackInd}.ProbOfExist = (sum(this.AssocWeightsMatrix(trackInd,:))*TrackList{trackInd}.ProbOfExist)/denom;
    %            TrackList{trackInd}.ProbOfExist = TrackList{trackInd}.ProbOfExist*Lj/(1-(1-Lj)*TrackList{trackInd}.ProbOfExist);
            end

            % Predict the PHD filter
            this.PHDFilter.predict();

            % Compute rhi as given by Eq. (16) in [2] 
            rhi = zeros(1, numMeas);
            for measInd = 1:numMeas
                rhi_tmp = 1;
                if(this.AssocWeightsMatrix>-1) % Check if beta exists
                    for trackInd = 1:numTracks
                        rhi_tmp = rhi_tmp*(1-this.AssocWeightsMatrix(trackInd,measInd+1));
                    end
                end
                rhi(measInd) = rhi_tmp;
            end

            % Set PHD Filter measurement and weights
            this.PHDFilter.MeasurementList = MeasurementList;
            this.PHDFilter.MeasWeights = rhi;

            % Rescale PHD likelihood matrix 
            measLikelihood = this.PHDFilter.MeasLikelihood;

            % Calculate p_i and w^{n,i} Eq. (21) of [2]
            w_ni = zeros(numMeas, this.PHDFilter.NumParticlesTotal);
            p_i = zeros(1,numMeas);
            for i = 1:numMeas
                w_ni(i,:) = (this.PHDFilter.ProbOfDetection*rhi(i)*measLikelihood(i,:)/ ...
                            (this.PHDFilter.ClutterRate + ... 
                                sum(this.PHDFilter.ProbOfDetection*measLikelihood(i,:).*this.PHDFilter.PredWeights,2))).*this.PHDFilter.PredWeights;
                p_i(i)= sum(w_ni(i,:),2);
            end

            % Select measurements to be used for spawning new tracks
            CritMeasurements = find(p_i>this.ProbOfConfirm);
            NewTracks = [];

            % Initiate new tracks
            for j = 1:size(CritMeasurements,2)
                measInd = CritMeasurements(j); % Index of measurement

                % Get particles and weights
                NewTrack = TrackX();
                NewTrack.addprop('Filter');
                NewTrack.addprop('TrackID');
                NewTrack.addprop('ProbOfExist');
                NewTrack.ProbOfExist = p_i(measInd);
                NewTrack.TrackID = randi(1000);
                NewTrack.Filter = copy(this.Filter);
                NewTrack.Filter.Particles = this.PHDFilter.PredParticles;
                NewTrack.Filter.Weights = w_ni(measInd,:)/sum(w_ni(measInd,:));
                % Resample particles to ensure they are correctly localised
                % around the measurement
                [NewTrack.Filter.Particles, NewTrack.Filter.Weights] = this.PHDFilter.Resampler.resample(NewTrack.Filter.Particles, NewTrack.Filter.Weights,NewTrack.Filter.NumParticles);
                NewTrack.ProbOfExist = p_i(measInd);
                TrackList{end+1} = NewTrack;
                %this.Params.NewTracks{end+1} = NewTrack; 
            end

            % 
            this.PHDFilter.MeasWeights(CritMeasurements) = 0;

            % 6) Update PHD filter
            this.PHDFilter.update();
            
            for trackInd = 1:numTracks
                if(TrackList{trackInd}.ProbOfExist<0.1)
                   TrackList{trackInd} = [];
                end
%                 TrackList{trackInd}.ProbOfExist = (1 - this.PHDFilter.ProbOfDeath)*TrackList{trackInd}.ProbOfExist;
%                 %Lj_total = (1-this.PHDFilter.ProbOfDetection*this.ProbOfGating) ...
%                 %            + this.PHDFilter.ProbOfDetection*this.ProbOfGating/TrackList{trackInd}.ClutterDensity*sum(TrackList{trackInd}.Filter.MeasLikelihood,2)';
%                 %Lj = this.AssocWeightsMatrix(trackInd,1)*Lj_total...
%                 %     + Lj_total*sum(this.AssocWeightsMatrix(trackInd,find(TrackList{trackInd}.ValidationMatrix)));
%                  denom = sum(this.AssocWeightsMatrix(trackInd,:))*TrackList{trackInd}.ProbOfExist + this.AssocWeightsMatrix(trackInd,1)*(1-TrackList{trackInd}.ProbOfExist)/((1-this.PHDFilter.ProbOfDetection*this.ProbOfGating));
%                  TrackList{trackInd}.ProbOfExist = (sum(this.AssocWeightsMatrix(trackInd,:))*TrackList{trackInd}.ProbOfExist)/denom;
    %            TrackList{trackInd}.ProbOfExist = TrackList{trackInd}.ProbOfExist*Lj/(1-(1-Lj)*TrackList{trackInd}.ProbOfExist);
            end
            if(numel(TrackList))
                TrackList = TrackList(~cellfun('isempty',TrackList));
            end
                %TrackList( :, ~any(TrackList,1) ) = [];

    %         % 7) Initiate tracks
    %         for j = 1:length(phd.Params.NewTracks)
    %             disp("Initating new track");
    %             if(isa(filter,'KalmanFilterX'))
    %                 filter.Params.x = sum(phd.Params.NewTracks{j}.particles.*phd.Params.NewTracks{j}.w,2);
    %                 filter.Params.P = weightedcov(phd.Params.NewTracks{j}.particles, phd.Params.NewTracks{j}.w);
    %             else
    %                 filter.Params.particles = phd.Params.NewTracks{j}.particles;
    %                 filter.Params.w = phd.Params.NewTracks{j}.w;
    %                 filter.Params.x = sum(phd.Params.NewTracks{j}.particles.*phd.Params.NewTracks{j}.w,2);
    %             end
    %             jpdaf.Params.TrackList{end+1}.TrackObj = copy(filter);
    %             jpdaf.Params.TrackList{end}.ExistProb = phd.Params.NewTracks{j}.ExistProb;
    %             newTrackId = 1;
    %             while(ismember(newTrackId,TrackIds))
    %                 newTrackId = newTrackId + 1;
    %             end
    %             jpdaf.Params.TrackList{end}.trackId = newTrackId;
    %             TrackIds = [TrackIds, newTrackId];
    %             jpdaf.Params.nTracks = jpdaf.Params.nTracks + 1;
        end 
    end
end

% function [jpdaf, phd, TrackIds] = ExistProbPHDSearchX(jpdaf, phd, filter, TrackIds)
%     %%ExistProbPHDSearchX Perform track management using existence
%     %                     probabilities and a PHD filter search track.
%  
%     % 3) Compute Existence Probabilities
%     for t = 1:jpdaf.Params.nTracks
%         jpdaf.Params.TrackList{t}.ExistProb = (1-phd.Params.pDeath)*jpdaf.Params.TrackList{t}.ExistProb;
%         c = sum(jpdaf.Params.AssocWeightsMatrix(t,:))*jpdaf.Params.TrackList{t}.ExistProb + jpdaf.Params.AssocWeightsMatrix(t,1)*(1-jpdaf.Params.TrackList{t}.ExistProb)/((1-jpdaf.Params.pDetect*jpdaf.Params.pGate));
%         jpdaf.Params.TrackList{t}.ExistProb = (sum(jpdaf.Params.AssocWeightsMatrix(t,:))*jpdaf.Params.TrackList{t}.ExistProb)/c;% 1-(pf.betta(1)*(1-pf.ExistProb))/c;
%     end
%     % 4) Predict the PHD filter
%     phd.Predict();
%     % 5) Compute rhi as given by Eq. (16) in [2] 
%     rhi = zeros(1, size(jpdaf.Params.DataList,2));
%     for m = 1:size(jpdaf.Params.DataList,2)
%         rhi_tmp = 1;
%         if(jpdaf.Params.AssocWeightsMatrix>-1) % Check if beta exists
%             for t = 1:jpdaf.Params.nTracks
%                 rhi_tmp = rhi_tmp*(1-jpdaf.Params.AssocWeightsMatrix(t,m+1));
%             end
%         end
%         rhi(m) = rhi_tmp;
%     end
%     phd.Params.rhi = rhi;%==1;
%     % 6) Update PHD filter
%     phd.Update();
%     % 7) Initiate tracks
%     for j = 1:length(phd.Params.NewTracks)
%         disp("Initating new track");
%         if(isa(filter,'KalmanFilterX'))
%             filter.Params.x = sum(phd.Params.NewTracks{j}.particles.*phd.Params.NewTracks{j}.w,2);
%             filter.Params.P = weightedcov(phd.Params.NewTracks{j}.particles, phd.Params.NewTracks{j}.w);
%         else
%             filter.Params.particles = phd.Params.NewTracks{j}.particles;
%             filter.Params.w = phd.Params.NewTracks{j}.w;
%             filter.Params.x = sum(phd.Params.NewTracks{j}.particles.*phd.Params.NewTracks{j}.w,2);
%         end
%         jpdaf.Params.TrackList{end+1}.TrackObj = copy(filter);
%         jpdaf.Params.TrackList{end}.ExistProb = phd.Params.NewTracks{j}.ExistProb;
%         newTrackId = 1;
%         while(ismember(newTrackId,TrackIds))
%             newTrackId = newTrackId + 1;
%         end
%         jpdaf.Params.TrackList{end}.trackId = newTrackId;
%         TrackIds = [TrackIds, newTrackId];
%         jpdaf.Params.nTracks = jpdaf.Params.nTracks + 1;
%     end
%     % 8) Delete tracks
%     del_tracks = 0;
%     del_flag = 0;
%     for t = 1:jpdaf.Params.nTracks
%         if(jpdaf.Params.TrackList{t}.ExistProb<0.1)
%             jpdaf.Params.TrackList{t} = [];
%             del_tracks = del_tracks + 1;
%             del_flag = 1;
%         end
%     end
%     if(del_flag)
%         jpdaf.Params.TrackList = jpdaf.Params.TrackList(~cellfun('isempty',jpdaf.Params.TrackList));
%         jpdaf.Params.nTracks = jpdaf.Params.nTracks - del_tracks;
%     end
% end