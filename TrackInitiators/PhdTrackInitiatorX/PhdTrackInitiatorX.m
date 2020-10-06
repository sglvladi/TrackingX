classdef PhdTrackInitiatorX < TrackInitiatorX 
% PhdTrackInitiatorX class
%
% Summary of PhdTrackInitiatorX:
% This is an implementation of a TrackInitiator that initiates tracks using
% a SMC-PHD filter, using the approach described in [2].
%
% PhdTrackInitiatorX Properties:
%   None
%
% PhdTrackInitiatorX Methods:
%   + PhdTrackInitiatorX - Constructor method
%
% (+) denotes public properties/methods
%
% [1]  P. Horridge and S. Maskell,  “Using a probabilistic hypothesis density 
%       filter to confirm tracks in a multi-target environment,” in 2011 
%       Jahrestagung der Gesellschaft fr Informatik, October 2011.
%
% April 2019 Lyudmil Vladimirov, University of Liverpool.
    
    properties
       TagGenerator
       PhdFilter
       InitFilter
       ConfirmThreshold
    end
        
    methods
        function this = PhdTrackInitiatorX(varargin)
        % PHDTRACKINITIATORX Constructor method
        %   
        % DESCRIPTION: 
        % * PhdTrackInitiatorX(config) 
        %   returns a PhdTrackInitiatorX object handle configured with the provided
        %   PHD filter instance PhdFilter and probability of gating ProbOfGating
            
            if(nargin==0)
                error('Not enough input arguments');
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.PhdFilter = config.PhdFilter;
                    this.TagGenerator = config.TagGenerator;
                    this.InitFilter = config.InitFilter;
                    this.ConfirmThreshold = config.ConfirmThreshold;
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
            this.PhdFilter = config.PhdFilter;
            this.TagGenerator = config.TagGenerator;
            this.InitFilter = config.InitFilter;
            this.DeleteThreshold = config.DeleteThreshold;
        end
        
        function [TrackList] = initiateTracks(this,TrackList, MeasurementList, AssocWeightsMatrix)
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
        
            numTracks = numel(TrackList);
            numMeas   = MeasurementList.NumMeasurements;
            
            % Predict the PHD filter
            this.PhdFilter.predict();

            % Compute rhi as given by Eq. (16) in [2] 
            rhi = zeros(1, numMeas);
            for measInd = 1:numMeas
                rhi_tmp = 1;
                if(AssocWeightsMatrix>-1) % Check if beta exists
                    for trackInd = 1:numTracks
                        rhi_tmp = rhi_tmp*(1-AssocWeightsMatrix(trackInd,measInd+1));
                    end
                end
                rhi(measInd) = rhi_tmp;
            end

            % Set PHD Filter measurement and weights
            this.PhdFilter.MeasurementList = MeasurementList;
            this.PhdFilter.MeasWeights = rhi;

            % Rescale PHD likelihood matrix 
            measLikelihood = this.PhdFilter.MeasurementLikelihoodsPerParticle;
            %measLikelihood = this.PhdFilter.MeasurementLikelihoodsPerComponent;

            % Calculate p_i and w^{n,i} Eq. (21) of [2]
            w_ni = zeros(numMeas, size(measLikelihood,2));
            p_i = zeros(1,numMeas);
            lambda = this.PhdFilter.Model.Clutter.pdf(this.PhdFilter.MeasurementList.Vectors);
            for i = 1:numMeas
                P_D = mean(this.PhdFilter.Model.Detection.pdf(this.PhdFilter.StatePrediction.Particles));
                %P_D = mean(this.PhdFilter.Model.Detection.pdf(this.PhdFilter.StatePrediction.Means));
                w_ni(i,:) = (P_D.*rhi(i).*measLikelihood(i,:)./ ...
                            (lambda(i) + sum(P_D.*measLikelihood(i,:).*this.PhdFilter.StatePrediction.Weights,2))).*this.PhdFilter.StatePrediction.Weights;
                p_i(i)= sum(w_ni(i,:),2);
            end
%             numNewTracks = floor(sum(p_i))
%             % Select measurements to be used for spawning new tracks
%             [p_i, ind] = sort(p_i, 'descend');
            CritMeasurements = find(p_i>this.ConfirmThreshold);
%             numCrit = numel(CritMeasurements);
%             if(numCrit && numNewTracks<numCrit)
%                 CritMeasurements = ind(CritMeasurements(1:numNewTracks));
%             else
%                 CritMeasurements = ind(CritMeasurements);
%             end
%             p_i = p_i(ind);
            % Initiate new tracks
            for j = 1:size(CritMeasurements,2)
                measInd = CritMeasurements(j); % Index of measurement
                
                measurement = this.PhdFilter.MeasurementList.Measurements(measInd);
                
                particles = this.PhdFilter.StatePrediction.Particles;
                weights = w_ni(measInd,:)/sum(w_ni(measInd,:));
                dist = ParticleDistributionX(particles, weights);
                dist.resample(5000);
                if(~isempty(measurement.Timestamp))
                    statePrior = ParticleStateX(dist,...
                                                measurement.Timestamp);
                else
                    statePrior = ParticleStateX(dist);
                end
                track_config = struct('Filter', copy(this.InitFilter), ...
                                      'Tag', this.TagGenerator.generate(), ...
                                      'ExistenceProbability', p_i(measInd));
                track_config.Filter.initialise('Model',this.InitFilter.Model,...
                                               'StatePrior', statePrior); 
                track_config.Trajectory = track_config.Filter.StatePosterior;
                newTrack = TrackX(track_config);
                TrackList{end+1} = newTrack; 
            end

% REMOVED: Should be handled by rhi  
%             % Adjust measurement weights to exclude measurements used for
%             % track initiation
%              this.PhdFilter.MeasWeights(CritMeasurements) = 0;

            % 6) Update PHD filter
            this.PhdFilter.update();
            
        end 
    end
end