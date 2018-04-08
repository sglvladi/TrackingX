classdef PhdTrackInitiatorX < TrackInitiatorX 
% PhdTrackInitiatorX Abstract class
%
% Summary of PhdTrackInitiatorX:
% This is an implementation of a Track Initiator which uses a PHD Filter
% to detect and initiate tracks.
%
% PhdTrackInitiatorX Properties:
%   PHDFilter
%   Filter
%   ConfirmThreshold
%   MeasurementList
%   UnassocWeights
%
% TrackInitiatorX Methods:
%   + TrackInitiatorX - Constructor method
%
% (+) denotes public properties/methods
%
% March 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
       MeasurementList
       UnassocWeights
       PHDFilter
       Filter
       ConfirmThreshold
       AssignedTrackIds = []
    end
        
    methods
        function this = PhdTrackInitiatorX(varargin)
        % PHDTRACKINITIATORX Constructor method
        %
        % Parameters
        % ----------
        % PHDFilter: SMC_PHDFilterX object
        %   A SMC_PHDFilterX object, which shall be used to perform
        %   filtering on the intensity of unconfirmed targets/tracks.
        % Filter: FilterX subclass object
        %   An instance of a subclass of a FilterX (e.g. KalmanFilterX,
        %   ParticleFilterX, etc) which shall be assigned newly initiated
        %   tracks.
        % ConfirmThreshold: scalar
        %   A scalar value in the range [0,1] which specifies the threshold
        %   probability of existance to be used to confirm tracks. 
        % 
        % Usage
        % -----
        % * PhdTrackInitiatorX(PHDFilter,Filter,ConfirmThreshold) returns 
        %   a PhdTrackInitiatorX object handle configured with the provided
        %   PHD filter instance PHDFilter, Filter instance Filter and
        %   confirmation threshold ConfirmThreshold.
        % * PhdTrackInitiatorX(config) can also be used, where config is a
        %   structure with fields config.PHDFilter, config.Filter and
        %   config.ConfirmThreshold.
        % * PhdTrackInitiatorX(___,Name,Value) instantiates an object 
        %   handle, configured with additional options specified by one or
        %   more Name,Value pair arguments.
            
            if(nargin==0)
                error('Not enough input arguments');
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.PHDFilter = config.PHDFilter;
                    this.Filter = config.Filter;
                    this.ConfirmThreshold = config.ConfirmThreshold;
                end
                return;
            elseif(nargin==3)
                this.PHDFilter = varargin{1};
                this.Filter = varargin{2};
                this.ConfirmThreshold = varargin{3};
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.PHDFilter = config.PHDFilter;
            this.Filter = config.Filter;
            this.ConfirmThreshold = config.ConfirmThreshold;
        end
        
        function NewTracks = initiateTracks(this,varargin)
        % INITIATETRACKS Perform track initiation
        % 
        % Parameters
        % ----------
        % MeasurementList: (NumObsDims x Nm) matrix, optional
        %   A matrix, whose columns correspond to individual measurements
        %   (default = this.MeasurementList)
        % UnassocWeights: (1 x Nm) row vector, optional
        %   A row vector, whose columns correspond to the probability that
        %   the respective measurements (columns) in MeasurementList have
        %   not been associated to any (already) confirmed tracks.
        %   (default = this. UnassocWeights)
        % ConfirmThreshold: scalar, optional
        %   A scalar value in the range [0,1] which specifies the threshold
        %   probability of existance to be used to confirm tracks.
        %
        % Returns
        % -------
        % NewTracks: (1 x Ni) cell array
        %   A cell array of TrackX objects, each corresponding to a newly
        %   inittiated track.
        %
        % Usage
        % -----
        % * NewTracks = initiateTracks(this,MeasurementList,UnassocWeights,ConfirmThreshold) 
        %   utilises the provided parameters to perform track initiation.
        %   Note that the provided parameters shall be used to update the
        %   instance's respective properties.
        % * NewTracks = initiateTracks(this) uses the instance's own
        
            if(nargin==4)
               this.ConfirmThreshold = ConfirmThreshold;
            end
            if(nargin==3)
               this.UnassocWeights = UnassocWeights;
            end
            if(nargin==2)
               this.MeasurementList = MeasurementList;
            end
            
            numMeas   = size(this.MeasurementList,2);

            % Predict the PHD filter
            this.PHDFilter.predict();

            % Set PHD Filter measurement and weights
            this.PHDFilter.MeasurementList = this.MeasurementList;
            this.PHDFilter.MeasWeights = this.UnassocWeights;

            % Rescale PHD likelihood matrix 
            measLikelihood = this.PHDFilter.MeasLikelihood;

            % Calculate p_i and w^{n,i} Eq. (21) of [2]
            w_ni = zeros(numMeas, this.PHDFilter.NumParticlesTotal);
            p_i = zeros(1,numMeas);
            for i = 1:numMeas
                w_ni(i,:) = (this.PHDFilter.ProbOfDetection*this.UnassocWeights(i)*measLikelihood(i,:)/ ...
                            (this.PHDFilter.ClutterRate + ... 
                                sum(this.PHDFilter.ProbOfDetection*measLikelihood(i,:).*this.PHDFilter.PredWeights,2))).*this.PHDFilter.PredWeights;
                p_i(i)= sum(w_ni(i,:),2);
            end

            % Select (critical) measurements to be used for initiating tracks
            CritMeasurements = find(p_i>this.ConfirmThreshold);

            % Initiate new tracks
            NumNewTracks = numel(CritMeasurements);
            if(NumNewTracks==0)
                NewTracks = {};
            else
                NewTracks = cell(size(CritMeasurements));
                for j = 1:NumNewTracks
                    measInd = CritMeasurements(j); % Index of measurement

                    % Get particles and weights
                    NewTrack = TrackX();
                    NewTrack.addprop('Filter');
                    NewTrack.addprop('TrackID');
                    NewTrack.addprop('ProbOfExist');
                    NewTrack.ProbOfExist = p_i(measInd);

                    % Generate track id
                    NewTrack.TrackID = randi(1000);
                    while(ismember(NewTrack.TrackID,this.AssignedTrackIds))
                        NewTrack.TrackID = randi(1000);
                    end
                    this.AssignedTrackIds = union(this.AssignedTrackIds,NewTrack.TrackID);

                    % Initiate filter
                    NewTrack.Filter = copy(this.Filter);
                    NewTrack.Filter.Particles = this.PHDFilter.PredParticles;
                    NewTrack.Filter.Weights = w_ni(measInd,:)/sum(w_ni(measInd,:));
                    % Resample particles to ensure they are correctly localised
                    % around the measurement
                    [NewTrack.Filter.Particles, NewTrack.Filter.Weights] = this.PHDFilter.Resampler.resample(NewTrack.Filter.Particles, NewTrack.Filter.Weights,NewTrack.Filter.NumParticles);
                    %NewTrack.ProbOfExist = p_i(measInd);

                    % append to new tracks
                    NewTracks{j} = NewTrack; 
                end
            end
            % Set weights of critical measurements to 0
            this.PHDFilter.MeasWeights(CritMeasurements) = 0;

            % Update PHD filter
            this.PHDFilter.update();
            
        end 
    end
end