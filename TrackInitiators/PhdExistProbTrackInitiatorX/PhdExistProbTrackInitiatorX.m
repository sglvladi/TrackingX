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
            config = parser.Results;
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
%             for trackInd = 1:numel(TrackList)
%                 TrackList{trackInd}.ProbOfExist = (1 - this.PHDFilter.ProbOfDeath)*TrackList{trackInd}.ProbOfExist;
%                 denom = sum(this.AssocWeightsMatrix(trackInd,:))*TrackList{trackInd}.ProbOfExist + this.AssocWeightsMatrix(trackInd,1)*(1-TrackList{trackInd}.ProbOfExist);%/((1-this.PHDFilter.ProbOfDetection*this.ProbOfGating));
%                 TrackList{trackInd}.ProbOfExist = (sum(this.AssocWeightsMatrix(trackInd,:))*TrackList{trackInd}.ProbOfExist)/denom;
%             end
            
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
            rhi = rhi==1;

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
            end

            % 
            this.PHDFilter.MeasWeights(CritMeasurements) = 0;

            % 6) Update PHD filter
            this.PHDFilter.update();
            
            for trackInd = 1:numTracks
                if(TrackList{trackInd}.ProbOfExist<0.1)
                   TrackList{trackInd} = [];
                end
            end
            if(numel(TrackList))
                TrackList = TrackList(~cellfun('isempty',TrackList));
            end
            
            % 3) Compute Existence Probabilities
            for trackInd = 1:numel(TrackList)
                TrackList{trackInd}.ProbOfExist = (1 - this.PHDFilter.ProbOfDeath)*TrackList{trackInd}.ProbOfExist;
%                 denom = sum(this.AssocWeightsMatrix(trackInd,:))*TrackList{trackInd}.ProbOfExist + this.AssocWeightsMatrix(trackInd,1)*(1-TrackList{trackInd}.ProbOfExist);%/((1-this.PHDFilter.ProbOfDetection*this.ProbOfGating));
%                 TrackList{trackInd}.ProbOfExist = (sum(this.AssocWeightsMatrix(trackInd,:))*TrackList{trackInd}.ProbOfExist)/denom;
            end
        end 
    end
end