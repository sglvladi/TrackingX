classdef EllipsoidalGaterX <GaterX
% ELLIPSOIDALGATERX Class 
%
% Summary of EllipsoidalGaterX:
% This is a class implementation of a multinomial resampler.
%
% EllipsoidalGaterX Properties:
%   + GateLevel - The desired number of standard deviations (away from the
%                 mean), which should be used as a threshold
%   + GatingProbability - The normalised percentage of the predicted measurement
%                    pdf, that should be incorporated in the validation
%                    region (i.e. th probability of gating)
%
% EllipsoidalGaterX Methods:
%   + EllipsoidalGaterX  - Constructor method
%   + gate - Perform gating and generate the validation matrix
% 
% See also SystematicResamplerX
    properties 
        GateLevel
        GatingProbability 
        Mapping
    end
    
    properties (SetAccess = immutable)
        NumMeasDims
        C = [2, pi, 4*pi/3, pi^2/2, 8*pi^2/15, pi^3/6, 16*pi^3/105, pi^4/24, 32*pi/945];
        % n-dim ellipsoid volumes up to 9-dims, as taken from
        % [http://oaji.net/articles/2014/1420-1415594291.pdf]
    end
    
    methods
        function this = EllipsoidalGaterX(varargin)
        % ELLIPSOIDALGATERX Constructor method
        % 
        % Parameters
        % ----------
        % NumMeasDims: scalar
        %   The number of observation dimensions, i.e. the dimensionality
        %   of the gate. For example, 1 = range, 2 = ellipse, 3 = ellipsoid
        % GateLevel: scalar, optional
        %   The desired number of standard deviations (away from the mean),
        %   which should be used as a gating threshold. If supplied, any
        %   value supplied for GatingProbability will be overriden and
        %   recomputed based on this value.
        % GatingProbability: normalised scalar, optional
        %   The normalised ([0,1]) percentage of the predicted measurement 
        %   pdf, that should be incorporated in the validation region 
        %   (i.e. the probability of gating)
        % Mapping: (NumMeasDimsx2) array specifying the mapping between
        %   state and measurement vectors, used to perform the gating.
        %   For example, we may have a state vector of the form
        %       x = [x_pos, x_vel, y_pos, y_vel]'
        %   and a measurement vector of the form
        %       y = [x_pos, x_vel, y_pos, y_vel]'
        %   where we want to perform gating based on x_pos,y_pos, then 
        %       mapping = [1,1;
        %                  3,3]
        %   i.e. y(1)->x(1) and y(2)->y(2)
        % Usage
        % -----
        % * eg = EllipsoidalGaterX(NumMeasDims, 'GateLevel', GateLevel) returns  
        %   an object configured to exclude any measurements falling further 
        %   that GateLevel standard deviations from a track's predicted measurement 
        %   mean. NumMeasDims specifies the dimensionality of the measurements.
        % * eg = EllipsoidalGaterX(NumMeasDims, 'GatingProbability',GatingProbability) returns an object
        %   configured to associate measurements which fall within the GatingProbability 
        %   percentile of the track's predicted measurement innovation covariance.
        %   NumMeasDims specifies the dimensionality of the measurements.
        % * eg = EllipsoidalGaterX('NumMeasDims', NumMeasDims, ___) is also valid
        %
        % NOTE: Either GateLevel or GatingProbability can be passed as input
        % arguments to the constructor, but NOT BOTH. This is done to
        % ensure 
        %
        % See also EllipsoidalGaterX/validate
            
            parser = inputParser;
            parser.addOptional('NumMeasDims',[]);
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            
            if(isfield(parser.Results,'NumMeasDims'))
                if(~isempty(parser.Results.NumMeasDims))
                    this.NumMeasDims = parser.Results.NumMeasDims;
                else
                    error('ELLIPSGATER:MISSINGNOBS','NumMeasDims must be a scalar (non-empty)');
                end
            else
                error('ELLIPSGATER:MISSINGNOBS','NumMeasDims must be a scalar (non-empty)');
            end
            
            config = parser.Unmatched;
            fields = fieldnames(config);
            if(numel(fields)>1)
                error('ELLIPSGATER:INVCONSTRARGS','Either GateLevel or GatingProbability can be passed as input arguments to the constructor, but NOT BOTH.');
            end
            fieldname = fields{1};
            if(strcmp(fieldname,'GateLevel'))
                this.GateLevel = config.GateLevel;
                this.GatingProbability = chi2cdf(this.GateLevel,this.NumMeasDims);
            elseif(strcmp(fieldname,'GatingProbability'))
                this.GatingProbability = config.GatingProbability;
                this.GateLevel = chi2inv(this.GatingProbability,this.NumMeasDims);
            end
            if(strcmp(fieldname,'Mapping'))
                this.Mapping = config.Mapping;
            else
                this.Mapping = [1:this.NumMeasDims;1:this.NumMeasDims]';
            end
        end
        
        function [ValidationMatrix, GateVolumes] = gate(this,varargin)
        % GATE Perform ellipsoidal gating to generate the validation
        % matrix indicating the valid associations between Tracks and
        % measurements
        %   
        % Parameters
        % ----------
        % TrackList: cell vector, optional
        %   A (1 x NumTracks) cell vector, where each cell contains a TrackX 
        %   object. It is required that each TrackX object should have a FilterX
        %   object as a property.
        % TrackPredMeasMeans: matrix, optional
        %   A (NumMeasDims x NumTracks) matrix, where each ith column represents
        %   the predicted measurement mean for the ith track. Can be used,
        %   in combination with TrackInnovErrCovars (see below), instead of
        %   a TrackList (see above).
        % TrackInnovErrCovars: matrix, optional
        %   A (NumMeasDims x NumMeasDims x NumTracks) matrix, where each
        %   (:,:,i) page represents the innovation covariance matrix for
        %   the ith track. Can be used, in combination with TrackInnovErrCovars
        %   (see above), instead of a TrackList (see above).
        % MeasurementList: matrix
        %   A (NumMeasDims x NumMeasurements) matrix, where each jth column
        %   corresponds to the jth measurement.
        %
        % Usage
        % -----
        % * [ValidationMatrix, GateVolumes] = gate(this,TrackList,MeasurementList) 
        %   returns the (NumTracks x NumMeasurements) validation matrix 
        %   ValidationMatrix indicating the valid associations between tracks  
        %   in TrackList (1 x NumTracks array of TrackX objects) and measurements 
        %   in MeasurmentList (NumMeasDims x NumMeasurements matrix), as well as a 
        %   (1 x NumTracks) vector GateVolumes, containing the volume of the computed
        %   gate for each track.
        % * [ValidationMatrix, GateVolumes] = ...
        %   gate(this,TrackPredMeasMeans,TrackInnovErrCovars,MeasurementList)
        %   returns the (NumTracks x NumMeasurements) validation matrix 
        %   ValidationMatrix indicating the valid associations between tracks  
        %   in with predicted measurement means given by the (NumMeasDims x NumTracks)
        %   TrackPredMeasMeans matrix, innovation covariances given by the 
        %   (NumMeasDims x NumMeasDims x NumTracks) TrackInnovErrCovars matrix, and  
        %   measurements in MeasurmentList (NumMeasDims x NumMeasurements matrix), as well as a 
        %   (1 x NumTracks) vector GateVolumes, containing the volume of the computed
        %   gate for each track.
        %
        % See also EllipsoidalGaterX/EllipsoidalGaterX
        
            if(nargin == 3)
                TrackList = varargin{1};
                MeasurementList = varargin{2};
                numTracks = numel(TrackList);
            else
                MeasurementPredictions = varargin{1};
                MeasurementList = varargin{2};
                numTracks = size(TrackPredMeasMeans,2);
            end
            
            numMeasurements = size(MeasurementList,2);
             
            % Pre-allocate memory
            GateVolumes = zeros(1,numTracks); 
            ValidationMatrix = zeros(numTracks,numMeasurements);    
            
            % Construct Validation and GateVolumes Matrices
            for trackInd = 1:numTracks
                
                if(nargin==3)
                    % Store pointer to track
                    Track = TrackList{trackInd};

                    % Extract predicted measurement and innovation covariance from filter
                    PredMeasMean = Track.Filter.MeasurementPrediction.Mean(this.Mapping(:,2));
                    InnovErrCovar = Track.Filter.MeasurementPrediction.Covar(this.Mapping(:,2),this.Mapping(:,2));
                else
                    % Extract predicted measurement and innovation covariance from filter
%                     PredMeasMean = TrackPredMeasMeans(this.Mapping(:,2),trackInd);
%                     InnovErrCovar = TrackInnovErrCovars(this.Mapping(:,2),this.Mapping(:,2),trackInd);
                end
                
                % Perform Gating
                c = this.C(this.NumMeasDims);
                GateVolumes(trackInd) = c*this.GateLevel^(this.NumMeasDims/2)*det(InnovErrCovar)^(1/2);    
                if(numel(MeasurementList)>0)
                    ValidationMatrix(trackInd,:) = this.mahalDist(MeasurementList(this.Mapping(:,2),:), PredMeasMean, InnovErrCovar, 2) < this.GateLevel;
                end
            end
            this.ValidationMatrix = ValidationMatrix;
            this.GateVolumes = GateVolumes;
        end
    end
    
    methods (Static)
        function D=mahalDist(x, m, C, use_log)
        % p=gaussian_prob(x, m, C, use_log)
        %
        % Evaluate the multi-variate density with mean vector m and covariance
        % matrix C for the input vector x.
        % Vectorized version: Here X is a matrix of column vectors, and p is 
        % a vector of probabilities for each vector.

            if nargin<4, use_log = 0; end

            d   = length(m);

            if size(x,1)~=d
               x=x';
            end
            N       = size(x,2);

            m       = m(:);
            M       = m*ones(1,N);
            denom   = (2*pi)^(d/2)*sqrt(abs(det(C)));
            %invC    = inv(C);
            mahal   = sum(((x-M)'/C).*(x-M)',2);   % Chris Bregler's trick

            switch use_log
            case 2
              D     = mahal;
            case 1
              D     = -0.5*mahal - log(denom);
            case 0
              numer = exp(-0.5*mahal);
              D     = numer/denom;
            otherwise
                error('Unsupported log type')
            end
        end
    end
end

