classdef EllipsoidalGaterX <GaterX
% ELLIPSOIDALGATERX Class 
%
% Summary of EllipsoidalGaterX:
% This is a class implementation of a multinomial resampler.
%
% EllipsoidalGaterX Properties:
%   None
%
% EllipsoidalGaterX Methods:
%    EllipsoidalGaterX  - Constructor method
%    gate - Perform gating and generate the validation matrix
% 
% See also SystematicResamplerX
    properties 
        GateLevel
        ProbOfGating 
    end
    
    properties (SetAccess = immutable)
        NumObsDims
    end
    
    methods
        function this = EllipsoidalGaterX(varargin)
        % ELLIPSOIDALGATERX Constructor method
        %   
        % DESCRIPTION: 
        % * eg = EllipsoidalGaterX(NumObsDims, 'GateLevel', GateLevel) returns  
        %   an object configured to exclude any measurements falling further 
        %   that GateLevel standard deviations from a track's predicted measurement 
        %   mean. NumObsDims specifies the dimensionality of the measurements.
        % * eg = EllipsoidalGaterX(NumObsDims, 'ProbOfGating',ProbOfGating) returns an object
        %   configured to associate measurements which fall within the ProbOfGating 
        %   percentile of the track's predicted measurement innovation covariance.
        %   NumObsDims specifies the dimensionality of the measurements.
        % * eg = EllipsoidalGaterX('NumObsDims', NumObsDims, ___) is also valid
        %
        % NOTE: Either GateLevel or ProbOfGating can be passed as input
        % arguments to the constructor, but NOT BOTH. This is done to
        % ensure 
        %
        % See also EllipsoidalGaterX/validate
            
            parser = inputParser;
            parser.addOptional('NumObsDims',[]);
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            
            if(isfield(parser.Results,'NumObsDims'))
                if(~isempty(parser.Results.NumObsDims))
                    this.NumObsDims = parser.Results.NumObsDims;
                else
                    error('ELLIPSGATER:MISSINGNOBS','NumObsDims must be a scalar (non-empty)');
                end
            else
                error('ELLIPSGATER:MISSINGNOBS','NumObsDims must be a scalar (non-empty)');
            end
            
            config = parser.Unmatched;
            fields = fieldnames(config);
            if(numel(fields)>1)
                error('ELLIPSGATER:INVCONSTRARGS','Either GateLevel or ProbOfGating can be passed as input arguments to the constructor, but NOT BOTH.');
            end
            fieldname = fields{1};
            if(strcmp(fieldname,'GateLevel'))
                this.GateLevel = config.GateLevel;
                this.ProbOfGating = chi2cdf(this.GateLevel,this.NumObsDims);
            elseif(strcmp(fieldname,'ProbOfGating'))
                this.ProbOfGating = config.ProbOfGating;
                this.GateLevel = chi2inv(this.ProbOfGating,this.NumObsDims);
            end
        end
        
        function [ValidationMatrix, GateVolumes] = gate(this,varargin)
        % GATE Perform ellipsoidal gating to generate the validation
        % matrix indicating the valid associations between Tracks and
        % measurements
        %   
        % DESCRIPTION: 
        % * [ValidationMatrix, GateVolumes] = gate(this,TrackList,MeasurementList) 
        %   returns the (NumTracks x NumMeasurements) validation matrix 
        %   ValidationMatrix indicating the valid associations between tracks  
        %   in TrackList (1 x NumTracks array of TrackX objects) and measurements 
        %   in MeasurmentList (NumObsDims x NumMeasurements matrix), as well as a 
        %   (1 x NumTracks) vector GateVolumes, containing the volume of the computed
        %   gate for each track.
        % * [ValidationMatrix, GateVolumes] = ...
        %   gate(this,TrackPredMeasMeans,TrackInnovErrCovars,MeasurementList)
        %   returns the (NumTracks x NumMeasurements) validation matrix 
        %   ValidationMatrix indicating the valid associations between tracks  
        %   in with predicted measurement means given by the (NumObsDims x NumTracks)
        %   TrackPredMeasMeans matrix, innovation covariances given by the 
        %   (NumObsDims x NumObsDims x NumTracks) TrackInnovErrCovars matrix, and  
        %   measurements in MeasurmentList (NumObsDims x NumMeasurements matrix), as well as a 
        %   (1 x NumTracks) vector GateVolumes, containing the volume of the computed
        %   gate for each track.
        %
        % See also EllipsoidalGaterX/EllipsoidalGaterX
        
            if(nargin == 3)
                TrackList = varargin{1};
                MeasurementList = varargin{2};
                NumTracks = numel(TrackList);
            else
                TrackPredMeasMeans = varargin{1};
                TrackInnovErrCovars = varargin{2};
                MeasurementList = varargin{3};
                NumTracks = size(TrackPredMeasMeans,2);
            end
            
            [NumObsDims, NumMeasurements] = size(MeasurementList);

            if(NumObsDims ~= this.NumObsDims)
                error('ELLIPSGATER:INVMEASDIMS','The number of rows in MeasurementList must be equal to this.NumObsDims');
            end
             
            % Pre-allocate memory
            GateVolumes = zeros(1,NumTracks); 
            ValidationMatrix = zeros(NumTracks,NumMeasurements);    
            
            % Construct Validation and GateVolumes Matrices
            for trackInd = 1:NumTracks
                
                if(nargin==3)
                    % Store pointer to track
                    Track = TrackList{trackInd};

                    % Extract predicted measurement and innovation covariance from filter
                    PredMeasMean = Track.Filter.PredMeasMean;
                    InnovErrCovar = Track.Filter.InnovErrCovar;
                else
                    % Extract predicted measurement and innovation covariance from filter
                    PredMeasMean = TrackPredMeasMeans(:,trackInd);
                    InnovErrCovar = TrackInnovErrCovars(:,:,trackInd);
                end
                
                % Perform Gating
                switch numel(PredMeasMean)
                    case 1 
                        C = 2;
                    case 2
                        C = pi;
                    case 3
                        C = 4*pi/3;
                    otherwise
                        error('ELLIPSGATER:INVMEASDIMS','Ellipsoidal Gating has only been implemented for observations of up to 3 dimensions!');
                end
                GateVolumes(trackInd) = C*this.GateLevel^(NumObsDims/2)*det(InnovErrCovar)^(1/2);    
                ValidationMatrix(trackInd,:) = this.mahalDist(MeasurementList, PredMeasMean, InnovErrCovar, 2) < this.GateLevel;
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

