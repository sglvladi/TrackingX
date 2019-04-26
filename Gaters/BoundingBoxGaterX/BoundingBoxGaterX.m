classdef BoundingBoxGaterX <GaterX
% BOUNDINGBOXGATERX Class 
%
% Summary of BoundingBoxGaterX:
% This is a class implementation of a bounding box gater, which gates
% measurements, given in the form of bounding boxes (x,y,w,h), depending on
% their proportional overlap with the estimated target bounding boxes.
%
% BoundingBoxGaterX Properties:
%   + OverlapThresh - The minimum overlap, proportional to the each
%                     target's bounding box.
%
% BoundingBoxGaterX Methods:
%   + BoundingBoxGaterX  - Constructor method
%   + gate - Perform gating and generate the validation matrix
% 
% See also SystematicResamplerX
    properties 
        OverlapThresh 
    end
    
    properties (SetAccess = immutable)
        NumObsDims
    end
    
    methods
        function this = BoundingBoxGaterX(varargin)
        % BoundingBoxGaterX Constructor method
        % 
        % Parameters
        % ----------
        % OverlapThresh: normalised scalar, optional
        %   The normalised ([0,1]) minimum percentage overlap between the estimated 
        %   target bounding box and any given measurement, such that a measurement
        %   should be gated.
        % 
        % Usage
        % -----
        % * eg = BoundingBoxGaterX(OverlapThresh) returns an object configured 
        %   to exclude any measurements that cover less than OverlapThresh % 
        %   of a given target state.
        % * eg = BoundingBoxGaterX('OverlapThresh', OverlapThresh) is also valid
        %
        % See also EllipsoidalGaterX/validate
            
            parser = inputParser;
            parser.addOptional('OverlapThresh',[]);
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            this.OverlapThresh = parser.Results.OverlapThresh;
          
        end
        
        function [ValidationMatrix, GateVolumes] = gate(this,varargin)
        % GATE Perform bbox gating to generate the validation
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
        %   A (NumObsDims x NumTracks) matrix, where each ith column represents
        %   the predicted measurement mean for the ith track. Can be used,
        %   in combination with TrackInnovErrCovars (see below), instead of
        %   a TrackList (see above).
        % TrackInnovErrCovars: matrix, optional
        %   A (NumObsDims x NumObsDims x NumTracks) matrix, where each
        %   (:,:,i) page represents the innovation covariance matrix for
        %   the ith track. Can be used, in combination with TrackInnovErrCovars
        %   (see above), instead of a TrackList (see above).
        % MeasurementList: matrix
        %   A (NumObsDims x NumMeasurements) matrix, where each jth column
        %   corresponds to the jth measurement.
        %
        % Usage
        % -----
        % * [ValidationMatrix] = gate(this,TrackList,MeasurementList) 
        %   returns the (NumTracks x NumMeasurements) validation matrix 
        %   ValidationMatrix indicating the valid associations between tracks  
        %   in TrackList (1 x NumTracks array of TrackX objects) and measurements 
        %   in MeasurmentList (NumObsDims x NumMeasurements matrix)
        % * [ValidationMatrix] = gate(this,TrackPredMeasMeans,TrackInnovErrCovars,MeasurementList)
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
                MeasurementList = varargin{3};
                NumTracks = size(TrackPredMeasMeans,2);
            end
            
            [~, NumMeasurements] = size(MeasurementList);
             
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
                else
                    % Extract predicted measurement and innovation covariance from filter
                    PredMeasMean = TrackPredMeasMeans(:,trackInd);
                end
                if(NumMeasurements>0)
                    intersect_1 = rectint(PredMeasMean', MeasurementList');
                    intersect_1 = intersect_1./(PredMeasMean(3)*PredMeasMean(4));
                    intersect_2 = rectint(MeasurementList', PredMeasMean')';
                    intersect_2 = intersect_2./(MeasurementList(3,:).*MeasurementList(4,:));
                    ValidationMatrix(trackInd,:) = (intersect_1>this.OverlapThresh | intersect_2>this.OverlapThresh);
                end
                GateVolumes(trackInd) = PredMeasMean(3)*PredMeasMean(4);
            end
            this.ValidationMatrix = ValidationMatrix;
            this.GateVolumes = GateVolumes;
        end
    end
end

