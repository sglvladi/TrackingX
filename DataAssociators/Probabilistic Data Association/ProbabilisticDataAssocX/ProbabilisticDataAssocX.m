classdef ProbabilisticDataAssocX < DataAssociatorX
% ProbabilisticDataAssocX class
%
% Summary of ProbabilisticDataAssocX:
% This is a class implementation of a Joint Probabilistic Data Association Filter.
%
% ProbabilisticDataAssocX Properties:
%   + TrackList - A (1 x Tracks) vector of TrackX objects
%   + MeasurementList - A MeasurementListX object, containing a (1 xTracks) 
%       array of MeasurementX objects
%   + Gater - A GaterX object used to perform gating
%   + Clusterer - A ClustererX object used to perform clustering of tracks
%   + ClutterModel - A ClutterModelX subclass used to model clutter
%
% ProbabilisticDataAssocX Methods:
%   + ProbabilisticDataAssocX  - Constructor method
%   + associate - Performs JPDAF association step
%   + updateTracks - Performs JPDAF update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% See also ParticleFilterX, KalmanFilerX.

    properties (Dependent)
        NumTracks
        NumMeasurements
    end
    
    properties
        ClutterModel
        DetectionModel
        
        ValidationMatrix
        GateVolumes
        LikelihoodMatrix
        AssocLikelihoodMatrix
        AssocWeightsMatrix
        ClusterList = []
        UnassocTrackInds = []
    end
    
    methods (Access=protected)
        function initialise_(this, config)
            if (isfield(config,'Gater'))
                this.Gater  = config.Gater;
            end
            if (isfield(config,'Clusterer'))
                this.Clusterer  = config.Clusterer;
            end
            if (isfield(config,'ClutterModel'))
                this.ClutterModel  = config.ClutterModel;
            end
            if (isfield(config,'DetectionModel'))
                this.DetectionModel  = config.DetectionModel;
            end
            if (isfield(config,'Hypothesiser'))
                this.Hypothesiser  = config.Hypothesiser;
            end
        end
    end
    
    methods
        function this = ProbabilisticDataAssocX(varargin)
        % ProbabilisticDataAssocX - Constructor method
        % 
        % Parameters
        % ----------
        % Gater: GaterX subclass, optional
        %   A gater object which should be used to perform gating of measurements. 
        %   (default = None, meaning that no gating is performed)
        % Clusterer: ClustererX subclass, optional
        %   A clusterer object which should be used to perform clustering of tracks. 
        %   (default = None, meaning that no clustering is performed)
        % Hypothesiser: HypothesiserX subclass, optional
        %   A hypothesiser object used generate and evaluate association
        %   hypotheses. (default = EfficientHypothesisManagementX());
        % DetectionModel: scalar, optional
        %   The target detection model
        % ClutterModel: ClutterModelX subclass
        %   A clutter model
        %
        % Usage
        % -----
        % * jpda = ProbabilisticDataAssocX(___,Name,Value) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        %  See also ProbabilisticDataAssocX/associate, ProbabilisticDataAssocX/updateTracks.   
                    
            this.Hypothesiser = EfficientHypothesisManagementX();
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise PDA with a certain set of parameters. 
        %   
        % DESCRIPTION: 
        % * initialise(jpda,gater,clusterer) initialises the ProbabilisticDataAssocX
        %   object jpda with the provided GaterX and ClustererX object handles gater and clusterer.
        % * initialise(jpda,___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        % INPUT ARGUMENTS:
        % * Gater               (GaterX) A gater object which should be used to  
        %                       perform gating of measurements. Default = None
        % * Clusterer           (ClustererX) A clusterer object which should be used to  
        %                       perform clustering of tracks. Default = None
        %
        %  See also ProbabilisticDataAssocX/associate, ProbabilisticDataAssocX/updateTracks.   
                    
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
        end
        
        function associate(this,TrackList,MeasurementList)
        % associate - Performs JPDAF association step
        %
        % Parameters
        % ----------
        % TrackList: TrackListX container, optional
        %   Array of TrackX objects to be associated. If skipped the function 
        %   will attempt to use an internally cached version of the
        %   variable.
        % MeasurementList: MeasurementListX container
        %   Array of MeasurementX objects. If skipped the function will attempt
        %   to use an internally cached version of the variable.
        %
        % See also ProbabilisticDataAssocX/initialise, ProbabilisticDataAssocX/updateTracks.
        
            if(nargin>1)
                this.TrackList = TrackList;
                this.MeasurementList = MeasurementList;
            end
            
            if(this.NumTracks>0)
                associate@DataAssociatorX(this);    
            else
                this.ValidationMatrix = zeros(1, this.NumMeasurements);
                this.ClusterList = [];
                this.UnassocTrackInds = [];
                this.AssocWeightsMatrix = -1; % Set betta to -1
            end
        end
        
        function predictTracks(this)
        % predictTracks - Performs JPDAF track prediction step
        %   
        % DESCRIPTION: 
        % * predictTracks(jpda) performs JPDAF update on all tracks contained
        %   in jpda.TrackList.
        %
        %   See also ProbabilisticDataAssocX/initialise, ProbabilisticDataAssocX/updateTracks.
            
            numTracks = this.NumTracks; 
            if(numTracks>0)
                % Compute weights and update each track
                for trackInd=1:numTracks
                    this.TrackList{trackInd}.Filter.predict();
                end    
            end
        end
        
        function updateTracks(this)
        % updateTracks - Performs JPDAF track update step
        %   
        % DESCRIPTION: 
        % * updateTracks(jpda) performs JPDAF update on all tracks contained
        %   in jpda.TrackList.
        %
        %   See also ProbabilisticDataAssocX/initialise, ProbabilisticDataAssocX/updateTracks.
            
            numTracks = this.NumTracks;
            if(numTracks>0)
                % Compute weights and update each track
                for trackInd=1:numTracks
                    ValidDataInd = find(this.ValidationMatrix(trackInd,:));    % Associated measurements
                    assocWeights = this.AssocWeightsMatrix(trackInd,[1 ValidDataInd+1]);
                    this.TrackList{trackInd}.Filter.updatePDA(assocWeights);
                    this.TrackList{trackInd}.Trajectory(end+1) = this.TrackList{trackInd}.Filter.StatePosterior;
                end    
            end
        end
        
        function numTracks = get.NumTracks(this)
            numTracks = numel(this.TrackList);
        end
        function numMeasurements = get.NumMeasurements(this)
            numMeasurements = this.MeasurementList.NumMeasurements;
        end
    end
    
    methods (Access = protected)
        
        function performGating(this)
            
            numTracks = size(this.TrackList,2);
            numMeasurements = this.MeasurementList.NumMeasurements;
            
            if(isa(this.Gater,'GaterX'))
                % Validation matix and volume
                [this.ValidationMatrix, this.GateVolumes] = this.Gater.gate(this.TrackList,this.MeasurementList.Vectors);
            else
                this.ValidationMatrix = ones(size(this.TrackList,2),size(this.MeasurementList.Vectors,2));
                this.GateVolumes = [];
            end
            
            for trackInd = 1:this.NumTracks
                if(~isprop(this.TrackList{trackInd},'GateVolume'))
                    this.TrackList{trackInd}.addprop('GateVolume');
                end
                this.TrackList{trackInd}.GateVolume = this.GateVolumes(trackInd);
                if(~isprop(this.TrackList{trackInd},'ValidationMatrix'))
                    this.TrackList{trackInd}.addprop('ValidationMatrix');
                end
                this.TrackList{trackInd}.ValidationMatrix = this.ValidationMatrix(trackInd,:);
            end
        end
        
        function computeLikelihoods(this)
            
            this.LikelihoodMatrix = zeros(this.NumTracks, this.NumMeasurements);
            
            % Compute Likelihood matrix
            for trackInd = 1:this.NumTracks          
                
                ValidMeasInds = find(this.TrackList{trackInd}.ValidationMatrix);
                if(~isempty(ValidMeasInds))
                    % Extract valid measurements
                    this.TrackList{trackInd}.Filter.MeasurementList = MeasurementListX(this.MeasurementList.Vectors(:,ValidMeasInds),...
                                                                                           this.MeasurementList.Timestamp);
                    % Evaluate the measurement likelihoods
                    this.LikelihoodMatrix(trackInd, ValidMeasInds) = this.TrackList{trackInd}.Filter.MeasurementLikelihoods;
                else
                    % Extract valid measurements
                    this.TrackList{trackInd}.Filter.MeasurementList = MeasurementListX();
                end
            end
        end
        
        function performClustering(this)
            % Do nothing. PDAF does not perform clustering
        end
        
        function evaluateAssociations(this)
            
            numMeasurements = this.NumMeasurements;
            numTracks = this.NumTracks;
            
            % Allocate memory for association weights 
            % (Dummy measurement at index 1)
            this.AssocLikelihoodMatrix = zeros(numTracks,numMeasurements+1);
            this.AssocWeightsMatrix = zeros(numTracks,numMeasurements+1);
            
            % Iterate over the tracks
            for trackInd=1:numTracks
                
                % Extract valid measurements
                measIndList = find(this.ValidationMatrix(trackInd,:));
                
                if(isempty(measIndList))
                % If no measurements associated to track
                    this.AssocLikelihoodMatrix(trackInd,:) = [1 zeros(1,numMeasurements)];
                    this.AssocWeightsMatrix(trackInd,:) = [1 zeros(1,numMeasurements)];
                    continue;
                else
                    
                    % Compute gating probability
                    Pg = 1;
                    if(isa(this.Gater,'EllipsoidalGaterX'))
                        Pg = this.Gater.GatingProbability;
                    end
                    
                    % Compute detection probability
                    Pd = 1;
                    if(isa(this.DetectionModel, 'DetectionModelX'))
                        Pd = this.DetectionModel.pdf(this.TrackList{trackInd}.Filter.StatePrediction.Mean);
                    end
                        
                    % Compute clutter density per unit volume
                    if(isempty(this.ClutterModel))
                        lambda = numel(measIndList)/this.GateVolumes(trackInd);
                    else
                        lambda = this.ClutterModel.pdf(this.TrackList{trackInd}.Filter.MeasurementPrediction.Mean);
                    end
                    lambda(lambda==0) = eps; % Ensure lambda is non-zero
                    
                    % Construct standard association likelihood matrix
                    this.AssocLikelihoodMatrix(trackInd,:) = ...
                        [lambda*(1-Pd*Pg), Pd*Pg*this.LikelihoodMatrix(trackInd,:)];
                    
                    % Hypothesise/Compute Association probabilities
                    this.AssocWeightsMatrix(trackInd,:) = this.Hypothesiser.hypothesise(this.AssocLikelihoodMatrix(trackInd,:));
                end
            end          
        end
    end
end