classdef GlobalNearestNeighbourDataAssocX < NearestNeighbourDataAssocX
% GlobalNearestNeighbourDataAssocX class
%
% Summary of GlobalNearestNeighbourDataAssocX:
% This is a class implementation of a Global Nearest Neighbour Data Association Filter.
%
% GlobalNearestNeighbourDataAssocX Properties:
%   + TrackList        A (1-by-NumTracks) vector of TrackX objects
%   + MeasurementList           A (1-by-NumMeas) vector of observations/measurements
%   + Gater             A GaterX object used to perform gating
%   + Clusterer         A ClustererX object used to perform clustering of
%                       tracks
%   + ClutterDensity    
%
% GlobalNearestNeighbourDataAssocX Methods:
%   + GlobalNearestNeighbourDataAssocX  - Constructor method
%   + associate                     - Performs JPDAF association step
%   + updateTracks                 - Performs JPDAF update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% See also ParticleFilterX, KalmanFilerX.
    
    methods
        function this = GlobalNearestNeighbourDataAssocX(varargin)
        % GlobalNearestNeighbourDataAssocX - Constructor method
        %   
        % DESCRIPTION: 
        % * gnn = GlobalNearestNeighbourDataAssocX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * gnn = GlobalNearestNeighbourDataAssocX(gater,clusterer) returns an 
        %   object handle, preconfigured with the provided GaterX and ClustererX 
        %   object handles gater and clusterer.
        % * gnn = GlobalNearestNeighbourDataAssocX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        % INPUT ARGUMENTS:
        % * Gater               (GaterX) A gater object which should be used to  
        %                       perform gating of measurements. Default = None
        % * Clusterer           (ClustererX) A clusterer object which should be used to  
        %                       perform clustering of tracks. Default = None
        % * DetectionProbability        (scalar) The target detection probability
        % * ClutterDensity      (scalar) The spatial density of clutter
        %
        %  See also GlobalNearestNeighbourDataAssocX/associate, GlobalNearestNeighbourDataAssocX/updateTracks.   
                    
            this@NearestNeighbourDataAssocX(varargin{:});
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise gnn with a certain set of parameters. 
        %   
        % DESCRIPTION: 
        % * initialise(gnn,gater,clusterer) initialises the GlobalNearestNeighbourDataAssocX
        %   object gnn with the provided GaterX and ClustererX object handles gater and clusterer.
        % * initialise(gnn,___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        % INPUT ARGUMENTS:
        % * Gater               (GaterX) A gater object which should be used to  
        %                       perform gating of measurements. Default = None
        % * Clusterer           (ClustererX) A clusterer object which should be used to  
        %                       perform clustering of tracks. Default = None
        %
        %  See also GlobalNearestNeighbourDataAssocX/associate, GlobalNearestNeighbourDataAssocX/updateTracks.   
                    
            initialise@NearestNeighbourDataAssocX(this,varargin{:});
        end
        
        function associate(this,TrackList,MeasurementList)
        % associate - Performs JPDAF association step
        %   
        % DESCRIPTION: 
        % * associate(gnn) performs data association on the object gnn, based 
        %   on its  internally stored gnn.TrackList and gnn.MeasurementList properties
        % * associate(gnn,TrackList,MeasurementList) performs data association on 
        %   the object gnn, based on the provided list of tracks TrackList 
        %   and list of observations MeasurementList  
        %   
        %   Usage:
        %       (gnn.Params.k = 1; % 1 sec)
        %       gnn.Predict();
        %
        %   See also GlobalNearestNeighbourDataAssocX/initialise, GlobalNearestNeighbourDataAssocX/updateTracks.
        
            % Call super class method
            switch(nargin)
                case 1
                    associate@NearestNeighbourDataAssocX(this);
                case 2
                    associate@NearestNeighbourDataAssocX(this,TrackList);
                case 3
                    associate@NearestNeighbourDataAssocX(this,TrackList,MeasurementList);
            end
        end   
    end
    
    methods (Access = protected)
        
        function performClustering(this)
            % Use the defined this.Clusterer to perform clustering
            if(~isempty(this.Clusterer))
                % Get all clusters
                this.ClusterList = this.Clusterer.cluster(this.ValidationMatrix);
            else
                this.ClusterList(1).TrackIndList = 1:this.NumTracks;
                this.ClusterList(1).MeasIndList = find(sum(this.ValidationMatrix,1));
            end
        end
        
        function evaluateAssociations(this)
            
            numTracks = this.NumTracks;
            numMeasurements = this.NumMeasurements;
            
            % Allocate memory for association likelihoods & weights 
            % (Dummy measurement at index 1)
            this.AssocLikelihoodMatrix = zeros(numTracks,numMeasurements+1);
            this.AssocWeightsMatrix = zeros(numTracks,numMeasurements+1);
                        
            % Loop through all clusters
            NumClusters = numel(this.ClusterList);
            for clusterInd=1:NumClusters
                % Create Hypothesis net for each cluster and 
                % populate association weights matrix

                % Extract track and measurement list for cluster
                MeasIndList = this.ClusterList(clusterInd).MeasIndList;
                TrackIndList = this.ClusterList(clusterInd).TrackIndList;
                numClusterMeasurements = length(MeasIndList);
                numClusterTracks = length(TrackIndList);
                
                % Compute gating probability
                Pg = this.Gater.GatingProbability;

                if(numClusterMeasurements==0)
                % If no measurements associated to cluster
                    % Simply set all cluster tracks to be missed
                    this.AssocLikelihoodMatrix(TrackIndList,:) = [1 zeros(1,numMeasurements)];
                    this.AssocWeightsMatrix(TrackIndList,:) = [1 zeros(1,numMeasurements)];
                    continue;
                else
                    
                    % Pre-allocate memory
                    this.ClusterList(clusterInd).AssocLikelihoodMatrix = ...
                        zeros(numClusterTracks, numClusterMeasurements+1);
                    this.ClusterList(clusterInd).AssocWeightsMatrix = ...
                        zeros(numClusterTracks, numClusterMeasurements+1);

                    % Compute New Track/False Alarm density for the cluster
                    for t = 1:numClusterTracks
                        trackInd = TrackIndList(t);

                        % Compute detection probability
                        Pd = this.DetectionModel.pdf(this.TrackList{trackInd}.Filter.StatePrediction.Mean);

                        % Compute clutter density per unit volume
                        if(isempty(this.ClutterModel))
                            lambda = numClusterMeasurements/sum(this.GateVolumes(TrackIndList));
                        else
                            lambda = this.ClutterModel.pdf(this.TrackList{trackInd}.Filter.MeasurementPrediction.Mean);
                        end
                        lambda(lambda==0) = eps; % Ensure lambda is non-zero

                        % Construct cluster association likelihood matrix
                        this.ClusterList(clusterInd).AssocLikelihoodMatrix(t,:) = ...
                            [lambda*(1-Pd*Pg), Pd*Pg*this.LikelihoodMatrix(trackInd,MeasIndList)];
                    end

                    % Compute likelihood matrix
                    this.AssocLikelihoodMatrix(TrackIndList,[1 MeasIndList+1]) = ...
                        this.ClusterList(clusterInd).AssocLikelihoodMatrix;

                    % Compute association weights
                    this.ClusterList(clusterInd).AssocWeightsMatrix =...
                        this.Hypothesiser.hypothesise(this.ClusterList(clusterInd).AssocLikelihoodMatrix);

                    % Store cluster association weights in global
                    % associatiopn weights matrix
                    this.AssocWeightsMatrix(TrackIndList, [1 MeasIndList+1]) = ...
                        this.ClusterList(clusterInd).AssocWeightsMatrix; %./sum(this.ClusterList(clusterInd).AssocWeightsMatrix,2);
                end 
            end
        end
    end
end