classdef JointProbabilisticDataAssocX < ProbabilisticDataAssocX
% JointProbabilisticDataAssocX class
%
% Summary of JointProbabilisticDataAssocX:
% This is a class implementation of a Joint Probabilistic Data Association Filter.
%
% JointProbabilisticDataAssocX Properties:
%   + TrackList        A (1-by-NumTracks) vector of TrackX objects
%   + MeasurementList           A (1-by-NumMeas) vector of observations/measurements
%   + Gater             A GaterX object used to perform gating
%   + Clusterer         A ClustererX object used to perform clustering of
%                       tracks
%   + ClutterDensity    
%
% JointProbabilisticDataAssocX Methods:
%   + JointProbabilisticDataAssocX  - Constructor method
%   + associate                     - Performs JPDAF association step
%   + updateTracks                 - Performs JPDAF update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% See also ParticleFilterX, KalmanFilerX.

    properties
    end
    
    methods
        function this = JointProbabilisticDataAssocX(varargin)
        % JointProbabilisticDataAssocX - Constructor method
        %   
        % DESCRIPTION: 
        % * jpda = JointProbabilisticDataAssocX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * jpda = JointProbabilisticDataAssocX(gater,clusterer) returns an 
        %   object handle, preconfigured with the provided GaterX and ClustererX 
        %   object handles gater and clusterer.
        % * jpda = JointProbabilisticDataAssocX(___,Name,Value,___) instantiates an  
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
        %  See also JointProbabilisticDataAssocX/associate, JointProbabilisticDataAssocX/updateTracks.   
                    
            this@ProbabilisticDataAssocX(varargin{:});
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise JPDA with a certain set of parameters. 
        %   
        % DESCRIPTION: 
        % * initialise(jpda,gater,clusterer) initialises the JointProbabilisticDataAssocX
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
        %  See also JointProbabilisticDataAssocX/associate, JointProbabilisticDataAssocX/updateTracks.   
                    
            initialise@ProbabilisticDataAssocX(this,varargin{:});
        end
        
        function associate(this,TrackList,MeasurementList)
        % associate - Performs JPDAF association step
        %   
        % DESCRIPTION: 
        % * associate(jpda) performs data association on the object jpda, based 
        %   on its  internally stored jpda.TrackList and jpda.MeasurementList properties
        % * associate(jpda,TrackList,MeasurementList) performs data association on 
        %   the object jpda, based on the provided list of tracks TrackList 
        %   and list of observations MeasurementList  
        %   
        %   Usage:
        %       (jpdaf.Params.k = 1; % 1 sec)
        %       jpdaf.Predict();
        %
        %   See also JointProbabilisticDataAssocX/initialise, JointProbabilisticDataAssocX/updateTracks.
        
            % Call super class method
            switch(nargin)
                case 1
                    associate@ProbabilisticDataAssocX(this);
                case 2
                    associate@ProbabilisticDataAssocX(this,TrackList);
                case 3
                    associate@ProbabilisticDataAssocX(this,TrackList,MeasurementList);
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
                Pg = 1;
                if(isa(this.Gater,'EllipsoidalGaterX'))
                    Pg = this.Gater.GatingProbability;
                end

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
                        Pd = 1;
                        if(isa(this.DetectionModel, 'DetectionModelX'))
                            Pd = this.DetectionModel.pdf(this.TrackList{trackInd}.Filter.StatePrediction.Mean);
                        end

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