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
        
        function predictTracks(this)
        % predictTracks - Performs JPDAF track prediction step
        %   
        % DESCRIPTION: 
        % * predictTracks(jpda) performs JPDAF update on all tracks contained
        %   in jpda.TrackList.
        %
        %   See also JointProbabilisticDataAssocX/initialise, JointProbabilisticDataAssocX/updateTracks.
        
            % Call super class method
            predictTracks@ProbabilisticDataAssocX(this);
        end
        
        function updateTracks(this)
        % updateTracks - Performs JPDAF track update step
        %   
        % DESCRIPTION: 
        % * updateTracks(jpda) performs JPDAF update on all tracks contained
        %   in jpda.TrackList.
        %
        %   See also JointProbabilisticDataAssocX/initialise, JointProbabilisticDataAssocX/updateTracks.
        
            % Call super class method
            updateTracks@ProbabilisticDataAssocX(this);
        end
        
    end
    
    methods (Access = protected)
        
        function performGating(this)
            % Call super class method
            performGating@ProbabilisticDataAssocX(this);
        end
        
        function computeLikelihoods(this)
            % Call super class method
            computeLikelihoods@ProbabilisticDataAssocX(this);
        end
        
        function performClustering(this)
            % Use the defined this.Clusterer to perform clustering
            if(isempty(this.Clusterer))
                % Do nothing
                return;
            else
                % Get all clusters
                [this.ClusterList, this.UnassocTrackInds] = this.Clusterer.cluster(this.ValidationMatrix);
            end
        end
        
        function evaluateAssociations(this)
            
            numTracks = this.NumTracks;
            numMeasurements = this.NumMeasurements;
            
            this.AssocLikelihoodMatrix = [ones(numTracks,1), zeros(numTracks,numMeasurements)];
            this.AssocWeightsMatrix = [ones(numTracks,1), zeros(numTracks,numMeasurements)];
            
            if(sum(sum(this.ValidationMatrix))==0)
                return;
            else
                if(isempty(this.Clusterer))
                                        
                    if(isa(this.Gater,'EllipsoidalGaterX'))
                        GatingProbability = this.Gater.GatingProbability;
                    else
                        GatingProbability = 1;
                    end
                    
                    for trackInd = 1:numTracks
                        % Compute New Track/False Alarm density for the cluster
                        if(isempty(this.ClutterModel))
                            if(isempty(this.GateVolumes))
                                error('Cannot perform non-parametric (J)PDA as no gater has been set and thus gate volumes cannot be computed.. TODO: Allow for search space volume (V) to be set in the future');
                            end
                            clutterDensity = GatingProbability*numMeasurements/sum(this.GateVolumes);
                        else
                            clutterDensity = this.ClutterModel.pdf(this.TrackList{trackInd}.Filter.MeasurementPrediction.Mean);
                        end
                        if(clutterDensity==0)
                            clutterDensity = eps;
                        end
                       this.AssocLikelihoodMatrix(trackInd,1) = clutterDensity*(1-this.DetectionProbability*GatingProbability);
                    end
                    
                    this.AssocLikelihoodMatrix(trackInd,2:end) = this.DetectionProbability*GatingProbability*this.LikelihoodMatrix;
                    
                    this.AssocWeightsMatrix = this.Hypothesiser.hypothesise(this.AssocLikelihoodMatrix);                                                           
                    
                else
                    % Get all clusters
                    [this.ClusterList, this.UnassocTrackInds] = this.Clusterer.cluster(this.ValidationMatrix);

                    % Allocate memory for association weights and fill in weights for unassociated tracks
                    this.AssocWeightsMatrix = zeros(numTracks, numMeasurements+1); % Dummy measurement weights at index 1
                    this.AssocWeightsMatrix(this.UnassocTrackInds,1) = 1;

                    % Create Hypothesis net for each cluster and populate association weights matrix
                    NumClusters = numel(this.ClusterList);
                    if NumClusters<3 && isempty(this.UnassocTrackInds)
                    end
                    for clusterInd=1:NumClusters
                        
                        % Extract track and measurement list for cluster
                        MeasIndList = this.ClusterList{clusterInd}.MeasIndList;
                        TrackIndList = this.ClusterList{clusterInd}.TrackIndList;
                        numClusterMeasurements = length(MeasIndList);
                        numClusterTracks = length(TrackIndList);
                        
                        % Compute probability of gating
                        if(isa(this.Gater,'EllipsoidalGaterX'))
                            GatingProbability = this.Gater.GatingProbability;
                        else
                            GatingProbability = 1;
                        end
                        
                        this.ClusterList{clusterInd}.AssocLikelihoodMatrix = ...
                            zeros(numClusterTracks, numClusterMeasurements+1);
                            
                        % Compute New Track/False Alarm density for the cluster
                        for t = 1:numClusterTracks
                            trackInd = TrackIndList(t);
                            if(~isprop(this.TrackList{trackInd},'ClutterDensity'))
                               this.TrackList{trackInd}.addprop('ClutterDensity');
                            end
                            if(isempty(this.ClutterModel))
                                if(isempty(this.GateVolumes))
                                    error('Cannot perform non-parametric (J)PDA as no gater has been set and thus gate volumes cannot be computed.. TODO: Allow for search space volume (V) to be set in the future');
                                end
                                clutterDensity = ...
                                    GatingProbability*numClusterMeasurements/sum(this.GateVolumes(TrackIndList));
                            else
                                clutterDensity = this.ClutterModel.pdf(this.TrackList{trackInd}.Filter.MeasurementPrediction.Mean);
                            end
                            this.ClusterList{clusterInd}.AssocLikelihoodMatrix(t,1) = clutterDensity*(1-this.DetectionProbability*GatingProbability);
                            this.TrackList{trackInd}.ClutterDensity = clutterDensity; 
                        end
                        
                        % Compute likelihood matrix
                        this.ClusterList{clusterInd}.AssocLikelihoodMatrix(:,2:end) = ...
                            this.DetectionProbability*GatingProbability*this.LikelihoodMatrix(TrackIndList,MeasIndList);
                        this.AssocLikelihoodMatrix(TrackIndList,[1 MeasIndList+1]) = ...
                            this.ClusterList{clusterInd}.AssocLikelihoodMatrix;
                        
                        % Compute association weights
                        this.ClusterList{clusterInd}.AssocWeightsMatrix =...
                            this.Hypothesiser.hypothesise(this.ClusterList{clusterInd}.AssocLikelihoodMatrix);
                        
                        % Store cluster association weights in global
                        % associatiopn weights matrix
                        this.AssocWeightsMatrix(TrackIndList, [1 MeasIndList+1]) = ...
                            this.ClusterList{clusterInd}.AssocWeightsMatrix./sum(this.ClusterList{clusterInd}.AssocWeightsMatrix,2);
                        
                    end 
                end
            end
        end
    end
end