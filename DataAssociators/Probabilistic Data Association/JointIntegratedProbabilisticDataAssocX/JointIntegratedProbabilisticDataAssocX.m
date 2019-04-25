classdef JointIntegratedProbabilisticDataAssocX < JointProbabilisticDataAssocX
% JointIntegratedProbabilisticDataAssocX class
%
% Summary of JointIntegratedProbabilisticDataAssocX:
% This is a class implementation of a Joint Integrated Probabilistic Data Association Filter.
%
% JointIntegratedProbabilisticDataAssocX Properties:
%   + TrackList        A (1-by-NumTracks) vector of TrackX objects
%   + MeasurementList           A (1-by-NumMeas) vector of observations/measurements
%   + Gater             A GaterX object used to perform gating
%   + Clusterer         A ClustererX object used to perform clustering of
%                       tracks
%   + ClutterDensity    
%
% JointIntegratedProbabilisticDataAssocX Methods:
%   + JointIntegratedProbabilisticDataAssocX  - Constructor method
%   + associate                     - Performs JPDAF association step
%   + updateTracks                 - Performs JPDAF update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% See also ParticleFilterX, KalmanFilerX.

    properties
        SurvivalProbability = 0.99
    end
    
    methods
        
        function this = JointIntegratedProbabilisticDataAssocX(varargin)
        % JointIntegratedProbabilisticDataAssocX - Constructor method
        %   
        % DESCRIPTION: 
        % * jpda = JointIntegratedProbabilisticDataAssocX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * jpda = JointIntegratedProbabilisticDataAssocX(gater,clusterer) returns an 
        %   object handle, preconfigured with the provided GaterX and ClustererX 
        %   object handles gater and clusterer.
        % * jpda = JointIntegratedProbabilisticDataAssocX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        % INPUT ARGUMENTS:
        % * Gater               (GaterX) A gater object which should be used to  
        %                       perform gating of measurements. Default = None
        % * Clusterer           (ClustererX) A clusterer object which should be used to  
        %                       perform clustering of tracks. Default = None
        % * DetectionProbability        (scalar) The target detection probability
        %
        %  See also JointIntegratedProbabilisticDataAssocX/associate, JointIntegratedProbabilisticDataAssocX/updateTracks.   
                    
            this@JointProbabilisticDataAssocX(varargin{:});
            % TODO:
            % Proper initiation!!
            % ==================>
        end
        
        function predictTracks(this)
        % predictTracks - Performs JPDAF track prediction step
        %   
        % DESCRIPTION: 
        % * predictTracks(jpda) performs JPDAF update on all tracks contained
        %   in jpda.TrackList.
        %
        %   See also JointIntegratedProbabilisticDataAssocX/initialise, JointIntegratedProbabilisticDataAssocX/updateTracks.
            
            for trackInd = 1:this.NumTracks
                this.TrackList{trackInd}.ExistenceProbability = this.SurvivalProbability*this.TrackList{trackInd}.ExistenceProbability;
            end
            % Call super class method
            predictTracks@JointProbabilisticDataAssocX(this);
        end
    end
    
    methods (Access = protected)
                
        function evaluateAssociations(this)
            
            numTracks = this.NumTracks;
            numMeasurements = this.NumMeasurements;
            
            % Allocate memory for association likelihoods & weights 
            % (Dummy measurement at index 1)
            this.AssocLikelihoodMatrix = zeros(numTracks, numMeasurements+1); 
            this.AssocWeightsMatrix = zeros(numTracks,numMeasurements+1);

            % Create Hypothesis net for each cluster and populate association weights matrix
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
                    % (No need to perform data-association)
                    this.AssocLikelihoodMatrix(TrackIndList,:) = [1 zeros(1,numMeasurements)];
                    this.AssocWeightsMatrix(TrackIndList,:) = [1 zeros(1,numMeasurements)];
                    
                    % Update existence probabilities
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
                        
                        % Extract existence probability
                        Pe = this.TrackList{trackInd}.ExistenceProbability;
                        
                        % Calculate the likelihood ratio
                        b = [lambda*(1 - Pe),...
                             lambda*(1 - Pd*Pg)*Pe];
                        w = b(1)/b(2);
                        
                        % Compute updated existence probability
                        this.TrackList{trackInd}.ExistenceProbability = 1/(1+w);
                    
                    end
                    continue;
                else
                    
                    % Pre-allocate memory
                    this.ClusterList(clusterInd).AssocLikelihoodMatrix = ...
                        zeros(numClusterTracks, numClusterMeasurements+1);
                    this.ClusterList(clusterInd).AssocWeightsMatrix = ...
                        zeros(numClusterTracks, numClusterMeasurements+1);
                        
                    % Compute New Track/False Alarm density for the cluster
                    w = zeros(1,numClusterTracks);
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
                        
                        % Extract existence probability
                        Pe = this.TrackList{trackInd}.ExistenceProbability;
                        
                        % Compute combined association likelihood
                        this.ClusterList(clusterInd).AssocLikelihoodMatrix(t,:) = ...
                            [lambda*((1-Pe)+(1-Pd*Pg)*Pe), Pd*Pg*Pe*this.LikelihoodMatrix(trackInd,MeasIndList)];
                        
                        % Calculate the likelihood ratio
                        b = [lambda*(1 - Pe),...
                             lambda*(1 - Pd*Pg)*Pe];
                        w(t) = b(1)/b(2); 
                    end
                    
                    % Compute likelihood matrix
                    this.AssocLikelihoodMatrix(TrackIndList,[1 MeasIndList+1]) = ...
                        this.ClusterList(clusterInd).AssocLikelihoodMatrix;
                    
                    % Generate joint association weights
                    this.ClusterList(clusterInd).AssocWeightsMatrix = ...
                        this.Hypothesiser.hypothesise(this.ClusterList(clusterInd).AssocLikelihoodMatrix);

                    % Update existence probabilities
                    for t = 1:numClusterTracks
                        trackInd = TrackIndList(t);
                        this.TrackList{trackInd}.ExistenceProbability = ...
                                this.ClusterList(clusterInd).AssocWeightsMatrix(t,1)/(1+w(t))...
                                + sum(this.ClusterList(clusterInd).AssocWeightsMatrix(t,2:end),2);
                    end

                    % Condition on existence
                    this.ClusterList(clusterInd).AssocWeightsMatrix(:,1) = ...
                        this.ClusterList(clusterInd).AssocWeightsMatrix(:,1)./(1+w)';
                    this.ClusterList(clusterInd).AssocWeightsMatrix = ...
                        this.ClusterList(clusterInd).AssocWeightsMatrix./sum(this.ClusterList(clusterInd).AssocWeightsMatrix,2);

                    % Store cluster association weights in global
                    % associatiopn weights matrix
                    this.AssocWeightsMatrix(TrackIndList, [1 MeasIndList+1]) = ...
                        this.ClusterList(clusterInd).AssocWeightsMatrix;
                end 
            end
        end
    end
end