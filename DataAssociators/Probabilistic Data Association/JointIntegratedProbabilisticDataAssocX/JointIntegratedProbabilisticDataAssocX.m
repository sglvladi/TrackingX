classdef JointIntegratedProbabilisticDataAssocX < JointProbabilisticDataAssocX
% JointIntegratedProbabilisticDataAssocX class
%
% Summary of JointProbabilisticDataAssocX:
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
%   + JointProbabilisticDataAssocX  - Constructor method
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
        function predictTracks(this)
        % predictTracks - Performs JPDAF track prediction step
        %   
        % DESCRIPTION: 
        % * predictTracks(jpda) performs JPDAF update on all tracks contained
        %   in jpda.TrackList.
        %
        %   See also JointProbabilisticDataAssocX/initialise, JointProbabilisticDataAssocX/updateTracks.
            
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
            
            if(isa(this.Gater,'EllipsoidalGaterX'))
                GatingProbability = this.Gater.GatingProbability;
            else
                GatingProbability = 1;
            end
            
            % If no valid measurement associations exist
            if(sum(sum(this.ValidationMatrix))==0)
                
                this.AssocLikelihoodMatrix = [ones(numTracks,1) zeros(numTracks,numMeasurements)];
                this.AssocWeightsMatrix = [ones(numTracks,1) zeros(numTracks,numMeasurements)];
                
                clutterDensity = eps;
                
                % Calculate the likelihood ratio
                a = [clutterDensity*(1-cellfun(@(x) x.ExistenceProbability, this.TrackList))',...
                     clutterDensity*(1-this.DetectionProbability*GatingProbability)*cellfun(@(x) x.ExistenceProbability, this.TrackList)'];
                w = a(:,1)./a(:,2);


                % Update existence probabilities
                for trackInd = 1:numTracks
                    this.TrackList{trackInd}.ExistenceProbability = this.AssocWeightsMatrix(trackInd,1)/(1+w(trackInd)) + sum(this.AssocWeightsMatrix(trackInd,2:end),2);
                end
                return;
            else
                if(isempty(this.Clusterer))
                    
                    % Compute combined association likelihood
                    this.AssocLikelihoodMatrix = zeros(numTracks, numMeasurements+1);
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
                    
                        this.AssocLikelihoodMatrix(trackInd,1) = ...
                            clutterDensity*((1-this.TrackList{trackInd}.ExistenceProbability)...
                            +(1-this.DetectionProbability*GatingProbability)*this.TrackList{trackInd}.ExistenceProbability);
                        this.AssocLikelihoodMatrix(trackInd,2:end) = ...
                            this.DetectionProbability*GatingProbability*this.TrackList{trackInd}.ExistenceProbability...
                                * this.LikelihoodMatrix(trackInd,:);
                    end
                    
                    % Generate joint association weights
                    this.AssocWeightsMatrix = this.Hypothesiser.hypothesise(this.AssocLikelihoodMatrix);                                       
                                        
                    
                    % Update existence probabilities
                    a = zeros(numTracks,2);
                    w = zeros(numTracks,1);
                    for trackInd = 1:numTracks
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
                        % Calculate the likelihood ratio
                        a(trackInd,:) = [clutterDensity*(1-this.TrackList{trackInd}.ExistenceProbability),...
                             clutterDensity*(1-this.DetectionProbability*GatingProbability)*this.TrackList{trackInd}.ExistenceProbability];
                        w(trackInd) = a(trackInd,1)./a(trackInd,2);
                        this.TrackList{trackInd}.ExistenceProbability = this.AssocWeightsMatrix(trackInd,1)/(1+w(trackInd)) + sum(this.AssocWeightsMatrix(trackInd,2:end),2);
                    end
                    
                    % Condition on existence
                    this.AssocWeightsMatrix(:,1) = this.AssocWeightsMatrix(:,1)./(1+w);
                    this.AssocWeightsMatrix = this.AssocWeightsMatrix./sum(this.AssocWeightsMatrix,2);
                           
                else
                    % Get all clusters
                    [this.ClusterList, this.UnassocTrackInds] = this.Clusterer.cluster(this.ValidationMatrix);

                    % Allocate memory for association weights and fill in weights for unassociated tracks
                    this.AssocWeightsMatrix = zeros(numTracks, numMeasurements+1); % Dummy measurement weights at index 1
                    this.AssocWeightsMatrix(this.UnassocTrackInds,1) = 1;
                    
                    % Process unassociated tracks
                    if(numel(this.UnassocTrackInds)>0)
                        clutterDensity = eps;

                        % Calculate the likelihood ratio
                        a = [clutterDensity*(1-cellfun(@(x) x.ExistenceProbability, this.TrackList(this.UnassocTrackInds)))',...
                             clutterDensity*(1-this.DetectionProbability*GatingProbability)*cellfun(@(x) x.ExistenceProbability, this.TrackList(this.UnassocTrackInds))'];
                        w = a(:,1)./a(:,2);

                        % Update existence probabilities
                        for i = 1:numel(this.UnassocTrackInds)
                            trackInd = this.UnassocTrackInds(i);
                            this.TrackList{trackInd}.ExistenceProbability = ...
                                this.AssocWeightsMatrix(trackInd,1)/(1+w(i)) + sum(this.AssocWeightsMatrix(trackInd,2:end),2);
                        end
                    end
                    
                    % Create Hypothesis net for each cluster and populate association weights matrix
                    NumClusters = numel(this.ClusterList);
                    for clusterInd=1:NumClusters
                        
                        % Extract track and measurement list for cluster
                        ObsIndList = this.ClusterList{clusterInd}.ObsIndList;
                        TrackIndList = this.ClusterList{clusterInd}.TrackIndList;
                        
                        % Compute New Track/False Alarm density for the cluster
                        if(isempty(this.ClutterDensity))
                            this.ClusterList{clusterInd}.ClutterDensity = ...
                                sum(sum(this.ValidationMatrix(TrackIndList,:)))...
                                /sum(this.GateVolumes(TrackIndList));
                        else
                            this.ClusterList{clusterInd}.ClutterDensity = this.ClutterDensity;
                        end
                        for t = 1:numel(TrackIndList)
                            trackInd = TrackIndList(t);
                            if(~isprop(this.TrackList{trackInd},'ClutterDensity'))
                               this.TrackList{trackInd}.addprop('ClutterDensity');
                            end
                            this.TrackList{trackInd}.ClutterDensity = this.ClusterList{clusterInd}.ClutterDensity; 
                        end
                        if(this.ClusterList{clusterInd}.ClutterDensity==0)
                            this.ClusterList{clusterInd}.ClutterDensity = 1;
                        end
                        
                        % Compute combined association likelihood
                        this.ClusterList{clusterInd}.AssocLikelihoodMatrix = ...
                                zeros(numel(TrackIndList), numel(ObsIndList)+1);
                        for i = 1:numel(TrackIndList)
                            trackInd = TrackIndList(i);
                            this.ClusterList{clusterInd}.AssocLikelihoodMatrix(i,1) = ...
                                this.ClusterList{clusterInd}.ClutterDensity*((1-this.TrackList{trackInd}.ExistenceProbability)...
                                +(1-this.DetectionProbability*GatingProbability)*this.TrackList{trackInd}.ExistenceProbability);
                            this.ClusterList{clusterInd}.AssocLikelihoodMatrix(i,2:end) = ...
                                this.DetectionProbability*GatingProbability*this.TrackList{trackInd}.ExistenceProbability...
                                    * this.LikelihoodMatrix(trackInd,ObsIndList);
                        end

                        % Generate joint association weights
                        this.ClusterList{clusterInd}.AssocWeightsMatrix = ...
                            this.Hypothesiser.hypothesise(this.ClusterList{clusterInd}.AssocLikelihoodMatrix);
                        
                        % Calculate the likelihood ratio
                        a = [this.ClusterList{clusterInd}.ClutterDensity*(1-cellfun(@(x) x.ExistenceProbability, this.TrackList(TrackIndList)))',...
                             this.ClusterList{clusterInd}.ClutterDensity*(1-this.DetectionProbability*GatingProbability)*cellfun(@(x) x.ExistenceProbability, this.TrackList(TrackIndList))'];
                        w = a(:,1)./a(:,2);

                        % Update existence probabilities
                        for i = 1:numel(TrackIndList)
                            trackInd = TrackIndList(i);
                            this.TrackList{trackInd}.ExistenceProbability = ...
                                this.ClusterList{clusterInd}.AssocWeightsMatrix(i,1)/(1+w(i))...
                                + sum(this.ClusterList{clusterInd}.AssocWeightsMatrix(i,2:end),2);
                        end
                        
                        % Condition on existence
                        this.ClusterList{clusterInd}.AssocWeightsMatrix(:,1) = ...
                            this.ClusterList{clusterInd}.AssocWeightsMatrix(:,1)./(1+w);
                        this.ClusterList{clusterInd}.AssocWeightsMatrix = ...
                            this.ClusterList{clusterInd}.AssocWeightsMatrix./sum(this.ClusterList{clusterInd}.AssocWeightsMatrix,2);
                        
                        % Store cluster association weights in global
                        % associatiopn weights matrix
                        this.AssocWeightsMatrix(TrackIndList, [1 ObsIndList+1]) = ...
                            this.ClusterList{clusterInd}.AssocWeightsMatrix;
                    end 
                end
            end
        end
    end
end