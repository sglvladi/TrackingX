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
    end
    
    methods (Access = protected)
                
        function evaluateAssociations(this)
            this.NumTracks = numel(this.TrackList);
            this.NumMeas = size(this.ValidationMatrix,2);
            
            if(isa(this.Gater,'EllipsoidalGaterX'))
                ProbOfGating = this.Gater.ProbOfGating;
            else
                ProbOfGating = 1;
            end
            
            % If no valid measurement associations exist
            if(sum(sum(this.ValidationMatrix))==0)
                
                this.AssocLikelihoodMatrix = [ones(this.NumTracks,1) zeros(this.NumTracks,this.NumMeas)];
                this.AssocWeightsMatrix = [ones(this.NumTracks,1) zeros(this.NumTracks,this.NumMeas)];
                
                clutterDensity = eps;
                
                % Calculate the likelihood ratio
                a = [clutterDensity*(1-cellfun(@(x) x.ProbOfExist, this.TrackList))',...
                     clutterDensity*(1-this.ProbOfDetect*ProbOfGating)*cellfun(@(x) x.ProbOfExist, this.TrackList)'];
                w = a(:,1)./a(:,2);


                % Update existence probabilities
                for trackInd = 1:this.NumTracks
                    this.TrackList{trackInd}.ProbOfExist = this.AssocWeightsMatrix(trackInd,1)/(1+w(trackInd)) + sum(this.AssocWeightsMatrix(trackInd,2:end),2);
                end
                return;
            else
                if(isempty(this.Clusterer))
                    % Compute New Track/False Alarm density for the cluster
                    if(isempty(this.ClutterDensity))
                        clutterDensity = sum(sum(this.ValidationMatrix))/sum(this.GateVolumes);
                    else
                        clutterDensity = this.ClutterDensity;
                    end
                    if(clutterDensity==0)
                        clutterDensity = eps;
                    end
                    
                    % Compute combined association likelihood
                    this.AssocLikelihoodMatrix = zeros(this.NumTracks, this.NumMeas+1);
                    for trackInd = 1:this.NumTracks
                        this.AssocLikelihoodMatrix(trackInd,1) = ...
                            clutterDensity*((1-this.TrackList{trackInd}.ProbOfExist)...
                            +(1-this.ProbOfDetect*ProbOfGating)*this.TrackList{trackInd}.ProbOfExist);
                        this.AssocLikelihoodMatrix(trackInd,2:end) = ...
                            this.ProbOfDetect*ProbOfGating*this.TrackList{trackInd}.ProbOfExist...
                                * this.LikelihoodMatrix(trackInd,:);
                    end
                    
                    % Generate joint association weights
                    this.AssocWeightsMatrix = this.Hypothesiser.hypothesise(this.AssocLikelihoodMatrix);                                       
                    
                    % Calculate the likelihood ratio
                    a = [clutterDensity*(1-cellfun(@(x) x.ProbOfExist, this.TrackList))',...
                         clutterDensity*(1-this.ProbOfDetect*ProbOfGating)*cellfun(@(x) x.ProbOfExist, this.TrackList)'];
                    w = a(:,1)./a(:,2);
                    
                    
                    % Update existence probabilities
                    for trackInd = 1:this.NumTracks
                        this.TrackList{trackInd}.ProbOfExist = this.AssocWeightsMatrix(trackInd,1)/(1+w(trackInd)) + sum(this.AssocWeightsMatrix(trackInd,2:end),2);
                    end
                    
                    % Condition on existence
                    this.AssocWeightsMatrix(:,1) = this.AssocWeightsMatrix(:,1)./(1+w);
                    this.AssocWeightsMatrix = this.AssocWeightsMatrix./sum(this.AssocWeightsMatrix,2);
                           
                else
                    % Get all clusters
                    [this.ClusterList, this.UnassocTrackInds] = this.Clusterer.cluster(this.ValidationMatrix);

                    % Allocate memory for association weights and fill in weights for unassociated tracks
                    this.AssocWeightsMatrix = zeros(this.NumTracks, this.NumMeas+1); % Dummy measurement weights at index 1
                    this.AssocWeightsMatrix(this.UnassocTrackInds,1) = 1;
                    
                    % Process unassociated tracks
                    if(numel(this.UnassocTrackInds)>0)
                        clutterDensity = eps;

                        % Calculate the likelihood ratio
                        a = [clutterDensity*(1-cellfun(@(x) x.ProbOfExist, this.TrackList(this.UnassocTrackInds)))',...
                             clutterDensity*(1-this.ProbOfDetect*ProbOfGating)*cellfun(@(x) x.ProbOfExist, this.TrackList(this.UnassocTrackInds))'];
                        w = a(:,1)./a(:,2);

                        % Update existence probabilities
                        for i = 1:numel(this.UnassocTrackInds)
                            trackInd = this.UnassocTrackInds(i);
                            this.TrackList{trackInd}.ProbOfExist = ...
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
                                this.ClusterList{clusterInd}.ClutterDensity*((1-this.TrackList{trackInd}.ProbOfExist)...
                                +(1-this.ProbOfDetect*ProbOfGating)*this.TrackList{trackInd}.ProbOfExist);
                            this.ClusterList{clusterInd}.AssocLikelihoodMatrix(i,2:end) = ...
                                this.ProbOfDetect*ProbOfGating*this.TrackList{trackInd}.ProbOfExist...
                                    * this.LikelihoodMatrix(trackInd,ObsIndList);
                        end

                        % Generate joint association weights
                        this.ClusterList{clusterInd}.AssocWeightsMatrix = ...
                            this.Hypothesiser.hypothesise(this.ClusterList{clusterInd}.AssocLikelihoodMatrix);
                        
                        % Calculate the likelihood ratio
                        a = [this.ClusterList{clusterInd}.ClutterDensity*(1-cellfun(@(x) x.ProbOfExist, this.TrackList(TrackIndList)))',...
                             this.ClusterList{clusterInd}.ClutterDensity*(1-this.ProbOfDetect*ProbOfGating)*cellfun(@(x) x.ProbOfExist, this.TrackList(TrackIndList))'];
                        w = a(:,1)./a(:,2);

                        % Update existence probabilities
                        for i = 1:numel(TrackIndList)
                            trackInd = TrackIndList(i);
                            this.TrackList{trackInd}.ProbOfExist = ...
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