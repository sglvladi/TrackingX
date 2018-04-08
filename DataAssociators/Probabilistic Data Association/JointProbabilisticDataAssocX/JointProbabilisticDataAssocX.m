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
        % * ProbOfDetect        (scalar) The target detection probability
        %
        %  See also JointProbabilisticDataAssocX/associate, JointProbabilisticDataAssocX/updateTracks.   
                    
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'Gater'))
                        this.Gater  = config.Gater;
                    end
                    if (isfield(config,'Clusterer'))
                        this.Clusterer  = config.Clusterer;
                    end
                    if (isfield(config,'ProbOfDetect'))
                        this.ProbOfDetect  = config.ProbOfDetect;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
            if (isfield(config,'Gater'))
                this.Gater  = config.Gater;
            end
            if (isfield(config,'Clusterer'))
                this.Clusterer  = config.Clusterer;
            end
            if (isfield(config,'ProbOfDetect'))
                this.ProbOfDetect  = config.ProbOfDetect;
            end
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
                    
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'Gater'))
                        this.Gater  = config.Gater;
                    end
                    if (isfield(config,'Clusterer'))
                        this.Clusterer  = config.Clusterer;
                    end
                    if (isfield(config,'ProbOfDetect'))
                        this.ProbOfDetect  = config.ProbOfDetect;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
            if (isfield(config,'Gater'))
                this.Gater  = config.Gater;
            end
            if (isfield(config,'Clusterer'))
                this.Clusterer  = config.Clusterer;
            end
            if (isfield(config,'ProbOfDetect'))
                this.ProbOfDetect  = config.ProbOfDetect;
            end
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
            if(sum(sum(this.ValidationMatrix))==0)
                this.AssocLikelihoodMatrix = [ones(this.NumTracks,1) zeros(this.NumTracks,this.NumMeas)];
                this.AssocWeightsMatrix = [ones(this.NumTracks,1) zeros(this.NumTracks,this.NumMeas)];
                return;
            else
                if(isempty(this.Clusterer))
                    % Compute New Track/False Alarm density for the cluster
                    if(isempty(this.ClutterDensity))
                        clutterDensity = sum(sum(this.ValidationMatrix))/sum(GateVolumes);
                    else
                        clutterDensity = this.ClutterDensity*sum(GateVolumes);
                    end
                    if(clutterDensity==0)
                        clutterDensity = eps;
                    end
                    this.AssocLikelihoodMatrix = ...
                        [ones(this.NumTracks,1)*clutterDensity*(1-this.ProbOfDetect*this.Gater.ProbOfGating), ...
                         this.ProbOfDetect*this.Gater.ProbOfGating*this.LikelihoodMatrix];
                    this.AssocWeightsMatrix = this.Hypothesiser.hypothesise(this.LikelihoodMatrix);
                else
                    % Get all clusters
                    [this.ClusterList, this.UnassocTrackInds] = this.Clusterer.cluster(this.ValidationMatrix);

                    % Allocate memory for association weights and fill in weights for unassociated tracks
                    this.AssocWeightsMatrix = zeros(this.NumTracks, this.NumMeas+1); % Dummy measurement weights at index 1
                    this.AssocWeightsMatrix(this.UnassocTrackInds,1) = 1;

                    % Create Hypothesis net for each cluster and populate association weights matrix
                    NumClusters = numel(this.ClusterList);
                    %this.NetList = cell(1,NumClusters);
                    for clusterInd=1:NumClusters
                        ObsIndList = this.ClusterList{clusterInd}.ObsIndList;
                        TrackIndList = this.ClusterList{clusterInd}.TrackIndList;
                        % Compute New Track/False Alarm density for the cluster
                        this.ClusterList{clusterInd}.ClutterDensity = ...
                            sum(sum(this.ValidationMatrix(TrackIndList,:)))...
                            /sum(this.GateVolumes(TrackIndList));
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
                        this.ClusterList{clusterInd}.AssocLikelihoodMatrix = ...
                            [ones(numel(TrackIndList),1)*this.ClusterList{clusterInd}.ClutterDensity*(1-this.ProbOfDetect*this.Gater.ProbOfGating), ...
                            this.ProbOfDetect*this.Gater.ProbOfGating*this.LikelihoodMatrix(TrackIndList,ObsIndList)];
                        this.AssocLikelihoodMatrix(TrackIndList,[1 ObsIndList+1]) = ...
                            this.ClusterList{clusterInd}.AssocLikelihoodMatrix;
                        this.ClusterList{clusterInd}.AssocWeightsMatrix =...
                            this.Hypothesiser.hypothesise(this.ClusterList{clusterInd}.AssocLikelihoodMatrix);
                        %this.ClusterList{clusterInd}.NetObj = this.Hypothesiser.NetObj;
                        this.AssocWeightsMatrix(TrackIndList, [1 ObsIndList+1]) = this.ClusterList{clusterInd}.AssocWeightsMatrix;
                    end               
                end
            end
        end
    end
end