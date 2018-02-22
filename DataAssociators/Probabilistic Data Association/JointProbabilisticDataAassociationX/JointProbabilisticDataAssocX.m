classdef JointProbabilisticDataAssocX < handle
% JointProbabilisticDataAssocX class
%
% Summary of JointProbabilisticDataAssocX:
% This is a class implementation of a Joint Probabilistic Data Association Filter.
%
% JointProbabilisticDataAssocX Properties:
%   + TrackList        A (1-by-NumTracks) vector of TrackX objects
%   + MeasurementList           A (1-by-NumObs) vector of observations/measurements
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
        NumTracks
        NumObs
        TrackList
        MeasurementList
        ValidationMatrix
        LikelihoodMatrix
        AssocWeightsMatrix
        Gater = []
        Clusterer = []
        ClutterDensity
        ProbOfDetect = 1
        ClusterList = []
        UnassocTrackInds = []
        NetList = []
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
        
            if(nargin>1)
                this.TrackList = TrackList;
                this.MeasurementList = MeasurementList;
            end
            % Get number of available tracks and observations
            this.NumTracks = size(this.TrackList,2);
            this.NumObs   = size(this.MeasurementList,2);
            
            if(~isempty(this.TrackList))
                % Validation matix and volume
                [this.ValidationMatrix, GateVolumes] = this.Gater.gate(this.TrackList,this.MeasurementList);
                this.LikelihoodMatrix = zeros(this.NumTracks, this.NumObs);

                % Compute Likelihood matrix
                for trackInd = 1:this.NumTracks          
                    % Extract valid measurements
                    ValidObsInds = find(this.ValidationMatrix(trackInd,:));
                    ValidObs = this.MeasurementList(:,ValidObsInds);
                    this.TrackList{trackInd}.Filter.Measurement = ValidObs;

                    if(isa(this.TrackList{trackInd}.Filter,'KalmanFilterX')...
                       ||isa(this.TrackList{trackInd}.Filter,'UnscentedParticleFilterX')...
                       ||isa(this.TrackList{trackInd}.Filter,'ExtendedParticleFilterX'))

                        % Extract predicted measurement and innovation covariance from filter
                        predMeasMean = this.TrackList{trackInd}.Filter.PredMeasMean;
                        innovCovar = this.TrackList{trackInd}.Filter.InnovErrCovar;

                        % Compute likelihood matrix
                        this.LikelihoodMatrix(trackInd, ValidObsInds) = ...
                            gauss_pdf(ValidObs, predMeasMean,innovCovar);

                    elseif(isa(this.TrackList{trackInd}.Filter,'ParticleFilterX'))    
                        % Compute likelihood matrix via expected likelihood
                        this.LikelihoodMatrix(trackInd, ValidObsInds) = ...
                            mean(this.TrackList{trackInd}.Filter.ObsModel.eval(...
                                ValidObs, this.TrackList{trackInd}.Filter.PredParticles),2)';
                    end
                end

                if(isempty(this.Clusterer))
                    % If no clusterer has been supplied
                    Cluster.TrackIndList = [];
                    Cluster.ObsIndList = [];
                    UnassocTrackInds = [];
                    for trackInd=1:this.NumTracks
                        ValidObs = find(this.ValidationMatrix(trackInd,:));
                        if(~isempty(ValidObs))
                            Cluster.TrackIndList = union(Cluster.TrackIndList,trackInd);
                            Cluster.ObsIndList = union(Cluster.ObsIndList,ValidObs);
                        else
                            UnassocTrackInds = union(UnassocTrackInds,trackInd); %#ok<*PROPLC>
                        end
                    end
                    this.UnassocTrackInds = UnassocTrackInds;
                    this.ClusterList(1) = Cluster;
                else
                    % Get all clusters
                    [this.ClusterList, this.UnassocTrackInds] = this.Clusterer.cluster(this.ValidationMatrix);
                end

                % Allocate memory for association weights and fill in weights for unassociated tracks
                this.AssocWeightsMatrix = zeros(this.NumTracks, this.NumObs+1); % Dummy measurement weights at index 1
                this.AssocWeightsMatrix(this.UnassocTrackInds,1) = 1;

                % Create Hypothesis net for each cluster and populate association weights matrix
                NumClusters = numel(this.ClusterList);
                this.NetList = cell(1,NumClusters);
                for clusterInd=1:NumClusters
                    Cluster = this.ClusterList{clusterInd};
                    ClustObsIndList = Cluster.ObsIndList;
                    ClustTrackIndList = Cluster.TrackIndList;
                    % Compute New Track/False Alarm density for the cluster
                    Cluster.ClutterDensity = ...
                        sum(sum(this.ValidationMatrix(ClustTrackIndList,:)))...
                        /sum(GateVolumes(ClustTrackIndList));
                    if(Cluster.ClutterDensity==0)
                        Cluster.ClutterDensity = 1;
                    end
                    ClustLi = [ones(numel(ClustTrackIndList),1)*Cluster.ClutterDensity*(1-this.ProbOfDetect*this.Gater.ProbOfGating), ...
                        this.ProbOfDetect*this.Gater.ProbOfGating*this.LikelihoodMatrix(ClustTrackIndList,ClustObsIndList)];
                    ClustVm = [ones(numel(ClustTrackIndList),1),  this.ValidationMatrix(ClustTrackIndList, ClustObsIndList)];
                    this.NetList{clusterInd} = buildEHMnet_trans(ClustVm, ClustLi);
                    this.AssocWeightsMatrix(ClustTrackIndList, [1, ClustObsIndList+1]) = this.NetList{clusterInd}.betta;
                end                   
            else
                fprintf('No tracks where found. Skipping JPDAF association step...\n');
                this.ValidationMatrix = zeros(1, size(this.Params.DataList,2));
                this.ClusterList = [];
                this.UnassocTrackInds = [];
                this.NetList = [];
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
        %   See also JointProbabilisticDataAssocX/initialise, JointProbabilisticDataAssocX/updateTracks.
        
            if(~isempty(this.TrackList))
                % Compute weights and update each track
                for trackInd=1:this.NumTracks
                    this.TrackList{trackInd}.Filter.predict();
                end    
            else
                fprintf('No tracks where found. Skipping JPDAF prediction step...\n');
            end
        end
        
        function updateTracks(this)
        % updateTracks - Performs JPDAF track update step
        %   
        % DESCRIPTION: 
        % * updateTracks(jpda) performs JPDAF update on all tracks contained
        %   in jpda.TrackList.
        %
        %   See also JointProbabilisticDataAssocX/initialise, JointProbabilisticDataAssocX/updateTracks.
        
            if(~isempty(this.TrackList))
                % Compute weights and update each track
                for trackInd=1:this.NumTracks
                    ValidDataInd = find(this.ValidationMatrix(trackInd,:));    % Associated measurements
                    assocWeights = this.AssocWeightsMatrix(trackInd,[1 ValidDataInd+1]);
                    this.TrackList{trackInd}.Filter.updatePDA(assocWeights);
                end    
            else
                fprintf('No tracks where found. Skipping JPDAF Update step...\n');
            end
        end
        
    end
end