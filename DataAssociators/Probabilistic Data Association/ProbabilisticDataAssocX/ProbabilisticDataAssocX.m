classdef ProbabilisticDataAssocX < DataAssociatorX
% ProbabilisticDataAssocX class
%
% Summary of ProbabilisticDataAssocX:
% This is a class implementation of a Joint Probabilistic Data Association Filter.
%
% ProbabilisticDataAssocX Properties:
%   + TrackList - A (1 x Tracks) vector of TrackX objects
%   + MeasurementList - A (1 x NumMeas) vector of observations/measurements
%   + Gater - A GaterX object used to perform gating
%   + Clusterer - A ClustererX object used to perform clustering of
%                       tracks
%   + ClutterDensity 
%   + LikelihoodMatrix - A (
%
% ProbabilisticDataAssocX Methods:
%   + ProbabilisticDataAssocX  - Constructor method
%   + associate - Performs JPDAF association step
%   + updateTracks                 - Performs JPDAF update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% See also ParticleFilterX, KalmanFilerX.

    properties
        NumTracks
        NumMeas
        ValidationMatrix
        GateVolumes
        LikelihoodMatrix
        AssocLikelihoodMatrix
        AssocWeightsMatrix
        ClutterDensity
        ProbOfDetect = 1
        ClusterList = []
        UnassocTrackInds = []
    end
    
    methods
        function this = ProbabilisticDataAssocX(varargin)
        % ProbabilisticDataAssocX - Constructor method
        %   
        % DESCRIPTION: 
        % * jpda = ProbabilisticDataAssocX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * jpda = ProbabilisticDataAssocX(gater,clusterer) returns an 
        %   object handle, preconfigured with the provided GaterX and ClustererX 
        %   object handles gater and clusterer.
        % * jpda = ProbabilisticDataAssocX(___,Name,Value,___) instantiates an  
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
        %  See also ProbabilisticDataAssocX/associate, ProbabilisticDataAssocX/updateTracks.   
                    
            this.Hypothesiser = EfficientHypothesisManagementX();
            
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
        %   See also ProbabilisticDataAssocX/initialise, ProbabilisticDataAssocX/updateTracks.
        
            if(nargin>1)
                this.TrackList = TrackList;
                this.MeasurementList = MeasurementList;
            end
            % Get number of available tracks and observations
            this.NumTracks = size(this.TrackList,2);
            this.NumMeas   = size(this.MeasurementList,2);
            
            if(~isempty(this.TrackList))
                
                associate@DataAssociatorX(this);
                                    
            else
                fprintf('No tracks where found. Skipping JPDAF association step...\n');
                this.ValidationMatrix = zeros(1, size(this.MeasurementList,2));
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
        %   See also ProbabilisticDataAssocX/initialise, ProbabilisticDataAssocX/updateTracks.
        
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
    
    methods (Access = protected)
        
        function performGating(this)
            % Validation matix and volume
            [this.ValidationMatrix, this.GateVolumes] = this.Gater.gate(this.TrackList,this.MeasurementList);
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
            this.LikelihoodMatrix = zeros(this.NumTracks, this.NumMeas);
            % Compute Likelihood matrix
            for trackInd = 1:this.NumTracks          
                % Extract valid measurements
                ValidObsInds = find(this.TrackList{trackInd}.ValidationMatrix);
                ValidObs = this.MeasurementList(:,ValidObsInds);
                this.TrackList{trackInd}.Filter.Measurement = this.MeasurementList(:,find(this.TrackList{trackInd}.ValidationMatrix));
                
                if(isempty(ValidObs))
                    this.LikelihoodMatrix(trackInd, ValidObsInds) = 0;
                else
                    if(isa(this.TrackList{trackInd}.Filter,'KalmanFilterX')...
                       ||isa(this.TrackList{trackInd}.Filter,'UnscentedParticleFilterX')...
                       ||isa(this.TrackList{trackInd}.Filter,'ExtendedParticleFilterX'))

                        % Extract predicted measurement and innovation covariance from filter
                        predMeasMean = this.TrackList{trackInd}.Filter.PredMeasMean;
                        innovCovar = this.TrackList{trackInd}.Filter.InnovErrCovar;

                        % Compute likelihood matrix
                        this.LikelihoodMatrix(trackInd, ValidObsInds) = ...
                            this.TrackList{trackInd}.Filter.Model.Obs.pdf(ValidObs,this.TrackList{trackInd}.Filter.PredStateMean);
                            %gauss_pdf(ValidObs, predMeasMean,innovCovar);

                    elseif(isa(this.TrackList{trackInd}.Filter,'ParticleFilterX'))    
                        % Compute likelihood matrix via expected likelihood
                        this.LikelihoodMatrix(trackInd, ValidObsInds) = ...
                            mean(this.TrackList{trackInd}.Filter.MeasLikelihood,2)';
    %                         mean(this.TrackList{trackInd}.Filter.ObsModel.eval(...
    %                             ValidObs, this.TrackList{trackInd}.Filter.PredParticles),2)';
                    end
                end
            end
        end
        
        function performClustering(this)
            % Do nothing. PDAF does not perform clustering
        end
        
        function evaluateAssociations(this)

            % Allocate memory for association weights and fill in weights for unassociated tracks
            this.AssocLikelihoodMatrix = zeros(this.NumTracks, this.NumMeas+1);
            this.AssocWeightsMatrix = zeros(this.NumTracks, this.NumMeas+1); % Dummy measurement weights at index 1

            %this.NetList = cell(1,NumClusters);
            for trackInd=1:this.NumTracks
                
                % Extract valid measurements
                obsIndList = find(this.ValidationMatrix(trackInd,:));
                
                if(isempty(obsIndList))
                    this.AssocLikelihoodMatrix(trackInd,:) = [1 zeros(1,this.NumMeas)];
                    this.AssocWeightsMatrix(trackInd,:) = [1 zeros(1,this.NumMeas)];
                    continue;
                else
                    % Compute New Track/False Alarm density for the cluster
                    if(isempty(this.ClutterDensity))
                        clutterDensity = numel(obsIndList)/this.GateVolumes(trackInd);
                    else
                        clutterDensity = this.ClutterDensity*this.GateVolumes(trackInd);
                    end
                    this.AssocLikelihoodMatrix(trackInd,:) = ...
                        [clutterDensity*(1-this.ProbOfDetect*this.Gater.ProbOfGating), ...
                        this.ProbOfDetect*this.Gater.ProbOfGating*this.LikelihoodMatrix(trackInd,:)];
                    this.AssocWeightsMatrix(trackInd,:) = this.Hypothesiser.hypothesise(this.AssocLikelihoodMatrix(trackInd,:));
                end
            end          
        end
    end
end