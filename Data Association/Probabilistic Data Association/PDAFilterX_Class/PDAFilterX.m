classdef PDAFilterX < handle
    % PDAFilterX class
    %
    % Summary of PDAFX:
    % This is a class implementation of a Probabilistic Data Association Filter.
    %
    % PDAFX Properties:
    %    - Params       = structure, with fields:
    %       .k                      = Time index. Can also act as a time interval (Dt), depending on the underlying models.
    %                                 The time index of each filter is synced with this value.
    %       .DataList               = All available observations at time k
    %       .TrackList              = A list of known tracks at time k
    %       .pDetect                = Probability of detection
    %       .performGating          = Set (true|false) to (enable|disable) gating. Default = true (*)
    %       .pGate                  = Probability of gating. Default = .993    | Only necessary if gating is enabled.
    %       .gateLevel              = Gate size as number of std. Default = 10 | .gateLevel can be computed by running 'chi2inv(.pGate, dim);', where dim is the dimensionality of the gate.
    %                                 (When no gating is performed, it is implied that all targets share all measurements, i.e. they form a single cluster)
    %       .adaptiveLambda (**)    = Set (true|false) to (enable|disable) adaptive computation of clutter density. Default = true
    %       .lambda                 = Clutter density (?). Clutter is assumed to be Poisson distributed with mean .lambda. Only necessary if .adaptiveLambda is disabled
    %
    %   (*) Gating and clustering provide great gains in terms of performance, and should ONLY be turned off for debugging/evaluation purposes.
    %  (**) Adaptive computation of lambda is done on a cluster-by-cluster basis. See Predict method comments for more info.
    %
    % PDAFX Methods:
    %    PDAFX          - Constructor method
    %    Predict         - Performs PDAF prediction step
    %    Update          - Performs PDAF update step
    %    Iterate         - Performs a complete PDAF iteration (Predict & Update)
    % 
    % PDAFX Example:


    properties
        Params
    end
    
    methods
        function this = PDAFilterX(Init)
        % PDAFX - Constructor method
        %   
        %   Inputs:
        %       Params    |=> Check class help for more details|
        %   
        %   Usage:
        %       jpdaf = PDAFX(Params); 
        %
        %   See also Predict, Update, Iterate, Smooth.
                    
            % Validate .pDetect ~~~~~~~~~~~~~~~~~~~~~>
            if ~isfield(Init, 'pDetect')
                error('[PDAF] Probability of detection (.pDetect) has not been set!');
            else
                this.Params.pDetect = Init.pDetect;
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            % Validate .performGating ~~~~~~~~~~~~~~~>
            if ~isfield(Init,'performGating')
                disp("[PDAF] Applying default settings '.performGating=true");
                this.Params.performGating = true;
            else
                this.Params.performGating = Init.performGating;
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
             % Validate .pGate, .gateLevel ~~~~~~~~~~>
            if (this.Params.performGating && (~isfield(Init,'pGate')||~isfield(Init,'gateLevel')))
               if (~isfield(Init,'pGate')&&~isfield(Init,'gateLevel'))
                   disp("[PDAF] Applying default settings '.pGate = .993' and '.gateLevel = chi2inv(.pGate,2)' = 10(approx)");
                   this.Params.pGate = .993;
                   this.Params.gateLevel = chi2inv(this.Params.pGate,2);
               elseif (~isfield(Init,'pGate'))
                   disp("[PDAF] Applying default setting '.pGate=0.993'");
                   this.Params.pGate = 0.993;
               else
                   disp("[PDAF] Applying default setting '.gateLevel = chi2inv(.pGate,2)' = 10(approx)");
                   this.Params.gateLevel = chi2inv(this.Params.pGate);
               end
            elseif (~this.Params.performGating)
                this.Params.pGate = 1;
                if (isfield(this.Params.pGate)||isfield(this.Params.gateLevel))
                    warning('[PDAF] provided values for .pGate and/or .gateLevel will be ignored because .performGating = false!');
                end
            else
                this.Params.pGate = Init.pGate;
                this.Params.gateLevel = Init.gateLevel;
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            % Validate .performClustering ~~~~~~~~~~~>
            if ~isfield(Init,'performClustering')
                if(~this.Params.performGating)
                    this.Params.performClustering = false;
                    disp("[PDAF] Gating turned off. Applying default setting '.performClustering = false'");
                else
                    this.Params.performClustering = true;
                    disp("[PDAF] Applying default setting '.performClustering = true'");
                end
            else
                this.Params.performClustering = Init.performClustering;
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            % Validate .adaptiveLambda ~~~~~~~~~~~~~~>
            if ~isfield(Init,'adaptiveLambda')
                this.Params.adaptiveLambda = true;
                disp("[PDAF] Applying default setting '.adaptiveLambda = true'");
            else
                this.Params.adaptiveLambda = Init.adaptiveLambda;
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            % Validate .lambda ~~~~~~~~~~~~~~>
            if ~isfield(Init,'lambda')
                if(~this.Params.adaptiveLambda)
                    error('[PDAF] A value for .lambda must be provided when .adaptiveLambda=false');
                end
            else
                this.Params.lambda = Init.lambda;
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
      
        end
        
        function Predict(this)
        % Predict - Performs PDAF prediction step
        %   
        %   Inputs:
        %       N/A
        %
        %   Properties required:
        %       Params:
        %           .TrackList : with Predict, getYpred, ObsModel.eval_likelihood methods
        %           .DataList       
        %
        %   Outputs:
        %       N/A
        %   
        %   Properties modified:
        %       Params:
        %           .ValidationMatrix
        %           .LikelihoodMatrix
        %           .AssocWeightMatrix
        %           .ClusterList
        %           .NetList
        %
        %   (NOTE: The time index/interval "this.Params.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (jpdaf.Params.k = 1; % 1 sec)
        %       jpdaf.Predict();
        %
        %   See also PDAFX, Update, Iterate, Smooth.
        
            
            % Get number of available tracks and observations
            this.Params.nTracks = size(this.Params.TrackList,2);
            this.Params.nData   = size(this.Params.DataList,2);
            
            % Validation matix and volume
            ValidationMatrix = zeros(this.Params.nTracks, this.Params.nData); % (Nt x Nm)
            LikelihoodMatrix = zeros(this.Params.nTracks, this.Params.nData); % (Nt x Nm)
            V_k = zeros(1, this.Params.nTracks); % Validation region volume per track (1 x Nt)
            
            if(~isempty(this.Params.TrackList))
                
                % Predict and construct Validation and Likelihood matrices
                for trackInd = 1:this.Params.nTracks
                    
                    TrackObj = this.Params.TrackList{trackInd}.TrackObj;
                    TrackObj.Params.k = this.Params.k;
                    
                    % Predict
                    TrackObj.Predict();
                    
                    % Extract predicted measurement and innovation covariance from filter
                    [yPred, S] = TrackObj.getYpred();
                    
                    % Perform Gating
                    switch numel(yPred)
                        case 1 
                            C = 2;
                        case 2
                            C = pi;
                        case 3
                            C = 4*pi/3;
                        otherwise
                            error('[PDAF] Gating has only been implemented for observations of up to 3 dimensions!');
                    end
                    V_k(trackInd) = C*this.Params.gateLevel^(numel(yPred)/2)*det(S)^(1/2);    
                    
                    % Perform gating (if enabled)
                    if(this.Params.performGating)
                        ValidationMatrix(trackInd,:) = this.mahalDist(this.Params.DataList, yPred, S, 2) < this.Params.gateLevel;
                    else
                         ValidationMatrix(trackInd,:) = ones(1,this.Params.nData);
                    end
                    
                    % Extract valid measurements
                    ValidDataInd = find(ValidationMatrix(trackInd,:));
                    ValidData = this.Params.DataList(:,ValidDataInd);
                    TrackObj.Params.y = ValidData;
                    
                    % Compute Likelihood matrix
                    if(isa(TrackObj,'KalmanFilterX')||isa(TrackObj,'UParticleFilterX')||isa(TrackObj,'EParticleFilterX'))
                        LikelihoodMatrix(trackInd, ValidDataInd) = TrackObj.ObsModel.eval(TrackObj.Params.k, ValidData, TrackObj.Params.xPred)';
                    elseif(isa(TrackObj,'ParticleFilterX'))
                        TrackObj.Params.LikelihoodMatrix = TrackObj.ObsModel.eval(TrackObj.Params.k, ValidData, TrackObj.Params.particles);
                        LikelihoodMatrix(trackInd, ValidDataInd) = sum(TrackObj.Params.LikelihoodMatrix,2)'/TrackObj.Params.Np;
                    end
                    this.Params.TrackList{trackInd}.TrackObj = TrackObj;
                end
                
                % Store validation and likelihood matrices
                this.Params.ValidationMatrix = ValidationMatrix;
                this.Params.LikelihoodMatrix = LikelihoodMatrix;
                
                % Get all clusters
                this.FormClusters(); % Each target forms a cluster 
                   
                % Allocate memory for association weights (?) and fill in weights for unassociated tracks
                AssocWeightsMatrix = zeros(this.Params.nTracks, this.Params.nData+1); % Dummy measurement weights at index 1

                % Create Hypothesis net for each cluster and populate association weights matrix
                nClusters = numel(this.Params.ClusterList);
                NetList = cell(1,nClusters);
                for clusterInd=1:nClusters
                    Cluster = this.Params.ClusterList{clusterInd};
                    ClustMeasIndList = Cluster.MeasIndList;
                    ClustTrackIndList = Cluster.TrackIndList;
                    % Compute New Track/False Alarm density for the cluster
                    Cluster.lambda = sum(sum(this.Params.ValidationMatrix(ClustTrackIndList,:)))/sum(V_k(ClustTrackIndList));
                    if(Cluster.lambda==0)
                        Cluster.lambda = 1;
                    end
                    ClustLi = [ones(numel(ClustTrackIndList),1)*Cluster.lambda*(1-this.Params.pDetect*this.Params.pGate), this.Params.pDetect*this.Params.pGate*this.Params.LikelihoodMatrix(ClustTrackIndList,ClustMeasIndList)];
                    ClustVm = [ones(numel(ClustTrackIndList),1),  this.Params.ValidationMatrix(ClustTrackIndList, ClustMeasIndList)];
                    this.Params.ClusterList{clusterInd} = Cluster;
                    NetList{clusterInd} = buildEHMnet_trans(ClustVm, ClustLi);
                    AssocWeightsMatrix(ClustTrackIndList, [1, ClustMeasIndList+1]) = NetList{clusterInd}.betta;
                end
                
                % Store computed EHM nets and association weights matrix
                this.Params.NetList = NetList;
                this.Params.AssocWeightsMatrix = AssocWeightsMatrix;
            else
                fprintf('No tracks where found. Skipping PDAF Predict step...\n');
                this.Params.ValidationMatrix = zeros(1, size(this.Params.DataList,2));
                this.Params.AssocWeightsMatrix = -1; % Set betta to -1
            end
        end
        
        function Update(this)
            if(~isempty(this.Params.TrackList))
                % Update all tracks
                for trackInd=1:this.Params.nTracks
                    ValidDataInd = find(this.Params.ValidationMatrix(trackInd,:));    % Associated measurements
                    assocWeights = this.Params.AssocWeightsMatrix(trackInd,[1 ValidDataInd+1]);
                    this.Params.TrackList{trackInd}.TrackObj.UpdatePDA(assocWeights);
                end    
            else
                fprintf('No tracks where found. Skipping PDAF Update step...\n');
            end
        end
        
        function FormClusters(this)
           
            % Form clusters of tracks sharing measurements
            ClusterList = cell(1,this.Params.nTracks);
            ClusterObj.MeasIndList = [];
            ClusterObj.TrackIndList = [];
            for t=1:this.Params.nTracks
                ClusterList{t} = ClusterObj;
                ClusterList{t}.MeasIndList = find(this.Params.ValidationMatrix(t,:));
                ClusterList{t}.TrackIndList = t;
            end
            this.Params.ClusterList = ClusterList;
        end
        
    end
    
    methods (Static)
        function D=mahalDist(x, m, C, use_log)
        % p=gaussian_prob(x, m, C, use_log)
        %
        % Evaluate the multi-variate density with mean vector m and covariance
        % matrix C for the input vector x.
        % Vectorized version: Here X is a matrix of column vectors, and p is 
        % a vector of probabilities for each vector.

            if nargin<4, use_log = 0; end

            d   = length(m);

            if size(x,1)~=d
               x=x';
            end
            N       = size(x,2);

            m       = m(:);
            M       = m*ones(1,N);
            denom   = (2*pi)^(d/2)*sqrt(abs(det(C)));
            invC    = inv(C);
            mahal   = sum(((x-M)'*invC).*(x-M)',2);   % Chris Bregler's trick

            switch use_log,
            case 2,
              D     = mahal;
            case 1,
              D     = -0.5*mahal - log(denom);
            case 0,
              numer = exp(-0.5*mahal);
              D     = numer/denom;
            otherwise
                error('Unsupported log type')
            end
        end
    end
end