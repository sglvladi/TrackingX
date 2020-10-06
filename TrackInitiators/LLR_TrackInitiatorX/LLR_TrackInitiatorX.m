 classdef LLR_TrackInitiatorX < TrackInitiatorX 
% LLR_TrackInitiatorX Abstract class
%
% Summary of LLR_TrackInitiatorX:
% This is a simple M out of N track initiator 
%
% LLR_TrackInitiatorX Properties:
%   TODO!!
%
% LLR_TrackInitiatorX Methods:
%   + LLR_TrackInitiatorX - Constructor method
%
% (+) denotes public properties/methods
%
% March 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
       TagGenerator
       TrackList
       MeasurementList
       AssocWeightsMatrix
       DataAssociator
       InitFilter
       CustomConfirmConditionFcn
       CustomDeleteConditionFcn
    end
        
    methods
        function this = LLR_TrackInitiatorX(varargin)
        % LLR_TrackInitiatorX Constructor method
        %   
        % DESCRIPTION: 
        % * LLR_TrackInitiatorX() 
        %   returns a LLR_TrackInitiatorX object handle configured with the provided
        %   PHD filter instance PHDFilter and probability of gating ProbOfGating
            
            if(nargin==0)
                error('Not enough input arguments');
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.TagGenerator = config.TagGenerator;
                    this.InitFilter = config.InitFilter;
                    this.DataAssociator = config.DataAssociator;
                    if(isfield(config,"CustomConfirmConditionFcn"))
                        this.CustomConfirmConditionFcn = config.CustomConfirmConditionFcn;
                    else
                       this.CustomConfirmConditionFcn=0;
                    end
                    if(isfield(config,"CustomDeleteConditionFcn"))
                        this.CustomDeleteConditionFcn = config.CustomDeleteConditionFcn;
                    else
                       this.CustomDeleteConditionFcn=0;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.TagGenerator = config.TagGenerator;
            this.InitFilter = config.InitFilter;
            this.DataAssociator = config.DataAssociator;
            if(isfield(config,"CustomConfirmConditionFcn"))
                this.CustomConfirmConditionFcn = config.CustomConfirmConditionFcn;
            else
               this.CustomConfirmConditionFcn=0;
            end
            if(isfield(config,"CustomDeleteConditionFcn"))
                this.CustomDeleteConditionFcn = config.CustomDeleteConditionFcn;
            else
               this.CustomDeleteConditionFcn=0;
            end
        end
        
        function [ConfirmedTrackList, TentativeTrackList, deletedTracks] = initiateTracks(this, ConfirmedTrackList, MeasurementList, AssocWeightsMatrix, LikelihoodMatrix)
        % INITIATETRACKS Perform track initiation
        %   
        % DESCRIPTION: 
        % * [TrackList] = initiateTracks(this,TrackList,MeasurementList) 
        %   utilises the provided lists of tracks (TrackList) and measurements
        %   (MeasurementList) 
        %
        % INPUTS:
        % * TrackList - (1-by-NumTracks) array of TrackX objects
        % * MeasurementList - (NumObsDims-by-NumMeas) matrix of
        %   measurements
        %
        % See also EllipsoidalGaterX/EllipsoidalGaterX
            
            THD = -8.9; %-5.9; 
            %% Confirmed tracks
            numConfirmedTracks = numel(ConfirmedTrackList);
            deletedTracks = {};
            if numConfirmedTracks>0
                for t=1:numConfirmedTracks
                    % Update LPR
                    Pd = this.DataAssociator.DetectionModel.pdf(ConfirmedTrackList{t}.Filter.StatePrediction.Mean);
                    Pg = this.DataAssociator.Gater.GatingProbability;
                    validInds = find(AssocWeightsMatrix(t,1:end));
                    DL = 0; 
                    for j=1:numel(validInds)
                        measInd = validInds(j);
                        if measInd == 1
                            DL = DL + AssocWeightsMatrix(t,1)*(1-Pd*Pg);
%                             ConfirmedTrackList{t}.LPR = ConfirmedTrackList{t}.LPR + AssocWeightsMatrix(t,1)*log((1-Pd*Pg));
                        else
                            lambda = this.DataAssociator.ClutterModel.pdf(ConfirmedTrackList{t}.Filter.MeasurementPrediction.Mean);
                            DL = DL + AssocWeightsMatrix(t,measInd)*Pd*LikelihoodMatrix(t,measInd-1)/lambda;
%                             ConfirmedTrackList{t}.LPR = ConfirmedTrackList{t}.LPR + AssocWeightsMatrix(t,measInd)*log(Pd*LikelihoodMatrix(t,measInd-1)/lambda);
                        end
                    end
                    ConfirmedTrackList{t}.LPR = ConfirmedTrackList{t}.LPR + log(DL);
                    % Update detection history
                    validInds = find(AssocWeightsMatrix(t,2:end));
                    validMeas = MeasurementListX(MeasurementList.getMeasurements(validInds));
                    ConfirmedTrackList{t}.Associations(end+1) = validMeas;
                    % Update maximum LPR
                    if (ConfirmedTrackList{t}.LPR>ConfirmedTrackList{t}.LPR_max)
                        ConfirmedTrackList{t}.LPR_max = ConfirmedTrackList{t}.LPR;
                    end
                end

                % Delete Tracks
                del_tracks = 0;
                del_flag = 0;
                for t = 1:numConfirmedTracks
                    % Check condition
                    del_cond = isnan(ConfirmedTrackList{t}.LPR) || ConfirmedTrackList{t}.LPR < ConfirmedTrackList{t}.LPR_max+THD;
                    if isa(this.CustomDeleteConditionFcn, 'function_handle')
                        c_del_cond = this.CustomDeleteConditionFcn(this,ConfirmedTrackList{t});
                        del_cond = del_cond || c_del_cond;
                    end
                    if(del_cond)
                        deletedTracks{end+1} = ConfirmedTrackList{t};
                        ConfirmedTrackList{t} = [];
                        del_tracks = del_tracks + 1;
                        del_flag = 1;
                    end
                end
                if(del_flag)
                    ConfirmedTrackList = ConfirmedTrackList(~cellfun('isempty',ConfirmedTrackList));
                    numConfirmedTracks = numConfirmedTracks - del_tracks;
                end
            end

            %% Tentative tracks 

            % Get indices of all measurements where there exist no possible association
            numMeas = MeasurementList.NumMeasurements;
            rhi = zeros(1, numMeas);
            for measInd = 1:numMeas
                rhi_tmp = 1;
                if(AssocWeightsMatrix>-1) % Check if beta exists
                    for trackInd = 1:numConfirmedTracks
                        rhi_tmp = rhi_tmp*(1-AssocWeightsMatrix(trackInd,measInd+1));
                    end
                end
                rhi(measInd) = rhi_tmp;
            end
            UnassocMeasInd = rhi==1; %any(sum(AssocWeightsMatrix,1)==0,1);

            this.DataAssociator.MeasurementList = MeasurementListX(MeasurementList.Measurements(UnassocMeasInd)); % Only use unassociated measurements
            this.DataAssociator.predictTracks();
            this.DataAssociator.associate();    
            this.DataAssociator.updateTracks();
            
            % Update Detection and Miss counters for all tracks
            for t=1:numel(this.DataAssociator.TrackList)
                % Update counts
                Pd = this.DataAssociator.DetectionModel.pdf(this.DataAssociator.TrackList{t}.Filter.StatePrediction.Mean);
                Pg = this.DataAssociator.Gater.GatingProbability;
                validInds = find(this.DataAssociator.AssocWeightsMatrix(t,1:end));
                DL = 0;
                for j=1:numel(validInds)
                    measInd = validInds(j);
                    if measInd == 1
                        DL = DL + this.DataAssociator.AssocWeightsMatrix(t,1)*(1-Pd*Pg);
%                         this.DataAssociator.TrackList{t}.LPR = this.DataAssociator.TrackList{t}.LPR + this.DataAssociator.AssocWeightsMatrix(t,1)*log((1-Pd*Pg));
                    else
                        lambda = this.DataAssociator.ClutterModel.pdf(this.DataAssociator.TrackList{t}.Filter.MeasurementPrediction.Mean);
                        DL = DL + this.DataAssociator.AssocWeightsMatrix(t,measInd)*Pd*this.DataAssociator.LikelihoodMatrix(t,measInd-1)/lambda;
%                         this.DataAssociator.TrackList{t}.LPR = this.DataAssociator.TrackList{t}.LPR + this.DataAssociator.AssocWeightsMatrix(t,measInd)*log(Pd*this.DataAssociator.LikelihoodMatrix(t,measInd-1)/lambda);
                    end
                end
                this.DataAssociator.TrackList{t}.LPR = this.DataAssociator.TrackList{t}.LPR + log(DL);
                % Update detection history
                validInds = find(this.DataAssociator.AssocWeightsMatrix(t,2:end));
                validMeas = MeasurementListX(this.DataAssociator.MeasurementList.getMeasurements(validInds));
                this.DataAssociator.TrackList{t}.Associations(end+1) = validMeas;
                % Update maximum LPR
                if (this.DataAssociator.TrackList{t}.LPR>this.DataAssociator.TrackList{t}.LPR_max)
                    this.DataAssociator.TrackList{t}.LPR_max = this.DataAssociator.TrackList{t}.LPR;
                end
            end

            % Delete or Confirm Tentative Tracks
            del_tracks = 0;
            del_flag = 0;
            for t = 1:numel(this.DataAssociator.TrackList)
                P_TM = 0.1; %0.1;       % Track miss probability (from Blackman & Popoli)
                P_FC = 10^-4; %10^-6;    % False confirm probability (from Blackman & Popoli)
                
%                 % High clutter region
%                 xk = this.DataAssociator.TrackList{t}.State.Mean;
%                 x = [-469.710602736631,-2572.92295458513,-2516.31234127506,-1320.74163107899,-835.307363807233, -469.710602736631];
%                 y = [-4915.41679013513,-4065.50122858976,-2858.60918927224,-1380.19032924249,-1261.08679763585, -4915.41679013513];
%                 b = find(inpolygon(xk(1,:),xk(3,:),x,y));
%                 if numel(b)>0
%                     P_FC = 10^-10;
%                 end
                % High and low thresholds for confirmation and deletion
                %  of tentative tracks
                gamma_low = log(P_TM/(1-P_FC));
                gamma_high = log((1-P_TM)/P_FC);
                
                P = this.DataAssociator.TrackList{t}.State.Covar;
                
                % Check condition
                del_cond = isnan(this.DataAssociator.TrackList{t}.LPR) || this.DataAssociator.TrackList{t}.LPR < gamma_low;
                if isa(this.CustomDeleteConditionFcn, 'function_handle')
                    c_del_cond = this.CustomDeleteConditionFcn(this,this.DataAssociator.TrackList{t});
                    del_cond = del_cond || c_del_cond;
                end
                if(del_cond)
                    this.DataAssociator.TrackList{t} = []; 
                    del_tracks = del_tracks + 1;
                    del_flag = 1;
                elseif(this.DataAssociator.TrackList{t}.LPR > gamma_high)
%                     disp("Confirming new track");
                    ConfirmedTrackList{end+1} = this.DataAssociator.TrackList{t};
                    this.DataAssociator.TrackList{t} = [];
                    del_tracks = del_tracks + 1;
                    del_flag = 1;
                end
            end
            if(del_flag)
                this.DataAssociator.TrackList = this.DataAssociator.TrackList(~cellfun('isempty',this.DataAssociator.TrackList));
            end

            % Initiate new tracks
            UnassocMeasInd = find(~any(this.DataAssociator.ValidationMatrix,1));
            for i=1:numel(UnassocMeasInd)
                
                measurement = this.DataAssociator.MeasurementList.Vectors(:,UnassocMeasInd(i));
                if(~isempty(this.DataAssociator.MeasurementList.Tags))
                    tag = this.DataAssociator.MeasurementList.Tags(UnassocMeasInd(i));
                else
                    tag = TagX(i);
                end
                stateMean = this.InitFilter.Model.Measurement.finv(measurement);
                dist = GaussianDistributionX(stateMean,... 
                                            this.InitFilter.StatePrior.Covar);
                if(~isempty(MeasurementList.Timestamp))
                    statePrior = GaussianStateX(dist,...
                                                MeasurementList.Timestamp);
                else
                    statePrior = GaussianStateX(dist);
                end
                track_config = struct('Filter', copy(this.InitFilter), ...
                                      'Tag', this.TagGenerator.generate(), ...
                                      'LPR', 1,...
                                      'LPR_max', 1,...
                                      'ExistenceProbability', 0.5,...
                                      'Associations', MeasurementListX(measurement,MeasurementList.Timestamp,tag));
                track_config.Filter.initialise('Model',this.InitFilter.Model,...
                                               'StatePrior', statePrior);       
                track_config.Trajectory = track_config.Filter.StatePosterior;
                newTrack = TrackX(track_config);

                this.DataAssociator.TrackList{end+1} = newTrack;
            end
            TentativeTrackList = this.DataAssociator.TrackList;
        end 
    end
end