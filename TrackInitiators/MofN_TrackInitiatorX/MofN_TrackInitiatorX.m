classdef MofN_TrackInitiatorX < TrackInitiatorX 
% MofN_TrackInitiatorX Abstract class
%
% Summary of MofN_TrackInitiatorX:
% This is a simple M out of N track initiator 
%
% MofN_TrackInitiatorX Properties:
%   None
%
% MofN_TrackInitiatorX Methods:
%   + MofN_TrackInitiatorX - Constructor method
%
% (+) denotes public properties/methods
%
% March 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
       TrackList
       MeasurementList
       AssocWeightsMatrix
       DataAssociator
       InitFilter
       ConfirmThreshold
       DeleteThreshold
       WindowSize
       CustomConfirmConditionFcn
       CustomDeleteConditionFcn
       TrackIdGenFcn
    end
        
    methods
        function this = MofN_TrackInitiatorX(varargin)
        % MofN_TrackInitiatorX Constructor method
        %   
        % DESCRIPTION: 
        % * MofN_TrackInitiatorX('PHDFilter',PHDFilter,'ProbOfGating',ProbOfGating) 
        %   returns a MofN_TrackInitiatorX object handle configured with the provided
        %   PHD filter instance PHDFilter and probability of gating ProbOfGating
            
            if(nargin==0)
                error('Not enough input arguments');
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.InitFilter = config.InitFilter;
                    this.DataAssociator = config.DataAssociator;
                    this.DeleteThreshold = config.DeleteThreshold;
                    this.ConfirmThreshold = config.ConfirmThreshold;
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
                    if(isfield(config,"TrackIdGenFcn"))
                        this.TrackIdGenFcn = config.TrackIdGenFcn;
                    else
                       this.TrackIdGenFcn=@(x) randi(10000);
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
            this.InitFilter = config.InitFilter;
            this.DataAssociator = config.DataAssociator;
            this.ConfirmThreshold = config.ConfirmThreshold;
            this.DeleteThreshold = config.DeleteThreshold;
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
            
            this.WindowSize = max(this.ConfirmThreshold,this.DeleteThreshold);
        end
        
        function [ConfirmedTrackList, TentativeTrackList] = initiateTracks(this,ConfirmedTrackList, MeasurementList, AssocWeightsMatrix)
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
        
            %% Confirmed tracks
            numConfirmedTracks = numel(ConfirmedTrackList);
            if numConfirmedTracks>0

                % Update detection history
                for t=1:numConfirmedTracks
                    track_detected = ~isempty(find(ConfirmedTrackList{t}.ValidationMatrix, 1));
                    ConfirmedTrackList{t}.DetectionHistory = [ConfirmedTrackList{t}.DetectionHistory(2:end),track_detected];
                end

                % Delete Tracks
                del_tracks = 0;
                del_flag = 0;
                for t = 1:numConfirmedTracks
                    M = this.DeleteThreshold(2);
                    % Check condition
                    del_cond = sum(ConfirmedTrackList{t}.DetectionHistory(end-(M-1):end)==0) >= this.DeleteThreshold(1);
                    if isa(this.CustomDeleteConditionFcn, 'function_handle')
                        del_cond = del_cond || this.CustomDeleteConditionFcn(this,ConfirmedTrackList{t});
                    end
                    if(del_cond)
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
            numMeas = size(MeasurementList,2);
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

            this.DataAssociator.MeasurementList = MeasurementList(:,UnassocMeasInd); % Only use unassociated measurements
            for j=1:numel(this.DataAssociator.TrackList)
                this.DataAssociator.TrackList{j}.Filter.predict();
            end
            this.DataAssociator.associate();    
            this.DataAssociator.updateTracks();
            
            % Update Detection and Miss counters for all tracks
            for t=1:numel(this.DataAssociator.TrackList)

                % Update counts
                track_detected = ~isempty(find(this.DataAssociator.TrackList{t}.ValidationMatrix, 1));
                this.DataAssociator.TrackList{t}.DetectionHistory = [this.DataAssociator.TrackList{t}.DetectionHistory(2:end),track_detected]; 
            end

            % Delete or Confirm Tentative Tracks
            del_tracks = 0;
            del_flag = 0;
            for t = 1:numel(this.DataAssociator.TrackList)
                P = this.DataAssociator.TrackList{t}.Filter.StateCovar;
                M_c = this.ConfirmThreshold(2);
                M_d = this.DeleteThreshold(2);
                % Check condition
                del_cond = sum(this.DataAssociator.TrackList{t}.DetectionHistory(end-(M_d-1):end)==0) >= this.DeleteThreshold(1);
                if isa(this.CustomDeleteConditionFcn, 'function_handle')
                    del_cond = del_cond || this.CustomDeleteConditionFcn(this,this.DataAssociator.TrackList{t});
                end
                if(del_cond)% || P(1,1)>sigmaDelThresh^2 || P(2,2)>sigmaDelThresh^2) 
                    this.DataAssociator.TrackList{t} = []; 
                    del_tracks = del_tracks + 1;
                    del_flag = 1;
                elseif(sum(this.DataAssociator.TrackList{t}.DetectionHistory(end-(M_c-1):end)==1)>= this.ConfirmThreshold(1))
                    disp("Confirming new track");
                    ConfirmedTrackList{end+1} = this.DataAssociator.TrackList{t};
                    this.DataAssociator.TrackList{t} = [];
                    del_tracks = del_tracks + 1;
                    del_flag = 1;
                end
            end
            if(del_flag)
                this.DataAssociator.TrackList = this.DataAssociator.TrackList(~cellfun('isempty',this.DataAssociator.TrackList));
                this.DataAssociator.NumTracks = this.DataAssociator.NumTracks - del_tracks;
            end

            UnassocMeasInd = find(~any(this.DataAssociator.ValidationMatrix,1));
            % Initiate new tracks
            for i=1:numel(UnassocMeasInd)
                              
                NewTrack = TrackX();
                NewTrack.addprop('Filter');
                NewTrack.addprop('TrackID');
                NewTrack.addprop('DetectionHistory');
                NewTrack.TrackID = this.TrackIdGenFcn(this);
                windowSize = max(this.ConfirmThreshold(2), this.DeleteThreshold(2));
                NewTrack.DetectionHistory = [-ones(1,windowSize-1),1];
                NewTrack.Filter = copy(this.InitFilter);
                stateMean = zeros(this.InitFilter.Model.Dyn.NumStateDims,1);
                for j=1:this.InitFilter.Model.Obs.NumObsDims
                    index = this.InitFilter.Model.Obs.Mapping(j);
                    if(index~=0)
                        stateMean(index) = this.DataAssociator.MeasurementList(j,UnassocMeasInd(i));
                    end
                end
                if(isa(this.InitFilter,'KalmanFilterX'))
                    NewTrack.Filter.StateMean = zeros(this.InitFilter.Model.Dyn.NumStateDims,1);
                    NewTrack.Filter.StateMean = stateMean;
                    %filter.Params.particles = filter.Params.gen_x0([jpdaf_init.Params.DataList(1,UnassocMeasInd(i)); jpdaf_init.Params.DataList(2,UnassocMeasInd(i));0;0],filter.Params.Np);
                    %filter.Params.x = sum(filter.Params.w.*filter.Params.particles,2);
                else
                    stateCovar =  NewTrack.Filter.StateCovar;
                    NewTrack.Filter.initialise(this.InitFilter.Model,'PriorStateMean',stateMean, 'PriorStateCovar',stateCovar);
                end 
                this.DataAssociator.TrackList{end+1} = NewTrack;
                this.DataAssociator.NumTracks = this.DataAssociator.NumTracks + 1;  
            end
            TentativeTrackList = this.DataAssociator.TrackList;
        end 
    end
end