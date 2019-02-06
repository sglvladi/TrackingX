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
       TagGenerator
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
                    this.TagGenerator = config.TagGenerator;
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
        
        function [ConfirmedTrackList, TentativeTrackList, deletedTracks] = initiateTracks(this,ConfirmedTrackList, MeasurementList, AssocWeightsMatrix)
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
            deletedTracks = {};
            if numConfirmedTracks>0

                % Update detection history
                for t=1:numConfirmedTracks
                    validInds = find(AssocWeightsMatrix(t,2:end));
                    if(numel(validInds)>1)
                    end
                    validMeas = MeasurementListX(MeasurementList.Vectors(validInds), MeasurementList.Timestamp);
                    track_detected = ~isempty(find(ConfirmedTrackList{t}.ValidationMatrix, 1));
                    ConfirmedTrackList{t}.DetectionHistory = [ConfirmedTrackList{t}.DetectionHistory(2:end),track_detected];
                    ConfirmedTrackList{t}.Trajectory(end+1) = ConfirmedTrackList{t}.Filter.StatePosterior;
                    ConfirmedTrackList{t}.Associations(end+1) = validMeas;
                end

                % Delete Tracks
                del_tracks = 0;
                del_flag = 0;
                for t = 1:numConfirmedTracks
                    M = this.DeleteThreshold(2);
                    % Check condition
                    del_cond = sum(ConfirmedTrackList{t}.DetectionHistory(end-(M-1):end)==0) >= this.DeleteThreshold(1);
                    if isa(this.CustomDeleteConditionFcn, 'function_handle')
                        c_del_cond = this.CustomDeleteConditionFcn(this,ConfirmedTrackList{t});
                        if(c_del_cond)
                            print('Deleting COVARIANCE based');
                        end
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

            this.DataAssociator.MeasurementList = MeasurementListX(MeasurementList.Vectors(:,UnassocMeasInd), MeasurementList.Timestamp); % Only use unassociated measurements
            this.DataAssociator.predictTracks();
            this.DataAssociator.associate();    
            this.DataAssociator.updateTracks();
            
            % Update Detection and Miss counters for all tracks
            for t=1:numel(this.DataAssociator.TrackList)
                % Update counts
                validInds = find(this.DataAssociator.ValidationMatrix(t,:));
                validMeas = MeasurementListX(MeasurementList.Vectors(:,validInds), MeasurementList.Timestamp);
                track_detected = ~isempty(validInds);
                this.DataAssociator.TrackList{t}.DetectionHistory = [this.DataAssociator.TrackList{t}.DetectionHistory(2:end),track_detected]; 
                this.DataAssociator.TrackList{t}.Trajectory(end+1) = this.DataAssociator.TrackList{t}.Filter.StatePosterior;
                this.DataAssociator.TrackList{t}.Associations(end+1) = validMeas;
            end

            % Delete or Confirm Tentative Tracks
            del_tracks = 0;
            del_flag = 0;
            for t = 1:numel(this.DataAssociator.TrackList)
                P = this.DataAssociator.TrackList{t}.State.Covar;
                M_c = this.ConfirmThreshold(2);
                M_d = this.DeleteThreshold(2);
                
                % Check condition
                del_cond = sum(this.DataAssociator.TrackList{t}.DetectionHistory(end-(M_d-1):end)==0) >= this.DeleteThreshold(1);
                if isa(this.CustomDeleteConditionFcn, 'function_handle')
                    c_del_cond = this.CustomDeleteConditionFcn(this,this.DataAssociator.TrackList{t});
                    if(del_cond)
                        print('Deleting COVARIANCE based');
                    end
                    del_cond = del_cond || c_del_cond;
                end
                if(del_cond)
                    disp("Deleting track");
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
            end

            % Initiate new tracks
            UnassocMeasInd = find(~any(this.DataAssociator.ValidationMatrix,1));
            for i=1:numel(UnassocMeasInd)
                
                measurement = this.DataAssociator.MeasurementList.Vectors(:,UnassocMeasInd(i));
                stateMean = this.InitFilter.Model.Measurement.finv(measurement);
                dist = GaussianDistributionX(stateMean,... 
                                            this.InitFilter.StatePrior.Covar);
                if(~isempty(MeasurementList.Timestamp))
                    statePrior = GaussianStateX(dist,...
                                                MeasurementList.Timestamp);
                else
                    statePrior = GaussianStateX(dist);
                end
                windowSize = max(this.ConfirmThreshold(2), this.DeleteThreshold(2));
                track_config = struct('Filter', copy(this.InitFilter), ...
                                      'Tag', this.TagGenerator.generate(), ...
                                      'DetectionHistory', [-ones(1,windowSize-1),1],...
                                      'ExistenceProbability', 0.5,...
                                      'Associations', MeasurementListX(measurement,MeasurementList.Timestamp));
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