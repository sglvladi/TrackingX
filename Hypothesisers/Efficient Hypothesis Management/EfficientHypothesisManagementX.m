classdef EfficientHypothesisManagementX < HypothesiserX
% EfficientHypothesisManagementX Abstract class
%
% Summary of EfficientHypothesisManagementX:
% This is an implementation of the Efficient Hypothesis Management (EHM)
% algorithm described in [1]
%
% HypothesiserX Properties:
%   None
%
% HypothesiserX Methods:
%   + HypothesiserX - Constructor method
%
% (+) denotes puplic properties/methods
%
% [1] Maskell, Simon et al. “Fast Mutual Exclusion.” (2004).
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        NetObj
        LikelihoodMatrix
        AssocWeightsMatrix
        Mode = 'Track-Oriented'
        timeout = 5
    end
    
    methods
        function this = EfficientHypothesisManagementX(varargin)
        % EFFICIENTHYPOTHESISMANAGEMENTX Constructor method
        % 
        % Parameters
        % ----------
        % Mode: string
        %   Specifies the mode of operation of the EHM algorithm, which can
        %   be set to either:
        %       1) 'Track-Oriented'
        %       2) 'Measurement-Oriented'
        %   (default='Track-Oriented')
        % 
        % Usage
        % -----
        % * EfficientHypothesisManagementX() creates a EHM object handle
        % * EfficientHypothesisManagementX('Mode',mode) creates an EHM
        %   object handle, configured to operate in Mode mode. Mode can be
        %   either 'Track-Oriented' (default) or 'Measurement-Oriented'.
        %
        % See also EfficientHypothesisManagementX/hypothesise
            % Return early if no arguments are supplied
            if(nargin==0)
                return;
            end
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'Mode'))
                        this.Mode  = config.Mode;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
            if (isfield(config,'Mode'))
                this.Mode  = config.Mode;
            end
        end
        
        function AssocWeightsMatrix = hypothesise(this,LikelihoodMatrix)
        % HYPOTHESISE Generate hypotheses and calculate association weights
        % 
        % Parameters
        % ----------
        % LikelihoodMatrix: (NumMeasurements+1 x NumTracks) matrix
        %   A matrix whose (i,j)th element is the likelihood of the (i-1)th
        %   measurement, given the measurement prediction of the jth track.
        %   The row (1,:) represents the dummy (null) measurement likelihood
        %   for each track.
        %
        % Returns
        % -------
        % AssocWeightsMatrix: (NumMeasurements+1 x NumTracks) matrix
        %   A matrix whose (i,j)th element is the joint probability that the 
        %   (i-1)th measurement is associated to the jth track. The row (1,:) 
        %   represents the dummy (null) measurement association probability
        %   for each target.
        %
        % Usage
        % -----
        % * hypothesise(ehm) generates all valid association hypotheses and  
        %   calculates the respective association weights, based on the
        %   likelihood matrix ehm.LikelihoodMatrix.
        % * hypothesise(ehm,lik) generates all valid association hypotheses and  
        %   calculates the respective association weights, based on the
        %   likelihood matrix lik. ehm.LikelihoodMatrix is also updated
        %   with lik.
        %
        % See also EfficientHypothesisManagementX/hypothesise
            % First check to see if a structure was received
            if(nargin==2)
                this.LikelihoodMatrix = LikelihoodMatrix;
            end
            ValidationMatrix = this.LikelihoodMatrix>0;
            [this.NetObj, timedOut] = this.buildNetTO(ValidationMatrix, this.timeout);
            if(timedOut)
                [NumTracks, NumMeasp1] = size(this.LikelihoodMatrix);
                AssocWeightsMatrix = [ones(NumTracks,1), zeros(NumTracks,NumMeasp1-1)];
                this.AssocWeightsMatrix = AssocWeightsMatrix;
            else
                [AssocWeightsMatrix, this.NetObj] = this.computeAssocWeights(this.NetObj,this.LikelihoodMatrix);
                this.AssocWeightsMatrix = AssocWeightsMatrix;
            end
        end
    end
    
    methods (Static)
        function [NetObj,timedOut] = buildNetTO(ValidationMatrix,timeOut)
        % BUILDNETTO Build Track-Oriented (TO) EHM net and compute respective 
        % association probabilities (betta).
        %
        % Parameters
        % ----------
        % ValidationMatrix: (T x M) matrix
        %   A matrix containing all possible measurement to track associations 
        %   (M: number of measurements, including dummy at index 1)
        %   (T: number of tracks)
        %
        % Returns
        % -------
        % NetObj: structure
        %   A net object with the following fields:
        %       + NetObj.NodeList - List of all Nodes (NodeObj) 
        %         contained in net.
        %       + NetObj.EdgeList - Cell Matrix (Nn-by-Nn), Nn being 
        %         the total number of nodes), where cell (i,j) contains 
        %         the tracks contained in edge from parent node i, to
        %         child node j.
        %
        % Usage
        % -----
        % * NetObj = buildNetTO(ValidationMatrix) generates and returns the
        %   EHM net NetObj.
        %
        % IMPORTANT REMINDER TO SELF:
        %   When computing the remainders for a node, we always look at the remaining
        %   association events in the !NEXT! layer and then get the difference
        %   between these and any entries in the node's MeasIndList.
        %
        % Author: Lyudmil Vladimirov
            
            % Set timedout to false
            timedOut = 0;
            tic;
            
            % Get number of tracks/layers
            TrackNum = size(ValidationMatrix,1); 
            LayerNum = TrackNum; % Layer 1 is root layer

            % Get number of Maasurements
            PointNum = size(ValidationMatrix,2);

            % Define Node Object
            NodeObj.TrackInd = 0; % Measurement index
            NodeObj.MeasIndList = []; % List of associated tracks
            NodeObj.ParentIndList = []; % List of Parent nodes
            NodeObj.ChildIndList = [];  % List of Child nodes
            NodeObj.Remainders = [ 1:PointNum ]'; % Remaining tracks

            % Define Net Object
            NetObj.NodeList = [];
            NetObj.EdgeList = [];
            NetObj.NodesPerTrack = cell(1,TrackNum);
            NetObj.ValidationMatrix = ValidationMatrix;

            % Create Root Node
            NetObj.NodeList{1} = NodeObj;

            % Build the net
            % ==============>
            % For every measurement/layer
            for j=1:LayerNum

                % Get index list (L_jm1_Ind) of nodes in previous layer (j-1)
                L_jm1_Ind = find(cell2mat(cellfun(@(sas)sas.TrackInd, NetObj.NodeList, 'uni', false))==j-1);

                % Get indices of all associated measurements for track j
                M_j = find(ValidationMatrix(j,:)); % i=0 for false alarm

                % For every node in L_jm1
                for i_jm1 = 1:size(L_jm1_Ind,2)

                    % Index of parent node
                    ParentInd = L_jm1_Ind(i_jm1);

                    % Get all measurements to consider
                    M_jm1 = union(1,setdiff(M_j,NetObj.NodeList{ParentInd}.MeasIndList)); 

                    % For every measurement in M_jm1
                    for i=1:size(M_jm1,2)
                        
                        if(toc > timeOut)
                            timedOut = 1;
                            return;
                        end
                        
                        % Get the track index
                        MeasInd = M_jm1(i);

                        % Init Match flag
                        match_flag = 0;

                        % Get index list L_j of nodes in current layer (j)
                        L_j_Ind = NetObj.NodesPerTrack{j};%find(dfg==j);

                        % if current Layer is empty
                        if isempty(L_j_Ind)
                            % Index of child node
                            ChildInd = size(NetObj.NodeList,2)+1;

                            % Create child node
                            NetObj.NodeList{ChildInd} = NodeObj;
                            NetObj.NodeList{ChildInd}.TrackInd = j;
                            NetObj.NodeList{ChildInd}.MeasIndList = union(NetObj.NodeList{L_jm1_Ind(i_jm1)}.MeasIndList(:,:), M_jm1(i));
                            NetObj.NodeList{ChildInd}.MeasIndList = NetObj.NodeList{ChildInd}.MeasIndList(:)';
                            NetObj.NodeList{ParentInd}.ChildIndList = union(NetObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                            NetObj.NodeList{ChildInd}.ParentIndList = union(NetObj.NodeList{ChildInd}.ParentIndList, ParentInd);

                            % Create edge from parent to child
                            NetObj.EdgeList{ParentInd, ChildInd} =[MeasInd];

                            % Compute remainders
                            M_rem_j = []; % T_j+1:mk[N^(j)_(i_j)]
                            for j_sub = j+1:LayerNum
                                M_rem_j = union(M_rem_j, find(ValidationMatrix(j_sub,:))');
                            end
                            NetObj.NodeList{ChildInd}.Remainders = setdiff(M_rem_j,setdiff(NetObj.NodeList{ChildInd}.MeasIndList,1)); 
                            NetObj.NodesPerTrack{j} = [NetObj.NodesPerTrack{j} ChildInd];
                        else

                            % Compute remainders (j-1)
                            M_rem_jm1_ti = []; %T_j+1:mk[N^(j-1)_(i_j-1),t_i]
                            for j_sub = j+1:LayerNum
                                M_rem_jm1_ti = union(M_rem_jm1_ti, find(ValidationMatrix(j_sub,:))');
                            end
                            R_jm1 = setdiff(M_rem_jm1_ti,setdiff(union(NetObj.NodeList{ParentInd}.MeasIndList, MeasInd),1));
                            R_jm1 = R_jm1(:); % Enforce that R_jm1 is a column vector

                            % For all nodes in current layer
                            for i_j=1:size(L_j_Ind,2)
                                ChildInd = L_j_Ind(i_j);

                                % Compute remainders (j)
                                R_j = NetObj.NodeList{ChildInd}.Remainders;

                                % If the node's list of remainders (R_j) is equal to R_jm1
                                if (isequal(R_jm1,R_j)) 
                                    % Simply add a new edge and update parent/child
                                    %   relationships
                                    NetObj.NodeList{ChildInd}.MeasIndList = union(NetObj.NodeList{ChildInd}.MeasIndList,MeasInd);
                                    if size(NetObj.EdgeList,1)<ParentInd
                                        NetObj.EdgeList{ParentInd, ChildInd} = MeasInd;
                                    else
                                        NetObj.EdgeList{ParentInd, ChildInd} = union(NetObj.EdgeList{ParentInd, ChildInd}, MeasInd);
                                    end
                                    NetObj.NodeList{ParentInd}.ChildIndList = union(NetObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                                    NetObj.NodeList{ChildInd}.ParentIndList = union(NetObj.NodeList{ChildInd}.ParentIndList, ParentInd);
                                    NetObj.NodeList{ChildInd}.MeasIndList = union(NetObj.NodeList{ChildInd}.MeasIndList, NetObj.NodeList{ParentInd}.MeasIndList);

                                    % set match_flag
                                    match_flag = 1;
                                    break;
                                end
                            end

                            if match_flag==0
                                % Index of child node
                                ChildInd = size(NetObj.NodeList,2)+1;

                                % Create child node
                                NetObj.NodeList{ChildInd} = NodeObj;
                                NetObj.NodeList{ChildInd}.TrackInd = j;
                                NetObj.NodeList{ChildInd}.MeasIndList = union(NetObj.NodeList{L_jm1_Ind(i_jm1)}.MeasIndList(:,:), M_jm1(i));
                                NetObj.NodeList{ChildInd}.MeasIndList = NetObj.NodeList{ChildInd}.MeasIndList(:)';
                                NetObj.NodeList{ParentInd}.ChildIndList = union(NetObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                                NetObj.NodeList{ChildInd}.ParentIndList = union(NetObj.NodeList{ChildInd}.ParentIndList, ParentInd);

                                % Create edge from parent to child
                                NetObj.EdgeList{ParentInd, ChildInd} =[MeasInd];

                                % Compute remainders
                                M_rem_j = []; % T_j+1:mk[N^(j)_(i_j)]
                                for j_sub = j+1:LayerNum
                                    M_rem_j = union(M_rem_j, find(ValidationMatrix(j_sub,:))');
                                end
                                R_j = setdiff(M_rem_j,setdiff(NetObj.NodeList{ChildInd}.MeasIndList,1));
                                NetObj.NodeList{ChildInd}.Remainders = R_j; 
                                NetObj.NodesPerTrack{j} = [NetObj.NodesPerTrack{j} ChildInd];
                            end
                        end
                    end
                end
            end
        end 
        
        function [AssocWeightsMatrix,NetObj] = computeAssocWeights(NetObj, LikelihoodMatrix)
        % COMPUTEASSOCWEIGHTS Computes the association weights
        %
        % DESCRIPTION:
        % * [AssocWeightsMatrix,NetObj]=computeAssocWeights(NetObj, LikelihoodMatrix)
        %
        % INPUTS:
        % * NetObj:             (Structure) with the following fields:
        %                       # NetObj.NodeList - List of all Nodes (NodeObj) 
        %                         contained in net.
        %                       # NetObj.EdgeList - Cell Matrix (Nn-by-Nn), Nn being 
        %                         the total number of nodes), where cell (i,j) contains 
        %                         the tracks contained in edge from parent node i, to
        %                         child node j.
        % * LikelihoodMatrix:   (T-by-M matrix) containing the measurement likelihood of
        %                       measurements, given the positions of all tracks.
        %                       (M: number of measurements, including dummy at index 1)
        %                       (T: number of tracks)
        %
        % OUTPUTS:
        % * AssocWeightsMatrix: (T-by-M matrix) containing the association weights
        %                       between tracks and measurements
        %                       (M: number of measurements, including dummy at index 1)
        %                       (T: number of tracks) 
        % * NetObj:             (Structure) with the following fields:
        %                       # NetObj.NodeList - List of all Nodes (NodeObj) 
        %                         contained in net.
        %                       # NetObj.EdgeList - Cell Matrix (Nn-by-Nn), Nn being 
        %                         the total number of nodes), where cell (i,j) contains 
        %                         the tracks contained in edge from parent node i, to
        %                         child node j.
        %
        %   IMPORTANT REMINDER TO SELF:
        %     When computing the remainders for a node, we always look at the remaining
        %     association events in the !NEXT! layer and then get the difference
        %     between these and any entries in the node's MeasIndList.
        %
        %   Author: Lyudmil Vladimirov
        
            % Get number of Tracks
            TrackNum = size(LikelihoodMatrix,1); 

            % Get number of Maasurements
            PointNum = size(LikelihoodMatrix,2);
           
            % Calculate the vector P_D
            p_D = zeros(size(NetObj.NodeList,2),1);
            p_D(1,1) = 1;
            % for every node
            for NodeInd=2:size(NetObj.NodeList,2)
                p_D_temp = 0;
                Node = NetObj.NodeList{NodeInd};

                % for every parent of the node
                for i = 1:size(Node.ParentIndList,2)

                    ParentInd = Node.ParentIndList(i);
                    p_D_m1 = p_D(ParentInd,1);

                    % for every measurement in the parent-child edge
                    MeasEdgeList = cell2mat(NetObj.EdgeList(ParentInd,NodeInd));
                    for t = 1:size(MeasEdgeList,2)
                        MeasInd = MeasEdgeList(t);
                        p_D_temp = p_D_temp + LikelihoodMatrix(Node.TrackInd, MeasInd)*p_D_m1;
                    end
                end
                p_D(NodeInd,1) = p_D_temp; 
            end
            NetObj.p_D = p_D; 

            % Calculate the vector P_U
            p_U = zeros(size(NetObj.NodeList,2),1);
            p_U(end,1) = 1;
            % for every node starting from the leaf
            for i=1:size(NetObj.NodeList,2)-1
                NodeInd = (size(NetObj.NodeList,2)) - i;
                p_U_temp = 0;
                Node = NetObj.NodeList{NodeInd};

                % for every child of the node
                for j = 1:size(Node.ChildIndList,2)
                    ChildInd = Node.ChildIndList(j);
                    ChildNode = NetObj.NodeList{ChildInd};
                    p_U_p1 = p_U(ChildInd,1);

                    % for every measurement in the parent-child edge
                    MeasEdgeList = cell2mat(NetObj.EdgeList(NodeInd, ChildInd));
                    for t = 1:size(MeasEdgeList,2)
                        MeasInd = MeasEdgeList(t);
                        p_U_temp = p_U_temp + LikelihoodMatrix(ChildNode.TrackInd, MeasInd)*p_U_p1;
                    end
                end
                p_U(NodeInd,1) = p_U_temp; 
            end
            NetObj.p_U = p_U;

            % Compute P_DT matrix
            p_DT = zeros(PointNum, size(NetObj.NodeList,2));
            for NodeInd=1:size(NetObj.NodeList,2)
                Node = NetObj.NodeList{NodeInd};
                for MeasInd = 1:PointNum

                    % Valid parents, where MeasInd belongs to edge
                    ValidParentIndList = [];
                    for j = 1:size(Node.ParentIndList,2)
                        ParentInd = Node.ParentIndList(j);
                        MeasEdgeList = cell2mat(NetObj.EdgeList(ParentInd, NodeInd));
                        if (ismember(MeasInd, MeasEdgeList)~=0)
                            ValidParentIndList = union(ValidParentIndList,ParentInd);
                        end
                    end 

                    for i = 1:size(ValidParentIndList,2)
                        ParentInd = ValidParentIndList(i);
                        p_D_m1 = p_D(ParentInd,1);
                        p_DT(MeasInd,NodeInd) = p_DT(MeasInd,NodeInd) + p_D_m1;
                    end
                end
            end
            NetObj.p_DT = p_DT;

            % Compute P_T matrix
            p_T = ones(PointNum, size(NetObj.NodeList,2));
            p_T(:,1) = zeros(PointNum,1);
            for NodeInd=2:size(NetObj.NodeList,2)
                Node = NetObj.NodeList{NodeInd};
                for MeasInd = 1:PointNum
                    p_T(MeasInd, NodeInd) = p_U(NodeInd)*LikelihoodMatrix(Node.TrackInd, MeasInd)*p_DT(MeasInd,NodeInd);
                end
            end
            NetObj.p_T = p_T;

            % Compute betta
            betta = zeros(TrackNum, PointNum);
            for TrackInd = 1:TrackNum
                % Get index list L_j of nodes in current Track layer (TrackInd)
                L_j_Ind = cell2mat(cellfun(@(sas)sas.TrackInd, NetObj.NodeList, 'uni', false ))==TrackInd;
                
                for MeasInd = 1:PointNum    
                    % Compute betta(j,i)
                    betta(TrackInd, MeasInd) = sum(p_T(MeasInd,L_j_Ind),2);
                end
            end

            % Normalise
            for j = 1:TrackNum
                betta(j,:) = betta(j,:)/sum(betta(j,:),2);
            end

            NetObj.betta = betta;
            AssocWeightsMatrix = betta;
        end
    end
 end
