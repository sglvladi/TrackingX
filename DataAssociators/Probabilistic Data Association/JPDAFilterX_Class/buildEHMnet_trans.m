%   buildEHMnet_trans.m                              Author: Lyudmil Vladimirov
%   ======================================================================>
%   Functionality: 
%       Build Track-Oriented (TO) EHM net and compute respective 
%       association probabilities (betta).
%       (normally executed for each cluster)
%
%   Input:
%       ValidationMatrix    - Matrix (T x M) containing all possible measurement
%                             to track associations (for cluster of interest).
%                             (Output from ObservationAssociation.m)
%       Li                  - Matrix (T x M) containing association likelihoods
%       (M: number of measurements, including dummy at index 1)
%       (T: number of tracks)
%
%   Output:
%       Structure NetObj:
%           #NetObj.NodeList - List of all Nodes (NodeObj) contained in net
%           #NetObj.EdgeList - Cell Matrix (n x n, n being the total number of
%                              nodes), where cell (i,j) contains the tracks
%                              contained in edge from parent node i, to
%                              child node j.
%           #NetObj.ValidationMatrix - Internal Validation matrix for
%                                      cluster.
%           #NetObj.p_D      - Computed p_D matrix
%           #NetObj.p_U      - Computed p_U matrix
%           #NetObj.p_DT     - Computed p_DT matrix
%           #NetObj.p_T      - Computed p_T matrix
%           #NetObj.betta    - Matrix (T x M) containing Association  
%                              probabilites computed for all T tracks and
%                              M measurements.
%
%   IMPORTANT REMINDER TO SELF:
%     When computing the remainders for a node, we always look at the remaining
%     association events in the !NEXT! layer and then get the difference
%     between these and any entries in the node's MeasIndList.
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


function NetObj = buildEHMnet_trans(ValidationMatrix, Li)

    % Get number of tracks/layers
    TrackNum = size(ValidationMatrix,1); 
    LayerNum = TrackNum; % Layer 1 is root layer

    % Augment ValidationMatrix
    %ValidationMatrix = [ones(1, TrackNum); ValidationMatrix']';

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
    NetObj.ValidationMatrix = ValidationMatrix;
    NetObj.Li = Li;

    % Create Root Node
    NetObj.NodeList{1} = NodeObj;

    % Build the net
    % ==============>
    % For every measurement/layer
    for j=1:LayerNum

        % Get index list (L_jm1_Ind) of nodes in previous layer (j-1)
        L_jm1_Ind = find(cell2mat(cellfun(@(sas)sas.TrackInd, NetObj.NodeList, 'uni', false))==j-1);

        % Get indeces of all associated measurements for track j
        M_j = [find(ValidationMatrix(j,:))]; % i=0 for false alarm

        % For every node in L_jm1
        for i_jm1 = 1:size(L_jm1_Ind,2)

            % Index of parent node
            ParentInd = L_jm1_Ind(i_jm1);

            % Get all measurements to consider
            M_jm1 = union(1,setdiff(M_j,NetObj.NodeList{ParentInd}.MeasIndList)); 

            % For every measurement in M_jm1
            for i=1:size(M_jm1,2)

                % Get the track index
                MeasInd = M_jm1(i);

                % Init Match flag
                match_flag = 0;

                dfg = cellfun(@(sas)sas.TrackInd, NetObj.NodeList, 'uni', false );
                dfg = cell2mat(dfg);

                % Get index list L_j of nodes in current layer (j)
                L_j_Ind = find(dfg==j);

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
                    end
                end

            end
        end
    end 

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
                p_D_temp = p_D_temp + Li(Node.TrackInd, MeasInd)*p_D_m1;
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
                p_U_temp = p_U_temp + Li(ChildNode.TrackInd, MeasInd)*p_U_p1;
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
            p_T(MeasInd, NodeInd) = p_U(NodeInd)*Li(Node.TrackInd, MeasInd)*p_DT(MeasInd,NodeInd);
        end
    end
    NetObj.p_T = p_T;

    % Compute betta
    betta = zeros(TrackNum, PointNum);
    for TrackInd = 1:TrackNum
        for MeasInd = 1:PointNum
            % Get index list L_j of nodes in current measurement layer (TrackInd)
            L_j_Ind = find(cell2mat(cellfun(@(sas)sas.TrackInd, NetObj.NodeList, 'uni', false ))==TrackInd);
            for j = 1:size(L_j_Ind, 2)
                NodeInd = L_j_Ind(j);
                betta(TrackInd, MeasInd) = betta(TrackInd, MeasInd) + p_T(MeasInd,NodeInd);%p_U(NodeInd)*Li(TrackInd, TrackInd)*p_DT(TrackInd, NodeInd);
            end
        end
    end
    
    %if(betta==

    % Normalise
    for j = 1:TrackNum
        betta(j,:) = betta(j,:)/sum(betta(j,:),2);
    end

    % Compute betta transpose
    betta_rescaled = betta./(betta(:,1)*ones(1,PointNum));
    betta_rescaled_reduced = betta_rescaled(:, 2:end);
    betta_modified = [ones(1, PointNum-1)*PointNum-1; betta_rescaled_reduced];
    betta_modified_transposed = betta_modified';
    betta_transpose = betta_modified_transposed./(sum(betta_modified_transposed,2)*ones(1,TrackNum+1));

    NetObj.betta = betta;
    NetObj.betta_trans = betta_transpose;
end