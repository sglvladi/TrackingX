function [DataList,newGroundTruth] = dataGen(ObsModel, GroundTruth, nTracks, lambda, V_bounds, ObscurRects, increment)
% DATAGEN Generate measurements from given a set of ground truth data and parameters
%
% INPUTS:
%   - 

    % Validate inputs
    if(nargin<7);   increment     = 1; end
    if(nargin<6);   ObscurRects   = []; end
    if(nargin<5);   V_bounds      = []; end
    if(nargin<4);   lambda        = 0; end
    
    
    % Initialise storage
    nTimesteps      = ceil(numel(GroundTruth)/increment);
    DataList        = cell(1,nTimesteps);
    newGroundTruth  = cell(1,nTimesteps);
    
    ki = 1;
    for k=1:increment:numel(GroundTruth)
        
        % Compute number of tracks and clutter measurements
        if(nargin<3); nTracks = size(GroundTruth(k),2); end
        nClutter = poissrnd(lambda);
        
        % Initialise measurement list for k-th timestep
        DataList{ki} = NaN(ObsModel.Params.yDim,nTracks+nClutter);
        newGroundTruth{ki} = NaN(ObsModel.Params.yDim,nTracks);
        
        % Generate true measurements
        for j=1:nTracks
            obscurred = false;
            for ij = 1:size(ObscurRects,1)
                if(GroundTruth{k}(1,j)>=ObscurRects(ij,1) && GroundTruth{k}(1,j)<=ObscurRects(ij,2) && GroundTruth{k}(2,j)>=ObscurRects(ij,3) && GroundTruth{k}(2,j)<=ObscurRects(ij,4))
                    obscurred = true;
                    break;
                end
            end
            if (~obscurred)
                DataList{ki}(:,j) = ObsModel.sample(1,[GroundTruth{k}(1,j); GroundTruth{k}(2,j)]);
                newGroundTruth{ki}(:,j) = GroundTruth{k}(:,j);
            end
        end
        
        % Generate clutter measurements
        for j=nTracks+1:nTracks+nClutter
            if(isa(ObsModel,'Polar2CartGaussModelX'))
                DataList{ki}(:,j) = [unifrnd(0,pi/2); unifrnd(0,14)];
            else
                DataList{ki}(:,j) = ObsModel.obs(0,[unifrnd(V_bounds(1),V_bounds(2)); unifrnd(V_bounds(3),V_bounds(4)); 0; 0]);
            end
        end
        ki = ki + 1;
    end
end