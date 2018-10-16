function [DataList, x, y] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, q , V_bounds, lambda, increment, Pd)
% gen_obs_cluttered_multi2 Generate measurements given ground truth data, measurements noise and clutter model. 
%
% INPUT ARGUMENTS:
% * TrackNum            (scalar) The number of true tracks present to
%                       detect. This is put in place to ensure that only a
%                       given number of targets is present in each scan.
% * x_true, y_true      x,y coordinates of true trajectories. each is a
%                       (NumTracks x NumIter) matrix
% * q                   Measurement noise standard deviation
% * V_bounds            The size of the surveillance region. Assumed to be
%                       a rectangle of the form [xmin,xmax,ymin,ymax].
%                       Gets used when sampling clutter.
% * lambda              The mean number of clutter measurements (Poisson)
%                       The actual number of clutter measurements is drawn
%                       from Pois(lambda)
% * increment           Can be used to skip measurements, such as to
%                       simulate a less responsive sensor. When set to 1,
%                       the measurements will be generated at the normal
%                       rate. Set to 2,3,4... to get every 2nd,3rd,4th...
%                       scan.
% * Pd                  Probability of detection. Assumed to be the same
%                       for all targets.
%   
%
%  See also JointProbabilisticDataAssocX/associate, JointProbabilisticDataAssocX/updateTracks.  
    DataList = {};
    k=1;
    x = [];
    y = [];
    for i=1:increment:size(x_true,1)
        tempDataList=[];
        for j=1:TrackNum
                p = rand;
                if(p<=Pd)
                    tempDataList(:,j) = mvnrnd([x_true(i,j), y_true(i,j)],[q^2,0;0,q^2]);
                    x(k,j) = x_true(i,j);
                    y(k,j) = y_true(i,j);
                else
                    j = j -1;
                end
        end
        No_of_clutter = poissrnd(lambda);
        for j=size(tempDataList,2)+1:No_of_clutter
            tempDataList(:,j) = [unifrnd(V_bounds(1),V_bounds(2)), unifrnd(V_bounds(3),V_bounds(4))];
        end
        tempDataList(:,randperm(size(tempDataList,2)));
        DataList{k} = tempDataList;
        k=k+1;
    end
end