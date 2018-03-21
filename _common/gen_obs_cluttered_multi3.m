function [DataList, x, y] = gen_obs_cluttered_multi3(TrackNum, x_true, y_true, q , lambda, increment)

    DataList = {};
    k=1;
    x = [];
    y = [];
    for i=1:increment:size(x_true,1)
        tempDataList=[];
        for j=1:TrackNum
            x(k,j) = x_true(i,j);
            y(k,j) = y_true(i,j);
            if(x(k,j)==0&&y(k,j)==0)
                tempDataList(:,j) = [ 0, 0];
            else
                tempDataList(:,j) = mvnrnd([x_true(i,j), y_true(i,j)],[q^2,0;0,q^2]);
            end
            %tempDataList(:,j) = mvnrnd([x_true(i,j), y_true(i,j)],[q^2,0;0,q^2]);
        end
        No_of_clutter = poissrnd(lambda);
        for j=TrackNum+1:No_of_clutter
            tempDataList(:,j) = [unifrnd(0,10), unifrnd(0,10)];
        end
        %tempDataList(:,randperm(size(tempDataList,2)));
        DataList{k} = tempDataList;
        k=k+1;
    end
end