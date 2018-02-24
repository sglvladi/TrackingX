function [DataList, x, y] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, q , V_bounds, lambda, increment)
%     q = 0.5
%     q_clutter = 8;
    DataList = {};
    k=1;
    x = [];
    y = [];
    for i=1:increment:size(x_true,1)
        tempDataList=[];
        for j=1:TrackNum
%             if( (x_true(i,j)>=1&&x_true(i,j)<=3&&y_true(i,j)>=3&&y_true(i,j)<=5) || (x_true(i,j)>=4&&x_true(i,j)<=6&&y_true(i,j)>=6&&y_true(i,j)<=7) || (x_true(i,j)>=8.75&&x_true(i,j)<=9.5&&y_true(i,j)>=4&&y_true(i,j)<=7))
%             %if(x(k,j)==0&&y(k,j)==0)
%                 tempDataList(:,j) = [ NaN, NaN];
%                 x(k,j) = NaN;
%                 y(k,j) = NaN;
%             else
                tempDataList(:,j) = mvnrnd([x_true(i,j), y_true(i,j)],[q^2,0;0,q^2]);
                x(k,j) = x_true(i,j);
                y(k,j) = y_true(i,j);
%             end
            %tempDataList(:,j) = mvnrnd([x_true(i,j), y_true(i,j)],[q^2,0;0,q^2]);
        end
        No_of_clutter = poissrnd(lambda);
        for j=TrackNum+1:No_of_clutter
            tempDataList(:,j) = [unifrnd(V_bounds(1),V_bounds(2)), unifrnd(V_bounds(3),V_bounds(4))];
        end
        tempDataList(:,randperm(size(tempDataList,2)));
        DataList{k} = tempDataList;
        k=k+1;
    end
end