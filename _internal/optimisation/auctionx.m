function [assigns, totvalue] = auction(valuematrix)

% [assigns, totvalue] = auction(valuematrix)
%
% First column is null assignment - any number of objects can be assigned
% this
%
% Author: Paul Horridge

bidincrement = 1e-6; % choose better later
[nobj, nhyps] = size(valuematrix); % nhyps includes null hypothesis
prices = zeros(1, nhyps);
assigns = ones(nobj, 1);
hyps2obj = zeros(1, nhyps);
happy = false(nobj, 1);

% Repeat while someone unhappy
while any(~happy)
    % Find unhappy person
    i = find(~happy, 1, 'first');
    % Get their net values
    netvals = valuematrix(i,:) - prices;
    [bestval, bestj] = max(netvals);
    if netvals(assigns(i)) >= bestval - bidincrement
        % Check if happy
        happy(i) = true;
    else
        % Get second best hypothesis
        nextval = max(netvals([1:bestj-1 bestj+1:end]));
        % Increase price if required
        if bestj > 1
            prices(bestj) = prices(bestj) + bestval - nextval + bidincrement;
            % De-assign original holder if required
            oldi = hyps2obj(bestj);
            if oldi > 0
                assigns(oldi) = 1;
                happy(oldi) = false;        
            end
        end
        % Assign best hypothesis to i
        assigns(i) = bestj;
        hyps2obj(bestj) = i;
        happy(i) = true;
    end
end
totvalue = sum(valuematrix((assigns'-1)*nobj + (1:nobj)));
