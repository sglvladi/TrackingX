function betta = LoopyBeliefPropagation(assocLikelihoods, convergenceThresh)
% LoopyBeliefPropagation - Implementation of the well known algorithm
%
% INPUTS:
%   assocLikelihoods    : A (T+1-by-M+1) table of association likelihoods for T
%                          targets and M measurements (index 1 stands for the dummy)%
%   convergenceThresh   : A convergence threshold (<< 1)

    [numTargets,numMeasurements] = size(assocLikelihoods); 
    numMeasurements = numMeasurements-1;
    
    muba = ones(numTargets,numMeasurements);
    muba_tilde = zeros(numTargets,numMeasurements);
    om = ones(1,numMeasurements);
    on = ones(1,numTargets);
    while max(max(abs(muba-muba_tilde)))>convergenceThresh        
        muba_tilde = muba;
        
        prodfact = muba .* assocLikelihoods(:,2:end);
        sumprod= assocLikelihoods(:,1)+ sum(prodfact,2);
        divider = sumprod(:,om) - prodfact;
        divider(divider==0) = eps;
        muab = assocLikelihoods(:,2:end) ./ divider;
        
        summuab= eps+sum(muab,1);
        divider =  summuab(on,:) - muab;
        divider(divider==0) = eps;
        muba = 1 ./ divider;
    end 
    betta = zeros(numTargets,numMeasurements+1);
    for i = 1:numTargets
        betta(i,1) = assocLikelihoods(i,1);
        betta(i,2:end) = assocLikelihoods(i,2:end).*muba(i,:);
    end
    betta = betta./sum(betta,2);
end