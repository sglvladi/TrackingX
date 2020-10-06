function [gamma_low, gamma_high] = llr_thresholds(P_TM,P_FC)
%LLR_THRESHOLDS Summary of this function goes here
%   Detailed explanation goes here
    gamma_low = log(P_TM/(1-P_FC));
    gamma_high = log((1-P_TM)/P_FC);
end

