P_TM = 0.1;
P_FC = 10^-9;

gamma_low = log(P_TM/(1-P_FC));
gamma_high = log((1-P_TM)/P_FC);
% [gamma_low, gamma_high] = llr_thresholds(P_TM, P_FC);