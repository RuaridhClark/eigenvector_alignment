% function to apply threshold filter for both subject group and random
% model comparisons. Outputs list of significant ROI pairs (pairs) and filtered p
% values from Welch's t-test (tt).
function [pairs,tt] = Threshold_Filter(tt,ea,cmprsnnames,roi_names,ad,mci,hc,sig)

[m_allow] = Reduce_By_Rand_Signif(cmprsnnames,ad,mci,hc);   % filter by significance against random models
tt(m_allow==0)=0;
tt(tt>sig)=0;

[row,col]=find(tt>0);
pairs=[roi_names(row),roi_names(col)];  % add significant ROI pair names
for i = 1 : length(row)
    pairs(i,3)={tt(row(i),col(i))};     % add p values for selected
    pairs(i,4)={ea(row(i),col(i))};     % add eig. align. values for selected
end