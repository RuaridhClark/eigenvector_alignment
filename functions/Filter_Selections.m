% function generates binary matrices of ROI pairs that produce significant
% eigenvector alignment distributions when compared to random models.
% Comparisons with random models are loaded using loadnames to define file
% name.
function [ad,mci,hc,prcnt_same] = Filter_Selections(loadnames,sig)
smA = zeros(length(loadnames),1); smM = smA; smH = smA;
for i = 1 : length(loadnames)
    load(loadnames{i},'save_ttAD','save_ttMCI','save_ttHC')             % load p values from random model comparison
    tt_A_R=save_ttAD{1}; tt_M_R=save_ttMCI{1}; tt_H_R=save_ttHC{1}; 
    tt_A_R(tt_A_R>sig)=0; tt_M_R(tt_M_R>sig)=0; tt_H_R(tt_H_R>sig)=0;   % remove values above threshold 'sig'
    tt_A_R(tt_A_R>0)=1; tt_M_R(tt_M_R>0)=1; tt_H_R(tt_H_R>0)=1;         % set nonzero entries to 1
    smA(i) = sum(sum(tt_A_R)); smM(i) = sum(sum(tt_M_R)); smH(i) = sum(sum(tt_H_R));% number of significant ROI pairs x2
    ad = tt_A_R; mci = tt_M_R; hc = tt_H_R;
end

prcnt_same(1)=(sum(sum(abs(ad))))/mean(smA);    % percent of consistent ROI pairs from different random model set comparisons 
prcnt_same(2)=(sum(sum(abs(mci))))/mean(smM);
prcnt_same(3)=(sum(sum(abs(hc))))/mean(smH);