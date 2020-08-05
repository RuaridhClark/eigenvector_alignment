function [M,t] = Process_Threshold(M)

    M(isnan(M))=0;      % set isnan entries to 0
    M=M(1:132,1:132);   % remove network nodes use only 132 ROIs
    [t] = CST(M,0);     % Cluster Span Threshold
    M(M<t)=0;
end