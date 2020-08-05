function [t] = CST(Adj,plot)
% [t] = CST(Adj,plot) is the Cluster Span Threshold, which generates an
% unbiased threshold for an adjacency matrix that is based on a clustering
% coefficient, C, that for a given network balances the number of triples
% that are clustered, forming loops, with those that are spanning, forming
% trees.

% Reference: Smith, Keith, et al. "Cluster-span threshold: An unbiased
% threshold for binarising weighted complete networks in functional
% connectivity analysis." 2015 37th Annual International Conference of the
% IEEE Engineering in Medicine and Biology Society (EMBC). IEEE, 2015.

c=[]; G=[]; C=[]; C_B=[];

t_all=0:0.01:1;
j=0;
for i = t_all        % Simple equation solver that checks between 0 an 1
    j=j+1;
    Adjhigh=Adj;
    Adjhigh(Adj<i)=0;
    Adjhigh(Adj>=i)=1;
    c(j)=num_conn_triples(Adjhigh);
    G(j) = loops3(Adjhigh);
    C(j) = G(j)/c(j);              
    C_B(j) = 4*C(j)*(1-C(j));    
end

[max_C,I] = max(C_B);   % Find t that best balances the number of clustered and spanning triples
t=t_all(I);                 % set threshold
