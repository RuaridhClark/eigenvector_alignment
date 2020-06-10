% Perform two-sample t test comparing eigenvector alignment between three subject
% groups: Alzheimer's Disease (1-9), amnestic Mild Cognitive Impairment
% (10-19), Health Control (20-29)
[tt_A_H,tt_A_M,tt_M_H] = ttest_fnctn();

function [tt_A_H,tt_A_M,tt_M_H] = ttest_fnctn()
    addpath([cd,'\functions'])  % add functions folder to path
    load('subjects.mat')        % contains 29 subjects: AD group ID 1-9; aMCI group ID 10-19; HC group ID 20-29
    AD_ID=1:9; aMCI_ID=10:19; HC_ID = 20:29;

    n_eigs = 3;         % number of dominant eigenvectors used for eigenvector alignment
    n_subjects = 29;    % number of subjects

    V_save = cell(n_subjects,1);
    for i = 1:n_subjects
        Adj=Z(:,:,i);
        [Adj] = process_threshold(Adj);     % Include only ROIs and apply Cluster Span Threshold
        [V] = ordered_eigvecs(Adj,n_eigs);  % Select dominant eigenvectors
        V_save{i}=V;
    end

    n_nodes = length(Adj);
    EA_all=zeros(n_nodes,n_subjects);	% Initialise variables
    tt_A_H = zeros(n_nodes,n_nodes); tt_A_M = tt_A_H; tt_M_H = tt_A_H;
    for i = 1:n_nodes
        for j=1:n_subjects
            V = V_save{j};
            [EA] = eig_align(V,i);      % Assess eigenvector alignment for node i w.r.t. all other nodes
            EA_all(:,j)=EA;
        end

        AD = EA_all(:,AD_ID);           % Eig. Align. values for node i
        aMCI = EA_all(:,aMCI_ID);
        HC = EA_all(:,HC_ID);

        for nn = 1 : n_nodes
            if i ~= nn
                [h,p,~,~] = ttest2(AD(nn,:),HC(nn,:));  % Matlab in-built ttest function
                if h == 1                               % Only record significant results
                    tt_A_H(nn,i)=p;                     % AD versus HC comparison
                end
                [h,p,~,~] = ttest2(AD(nn,:),aMCI(nn,:));
                if h == 1
                    tt_A_M(nn,i)=p;                     % AD versus aMCI comparison
                end
                [h,p,~,~] = ttest2(aMCI(nn,:),HC(nn,:));
                if h == 1
                    tt_M_H(nn,i)=p;                     % aMCI versus HC comparison
                end
            end
        end
    end
end
