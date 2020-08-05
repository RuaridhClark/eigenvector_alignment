% Perform welch's two-sample t test comparing eigenvector alignment between three subject
% groups: Alzheimer's Disease (1-9), amnestic Mild Cognitive Impairment
% (10-19), Health Control (20-29)
% Output:
% AH_pairs :- significant ROI pairs alongside their p value and mean align. change for AD vs HC
% MH_pairs :- significant ROI pairs alongside their p value and mean align. change for aMCI vs HC
% AM_pairs :- significant ROI pairs alongside their p value and mean align. change for AD vs aMCI
% tt_A_H :- p values for eig. align. changes in AD versus HC
% tt_A_M :- p values for eig. align. changes in AD versus aMCI
% tt_M_H :- p values for eig. align. changes in aMCI versus HC 
% roi_signif_sum :- number of significant eig. align. changes for each ROI in each comparison (AH,MH,AM) 

addpath([cd,'\functions'],[cd,'\data'])          % add functions folder to path

sig =0.05;          % define significance threshold
[AH_pairs,MH_pairs,AM_pairs,tt_A_H,tt_M_H,tt_A_M,roi_signif_sum] = Ttest_Filtered(sig);

function [AH_pairs,MH_pairs,AM_pairs,tt_A_H,tt_M_H,tt_A_M,roi_signif_sum] = Ttest_Filtered(sig)
    
    [tt_A_H,tt_A_M,tt_M_H,ea_A_H,ea_A_M,ea_M_H] = Ttest_Fnctn();

    loadnames={'rand1000_1.mat','rand1000_2.mat','rand1000_3.mat'}; % three sets of 1000 random model comparisons
    [ad,mci,hc,~] = Filter_Selections(loadnames,sig);               % filter ROI pair selections 

    load('subjects.mat','roi_names');
    cmprsnnames=[{'AD'},{'HC'}];
    [AH_pairs,tt_A_H] = Threshold_Filter(tt_A_H,ea_A_H,cmprsnnames,roi_names,ad,mci,hc,sig);

    cmprsnnames=[{'MC'},{'HC'}];
    [MH_pairs,tt_M_H] = Threshold_Filter(tt_M_H,ea_M_H,cmprsnnames,roi_names,ad,mci,hc,sig);

    cmprsnnames=[{'AD'},{'MC'}];
    [AM_pairs,tt_A_M] = Threshold_Filter(tt_A_M,ea_A_M,cmprsnnames,roi_names,ad,mci,hc,sig);

    roi_signif_sum=[(1:132)',sum(logical(tt_A_H))',sum(logical(tt_M_H))',sum(logical(tt_A_M))'];
end

function [tt_A_H,tt_A_M,tt_M_H,ea_A_H,ea_A_M,ea_M_H] = Ttest_Fnctn()

    load('subjects.mat','Z')        % contains 29 subjects: AD group ID 1-9; aMCI group ID 10-19; HC group ID 20-29
    AD_id=1:9; aMCI_id=10:19; HC_id = 20:29;

    n_eigs = 3;         % number of dominant eigenvectors used for eigenvector alignment
    n_subjects = 29;    % number of subjects

    V_save = cell(n_subjects,1);
    for i = 1:n_subjects
        adj=Z(:,:,i);
        [adj,~] = Process_Threshold(adj);   % Include only ROIs and apply Cluster Span Threshold
        [V] = Ordered_Eigvecs(adj,n_eigs);  % Select dominant eigenvectors
        V_save{i}=V;
    end

    n_nodes = length(adj);
    EA_all=zeros(n_nodes,n_subjects);       % Initialise variables
    tt_A_H = zeros(n_nodes,n_nodes); tt_A_M = tt_A_H; tt_M_H = tt_A_H;
    ea_A_H = zeros(n_nodes,n_nodes); ea_A_M = ea_A_H; ea_M_H = ea_A_H;
    for i = 1:n_nodes
        for j=1:n_subjects
            V = V_save{j};
            [EA] = Eig_Align(V,i);      % Assess eigenvector alignment for node i w.r.t. all other nodes
            EA_all(:,j)=EA;
        end

        AD = EA_all(:,AD_id);           % Eig. Align. for AD subjects
        aMCI = EA_all(:,aMCI_id);       % Eig. Align. for aMCI subjects
        HC = EA_all(:,HC_id);           % Eig. Align. for HC subjects

        for nn = 1 : n_nodes
            if i ~= nn
                [~,p,~,~] = ttest2(AD(nn,:),HC(nn,:),'Vartype','unequal');  % Matlab in-built ttest function
                tt_A_H(nn,i)=p;                             % p values AD versus HC comparison
                ea_A_H(nn,i)=mean(AD(nn,:))-mean(HC(nn,:)); % mean eig. align. change

                [~,p,~,~] = ttest2(aMCI(nn,:),HC(nn,:),'Vartype','unequal');
                 tt_M_H(nn,i)=p;                                % p values aMCI versus HC comparison
                 ea_M_H(nn,i)=mean(aMCI(nn,:))-mean(HC(nn,:));  % mean eig. align. change

                [~,p,~,~] = ttest2(AD(nn,:),aMCI(nn,:),'Vartype','unequal');
                tt_A_M(nn,i)=p;                                 % p values AD versus aMCI comparison
                ea_A_M(nn,i)=mean(AD(nn,:))-mean(aMCI(nn,:));   % mean eig. align. change
            end
        end
    end
end
