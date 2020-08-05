% Perform two-sample t test comparing eigenvector alignment between three subject
% groups: Alzheimer's Disease (1-9), amnestic Mild Cognitive Impairment
% (10-19), Health Control (20-29)
% Output:
% tt_A_R :- p values for eig. align. changes in AD versus random
% tt_M_R :- p values for eig. align. changes in aMCI versus random
% tt_H_R :- p values for eig. align. changes in HC versus random 
% ea_A_R :- mean eig. align. changes in AD versus random
% ea_M_R :- mean eig. align. changes in aMCI versus random
% ea_H_R :- mean eig. align. changes in HC versus random 

addpath([cd,'\functions'],[cd,'\data'])          % add functions folder to path

n_rand = 10;  % number of randomly generated connectivity matrices used in comparison
[tt_A_R,tt_M_R,tt_H_R,ea_A_R,ea_M_R,ea_H_R] = Ttest_Fnctn(n_rand);  % perform t-test

function [tt_A_R,tt_H_R,tt_M_R,ea_A_R,ea_H_R,ea_M_R] = Ttest_Fnctn(n_rand)

    [rand_Z] = Create_Rand_Adj(n_rand); % create n_rand random connectivity matrices
    load('Data\subjects.mat')           % contains 29 subjects: AD group ID 1-9; aMCI group ID 10-19; HC group ID 20-29
    AD_id=1:9; MCI_id=10:19; HC_id = 20:29; rand_id = 30:(29+n_rand);
    
    Z=Z(1:132,1:132,:);
    Z(:,:,rand_id)=rand_Z;         % add random dataset

    n_eigs = 3;                         % number of dominant eigenvectors used for eigenvector alignment
    n_subjects = (29+n_rand);           % number of subjects

    V_save = cell(n_subjects,1);
    for i = 1:n_subjects
        adj=Z(:,:,i);
        [adj,~] = Process_Threshold(adj);	% Include only ROIs and apply Cluster Span Threshold
        [V] = Ordered_Eigvecs(adj,n_eigs);  % Select dominant eigenvectors
        V_save{i}=V;
    end

    n_nodes = length(adj);
    EA_all=zeros(n_nodes,n_subjects);       % Initialise variables
    tt_A_R = zeros(n_nodes,n_nodes); tt_H_R = tt_A_R; tt_M_R = tt_A_R;
    ea_A_R = zeros(n_nodes,n_nodes); ea_H_R = ea_A_R; ea_M_R = ea_A_R;
    for i = 1:n_nodes
        for j=1:n_subjects
            V = V_save{j};
            [EA] = Eig_Align(V,i);          % Assess eigenvector alignment for node i w.r.t. all other nodes
            EA_all(:,j)=EA;
        end

        AD = EA_all(:,AD_id);               % Eig. Align. for AD subjects
        aMCI = EA_all(:,MCI_id);            % Eig. Align. for aMCI subjects
        HC = EA_all(:,HC_id);               % Eig. Align. for HC subjects
        Rand = EA_all(:,rand_id);           % Eig. Align. for random matrices

        for nn = 1 : n_nodes
            if i ~= nn
                [~,p,~,~] = ttest2(AD(nn,:),Rand(nn,:),'Vartype','unequal');  % Matlab in-built Welch's t-test function
                tt_A_R(nn,i)=p;                                 % p values AD versus random comparison
                ea_A_R(nn,i)=mean(AD(nn,:))-mean(Rand(nn,:)); 	% mean eig. align. change 
                
                [~,p,~,~] = ttest2(aMCI(nn,:),Rand(nn,:),'Vartype','unequal');
                tt_M_R(nn,i)=p;                                 % p values aMCI versus random comparison
                ea_M_R(nn,i)=mean(aMCI(nn,:))-mean(Rand(nn,:)); % mean eig. align. change
                
                [~,p,~,~] = ttest2(HC(nn,:),Rand(nn,:),'Vartype','unequal');
                tt_H_R(nn,i)=p;                                 % p values HC versus random comparison
                ea_H_R(nn,i)=mean(HC(nn,:))-mean(Rand(nn,:));  	% mean eig. align. change
            end
        end
    end
end
