# eigenvector_alignment

Eigenvector alignment uses the dominant system eigenvectors to assess structural changes in a network. This approach is applied to a dataset, included here as subjects.mat, of 29 subjects from three groups: Alzheimer's Disease, amnestic Mild Cognitive Impairment and Health Control. A two-sample t test is carried out to determine significant structural changes as described in the article "Eigenvector Alignment: assessing functional network changes in amnestic mild cognitive impairment and Alzheimer's disease".

**eig_align.m** - assess eigenvector alignment between a selected node n and all other network nodes, using the system's dominant eigenvectors V.

**ttest_eig_align.m** - perform two-sample t test comparing eigenvector alignment between three subject groups.

**subjects.mat** - the fMRI resting state data that is analysed in this work is from the 'Resting-state fMRI in Dementia Patients' dataset [2] (Harvard Dataverse).

[1] Clark, R.A., Nikolova, N., McGeown, W.J. and Macdonald, M., 2020. Eigenvector alignment: assessing functional network changes in amnestic mild cognitive impairment and Alzheimer's disease. bioRxiv.

[2] Mascali D, DiNuzzo M, Gili T, Moraschi M, Fratini M, Maraviglia B, et al., 2015. Resting-state fMRI in dementia patients.  Harvard Dataverse.  Availablefrom: http://doi.org/10.7910/DVN/29352
