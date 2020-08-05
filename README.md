# Eigenvector_Alignment

Eigenvector alignment uses the dominant system eigenvectors to assess structural changes in a network. This approach is applied to a dataset, included here as subjects.mat, of 29 subjects from three groups: Alzheimer's Disease, amnestic Mild Cognitive Impairment and Health Control. A two-sample t test is carried out to determine significant structural changes as described in the article "Eigenvector Alignment: assessing functional network changes in amnestic mild cognitive impairment and Alzheimer's disease" [1].

**Eig_Align.m** - assess eigenvector alignment between a selected node n and all other network nodes, using the system's dominant eigenvectors V.

**Ttest_Eig_Align_AMH.m** - perform Welch's t test comparing eigenvector alignment between three subject groups, where results are filtered by comparison with random models.

**Ttest_Eig_Align_Rand.m** - perform Welch's t test comparing eigenvector alignment for each subject group with a set number of random models of functional connectivity matrices.

### Data folder:
**subjects.mat** - the fMRI resting state data that is analysed in this work is from the 'Resting-state fMRI in Dementia Patients' dataset [2] (Harvard Dataverse).

**rand1000_(set).mat** data from three sets of comparisons with 1000 random models. Data used to filter t test results for subject groups comparison.

### References:
[1] Clark, R.A., Nikolova, N., McGeown, W.J. and Macdonald, M., 2020. Eigenvector alignment: assessing functional network changes in amnestic mild cognitive impairment and Alzheimer's disease. bioRxiv.

[2] Mascali D, DiNuzzo M, Gili T, Moraschi M, Fratini M, Maraviglia B, et al., 2015. Resting-state fMRI in dementia patients.  Harvard Dataverse.  Availablefrom: http://doi.org/10.7910/DVN/29352

[![DOI](https://zenodo.org/badge/271226567.svg)](https://zenodo.org/badge/latestdoi/271226567)
