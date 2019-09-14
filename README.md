# High-Dimensional-Source-MEG-FC

These scripts and functions are a template. Some code is designed to run using a parallel cluster in Matlab (parfor). Frequently, analysis are run for all the subjects automatically. However, variables should be modify to reflect your directory path accordingly, and also your other information associated with your particular analyses.

NOTE: Requires SPM12 toolbox. We are not responsible for the use of this code. Users will have to modify it to their particular interests.

Keep tuned! We will make available a friendly GUI in a future research work!

CODE **** scADSpain_preproc_parfor.m - data preprocessing ****

1. Read the ELEKTA-Neuromag "fif" file and read it using SPM toobox, which generates an object containing the data.
2. Data is downsample to 200 Hz using an SPM function: "spm_eeg_downsample".
3. Individual mri image is read and meshes are extracted using SPM tools.
4. Anatomical and MEG spaces are coregistered using the semi-automatic tool "spm_coreg", which was programmed by the authors with the purpose of providing a more flexible coregistration.
5. The SPM MEG data object is saved with the fields updated for the coregistered case.

NOTE: Some variables have to be setup accordingly to read your data, for all the subjects

CODE **** spm_coreg.m - semi-automatic MRI-MEG coregistration tool ****

NOTE: It is used only inside "scADSpain_preproc_parfor" as part of the coregistration process.

CODE **** scADSpain_invsol_parfor.m - estimate the source activity ****

1. Reads the SPM12 object containing participants data.
2. Pass-band data filtering using SPM12 function "spm_eeg_filter".
3. Runs the Bayesian minimum norm method as implemented in SPM12 (COH) using the function "spm_eeg_invert", for 2s segmented data with applied Hanning tapering per segment.
4. FFT analysis per segment.
5. FFT coefficients saved to hard-disk for post-hoc FC analysis

CODE **** scADSpain_sourceFC_parfor.m - Run the source FC analysis ****

1. FC analysis per frequency or frequency band (uncomment code as suits your needs)
2. Estimate imaginary coherence (iCOH) and the EIC method, other methods can be inserted in the code as needed.
3. Save the outcome to hard-disk

CODE **** scADSpain_clustperm_indsel_parfor.m - Run the cluster-permutation analysis ****

This script is controlled by several flag variables to control the running of particular cell code. The flag variables are:

1. fcompute_block - to estimate the sparse supra-threshold FC indices for each FC block and frequency band.
2. itype - currently take values from 1 to 7 to compute above analysis for different methods (Wilcoxon, Spearman, standard correlation, etc.) and using different psychological tests (MMSE, IRM, DRM).
3. fcompute_stat - assemble the FC indices for all the blocks, separately for each frequency band, and save the outcome.
4. fvisualize_stat - it is settled to true in the code because we always show some part of the outcome at the end as an inspection of the results. This last part use the function "plot_bsurfconn_spm" to provide a visualization of the cortico-cortical FC map in the SPM12 template surface.

NOTE: Users have to specify the particular thresholds for the upper and lower distribution tails, correspondingly to the rank-sum, or correlation analyses, accordingly with its particular dataset and analysis.

CODE **** plot_bsurfconn_spm.m - Provide visualization of FC map for spm template surface ****

NOTE: It is used only inside "scADSpain_clustperm_indsel_parfor".

CODE **** plot_bsurfconn.m - Provide visualization of FC map for generic surfaces ****

NOTE: It is used only inside "plot_bsurfconn_spm".
