# clinical_data_grouping

`evaluateUniformity.ipynb` used to identify T1w scan groups in CHOP clinical data for 3T. Process is done for T2w and FLAIR data with the same process as well as for 1.5T T1w data. The lab's [wiki](https://bgdlab.github.io/research/scan_identification_wiki.html) explains the process more in depth.

`identifyFilesFromCuBIDS.py` used to edit CuBIDS output so that files are renamed like: sub-HM3WPPESR_ses-439563557339procId006437ageDays_**acq-TSEStandardized3TScannerId45886Slices98Voxel0p898mm**_rec-NORM_run-001_T2w.nii.gz where scan group (e.g. TSE, MPRAGE, TIRM), magnetic field strength, number of slices, and voxel size are listed in the *acq* BIDS key
