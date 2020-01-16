## Eye movement and blink-related artifact correction in EEG and MEG data

This repository contains reference implementations of 
* the sparse generalized eye artifact subspace subtraction (SGEYESUB) algorithm presented in [1].
* 4 other eye artifact correction algorithms presented in [2-5].

After the algorithm parameters are fitted to calibration data, eye movement and blink-related eye artifacts can be corrected offline and online.
In [1,2] calibration data was recorded using a visually guided paradigm. A reference implementation using [Psychtoolbox](http://psychtoolbox.org/) and [labstreaminglayer](https://github.com/sccn/labstreaminglayer) is provided in the subfolder `paradigm`.

This repository comes also with a demonstration dataset containing electroencephalographic (EEG) and electrooculographic (EOG) activity of one person. The data is stored in [eeglab](https://sccn.ucsd.edu/eeglab/index.php) format.

#### Public EEG dataset

The pre-processed EEG dataset investigated in [1] is [publicly available](https://osf.io/2qgrd/) on OSF [6].

### Getting started
* Download or clone the repository
* Startup the [eeglab](https://sccn.ucsd.edu/wiki/Chapter_01:_Loading_Data_in_EEGLAB#Installing_EEGLAB_and_tutorial_files) toolbox
* Open the `demo_main.m` script. The script loads a calibration dataset (`demo_trainset.set`) and an evaluation dataset (`demo_testset.set`). 
* Both demo datasets contain continuous recordings. Before the algorithms are fitted and evaluated, the datasets are pre-processed in the script `demo_preprocessing.m`
* The detailed pre-processing steps are presented in [1]. 
* Next an object of the algorithm is created with `algo = sgeyesub()`
* The object is fitted to the calibration data `algo.fit(X_trn, y_trn, eeg_chan_idxs)` where `X_trn` and `y_trn` contain the M/EEG (and EOG) signals and labels.
* New samples (data) are corrected with `x_corrected = algo.apply(x)`.

### References

[1] Kobler, R. J., Sburlea, A. I., Lopes-Dias, C., Schwarz, A., Hirata, M. and Müller-Putz, G. R. "Corneo-retinal-dipole and eyelid-related eye artifacts can be corrected offline and online in electroencephalographic and magnetoencephalographic signals.", submitted

[2] Kobler, R. J., Sburlea, A. I., and Müller-Putz G.R., "A Comparison of Ocular Artifact Removal Methods for Block Design Based Electroencephalography Experiments." In Proceedings of the 7th Graz Brain-Computer Interface Conference, 236–41, 2017. https://doi.org/10.3217/978-3-85125-533-1-44

[3] Schlögl, A., Keinrath, C., Zimmermann, D., Scherer, R., Leeb, R., and Pfurtscheller, R. "A Fully Automated Correction Method of EOG Artifacts in EEG Recordings." Clinical Neurophysiology 118, no. 1 (2007): 98–104. https://doi.org/10.1016/j.clinph.2006.09.003

[4] Plöchl, M., Ossandón, J. P., and König P. "Combining EEG and Eye Tracking: Identification, Characterization, and Correction of Eye Movement Artifacts in Electroencephalographic Data." Frontiers in Human Neuroscience 6, (2012): 1–23. https://doi.org/10.3389/fnhum.2012.00278

[5] Zhou, X., Gerson, A. D., Lucas C Parra, L. C., and Paul Sajda, P. "EEGLAB Plugin EYESUBTRACT," (2005). Retrieved from http://sccn.ucsd.edu/eeglab/plugins/eyesubtract1.0.zip

[6] Kobler, R. J., Sburlea, A. I., Lopes-Dias, C., Schwarz, A., Mondini, V., and Müller-Putz, G. R. "EEG eye artifact dataset." (2020) Retrieved from https://osf.io/2qgrd 

### Acknowledgements
This work was supported by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (Consolidator Grant 681231 'Feel Your Reach').
