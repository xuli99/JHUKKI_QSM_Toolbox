# JHU/KKI_QSM_Toolbox

## Version 3.3

*************************************************************
Previous public version:
JHUKKI_QSMToolbox_v3.0 is available @ http://godzilla.kennedykrieger.org/QSM/
with some example 3T data and a basic manual


JHU/KKI_QSM_Toolbox is a MATLAB software with GUI for doing Quantitative Susceptibility Mapping (QSM) processing from MR phase data acquired with gradient echo sequence (GRE). Command line version without GUI was added after version 3.3.
*************************************************************

## Installation
1.This is a MATLAB based toolbox, need MATLAB installation
2.This toolbox also uses processing routine from FSL, e.g. bet, thus need to have FSL installed.

## Usage 
1. To start, open MATLAB and go to the folder where this toolbox is located.
2. (a) run QSM.m for the GUI
   (b) To use the command line version, run QSM_cluster(ParamSetFile.m, LogFile.txt)
   where ParamsSetFile.m can be edited to select QSM reconstruction parameters
   LogFile.txt is the log file saved.check the ParamsSetting_cluster.m for a ParamSetFile template
check makeTable.m for a Table making template

## Supported input/output format 
Supported input file format
a. PAR/REC (Philips v4.2 and above)

b. DICOM (It may work for Philips/GE/Siemens/Canon data, but not fully tested. For DICOM files, the GRE Magnitude and Phase data should be sorted first into a seperate folder. This can be done using 3rd party DICOM packages, e.g. Horos) 

c. EnhancedDICOM (Philips) 

d. .mat (MATLAB, customized, check "readerwrapper.m" under /QSM_Utilily/fileIO/)

e. Bruker 2dseq data (not fully tested)

Output file format
a. internal/final results are saved in MATLAB native format (.mat)

b. some results are also saved with a copy of NIFTI format (.nii.gz)


## Documentation
For detailed use of the Toolbox, read the manual (v3.0, will update later)

## A brief description of the toolbox output is given below:
*_PhaseUnwrapped_echo2-5[Path/Lap/NLFPath].mat: Unwrapped phase using [path/laplacian/nonlinear fitting+path] method, echo 2-5 (user-defined) was selected for QSM processing

*GREMag#.nii.gz 			: GRE Magnitude Image of selected echo for brain masking

*GREMag#_brain.nii.gz		: GRE Magnitude Image of brain only, using FSL BET

*GREMag#_brain_mask.nii.gz	: brain mask from FSL BET

*_brain_mask_echo#_r#_t##.mat	: final brain mask used for QSM

*[SMV/PDF/LBVSMV/iRSHARP]_[avg/LSlope/NLSlope]_echo2-5.mat: Frequency map (in Hz) or local field map generated using [VSHARP/PDF/LBV+VSHARP/iRSHARP] method, taking [echo averaging/linear fitting/nonlinear fitting], with echo2-5(user-defined). Have corresponding .nii.gz file too.

*R2star.mat/R2star.nii.gz	: R2 star map

*chi_[iLSQR/TKD/iTKD/MEDI/SFCR/SFCR+0/nSFCR/nSFCR+0/FANSI/NDI/TFI]: QSM map calculated using different methods. Have corresponding .nii.gz file.

--------------------------------------------------------------

## Update-Note
2019-09-03: Add back in MS-SFCR+0 for ARIC processing

2021-05-01: Added initial phase offset removal and weighted echo averaging

2021-06-01: Added FANSI and nSFCR with L1-nonlinear data-fidelity

2021-06-29: Added cluster version without GUI

2021-10-23: Added mcpc-3Ds/ASPIRE for coil-combination

## Disclaimer
This Matlab package is available publicly, in the hope that it will be useful, but without any warranty and without even the implied warranty of merchantability or fitness for a particular purpose. It is designed to help researchers in the field of Magnetic Resonance Imaging who may have interests in doing QSM. It is not intended for use in a clinical setting.

## Contributions
Major GUI and pipeline Authors: 
 Jiri van Bergen and Xu Li
 Affiliation: Radiology @ JHU & Kirby Center @ Kennedy Krieger Institute (KKI)
 Contact via xuli@mri.jhu.edu
 
 Other contributors include:
 Joseph Gillen, Craig Jones, Jonathan Farrell, Issel Lim and Jiadi Xu (Kirby Center @ KKI)
 Lijun Bao, Jinsheng Fang (Xiamen University)
 
 Some Codes are from or modifed from other publicly availabe codes including
 1. MEDI Toolbox:		http://weill.cornell.edu/mri/pages/qsm.html
 2. Berkin Bilgic Software: 	http://martinos.org/~berkin/software.html
 3. FANSI Toolbox: 		https://gitlab.com/cmilovic/FANSI-toolbox

