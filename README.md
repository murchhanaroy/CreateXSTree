Program to calculate rate using XSModel.
XSModel return DiffXS in ub/MeV-Sr for inelas or ub for elas.
 
To be improved:
1) HMS and SHMS acceptance.
   One need to run single-arm-simulation program to get the angle acceptance. 
   This angle acceptance should have a z dependence.  I am not sure if the
   single-arm-simulation program has already taken z dependence into account.
   If not, one need to develop a transportation packing for single-arm-simulation
   program first.

2) Make A1NRates.cc as a class such that it can be used in root cint.

First release version v1.0.0.

1. Hard-coded HMS and SHMS collimator size and angle acceptance, need to be improved.
2. Hard-coded target material, it will be better to separate them.

Notes on running this code to get rate weighted counts:

1. Run mc-single-arm : Run mc-single-arm.f with specific input file for a particular kinematic settings. In the input file, the spectrometer angle is entered as an absolute value, no need to specify negative sign for HMS.
2. mc-single-arm will generate root files in worksim/ directory with required branches.
3. Once mc-single-arm outputs rootfiles, go to CreateXSTree repository to run the code to get cross sections.
4. To run XSTree code, do ./a1n xstree <pElasOnly> <pBeam_GeV> <pDetectorAngle_deg> <pDetectorMomentum_GeV> <pDetector=1,10 HMS|2,20 SHMS> <rootfile>
The pDetectorAngle_Deg should have a negative sign for HMS.
5. This code will add new branches with rates to the rootfile specified (generated from mc-single-arm).