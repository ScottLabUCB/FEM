# FEM

These MATLAB codes are used for processing and analyzing fluctuation electron microscopy (FEM) data collected on transmission electron microscopes. They were developed at the Molecular Foundry and Lawrence Berkeley National Laboratory. They were developed by Colin Ophus, Ellis Kennedy, Mary Scott, and Peter Ercius. The codes provide information on structural variance and relative strains within amorphous thin films. Contact E.K. at ellisrae@berkeley.edu for assistance.

The workflow is as follows:

1 - use dm4Reader to read in and name a stack of .dm4 files <br/>
  Ex: test = dm4Reader();

2 - create a full mask and ring mask from the mean stack image with RealspaceLattice01; see mask.png for full mask example
  Ex: test_fullmask = RealspaceLattice01(mean(test.cube,3));
  
3* - run FEM01_MaskCoords
  Ex: sFEM_test = FEM01_MaskCoords(test.cube, test_fullmask, test_ringmask);

4* - run FEM02_FitFunction; coefs input is optinal; approximate the fitting parameters within body of code
  Ex: sFEM_test = FEM02_FitFunction(sFEM_test);
  
5* - run FEM03_Variance
  Ex: sFEM_test = FEM03_Variance(sFEM-test, test.cube);
  
6* - run FEM04_Strain
  Ex: sFEM_test = FEM01_MaskCoords(sFEM_test);
  
Altenately, steps 3-6 can be performed together using loopFEM; coefs input is optional
  Ex: sFEM_test = loopFEM(test.cube, test_fullmask, coefs);
  
For variance plots, plot sFEM_test.yy on th y-axis and sFEM_test.polarRadius on th x-axis.
For percent relative strain, refer to sFEM.strainMeas.
