# FEM

NOTE: The FEM4X series and loopFEM04 are updated versions of the original FEM processing scripts. The README has been updated to reflect the newer versions.

sampleStack.mat is a stack of 10 FEM images. It is already a data cube. The scale for the stack should be set to 0.022 in FEM41. Scale is in inverse nm and is dependent on camera length.

These MATLAB codes are used for processing and analyzing fluctuation electron microscopy (FEM) data collected on transmission electron microscopes. They take in 3D stacks of FEM data, perform elliptical correction (when desired), and provide information on intensity variance (medium-range ordering) and relative strain. They were developed at the Molecular Foundry and Lawrence Berkeley National Laboratory. They were developed by Colin Ophus, Ellis Kennedy, Mary Scott, and Peter Ercius. The codes provide information on structural variance and relative strains within amorphous thin films. Contact E.K. at ellisrae@berkeley.edu for assistance.

The workflow is as follows: <br/>

__1__ - Use dm4Reader.m to read in and name a stack of .dm4 files for experimental data. Select the file from a directory or type in file pathway. For simulated data, you can skip this step as long as data is a 3D stack of images. <br/>
  
  Command: stack = dm4Reader();

__2__ - Create a full mask and ring mask from the mean stack image with RealspaceLattice01.m. Requires RealspaceLattice01.fig to be in same folder. Do not keep beam stop or central bright spot. <br/>
  
  Command: stack_fullmask = RealspaceLattice01(mean(stack.cube,3));
  Command: stack_ringmask = RealspaceLattice01(mean(stack.cube,3));
  
   <img src="https://github.com/ScottLabUCB/FEM/blob/master/exampleFullMask.PNG" height="300"> <img src="https://github.com/ScottLabUCB/FEM/blob/master/exampleRingMask.PNG" height="300">
  
__3__ - Run FEM41.m. Looking at your mean FEM image, the ring pixel intensity should be between 500 and 2000; if it is not, multiply your stack by an appropriate factor before running FEM41.m. Additionally, set your pixel size in the body of the code if it is not included in the metadata of your stack. <br/>
  
  Command: sFEM_stack = FEM41(stack.cube, stack_fullmask, stack_ringmask);
  
  If intensity is too low: <br/>
  Command: sFEM_stack = FEM41(stack.cube `*`100, stack_fullmask, stack_ringmask);

__4__ - Run FEM42.m to normalize images within mask areas. <br/>
  Command: sFEM_stack = FEM42(sFEM_stack);
  
__5__ - Run FEM43.m. The second input is optional. Running the command with the second input will fit your FEM patterns to the input parameters. This should be done when looking at relative strain, but not for medium-range order variance. Runnng the command without the second input will fit the FEM pattern based on initial guesses and output the fit parameters. The initial guess parameters must be added to the body of the code directly. See image below for example output from FEM43. The ring and center should have the greates intensity and the fitting and actual image slices should be aligned.

Initial guess values correspond to: 1) y-coordinate of center, 2) x-coordinate of center, 3) C ellipse fit parameter (should be 1.0 for guess), 4) B ellipse fit parameter (should be 0.0 for guess), 5) intensity beyond first ring or background intensity, 6) intensity between central bright spot and first ring, 7) first ring width, 8) first ring radius, 9) second ring radius, 10) inner Gaussian tail of ring, 11) outer Gaussian tail. The most important parameters for fitting are 1) through 8). Paramters 9) through 11) can be faily off and the FEM pattern will still be fit. If there is trouble with fitting, check that the intensity of your mean ring is between 500 and 2000 or scale your parameters accordingly. <br/>
  
  Command: sFEM_stack = FEM43(sFEM_stack);
  
  Sample parameters for experimental data: <br/>
   sFEM.coefsInit = [270 260 1.00 0.00 200 10000 20 150 200 15 15];
   
   <img src="https://github.com/ScottLabUCB/FEM/blob/master/exampleFEM43.PNG" height="300">
  
__6__ - Run FEM44.m for a polar transformation of all images. This code has dependency on convolve02.m and exindex.m. These are NOT my or my colleagues' code. They were developed by D. Young and can be found on the Matlab File Exchange. <br/>
  
  Command: sFEM_stack = FEM44(sFEM_stack);
  
__7__ - Run FEM45.m to compute mean and variance. A plot of the variance will be displayed. sFEM_stack.radialVarNorm are the y-axis variance values and sFEM_stack.polarRadius are the x-axis k scattering vector (in inverse nm) values.  <br/>
  
  Command: sFEM_stack = FEM45(sFEM_stack);
  
  <img src="https://github.com/ScottLabUCB/FEM/blob/master/exampleFEM45.PNG" height="300">
  
__8__ - Run FEM46.m to calculate relative [exx eyy exy] strains. Multiply sFEM_stack.strainMeas by 100 for percent relative strain. <br/>
  
  Command: sFEM_test = FEM01_MaskCoords(sFEM_stack);
  
  
__Note:__ Altenately, steps 3-8 can be performed together using loopFEM04. The coefs input is optional. <br/>
  
  Command: sFEM_stack = loopFEM04(stack.cube, stack_fullmask, stack_ringmask, coefs);
  
For variance plots, plot sFEM_stack.radialVarNorm on the y-axis and sFEM_stack.polarRadius on th x-axis. <br/>
For percent relative strain, refer to sFEM_stack.strainMeas.


The script dm4Reader.m was developed by Peter Ercius at the Molecular Foundry, Lawrence Berkeley National Laboratory. exindex.m and convolve02.m were developed by David Young.
