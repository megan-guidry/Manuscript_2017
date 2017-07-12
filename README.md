# Manuscript_2017
Models and code necessary for exploring the affect [Ca2+]i has on the shape of the isometric and work-loop end-systolic (ES) curve 

0) The models required to replicate Kenneth's Matlab Isometric Force vs. time transients are:
  
**As of 12-7-2017 the OpenCMISS combined model and Kenneth's Matlab model do not produce identical results.  It is possible to replicate Kenneth's results if the output time-step of the OpenCMISS model is set to 1... however, testing the model in OpenCOR with a Euler forward solver revealed that using time-step of 1 results in an unconverged solution.  Andre is investigating OpenCMISS to see why the output time-step is an issue...
