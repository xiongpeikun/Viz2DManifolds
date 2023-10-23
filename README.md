# Viz2DManifolds
Matlab code for extracting 2D manifold of null points in 3D space.
Briefly introduction and result: https://youtu.be/bHYBmKZsoGY

The code saves the geodesic band for every recursion in the .mat (saved in '/bands') and Tecplot binary files (saved in '/plotfile').
This code works for the regular cubic grid data with a spheric inner boundary. 
The inner boundary is located at [0,0,0] with a radius of 4 Re.

Instruction for running code with Mex :

Strongly recommend compiling the Mex file of rk4.m in Matlab Coder.
Before compiling the Mex file, you must run "initcircles.m" to generate the necessary variables.
You can use flowline.m for Autodefine Input Types and checking issues.
The name of the Mex function should be "rk4_mex" (default name).
Be aware you need to recompile and overwrite the mex file when you change 
the size of cubic data (Nxyz).

If you don't want to compile the Mex file. Please set mex=0 in GLSmain.m.




Email: crazyxpk[at]gmail.com 
Please change [at] to @

