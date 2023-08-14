# Viz2DManifolds
Matlab code for extracting 2D manifold of null points in 3D space

Instruction for running code with Mex :
Strongly recommend compiling the Mex file of rk4.m in Matlab Coder.
Before compiling the Mex file, you must run "initcircles.m" to generate the necessary variables.
You can use flowline.m for Autodefine Input Types and checking issues.
The name of the Mex function should be rk4_mex (default name).
Be aware you need to recompile and overwrite the mex file when you change 
the size of cubic data (Nxyz).

If you don't want to compile the Mex file. Please change the "rk4_mex" function to "rk4" in files "timestp.m" and "intersectfr.m".

This code works for the regular cubic grid data with a spheric inner boundary. 
The inner boundary is located at [0,0,0] with a radius of 4 Re.



For any questions:
Email: crazyxpk[at]gmail.com 
Please change [at] to @

