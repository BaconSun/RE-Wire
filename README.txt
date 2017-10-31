Implementation for the wire art paper:
Project: Image-based Reconstruction of Wire Art
Venue: SIGGRAPH 2017 

Introduction:
	We provide the inputs of all the examples in the paper (in the folder "Models"), the source code of our main algorithm (in the folder "Reconstruct3DCurves" and "Final"), and also the main code for preprocessing (in the folder "Thinning").
	The source code of our main algorithm consists of two parts. The 3D Candidate Curve Generation (Section 4.1 of our paper) was developed in C++, and the source code of this part is included in the folder "Reconstruct3DCurves". It was compiled as CorresUI.exe by Visual Studio 2015 x64 Release (in the folder "matlabExe" under the folder "Final"), and was called in Matlab. The 3D Curve Selection (Section 4.2 of our paper) and the 3D wire decomposition (Section 5 of our paper) were developed in Matlab, the source code is included in the folder "Final".
 
Dependencies:
	Gurobi Optimizer 7.0 with Matlab interface

Platform:
	The code was developed and tested on Windows 10 x64 with Matlab R2016a x64. (noted: The C++ code was developed and compiled with Visual Studio 2015 x64.)

Run examples:
    In the folder "Models", we provide the inputs of all the examples in the paper. 
    In MATLAB, change the model_name in end2end.m, and directly run it. 

Preprocessing:
	For all the examples in the paper, the inputs for our main algorithm are proveided so that users can directly run end2end.m (in the folder "Final").  
	For a new example, users need to run VisualSFM (http://ccwu.me/vsfm/) to obtain the camera parameters (noted: please use the function "Use shared calibration"), extract image curves like "viewx_bd0.jpg", and process "viewx_bd0.jpg" for the thinning and segmentation of image curves by running main.m in the folder "Thinning".
	
Parameters:
    The parameters of all the models are in setup_NL.m (for the first part of the algorithm) and setup.m (for the second part of the algorithm).
	
For any questions or comments, please contact liulingjie0206@gmail.com.