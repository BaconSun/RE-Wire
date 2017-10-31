# RE-Wire
## What is RE-Wire  
This repo is the linux version of the WireArt project. I changed the VS project into a Cmake project, which is workable on Linux now.  
The project works fine on newest OpenCV, Eigen, Gurobi and Matlab thanks to their compability to the older version of themselves.  
Make sure the OpenCv, Eigen has installed on the computer and the Gurobi is set up for the matlab. The ANN ensential part of ANN library is included in the repo.  

## Changes between orginal repo
There are few changes between this one and the original one. I only changed the shell used in main.cpp from windows version to linux version,  '\' into '/' for directory path, and the usage of stat_code of file/dir. Also, in the ProjectImage.cpp, I changed a little bit in the Func() *FindClosestPointOnImgEachPoint* to return the correct value. Now it works fine on Linux.

## Further Work
The current version of WireArt is powerful but slow. I decide to change it into a weaker but much quicker one, to deal with the real-time construction of simple topology of wires. The combination of ROS will make it capbale of being applied with robots.