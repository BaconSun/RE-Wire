The format of the file of camera parameters

For all the examples in the paper, the inputs for our main algorithm are proveided so that users can directly run end2end.m (in the folder "Final").  
For a new example, users need to run VisualSFM (http://ccwu.me/vsfm/) to obtain the camera parameters (noted: please use the function "Use shared calibration"), extract image curves like "viewx_bd0.jpg", and process "viewx_bd0.jpg" for the thinning and segmentation of image curves by running main.m in the folder "Thinning".

The format of the file of camera parameters is defined as follows:
# Focal Length (of the undistorted image)
# 2-vec Principal Point (image center)
# 3-vec Translation T (as in P = K[R T])
# 3x3 Matrix format of R