//
// Created by bacon on 1/2/18.
//

#ifndef REWIRE_PROJECTIMAGE_H
#define REWIRE_PROJECTIMAGE_H

#include <Eigen/Dense>
#include "image_process_helper.h"

class ProjectImage {
public:
    ProjectImage() : ProjectImage(100){};
    explicit ProjectImage(int threshold_value);
    ~ProjectImage();
    void read_data(const std::string& filepath, int frame);
    void process_image();

    // some parameters used in the processing
    int threshold_value;

    Eigen::Matrix3d K;     // the camera intrinsic matrix
    Eigen::Matrix3d K_i;   // inverse of K
    Eigen::Matrix3d R;     // the rotation matrix to the world
    Eigen::Vector3d T;     // the translation matrix to the world
    Eigen::MatrixXd P;     // the Projection matrix of the camera model
    cv::Mat original_image;
    cv::Mat rectified_image;    // the rectified image computed by cv built-in function
    cv::Mat computed_image;     // the rectified image computed by matrix
//    cv::Mat conencted_computed_image;   // based on computed_image, we filled the gaps in it
    cv::Mat skeleton;
    std::vector<std::vector<cv::Point>> curve_segments;
};

#endif //REWIRE_PROJECTIMAGE_H
