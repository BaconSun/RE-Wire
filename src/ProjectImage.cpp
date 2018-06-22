//
// Created by bacon on 1/2/18.
//

#include "ProjectImage.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <fstream>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;
using namespace std;

void ProjectImage::read_data(const string& filepath, int frame)
{
    // Read the image
    original_image = imread(filepath+"_"+to_string(frame)+".png", 0);

    // For test, inverse the image gray value
//    for (int i=0; i<original_image.rows; i++)
//        for (int j=0; j<original_image.cols; j++)
//        {
//            original_image.at<uchar>(i,j) = 255 - int(original_image.at<uchar>(i,j));
//        }

    // Read the projection matrix and other matrix
    ifstream camera_config(filepath+".txt", fstream::in);
    for (int i = 0; i < 9; i++)
        camera_config >> K(i/3, i%3);
    for (int i = 0; i < 12; i++)
        camera_config >> P(i/4, i%4);
    for (int i=0; i < 9; i++)
        camera_config >> R(i/3, i%3);
    for (int i=0; i < 3; i++)
        camera_config >> T(i);
    K_i = K.inverse();

    process_image();
}


void ProjectImage::process_image()
{
    threshold(original_image, skeleton, threshold_value, 255, THRESH_BINARY);
    skel(skeleton, skeleton);
    curve_segment(skeleton, curve_segments);
    total_sample_points = 0;
    for (const auto & i : curve_segments)
    {
        total_sample_points += i.size();
    }
}

ProjectImage::ProjectImage(int threshold_value) : threshold_value(threshold_value)
{
    this->P.resize(3,4);
}

ProjectImage::~ProjectImage() = default;