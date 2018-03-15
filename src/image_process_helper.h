//
// Created by bacon on 1/2/18.
//

#ifndef REWIRE_IMAGE_PROCESS_HELPER_H
#define REWIRE_IMAGE_PROCESS_HELPER_H

#include <opencv2/core/core.hpp>

void GetLutSkel(cv::Mat& Lut);
void applylut_1(cv::Mat &src,cv::Mat &dst);
void applylut_joint(cv::Mat &src,cv::Mat &dst);
void applylut_8(cv::Mat &src,cv::Mat &dst,cv::Mat& lut);
void skel(cv::Mat &src,cv::Mat &dst);
void endp(cv::Mat &src,cv::Mat &dst);
void jointp(cv::Mat &src,cv::Mat &dst);
void curve_segment(cv::Mat &src, std::vector<std::vector<cv::Point>>& curve_segments);

#endif //REWIRE_IMAGE_PROCESS_HELPER_H
