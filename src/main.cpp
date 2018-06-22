//
// Created by bacon on 12/16/17.
//

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "ProjectImage.h"
#include "ReconPairs.h"

using namespace std;
using namespace cv;
using namespace Eigen;

int main (int argc, char* argv[])
{
    string workspace = "/home/bacon/catkin_ws/src/rewire/data/Curve/";

    const int camera_num = 5;

    shared_ptr<ProjectImage> camera[camera_num];
    for (int i  = 0; i < camera_num; i++)
    {
        camera[i] = make_shared<ProjectImage>();
        camera[i]->read_data(workspace+"Camera"+to_string(i+1), 1);
        cout << camera[i]->curve_segments.size() << " ";
//        imshow(to_string(i), camera[i]->skeleton);
    }
    cout << endl;
//    waitKey(0);

    shared_ptr<ReconPairs> recon = make_shared<ReconPairs>(camera[0], camera[3], camera[4]);

/*
//    Vector4d center;
//    center[0] = -0.29076; center[1] = -0.49210; center[2] = 0.45002; center[3] = 1.0;
//    for (int i = 0; i < 2; i++)
//    {
//        recon->center[i] = recon->view[i]->P * center;
//        recon->center[i] /= recon->center[i][2];
//        for (int j = 0; j < 3; j++)
//            recon->center[i](j) = round(recon->center[i](j));
//    }
*/

    recon->rectify();
    recon->find_pairs();
    recon->compute_3d();

    /*
//    imshow("ImageL Original", recon->view[1]->original_image);
//    imshow("ImageR Original", recon->view[0]->original_image);
//    imshow("ImageL", recon->view[1]->skeleton);
//    imshow("ImageR", recon->view[0]->skeleton);
//    imshow("ImageL After Rectify", recon->view[1]->rectified_image);
//    imshow("ImageR After Rectify", recon->view[0]->rectified_image);

//    Mat testImage[2];
//    for (int idx=0; idx<2; idx++)
//    {
//        testImage[idx] = Mat::zeros(recon->view[idx]->original_image.rows, recon->view[idx]->original_image.cols, CV_8UC1);
//        for (auto c: recon->rect_segments[idx])
//            for (auto p: c)
//            {
////                cout << idx << " " << p.anchor << endl;
//                testImage[idx].at<uchar>(p.anchor[1],p.anchor[0]) = 255;
//            }
//    }
//    cout << recon->view[0]->original_image.rows << " " << recon->view[0]->original_image.cols << endl;
//    imshow("Image compute by hand R", testImage[0]);
//    imshow("Image compute by hand L", testImage[1]);

//    waitKey(0);
//    return 0;

//    recon->compute_fundamental();
//    double sum_F = 0;
//    double result = 0;
//    double dist = 0;
//    double distorted_result = 0;
//    double distorted_dist = 0;
//    int count = 0;
//    for (double i = -5; i <= 5; i+=1.0)
//        for (double j = -5; j <= 5; j+=1.0)
//            {
//                Vector4d homo_p(i,j,0,1);
//                Vector3d main_p, neigh_p, para;
//                sum_F += (recon->neigh_view->P * homo_p).transpose() * recon->fundamental * recon->main_view->P * homo_p;
//                main_p = recon->main_view->P * homo_p;
//                neigh_p = recon->neigh_view->P * homo_p;
////                cout << main_p[2] << " " << neigh_p[2] << endl;
//                main_p[0] = round(main_p[0]/main_p[2]); main_p[1] = round(main_p[1]/main_p[2]); main_p[2] = 1.0;
//                neigh_p[0] = round(neigh_p[0]/neigh_p[2]); neigh_p[1] = round(neigh_p[1]/neigh_p[2]); neigh_p[2] = 1.0;
//                result += fabs(neigh_p.transpose() * recon->fundamental * main_p);
//                para = recon->fundamental * main_p;
//                dist += result / sqrt(para[0] * para[0] + para[1] * para[1]);
//                neigh_p[0] += 1; neigh_p[1] += 1;
//                distorted_result += fabs(neigh_p.transpose() * recon->fundamental * main_p);
//                distorted_dist += distorted_result / sqrt(para[0] * para[0] + para[1] * para[1]);
//                count ++;
//            }
//
//    cout << "Result: " << result/count << endl
//         << "Distance:" << dist/count << endl
//         << "Distorted result: " << distorted_result/count << endl
//         << "Distorted distance: " << distorted_dist/count << endl
//         << "Total error of F:" << sum_F/count << endl;
//
//    return 0;
     */

    int total_points = 0;
    for (auto &i:camera[0]->curve_segments)
        total_points += i.size();
    cout << "Main Total points: " << total_points << endl;

    total_points = 0;
    for (auto &i:camera[1]->curve_segments)
        total_points += i.size();
    cout << "Neigh Total points: " << total_points << endl;

    total_points = 0;
    for (auto &i:recon->rect_segments[0])
        total_points += i.size();
    cout << "Rectified Main Total points: " << total_points << " Rectified Main segments: " << recon->rect_segments[0].size() << endl;

    total_points = 0;
    for (auto &i:recon->rect_segments[1])
        total_points += i.size();
    cout << "Rectified Neigh Total points: " << total_points << " Rectified Neigh segments: " << recon->rect_segments[1].size() << endl;

    cout << "Total curve num: " << recon->total_curve_num << endl;
    cout << "Candidate size: " << recon->candidates.size() << endl;



    Mat combine[2];

    /*

//    for (int i = 0; i < 2; i++)
//    {
//        combine[i] = recon->view[i]->computed_image.clone();
//    }
//
//    imshow("main", combine[0]);
//    imshow("neigh", combine[1]);
//    waitKey(0);

//    int mature_segment = 0;
//    for (int idx = 0; idx < recon->candidates.size(); idx++)
//    {
//        if (recon->candidates[idx].curve.size() > 5)
//            mature_segment ++;
//        cout << endl << "Curve candidate Num: " << idx+1 << endl;
//        Mat combine[2];
//        for (int i = 0; i < 2; i++)
//        {
//            combine[i] = recon->view[i]->computed_image.clone();
//            cvtColor(combine[i], combine[i], cv::COLOR_GRAY2BGR);
//            for (int j = 0; j < recon->candidates[idx].segment[i].size(); j++)
//            {
////            cout << "Point " << j+1 << " : " << recon->candidates[i].curve[j] << endl;
//                Vector3d temp = recon->candidates[idx].segment[i][j];
//                temp = temp / temp[2];
//                if (temp[0] >= 0 and temp[0] <= 640 and temp[1] >= 0 and temp[1] <= 480)
//                {
//                    cout << temp << endl << endl;
//                    combine[i].at<Vec3b>(Point(int(temp[0]), int(temp[1])))[0] = 255;
//                    combine[i].at<Vec3b>(Point(int(temp[0]), int(temp[1])))[1] = 255;
//                    combine[i].at<Vec3b>(Point(int(temp[0]), int(temp[1])))[2] = 0;
//                }
//            }
//        }
//        imshow("main", combine[0]);
//        imshow("neigh", combine[1]);
//        waitKey(0);
//    }
//    cout << "No. of segments with enough candidates: " << mature_segment << endl;

//    for (auto &i: recon->candidates)
//    {
//        cout << i.main_idx << " ";
//    }
//
//    cout << endl;
//    for (auto &i: recon->candidates)
//    {
//        cout << i.neigh_idx << " ";
//    }
//
//    cout << endl;
//
//    total_points = 0;
//    for (auto &i: recon->candidates)
//    {
//        cout << i.main_idx << " " << i.neigh_idx << endl;
//        for (auto &j: i.main_curve)
//        {
//
//            cout << j << endl;
//        }
//        for (auto &j: i.neigh_curve)
//        {
//
//            cout << j << endl;
//        }
//        cout << endl;
//        total_points += i.main_curve.size();
//    }
//
//
//
//    cout << "Total points: " << total_points << endl;
//

//    Mat combine[2];
    */

    ofstream fout;
    fout.open(workspace+"3dpoints.txt");
    for (const auto &c: recon->candidates)
    {
        fout << "Curve Number:" << c.idx[0] << " " << c.idx[1] << endl;
        int counter = 1;
        for (const auto &p: c.curve)
        {
            fout << counter << " Point: " << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << endl;
            counter ++;
        }
        fout << endl;
    }

    for (auto &j: recon->candidates)
    {
//        combine[i] = recon->view[i]->skeleton.clone();
//
//        cvtColor(combine[i], combine[i], cv::COLOR_GRAY2BGR);
//        Vector3d temp_center = recon->view[i]->P * center;
//        temp_center = temp_center / temp_center[2];
//
//        cout << "Image: " << i << "cneter: "  << temp_center << endl;
//
//        combine[i].at<Vec3b>(Point(int(temp_center[0]), int(temp_center[1])))[0] = 255;
//        combine[i].at<Vec3b>(Point(int(temp_center[0]), int(temp_center[1])))[1] = 0;
//        combine[i].at<Vec3b>(Point(int(temp_center[0]), int(temp_center[1])))[2] = 255;

        for (int i = 0; i < 2; i++)
        {
            combine[i] = recon->view[i]->skeleton.clone();
            cvtColor(combine[i], combine[i], cv::COLOR_GRAY2BGR);

            for (auto &p: j.curve)
            {
                Vector3d temp = recon->view[i]->P * p;
                temp = temp / temp[2];
                if (temp[0] >= 0 and temp[0] <= 640 and temp[1] >= 0 and temp[1] <= 480)
                {
                    combine[i].at<Vec3b>(Point(int(temp[0]), int(temp[1])))[0] = 255;
                    combine[i].at<Vec3b>(Point(int(temp[0]), int(temp[1])))[1] = 255;
                    combine[i].at<Vec3b>(Point(int(temp[0]), int(temp[1])))[2] = 0;
                }
                else
                {
                    cout << "Point not in image! Image: " << i << " Point: " << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << endl;
                }
            }

        }
        imshow("main", combine[0]);
        imshow("neigh", combine[1]);
        waitKey(0);
    }

/*
//    for (int i=0; i<recon->candidates.size(); i++)
//    {
//        if (i%2==0)
//        {
//            for (auto &n: recon->candidates[i].main_curve)
//            {
//                main.at<Vec3b>(Point(int(n[0]), int(n[1])))[1]=0;
//            }
//            for (auto &n: recon->candidates[i].neigh_curve)
//            {
//                neigh.at<Vec3b>(Point(int(n[0]), int(n[1])))[1]=0;
//            }
//        }
//        else
//        {
//            for (auto &n: recon->candidates[i].main_curve)
//            {
//                main.at<Vec3b>(Point(int(n[0]), int(n[1])))[2]=0;
//            }
//            for (auto &n: recon->candidates[i].neigh_curve)
//            {
//                neigh.at<Vec3b>(Point(int(n[0]), int(n[1])))[2]=0;
//            }
//        }
//    }
 */

//    imshow("main", combine[0]);
//    imshow("neigh", combine[1]);
////    imshow("ImageL", recon->view[1]->skeleton);
////    imshow("ImageR", recon->view[0]->skeleton);
//    waitKey(0);

    // To release the memory used by ANN
    annClose();

    return 0;
}