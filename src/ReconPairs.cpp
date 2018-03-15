//
// Created by bacon on 1/2/18.
//

#include "ReconPairs.h"
#include <set>

#define DEBUG 1

#if DEBUG
#include <iostream>
#include <opencv2/core/eigen.hpp>
#include <opencv2/calib3d.hpp>
#include <opencv2/imgproc.hpp>

#endif

using namespace Eigen;
using namespace cv;
using namespace std;

HomoPoint::HomoPoint(const Point& p, int i) :idx(i)
{
    point[0] = static_cast<int>(round(p.x));
    point[1] = static_cast<int>(round(p.y));
    point[2] = 1;
}

HomoPoint::HomoPoint(Vector3d p, int i) :point(move(p)), idx(i) {}

//void Curve::compute_3d(Matrix3d Q_inv)
//{
////    disparity_segment.resize(segment[0].size());
//    curve.resize(segment[0].size());
//    for (int i = 0; i < disparity_segment.size(); i++)
//    {
////        disparity_segment[i][0] = segment[0][i][0];
////        disparity_segment[i][1] = segment[0][i][1];
////        disparity_segment[i][2] = segment[1][i][0] - segment[0][i][0];
////        disparity_segment[i][3] = 1.0;
//        curve[i] = Q_inv * segment[1][i];
//        curve[i] /= (segment[1][i][0] - segment[0][i][0])/scale;
////        curve[i] = M * disparity_segment[i];
////        curve[i] = curve[i] / curve[i][3];
//    }
//}

bool RectifiedPoint::operator==(const RectifiedPoint &p) const {
    return y == p.y and xmin-1 <= p.xmin and xmax+1 >= p.xmax;
}

void RectifiedPoint::add(Eigen::Vector3d p) {
    int new_x = static_cast<int>(round(p[0]));
    if (new_x < xmin)
    {
        xmin = new_x;
    }
    if (new_x > xmax)
    {
        xmax = new_x;
    }
    avatars.emplace_back(std::move(p));
}

RectifiedPoint::RectifiedPoint(Eigen::Vector3d p) {
    anchor[0] = static_cast<int>(round(p[0]));
    anchor[1] = static_cast<int>(round(p[1]));
    anchor[2] = static_cast<int>(round(p[2]));
    xmin = xmax = anchor[0];
    y = anchor[1];
}

ReconPairs::ReconPairs(shared_ptr<ProjectImage> mv, shared_ptr<ProjectImage> nv)
{
    view[0] = move(mv);
    view[1] = move(nv);
//    compute_fundamental();

    for (int idx = 0; idx < 2; idx++)
    {
        curve_segments[idx].resize(view[idx]->curve_segments.size());
        for (int i = 0; i < view[idx]->curve_segments.size(); i++)
        {
            curve_segments[idx][i].resize(view[idx]->curve_segments[i].size());
            for (int j = 0; j < view[idx]->curve_segments[i].size(); j++)
            {
                curve_segments[idx][i][j] = Vector3d(view[idx]->curve_segments[i][j].x, view[idx]->curve_segments[i][j].y, 1.0);
            }
        }
    }
}

ReconPairs::~ReconPairs() = default;

void ReconPairs::compute_fundamental()
{
    Matrix3d R;
    R = view[1]->R * view[0]->R.transpose();
    MatrixXd T;
    T = view[1]->T - R * view[0]->T;
    Matrix3d tx;
    tx << 0.0, -T(2,0), T(1,0),
            T(2,0), 0.0, -T(0,0),
            -T(1,0), T(0,0), 0.0;
    fundamental = view[1]->K.inverse().transpose() * tx * R * view[0]->K.inverse();
}

void ReconPairs::rectify() {
    Matrix3d R_;
    MatrixXd T_, D_;
    Mat P[2];
    Mat map_x[2], map_y[2];
    Mat R_after[2], P_after[2], R, T, D;
    Rect validROIL, validROIR;
    Size imageSize = view[0]->original_image.size();

    R_ = view[0]->R * view[1]->R.transpose();
    T_ = view[0]->T - R_ * view[1]->T;
    D_ = MatrixXd::Zero(5, 1);

    eigen2cv(R_, R);
    eigen2cv(T_, T);
    eigen2cv(D_, D);
    eigen2cv(view[0]->K, P[0]);
    eigen2cv(view[1]->K, P[1]);

    stereoRectify(P[1], D, P[0], D, imageSize, R, T, R_after[1], R_after[0], P_after[1], P_after[0], Q,
                  CV_CALIB_ZERO_DISPARITY,
                  1, imageSize, &validROIL, &validROIR);

    Matrix3d R_after_M[2], P_after_M[2];

    for (int idx = 0; idx < 2; idx++) {
        view[idx]->computed_image = Mat::zeros(view[idx]->original_image.rows, view[idx]->original_image.cols, CV_8UC1);
        initUndistortRectifyMap(P[idx], D, R_after[idx], P_after[idx], imageSize, CV_32FC1, map_x[idx], map_y[idx]);
        remap(view[idx]->skeleton, view[idx]->rectified_image, map_x[idx], map_y[idx], INTER_LINEAR);

        cv2eigen(R_after[idx], R_after_M[idx]);
        cv2eigen(P_after[idx], P_after_M[idx]);
        P_after_M[idx] = P_after_M[idx].block(0, 0, 3, 3);

        rect_segments[idx].resize(curve_segments[idx].size());
        for (int i = 0; i < curve_segments[idx].size(); i++) {
            for (int j = 0; j < curve_segments[idx][i].size(); j++) {
                Vector3d temp_v = P_after_M[idx] * R_after_M[idx] * view[idx]->K_i * curve_segments[idx][i][j];
                temp_v = temp_v / temp_v[2];
                RectifiedPoint temp_p(temp_v);
                view[idx]->computed_image.at<uchar>(temp_p.anchor[1], temp_p.anchor[0]) = 255;
                auto p = find(rect_segments[idx][i].begin(), rect_segments[idx][i].end(), temp_p);
                if (p != rect_segments[idx][i].end()) {
                    p->add(temp_v);
                } else {
                    temp_p.add(temp_v);
                    rect_segments[idx][i].emplace_back(move(temp_p));
                }
            }
        }
//        rectified_center[idx] = P_after_M[idx] * R_after_M[idx] * view[idx]->K_i * center[idx];
//        rectified_center[idx] /= rectified_center[idx](2);
    }
    KR_inv = (P_after_M[0] * R_after_M[0] * view[0]->R).inverse();
    KT = P_after_M[1] * R_after_M[1] * view[1]->T;
    scale = 1 / P_after[0].at<double>(0, 3);
//    cout << P_after[0] << endl << P_after[1] << endl;
//    cv2eigen(P_after[0], P_matrix[0]);
//    cv2eigen(P_after[1], P_matrix[1]);
//    cout << R_after_M[0] << endl << R_after_M[1] << endl << R_after_M[0] * R_after_M[1] << endl;
//
//    cout << endl << "Magic equation:" << endl << P_after_M[0] * R_after_M[0] * view[0]->K_i * view[0]->P << endl
//         << P_after_M[1] * R_after_M[1] * view[1]->K_i * view[1]->P << endl;

    //// 1. connect a 8-adjacent path 2. thinning it 3. repair
    //// To thin it we should map both rect-segment onto image (whose size may not be equal to the original one)
    //// Then, use the thin function to thin the segments

}

void ReconPairs::find_pairs()
{
    int current = 0;    // current group of the segments
    int hit = 0;

    for (int i = 0; i < rect_segments[0].size(); i++)
    {
        multiset<int> inter;     // intersection group for the current point
        int inter_num = 0;
        multiset<int> pre_inter; // intersection group for the previous point
        int pre_inter_num = 0;
        vector<Curve> temp_segment;
        for (int j = 0; j < rect_segments[0][i].size(); j++)
        {
            HomoPoint temp_main(rect_segments[0][i][j].anchor, i);
            vector<HomoPoint> neigh_candidates;
            for (int k = 0; k < rect_segments[1].size(); k++)
            {
                for (int l = 0; l < rect_segments[1][k].size(); l++)
                {
                    if (rect_segments[0][i][j].anchor[1] == rect_segments[1][k][l].anchor[1])
                    {
                        HomoPoint temp_neigh(rect_segments[1][k][l].anchor, k);
                        inter.emplace(k);
                        neigh_candidates.emplace_back(move(temp_neigh));
                        inter_num += 1;
                        hit++;
                    }
                }
            }
#if DEBUG

//            cout << "inter_num: " << inter_num << endl;
#endif

            // if the past index set is the same as the current one, we just need to push these points into the right vector
            if (pre_inter == inter && pre_inter_num == inter_num)
            {
                if (inter_num !=0 )
                {
                    for (auto &temp_s:temp_segment) {
                        double min_distance = 1e8;
                        vector<HomoPoint>::iterator temp;
                        for (auto hp = neigh_candidates.begin(); hp != neigh_candidates.end(); hp++) {
                            if (hp->idx == temp_s.idx[1]) {
//                            Vector3f temp_vector = temp_s.neigh_curve.back()[0] + hp->point;
                                double dist = (temp_s.segment[1].back() - hp->point).norm();
                                if (dist < min_distance) {
                                    temp = hp;
                                    min_distance = dist;
                                }
                            }
                        }
                        temp_s.segment[0].emplace_back(temp_main.point);
                        temp_s.segment[1].emplace_back(move(temp->point));
                        neigh_candidates.erase(temp);
                    }
                }
            }
            // if the past index is different form the current one, it means that we encounter one branching point
            // we have to cut the curve in the main_view into a new one, push back the temp vector and clear them (use std::move())
            else
            {
                if (!temp_segment.empty())
                {
                    for (auto &m: temp_segment)
                    {
                        candidates.emplace_back(move(m));
                    }
                }
                temp_segment.clear();
                if (inter_num !=0 )
                {
                    for (auto &m: neigh_candidates)
                    {
                        Curve temp;
                        temp.idx[0] = current;
                        temp.segment[0].emplace_back(move(temp_main.point));
                        temp.idx[1] = m.idx;
                        temp.segment[1].emplace_back(move(m.point));
                        temp_segment.emplace_back(move(temp));
                    }
                    current += 1;
                }
            }

            neigh_candidates.clear();
            swap(pre_inter, inter);
            inter.clear();
            pre_inter_num = inter_num;
            inter_num = 0;
        }
        for (auto &m: temp_segment)
        {
            candidates.emplace_back(move(m));
        }
    }
    total_curve_num = current;
#if DEBUG
    cout << "hit " << hit << endl;
#endif
}

//void ReconPairs::find_pairs2()
//{
//    /* For each point in the main view, find the corresponding point in the neighbour view
//     * If there is only one match point in the neighbour view, then it is
//     * Else, both candidates should be taken into account. We should record which segment each candidate belongs to
//     * in the neighbour view, and record it as the intersection set.
//     * If the number of intersections changes, or the intersection set changes, it means that the intersection situation changed
//     * Thus, we should cut off the segment in the main view into a new one, add it to @main_segments,
//     * and add those matching points in the neighbour view to @neigh_segments
//     * @neigh_segments may need linear interpolation to get a more precise result for further 3D construction
//     */
//    Vector3d para;
//
//    int current = 0;    // current group of the segments
//    int hit = 0;
//
//    for (auto &main_curve: main_view->curve_segments)
//    {
//        multiset<int> inter;     // intersection group for the current point
//        int inter_num = 0;
//        multiset<int> pre_inter; // intersection group for the previous point
//        int pre_inter_num = 0;
//        auto main_idx = int(&main_curve - &(neigh_view->curve_segments[0]));
//        vector<Curve> temp_segment;
//        for (auto &main_p: main_curve)
//        {
//            HomoPoint temp_main(main_p, main_idx);
//            // compute the parameter of the epipolar line
//            para = fundamental * temp_main.point;
//            para = para / (sqrt(para[0] * para[0] + para[1] * para[1]));
//
//            vector<HomoPoint> neigh_candidates;
//
//            for (auto &neigh_curve: neigh_view->curve_segments)
//            {
//                for (auto &neigh_p: neigh_curve)
//                {
//                    double dist = neigh_p.x * para[0] + neigh_p.y * para[1] + para[2];
//                    if (fabs(dist) < MAX_DISTANCE_TOLERANCE)
//                    {
//                        HomoPoint temp_neigh(neigh_p, int(&neigh_curve - &(neigh_view->curve_segments[0])));
//                        inter.emplace(temp_neigh.idx);
//                        neigh_candidates.emplace_back(move(temp_neigh));
//                        inter_num += 1;
//                        hit++;
//                    }
//                }
//            }
//#if DEBUG
//
////            cout << "hit " << hit << endl;
//#endif
//
//            // if the past index set is the same as the current one, we just need to push these points into the right vector
//            if (pre_inter == inter && pre_inter_num == inter_num)
//            {
//                if (inter_num !=0 )
//                {
//                    for (auto &temp_s:temp_segment) {
//                        double min_distance = 1e8;
//                        vector<HomoPoint>::iterator temp;
//                        for (auto hp = neigh_candidates.begin(); hp != neigh_candidates.end(); hp++) {
//                            if (hp->idx == temp_s.neigh_idx) {
////                            Vector3f temp_vector = temp_s.neigh_curve.back()[0] + hp->point;
//                                double dist = (temp_s.neigh_curve.back() + hp->point).norm();
//                                if (dist < min_distance) {
//                                    temp = hp;
//                                    min_distance = dist;
//                                }
//                            }
//                        }
//                        temp_s.main_curve.emplace_back(temp_main.point);
//                        temp_s.neigh_curve.emplace_back(move(temp->point));
//                        neigh_candidates.erase(temp);
//                    }
//                }
//            }
//            // if the past index is different form the current one, it means that we encounter one branching point
//            // we have to cut the curve in the main_view into a new one, push back the temp vector and clear them (use std::move())
//            else
//            {
//                if (!temp_segment.empty())
//                {
//                    for (auto &i: temp_segment)
//                    {
//                        candidates.emplace_back(move(i));
//                    }
//                }
//                temp_segment.clear();
//                if (inter_num !=0 )
//                {
//                    for (auto &i: neigh_candidates)
//                    {
//                        Curve temp;
//                        temp.main_idx = current;
//                        temp.main_curve.emplace_back(move(temp_main.point));
//                        temp.neigh_idx = i.idx;
//                        temp.neigh_curve.emplace_back(move(i.point));
//                        temp_segment.emplace_back(move(temp));
//                    }
//                    current += 1;
//                }
//            }
//
//            neigh_candidates.clear();
//            swap(pre_inter, inter);
//            inter.clear();
//            pre_inter_num = inter_num;
//            inter_num = 0;
//        }
//        for (auto &i: temp_segment)
//        {
//            candidates.emplace_back(move(i));
//        }
//    }
//    total_curve_num = current;
//#if DEBUG
//    cout << "hit " << hit << endl;
//#endif
//}

void ReconPairs::compute_3d()
{
    double w;
    for (auto &c : candidates)
    {
        c.curve.resize(c.segment[0].size());
        for (int i = 0; i < c.curve.size(); i++)
        {
            w = (c.segment[0][i][0] - c.segment[1][i][0]) * scale;
            c.curve[i] = KR_inv * (c.segment[1][i] - w * KT);
//            cout << c.curve[i] << endl;
            c.curve[i].conservativeResize(4);
//            cout << c.curve[i] << endl << endl;
            c.curve[i](3) = w;
            c.curve[i] /= c.curve[i](3);
        }
    }
//    VectorXd temp;
//    w = (rectified_center[0][0] - rectified_center[1][0]) * scale;
//    temp = KR_inv * (rectified_center[1] - w * KT);
//    temp.conservativeResize(4);
//    temp(3) = w;
//    temp /= temp(3);
//    cout << endl << "Compute Center: " << endl << "In main: " << center[0] << endl << "In neigh: " << center[1]
//         << endl << "In rectified main:" << rectified_center[0] << endl << "In rectified neigh: " << rectified_center[1]
//         << endl << "3D position:" << temp << endl;
//    VectorXd temp2;
//    for (int i = 0; i < 2; i++)
//    {
//        temp2 = P_matrix[i] * temp;
//        temp2 /= temp2[2];
//        cout << temp2 << endl;
//    }
}