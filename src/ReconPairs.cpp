//
// Created by bacon on 1/2/18.
//

#include "ReconPairs.h"
#include <set>

#define DEBUG 1

#if DEBUG
#include <iostream>
#include <algorithm>
#include <opencv2/core/eigen.hpp>
#include <opencv2/calib3d.hpp>
#include <opencv2/imgproc.hpp>
#include <gurobi_c++.h>

#endif

using namespace Eigen;
using namespace cv;
using namespace std;

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

// Just used to check whether the latter point is inside the range of previous point.
bool RectifiedPoint::operator == (const RectifiedPoint &p) const {
    return y == p.y and xmin-1 <= p.xmin and xmax+1 >= p.xmax;
}

void RectifiedPoint::add(const Vector3d& p) {
    auto new_x = static_cast<int>(round(p[0]));
    if (new_x < xmin)
    {
        xmin = new_x;
    }
    if (new_x > xmax)
    {
        xmax = new_x;
    }
//    avatars.emplace_back(std::move(p));
    set_center();
}

void RectifiedPoint::set_center()
{
    center = Vector3d ((xmax+xmin)/2.0, y, 1.0);
}

RectifiedPoint::RectifiedPoint(const Vector3d& p) {
    anchor[0] = p[0];
    anchor[1] = p[1];
    anchor[2] = p[2];
    xmin = xmax = static_cast<int>(round(anchor[0]));
    y = static_cast<int>(round(anchor[1]));
    set_center();
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

ReconPairs::ReconPairs(shared_ptr<ProjectImage> mv, shared_ptr<ProjectImage> nv, shared_ptr<ProjectImage> tv): ReconPairs(move(mv), move(nv))
{
    set_third_view(move(tv));
}

void ReconPairs::set_third_view( shared_ptr<ProjectImage> tv)
{
    third_view = move(tv);
    build_kd_tree();
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
                temp_p.idx = i;
                view[idx]->computed_image.at<uchar>(static_cast<int>(round(temp_p.anchor[1])), static_cast<int>(round(temp_p.anchor[0]))) = 255;
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

    /*
    // Try to connect the gaps between the points in the rect-segment
    for (auto & rect_segment : rect_segments) {
        for (auto & segment : rect_segment)
        {
            int i = 0;
            while (i != (segment.size()-1))
            {
                // for two points, if their y are continuous, just change there xmax and xmin; else, insert another points
                // to make sure the line is continuous.

                if (segment[i].y - segment[i+1].y <= 1 && segment[i].y - segment[i+1].y >= -1)
                {
                    auto temp = static_cast<int>(round((segment[i].xmax + segment[i].xmin + segment[i+1].xmax + segment[i+1].xmin) / 4.0));
                    if (segment[i].xmax + segment[i].xmin < segment[i+1].xmax + segment[i+1].xmin)
                    {
                        segment[i].xmax = temp; segment[i+1].xmin = temp;
                    }
                    else
                    {
                        segment[i].xmin = temp; segment[i+1].xmax = temp;
                    }
                }
                else
                {
                    vector<RectifiedPoint> temp_segment;
                    auto prev_center = static_cast<int>(round((segment[i].xmin + segment[i].xmax) / 2.0));
                    auto next_center = static_cast<int>(round((segment[i+1].xmin + segment[i+1].xmax) / 2.0));
                    int d = segment[i].y - segment[i+1].y;
                    d = (d > 0 ? d : -d);
                    double step = (prev_center - next_center) / 2.0 / d;
                    step = (step > 0 ? step : -step);

                    for (int j = 1; j < d; j++)
                    {
                        Vector3d temp_vector;
                        temp_vector[0] = ((d-j)*prev_center + j*next_center) * 1.0 / d;
                        temp_vector[1] = ((d-j)*segment[i].y + j*segment[i+1].y) * 1.0 / d;
                        temp_vector[2] = 1;
                        RectifiedPoint temp(temp_vector);
                        temp.xmax = static_cast<int>(round(temp_vector[0] + step));
                        temp.xmin = static_cast<int>(round(temp_vector[0] - step));
                        temp.add(temp_vector);
                        temp_segment.emplace_back(move(temp));
                    }

                    if (segment[i].xmax + segment[i].xmin < segment[i+1].xmax + segment[i+1].xmin)
                    {
                        segment[i].xmax = static_cast<int>(round(prev_center + step));
                        segment[i+1].xmin = static_cast<int>(round(next_center - step));
                    }
                    else
                    {
                        segment[i].xmin = static_cast<int>(round(prev_center - step));
                        segment[i+1].xmax = static_cast<int>(round(next_center + step));
                    }

                    segment.insert(segment.begin()+i+1, temp_segment.begin(), temp_segment.end());
                    i += d-1;
                }
                i++;
            }
        }
    }


    // Draw the rect_segments on the connected_computed_image
    for (int idx = 0; idx < 2; idx++) {
        view[idx]->conencted_computed_image = Mat::zeros(view[idx]->original_image.rows, view[idx]->original_image.cols, CV_8UC1);
        for (auto segment : rect_segments[idx])
        {
            for (auto rect_point: segment)
            {
                for (int i = rect_point.xmin; i <= rect_point.xmax; i++)
                {
                    view[idx]->conencted_computed_image.at<uchar>(rect_point.y, i) = 255;
                }
            }
        }
    }
     */

}

void ReconPairs::find_pairs()
{
    int current = 0;    // current group of the segments
    int hit = 0;

    for (const auto& main_segment : rect_segments[0])
    {
        multiset<int> inter;     // intersection group for the current point
        int inter_num = 0;
        multiset<int> pre_inter; // intersection group for the previous point
        int pre_inter_num = 0;
        vector<Curve> temp_segment;
        for (const auto& main_point : main_segment)
        {
            vector<RectifiedPoint> neigh_candidates;
            for (const auto& neigh_segment : rect_segments[1])
            {
                for (const auto& neigh_point : neigh_segment)
                {
                    if (main_point.y == neigh_point.y)
                    {
                        inter.emplace(neigh_point.idx);
                        neigh_candidates.emplace_back(neigh_point);
                        inter_num += 1;
                        hit++;
                    }
                }
            }
#if DEBUG

//            cout << "inter_num: " << inter_num << endl;
#endif

            // if the past index set is the same as the current one, we just need to push these points into the correct vector
            if (pre_inter == inter and pre_inter_num == inter_num)
            {
                if (inter_num !=0 )
                {
                    for (auto &temp_s:temp_segment) {
                        double min_distance = 1e8;
                        vector<RectifiedPoint>::iterator temp;
                        for (auto hp = neigh_candidates.begin(); hp != neigh_candidates.end(); hp++) {
                            if (hp->idx == temp_s.idx[1]) {
                                double dist = (temp_s.segment[1].back().center - hp->center).norm();
                                if (dist < min_distance) {
                                    temp = hp;
                                    min_distance = dist;
                                }
                            }
                        }
                        temp_s.segment[0].emplace_back(main_point);
                        temp_s.segment[1].emplace_back(*temp);
                        neigh_candidates.erase(temp);
                    }
                }
            }
                // if the past index is different form the current one, it means that we encounter one branching point
                // we have to cut the curve in the main_view into a new one, push back the temp vector and clear them
            else
            {
                // Push current temp segments into candidates. Clear the vector
                candidates.insert(candidates.end(), temp_segment.begin(), temp_segment.end());
                temp_segment.clear();

                if (inter_num !=0 )
                {
                    for (const auto &neigh_point: neigh_candidates)
                    {
                        Curve temp;
                        temp.idx[0] = current;
                        temp.segment[0].emplace_back(main_point);
                        temp.idx[1] = neigh_point.idx;
                        temp.segment[1].emplace_back(neigh_point);
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
        candidates.insert(candidates.end(), temp_segment.begin(), temp_segment.end());
    }

    total_curve_num = current;
#if DEBUG
    cout << "hit " << hit << endl;
#endif
}

void ReconPairs::compute_3d()
{
    // For sparse curve
    double w;
    for (auto &c : candidates)
    {
        c.sparse_curve.resize(c.segment[0].size());
        for (int i = 0; i < c.sparse_curve.size(); i++)
        {
            w = (c.segment[0][i].center[0] - c.segment[1][i].center[0]) * scale;
            VectorXd temp_point;
            temp_point = KR_inv * (c.segment[1][i].center - w * KT);
            temp_point.conservativeResize(4);
            temp_point(3) = w;
            temp_point /= temp_point(3);
            c.sparse_curve[i] = temp_point.head(4);
        }
    }

    // For dense curve
    for (auto &c : candidates)
    {
        for (int i = 0; i < c.segment[0].size(); i++)
        {
            int main_step = c.segment[0][i].xmax - c.segment[0][i].xmin;
            int neigh_step = c.segment[1][i].xmax - c.segment[1][i].xmin;
            if (main_step == 0)
            {
                c.curve.emplace_back(c.sparse_curve[i]);
            }
            else
            {
                for (int x = 0; x <= main_step; x++)
                {
                    Vector3d main_point(c.segment[0][i].xmin + x, c.segment[0][i].y, 1.0);
                    Vector3d neigh_point(c.segment[1][i].xmin + neigh_step * 1.0 / main_step * x, c.segment[0][i].y, 1.0);
                    w = (main_point[0] - neigh_point[0]) * scale;
                    VectorXd temp_point;
                    temp_point = KR_inv * (neigh_point - w * KT);
                    temp_point.conservativeResize(4);
                    temp_point(3) = w;
                    temp_point /= temp_point(3);
                    c.curve.emplace_back(move(temp_point.head(4)));
                }
            }
        }
    }

    /*
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
     */
}

void ReconPairs::build_kd_tree()
{
    // build the ANN tree for further use
    ANNpointArray points;
    int counter = 0;

    points = annAllocPts(third_view->total_sample_points, 2);
    idx_array.clear();
    idx_array.resize(static_cast<unsigned long>(third_view->total_sample_points));

    for (int i = 0; i < third_view->curve_segments.size(); i++)
    {
        for (int j = 0; j < third_view->curve_segments[i].size(); j++)
        {
            idx_array[counter].s = i;
            idx_array[counter].p = j;
            points[counter][0] = static_cast<double>(third_view->curve_segments[i][j].x);
            points[counter][1] = static_cast<double>(third_view->curve_segments[i][j].y);
            counter ++;
        }
    }

    kd_tree = std::make_shared<ANNkd_tree>(points, third_view->total_sample_points, 2);
}

double ReconPairs::find_nearest_point(const Vector4d& target_point, int& idx)
{
    ANNpoint query_point;
    ANNidxArray idxs;
    ANNdistArray dists;
    double eps = 0.0;

    query_point = annAllocPt(2);
    idxs = new ANNidx[1];
    dists = new ANNdist[1];

    Vector3d temp = third_view->P * target_point;
    temp /= temp[2];

    query_point[0] = temp[0];
    query_point[1] = temp[1];

    kd_tree->annkSearch(query_point, 1, idxs, dists, eps);

    idx = idxs[0];
    double distance = dists[0];

    delete[] idxs;
    delete[] dists;

    return distance;

}

void ReconPairs::filter_curves()
{
    //// Fisrt, seive out those segments whose points are too few
    auto itr = candidates.begin();
    while (itr < candidates.end())
    {
        if (itr->curve.size() < MINIMUM_CURVE_SIZE)
        {
            candidates.erase(itr);
        }
        else
        {
            itr++;
        }
    }

    //// Then, calculate the unary confidence score for each curve and filter out those whose cost is greater than threshold

    // The diagonal length of the third view.
    double diag = sqrt(pow(third_view->original_image.rows, 2) + pow(third_view->original_image.cols, 2));
    for (auto &c: candidates)
    {
        double distances = 0;
        int idx;
        for (const auto &p: c.curve)
        {
            distances += find_nearest_point(p, idx);
        }
        c.score = distances / c.curve.size() / diag;
    }

    itr = candidates.begin();
    while (itr < candidates.end())
    {
        if (itr->score > UNARY_THRESHOLD)
        {
            candidates.erase(itr);
        }
        else
        {
            itr++;
        }
    }

    // Return if candidate size are too few.
    if (candidates.size() <= 1)
    {
        return;
    }

    // calculate all start_point end_point and directions for all curve candidates for further use
    for (auto &c : candidates)
    {
        c.s_pt = c.curve.front();
        c.e_pt = c.curve.back();
        c.s_dire = (*(c.curve.begin()+1) - c.s_pt).normalized();
        c.e_dire = (*(c.curve.end() - 2) - c.e_pt).normalized();
    }

    //// Afterwards, compute the binary pairwise cost for each possible pairs and find the unique set with optimization
    int num_of_var;         // number of all variables. denote as n
    vector<double> unary_cost;    // list of unary_cost
    vector<double> binary_cost;   // list of binary_cost, size:n*(n-1), n rows, (n-1) cols. Omit the 0 score for oneself
    vector<vector<int>> unique_list;    // one list can only have one valid candidate. their idx[0]s should be same

    num_of_var = static_cast<int>(candidates.size());
    unary_cost.resize(static_cast<size_t>(num_of_var));
    binary_cost.resize((num_of_var * (num_of_var - 1))/2);
    for (int i = 0; i < num_of_var; i++)
    {
        bool found = false;
        for (auto &l: unique_list)
        {
            if (candidates[l.front()].idx[0] == candidates[i].idx[0])
            {
                found = true;
                l.emplace_back(i);
                break;
            }
        }
        if (not found)
        {
            vector<int> temp = {i};
            unique_list.emplace_back(move(temp));
        }

        unary_cost[i] = candidates[i].score * UNARY_SCALE; // multiply the UNARY_SCALE as compemsation for lambda in the paper

        for (int j = 0; j < num_of_var - i - 1; j++)
        {
            int idx = (i * (2*num_of_var - 1 - i)) / 2 + j;
            double dists[4], angles[4];

            dists[0] = (candidates[i].s_pt - candidates[j].s_pt).norm();
            dists[1] = (candidates[i].s_pt - candidates[j].e_pt).norm();
            dists[2] = (candidates[i].e_pt - candidates[j].s_pt).norm();
            dists[3] = (candidates[i].e_pt - candidates[j].e_pt).norm();

            angles[0] = candidates[i].s_dire.adjoint() * candidates[j].s_dire;
            angles[1] = candidates[i].s_dire.adjoint() * candidates[j].e_dire;
            angles[2] = candidates[i].e_dire.adjoint() * candidates[j].s_dire;
            angles[3] = candidates[i].e_dire.adjoint() * candidates[j].e_dire;

            auto min_idx = min_element(dists, dists+4);
            binary_cost[idx] = *min_idx + MU * (angles[min_idx-dists] + 1) / 2;
        }
    }

    GRBEnv env = GRBEnv();
    vector<GRBVar> var_list;
    vector<GRBVar> vars1, vars2;
    GRBModel MIPmodel = GRBModel(env);

    var_list.resize(num_of_var);
    for (auto &var: var_list)
    {
        var = MIPmodel.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
    vars1.resize(binary_cost.size());
    vars2.resize(binary_cost.size());
    for (int i = 0; i < num_of_var; i++)
    {
        for (int j = 0; j < num_of_var - i - 1; j++)
        {
            int idx = (i * (2 * num_of_var - 1 - i)) / 2 + j;
            vars1[idx] = var_list[i];
            vars2[idx] = var_list[i+j+1];
        }
    }


    // add objective
    GRBQuadExpr quad_obj;

    quad_obj.addTerms(&unary_cost[0], &var_list[0], num_of_var);
    quad_obj.addTerms(&binary_cost[0], &vars1[0], &vars2[0], binary_cost.size());

    MIPmodel.setObjective(quad_obj, GRB_MINIMIZE);

    // add constraints for unique_list
    for (auto &l: unique_list)
    {
        cout << "Constrain No.: " << &l - &unique_list.front() << endl;
        GRBLinExpr constraint;
        for (auto &t: l)
        {
            constraint += var_list[t];
            cout << t << " ";
        }
        cout << endl << endl;
        MIPmodel.addConstr(constraint, GRB_EQUAL, 1.0);
    }

    // Optimize model
    MIPmodel.optimize();

    for(int i = 0; i < num_of_var; i++)
    {
        cout << "Result: " << var_list[i].get(GRB_DoubleAttr_X) << " Index: " << candidates[i].idx[0] << " " << candidates[i].idx[1] << endl;
    }

    //// Finally, use mTSP to solve the problem.

}
