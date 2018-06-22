//
// Created by bacon on 1/2/18.
//

#ifndef REWIRE_RECONPAIRS_H
#define REWIRE_RECONPAIRS_H

#include <memory>
#include "ProjectImage.h"
#include <ANN/ANN.h>

struct RectifiedPoint
{
    Eigen::Vector3d anchor; // Vector3i
    Eigen::Vector3d center;
    int xmin, xmax, y;
    int idx;
//    std::vector<Eigen::Vector3d> avatars;

    bool operator==(const RectifiedPoint &) const;
    void add(const Eigen::Vector3d&);
    void set_center();
    explicit RectifiedPoint(const Eigen::Vector3d&);
};

struct Curve
{
    std::vector<Eigen::Vector4d> curve;   // 3D coordinates of points on the curve
    std::vector<Eigen::Vector4d> sparse_curve; // 3D coordinates of points, each RectifiedPoint will have only one corresponding point
    std::vector<RectifiedPoint> segment[2];
//    std::vector<Eigen::Vector4d> disparity_segment;
    int idx[2];   // Unique index to mark that there should only one valid candidate for one main_idx
//    void compute_3d(Eigen::Matrix3d Q_inv);
};

struct IndexMap
{
    int s;      // The no. of which segment it belongs to
    int p;      // The no. of the point inside this segment
};

// In this class, 0 means the main view, which is the right view.
// 1 means the neighbouring view, which is the left view.
class ReconPairs
{
public:
    ReconPairs(std::shared_ptr<ProjectImage>, std::shared_ptr<ProjectImage>);
    ReconPairs(std::shared_ptr<ProjectImage>, std::shared_ptr<ProjectImage>, std::shared_ptr<ProjectImage>);
    ~ReconPairs();

    double MAX_DISTANCE_TOLERANCE = 1e-1;

    std::shared_ptr<ProjectImage> view[2];      // the Main view and the Neighbour view
    std::shared_ptr<ProjectImage> third_view;   // the view for validation and selection of 3D curves
    std::shared_ptr<ANNkd_tree> kd_tree;         // the kd-tree is built based on the third view for ANN search
    std::vector<IndexMap> idx_array;                // the index array to map the idx from searching result to third_view->curve_segments
    int total_curve_num;
    std::vector<Curve> candidates;      // Candidates of possible 3D curves

    std::vector<std::vector<Eigen::Vector3d>> curve_segments[2];

    std::vector<std::vector<RectifiedPoint>> rect_segments[2];

    Eigen::Matrix3d fundamental;        // the fundamental matrix of the system with these two cameras
    cv::Mat Q;                          // re-projection matrix. used with disparity and (x,y) pixel coordinates
    Eigen::Matrix3d KR_inv;             // the matrix used to compute the 3d coordinates
    Eigen::Vector3d KT;                 // W{world} = KR_inv * (A{picture}-w * KT). This function used to compute the 3D coordinates
                                        // of the point. W is 3d vector, w is the 4th value of 4d homogeneous coordinate.
    double scale;                       // the scale of the 3d coordinates

    void find_pairs();                  // compute the candidates
    void compute_fundamental();
    void compute_3d();
    void rectify();                     // rectify two images given camera matrix
    void filter_curves();               // filter out the undesired 3D curve candidates
    void set_third_view(std::shared_ptr<ProjectImage>);
    void build_kd_tree();
    double find_nearest_point(const Eigen::Vector3d&, int&);

//    Eigen::Vector3d center[2];
//    Eigen::Vector3d rectified_center[2];
//    Eigen::MatrixXd P_matrix[2];
};


#endif //REWIRE_RECONPAIRS_H
