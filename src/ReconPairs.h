//
// Created by bacon on 1/2/18.
//

#ifndef REWIRE_RECONPAIRS_H
#define REWIRE_RECONPAIRS_H

#include <memory>
#include "ProjectImage.h"

struct Curve
{
    std::vector<Eigen::VectorXd> curve;   // 3D coordinates of points on the curve
    std::vector<Eigen::Vector3d> segment[2];
    std::vector<Eigen::Vector4d> disparity_segment;
    int idx[2];   // Unique index to mark that there should only one valid candidate for one main_idx
//    void compute_3d(Eigen::Matrix3d Q_inv);
};

struct HomoPoint
{
    HomoPoint(const cv::Point&, int);
    HomoPoint(Eigen::Vector3d, int);
    Eigen::Vector3d point;  // Vector3i
    int idx;
};

struct RectifiedPoint
{
    Eigen::Vector3d anchor; // Vector3i
    Eigen::Vector3d center;
    double xmin, xmax, y;
    std::vector<Eigen::Vector3d> avatars;
    bool operator==(const RectifiedPoint &) const;
    void add(Eigen::Vector3d);

    explicit RectifiedPoint(Eigen::Vector3d);
};

// In this class, 0 means the main view, which is the right view.
// 1 means the neighbouring view, which is the left view.
class ReconPairs
{
public:
    ReconPairs(std::shared_ptr<ProjectImage>, std::shared_ptr<ProjectImage>);
    ~ReconPairs();

    double MAX_DISTANCE_TOLERANCE = 1e-1;

    std::shared_ptr<ProjectImage> view[2];     // the Main view
    int total_curve_num;
    std::vector<Curve> candidates;      // Candidates of possible 3D curves

    std::vector<std::vector<Eigen::Vector3d>> curve_segments[2];

    std::vector<std::vector<RectifiedPoint>> rect_segments[2];

    Eigen::Matrix3d fundamental;       // the fundamental matrix of the system with these two cameras
    cv::Mat Q;                         // reprojection matrix. used with disparity and (x,y) pixel coordiantes
    Eigen::Matrix3d KR_inv;             // the matrix used to compute the 3d coordinates
    Eigen::Vector3d KT;                 // W{world} = KR_inv * (A{picture}-w * KT). This function used to compute the 3D coordinates
                                        // of the point. W is 3d vector, w is the 4th value of 4d homogeneous coordinate.
    double scale;                      // the scale of the 3d coordinates

    void find_pairs();                 // compute the candidates
    void compute_fundamental();
    void compute_3d();
    void rectify();                    // rectify two images given camera matrix
    void filter_curves(std::shared_ptr<ProjectImage>);      // Given a third view, filter out the undesired 3D curve candidates

//    Eigen::Vector3d center[2];
//    Eigen::Vector3d rectified_center[2];
//    Eigen::MatrixXd P_matrix[2];
};


#endif //REWIRE_RECONPAIRS_H
