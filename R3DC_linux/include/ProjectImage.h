#pragma once
#include "LibParty.h"
//static int height = 3648, width = 5472;
//static int height = 540, width = 960;

class ProjectImage {
	float focalLength;
	Vector2f cameraCenter;

public:
	std::shared_ptr<ANNkd_tree> kdTree;
	Matrix3f K;
	ProjectImage() {};
	~ProjectImage() { 
	};
	Matrix3f rot;
	Vector3f trans;
	MatrixXf projectionMatrix;
	vector<Vector2f> img2DPoints;
	vector<Vector2f> img2DDires;
	vector<Vector2f> neiImg2DPoints;
	vector<int> pointsNumInGroup;
	vector<int> groupIdxStart;
	vector<int> neiPointsNumInGroup;
	vector<int> neiGroupIdxStart;
	Matrix3f FMtoRefView;
	Vector2f epipole;
	vector<Vector2f> tangentPoints;
	vector<int> tangentSign;
	vector<int> tangentIdx;
	vector<Vector3f> tangentPara;
	vector<Vector2i> tangentLeftRight;
	vector<int> reconNumInGroup;
	vector<int> reconGroupIdxStart;
	vector<int> epiAreaTag; 
	vector<vector<int>> epiAreaIdx;

	vector<Vector2f> intersect;
	vector<int> intersectNum;
	vector<int> intersectIdxStart;

	int removeOneIntersect(int part, int k);
	void LoadView(string path, int Idxnum, char type);
	void LoadCamera(string filename, char type);
	void LoadImagePoints(string filename);
	void LoadImagePoints2(string filename);
	void LoadImagePoints3(string filename);
	void LoadImagePoints4(string filename);
	//void DrawProjectImage(vector<Vector3f> points3D, string name);
	void BuildKdTree(vector<Vector2f> imgPoints);
	float FindClosestPointOnImgEachPoint(Vector3f point3D, Vector2f& nearestPoint, int &idx);
	float FindClosestPointOnImgEachPoint2(Vector3f point3D, Vector2f& nearestPoint, int &idx);
	void FindN_NearestPointOnImgEachPoint(Vector3f point3D, vector<Vector2f>& nearestPoint, vector<int> &idx, int numNearest);
	static void CalculateFirstMissingSE(vector<ProjectImage> views, vector<vector<Vector2i>>& firstMissingSE);
	void RcCalculateB(int num, vector<Vector2f> nearestPoint, MatrixXf& para_B, VectorXf& para_b);
	void RcCalculateEachPointB(Vector2f nP, int i, int num, MatrixXf& para_B, Vector2f& para_b);
	void refresh_eachPart();
	void refresh_eachIter();
};