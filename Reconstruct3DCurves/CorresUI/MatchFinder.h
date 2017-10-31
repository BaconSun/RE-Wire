#pragma once

#include "LibParty.h"
#include "ProjectImage.h"
#include "QuickSort.h"

class MatchFinder
{
public:
	static Vector3f computeClosest3DPoint(Vector2f refPixel, Vector2f neighborPixel, MatrixXf refCam, MatrixXf neighborCam);
	static Matrix3f ComputeF(ProjectImage &refView, ProjectImage &neighborView);
	static int FindIntersectPoints(Vector3f para, vector<Vector2f> &pointsOnCurve2, vector<Vector2f> &intersect, int intersectSum, float dist_thre, float points_thre);
	static int FindIntersectPoints2(Vector3f para, ProjectImage& view);
	static int FindIntersectPoints3(Vector3f para, ProjectImage& view, float p2ldist);
	static bool CheckIntersection(Vector3f para, Vector2f intersect, Vector2f norm, float p2ldist, vector<Vector2f> tangentPoints);
	static void FindIntersectForView(int count, ProjectImage &refView, ProjectImage &neighborView, vector<ProjectImage> views, vector<int> t, float project_mindist);
	static void FindIntersectForView_epi(int count, vector<ProjectImage>& views, vector<int> t, float epiAreaThre, vector<Vector3f> &reconP);
	static int FindIntersectForOnePoint_epi(Vector3f para, ProjectImage view, vector<Vector2f>& intersect);
	static VectorXf SolveForCP(MatrixXf paraBAll, VectorXf parabAll, vector<Vector3f>& updateCP);
	static void ShowTangentCut(string path, ProjectImage view, int vidx, float scale, Vector3f paraStart, Vector3f paraEnd, int vRefIdx, Vector3f refStart, Vector3f refEnd);
};
