#pragma once

#include "LibParty.h"
#include "ProjectImage.h"
#include "ReconCurve.h"

class UI
{
protected:
	static Mat img_a_, img_b_;
	static vector<Point2i> handles_a_, handles_b_;
	static float scale_;
	static vector<Point2i> points_a_, points_b_;
	static vector<int> corres_num_;

public:
	static Matrix3f em_a_;
	static vector<Vector2f> corres_a_, corres_b_;
	static void SetImages(Mat img_a, Mat img_b, float scale, Matrix3f& em_a, vector<Vector2f> &pointsOnCurve, vector<Vector2f> &intersection, vector<int> &intersectStart, vector<int> interNum);
	static void SetImages(Mat img_a, Mat img_b, float scale, Matrix3f& em_a, vector<Vector2f>& corres_a, vector<Vector2f>& corres_b);
	static void Show();
	static void OnlyShow();
	static void OnlyShowEndPoints();
	static void OnlyShowTangentPoints();
	static void MouseEventb(int event, int x, int y, int flags, void* param);
	static void EndPointsOnLine(Vector2f& refpoint, vector<Vector2f>& endpoints, char which);
	static void ShowProcess(int count, string path, int refIdxnum, int neiIdxnum, ProjectImage refView, ProjectImage neighborView, float scaleinv);
	static void ShowLongMissingCurve(int iter, string path, vector<ProjectImage> views, vector<vector<Vector2i>> missingStartEndList, float scale);
	static void ShowFullMissingCurve(int iter, string path, vector<ProjectImage> views, vector<vector<int>> takenList, float scale);
	static void ShowReconCurveEachIter(int iter, string path, vector<ProjectImage> views, vector<ReconCurve> reconCurveList, float scale);
	static void ShowTakenList(string path, int vidx, ProjectImage view, vector<int> gidx, vector<int> curTakenList, float scale);
	static void ZoomImage(string path, vector<ProjectImage> views, float scale);
	static void ShowTangentOnPoints(string path, vector<ProjectImage> views, float scale);
	static void OnlyShowTangentPoints(string path, vector<ProjectImage> views, float scale);
	static void OnlyShowTangentPoints2(string path, ProjectImage view, int vidx, int offstep, float scale);
};