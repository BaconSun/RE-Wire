#pragma once
#include "LibParty.h"
#include "ProjectImage.h"
//#include "ui.h"
#include "MatchFinder.h"

class ReconCurve {
public:
	vector<Vector3f> curve;
	int iter;
	int partIdx;
	int clusterIdx;
	int idxInCluster;
	int refIdx;
	int startRef;
	int endRef;
	string curveName;
	Vector3f color;

	ReconCurve() {};
	~ReconCurve() {};

	static void SelectReconCurve(vector<ProjectImage> Views, vector<ReconCurve> &reconCurveList, float select_thre);
	static void SelectReconCurve2(vector<ProjectImage> Views, vector<ReconCurve> &reconCurveList, float select_thre, float dire_thre, string path, int iter);
	static void SelectReconCurve3(vector<ProjectImage> Views, vector<ReconCurve> &reconCurveList, float select_thre, float dire_thre, string path, int iter);
	static void SelectReconCurve4(vector<ProjectImage> Views, vector<ReconCurve> &reconCurveList, float select_thre, float dire_thre, string path, int iter);
	static void SaveSelectCurve(vector<ReconCurve> &reconCurveList, string path3, vector<ReconCurve> &totalReconCurveList, vector<int> &partCurveNum);
	static void SavePoints(vector<Vector3f> points, string path);
	static void SaveColorPoints(vector<Vector3f> points, string filename, Vector3f color);
	static void SaveEndPoints(vector<ReconCurve> reconCurves, string filename);
	static void FindMissingCurve(int refIdx, vector<ProjectImage> Views, vector<ReconCurve> reconCurveList, vector<vector<int>> &takenList, vector<vector<Vector2i>> &missingStartEndList, int offsize);
	static void GenerateCurve(vector<ProjectImage> Views, vector<vector<int>> takenList, vector<vector<Vector2i>> &missingStartEndList, int offsize);
	static void MapToViews(vector<ProjectImage>& Views, vector<vector<Vector2i>> missingStartEndList);
	static vector<int> CombineReconClusterOneIter(ProjectImage view, vector<int> finalTaken, vector<int> curTaken, vector<ReconCurve> reconIter, vector<int>& combineStartIdx);
	static vector<int> CombineReconClusterOneIter2(ProjectImage view, vector<int> finalTaken, vector<int> curTaken, vector<ReconCurve> reconIter, int combineStart, vector<Vector2i> missingSE, vector<int>& combineStartIdx);
	static void ProjPointToTakenList(vector<int> projIdx, ProjectImage view, int viewIdx, vector<vector<int>> &takenList);
	static void CalculateCDE(vector<ReconCurve> subcurves, vector<ProjectImage> Views, vector<int> partCurveNum, vector<string> curvesName, string outpath, string solver, string model);
	static void CalculateCDE2(vector<int> combineStartIdx, vector<ReconCurve> subcurves, vector<ProjectImage> Views, vector<int> partCurveNum, string outpath, string solver, string model, string suffix);
	static void CalculateCDE3(vector<int> combineStartIdx, vector<ReconCurve> subcurves, vector<ProjectImage> Views, vector<int> partCurveNum, string outpath, string model);
	void EndDirection(int pointNum, Vector3f& end1, Vector3f& end2);
	static float FindNearestPair(vector<Vector3f> curv1ends, vector<Vector3f> curv2ends, vector<int>& pair);
	static void FindIsolateCurves(vector<ReconCurve> &totalReconCurveList, vector<int> &totalReconStartIdx, float isoD, vector<int>& isoList);
	static void FindIsolateCurves2(vector<vector<Vector3f>> &dummyReconCurveList, vector<vector<int>> dummyIdx, float isoD, vector<int>& isoList, vector<int>& notIsoList);
	static void FindIsolateCurves3(vector<vector<Vector3f>> &dummyReconCurveList, vector<vector<int>> dummyIdx, float isoDscale, vector<int>& isoList, vector<int>& notIsoList);
	static void DummyCombinedCurve(vector<ReconCurve> totalReconCurveList, vector<vector<Vector3f>> &dummyReconCurveList, vector<vector<int>> &dummyIdx, string filename);
	static void IsoMapToInfo(vector<ReconCurve> &totalReconCurveList, vector<int> isoList, vector<int>& totalReconStartIdx, vector<int>& reconNumIter, vector<int>& partCurveNum);
	static void DeleteFileInDire(string path);
	static vector<Vector3f> FittingEpiArea(string path, int iter, int which, int nn, Vector3f startPoint, Vector3f endPoint, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color);
	static vector<Vector3f> FittingEpiArea2(string path, int iter, int which, int nn, Vector3f startPoint, Vector3f endPoint, vector<ProjectImage> views, vector<int> t, int numori, float lamda, int fittingIterNum, Vector3f color);
	static bool FittingEpiArea3(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, Vector3f startPoint, Vector3f endPoint, int tangentIdx, Vector3f para, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color);
	static bool FittingEpiArea32(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, Vector3f startPoint, Vector3f endPoint, int tangentIdx, Vector3f para, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color, float errorMax);
	static bool FittingEpiArea5(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, int startV1Idx, int endV1Idx, Vector3f startPoint, Vector3f endPoint, Vector3f para, int epiAreaIdx, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color);
	static bool FittingEpiArea6(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, int startV1Idx, int endV1Idx, Vector3f startPoint, Vector3f endPoint, Vector3f para, int epiAreaIdx, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color);
	static bool FittingEpiArea7(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, int startV1Idx, int endV1Idx, Vector3f startPoint, Vector3f endPoint, Vector3f para, int epiAreaIdx, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color);
	static vector<vector<Vector3f>> ReconEpiArea(int start, int startPointIdx, int endPointIdx, vector<ProjectImage> &views, vector<int> t, float considerArea, float epiAreaThre, string path, int iter, int count, int which, int num, float lamda, int fittingIterNum, Vector3f color, float errorMax);
};