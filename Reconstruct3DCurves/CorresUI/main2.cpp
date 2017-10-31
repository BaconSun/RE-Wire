#include "ui.h"
#include "MatchFinder.h"
#include "ProjectImage.h"
#include "QuickSort.h"
#include <algorithm>
#include "CurveOperation.h"
#include "ReconCurve.h"

/* Parameters for finding intersection points */
float dist_thre = 5.0, points_thre = 60.0, recondist_thre = 0.2;
/* Parameter for connecting intersection points */
float connect_thre = 300.0;
float contiThre = 100.0;
int select_subcurve = -1;
float project_mindist = 100.0;
int offsize = 10;
float onecurve_dist = 60.0;
float select_thre = 5.0;

int main()
{
	vector<int> viewIdx = { 1, 2, 3 };
	char type = 'm';
	string path = "E:/SfmRecon/EdgeDetection/";
	string group_name = "elephant";
	
	path += group_name;
	int viewNum = viewIdx.size();
	/* Load Image camera parameters and pixel information*/
	int offstep = 100;
	vector<ProjectImage> Views(viewNum);
	for (int i = 0; i < viewNum; i++) {
		Views[i].LoadView(path, viewIdx[i], type);
		/* Compute the Fundamental Matrix and tangent points*/
		if (i != 0) {
			Views[i].FMtoRefView = MatchFinder::ComputeF(Views[0], Views[i]);
			Views[i].tangentPoints = CurveOperation::FindTangentPoints(Views[0], Views[i], offstep);
		}
	}

	vector<ReconCurve> reconCurveList;
	for (int count = 0; count < Views[0].pointsNumInGroup.size(); count++) {
		MatchFinder::FindIntersectForView(count, Views[0], Views[1], Views, project_mindist);

		vector<int> startPointIdx = CurveOperation::SeperateCurve(Views[1], contiThre);

		vector<vector<vector<int>>> totalPairString(startPointIdx.size() - 1);
		vector<vector<int>> totalLength(startPointIdx.size() - 1);
		int numCurve = CurveOperation::GenerateTrail(startPointIdx, Views[1], totalPairString, totalLength);

		string path1 = path + "/recon/string";
		CurveOperation::SaveTrail(numCurve, startPointIdx, totalPairString, Views[1], Views[0], reconCurveList, path1, count, offsize, viewIdx[0]);

		//vector<vector<Vector3f>> newReconList;
		//int numCurve2 = CurveOperation::GenerateLongTrail(startPointIdx, Views, refView, reconCurveList, onecurve_dist, newReconList);
		//CurveOperation::SaveTrail2(newReconList, path, count);

		string path2 = path + "/recon_all/string" + std::to_string(count) + ".off";
		CurveOperation::SaveAllReconPoints(count, Views[0], Views[1], path2);
		
		float scale = 0.5; 
		UI::ShowProcess(count, path, viewIdx[0], viewIdx[1], Views[0], Views[1], scale);

		for (int i = 1; i < viewNum; i++) {
			Views[i].intersect.clear();
			Views[i].intersectNum.clear();
			Views[i].intersectIdxStart.clear();
			Views[i].epipolarLine.clear();
		}
	}
	cout<<reconCurveList.size() << endl;
	ReconCurve::SelectReconCurve(Views, reconCurveList, select_thre);
	cout << reconCurveList.size() << endl;

	//vector<vector<int>> takenList(Views.size());
	//for (int i = 0; i < Views.size(); i++) {
	//	vector<int> taken(Views[i].img2DPoints.size(), 0);
	//	takenList[i] = taken;
	//}
	//vector<vector<Vector2i>> missingStartEndList(Views.size());
	//ReconCurve::FindMissingCurve(Views, reconCurveList, takenList, missingStartEndList, offsize);
	///*Debug*/
	//string path3 = path + "/x1";
	//ReconCurve::SaveSelectCurve(reconCurveList, path3);
	
	return 0;
}
