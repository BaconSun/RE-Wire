#include "ui.h"
#include "MatchFinder.h"
#include "ProjectImage.h"
#include "QuickSort.h"
#include <algorithm> 
#include "CurveOperation.h"
#include "ReconCurve.h"

int main(int argc, char *argv[])
{
	if (argc < 10)
	{
		std::cerr << "usage" << std::endl;
		std::cerr << "param1: path" << std::endl;
		std::cerr << "param2: outpath" << std::endl;
		std::cerr << "param3: group_name" << std::endl;
		std::cerr << "param4: select_thre" << std::endl;
		std::cerr << "param5: epiAreaThre" << std::endl;
		std::cerr << "param6: reguLamda" << std::endl;
		std::cerr << "param7: isoDScale" << std::endl;
		std::cerr << "param8: dire_thre" << std::endl;
		std::cerr << "param9: tangentThre" << std::endl;
		std::cerr << "param10: offsize" << std::endl;
		std::cerr << "param11: considerArea" << std::endl;
		std::cerr << "param12: comb" << std::endl;
		return 0;
	}

	int viewNum = 3;
	char type = 'r'; 
	string path = argv[1];//"../../Models/";
	std::cerr << "path " << path << std::endl;
	string outpath = argv[2];//"../../Models/";
	std::cerr << "outpath " << outpath << std::endl;
	string group_name = argv[3];// "turtle_real";
	std::cerr << "group_name " << group_name << std::endl;
	float select_thre = std::stod(argv[4]);
	float epiAreaThre= std::stod(argv[5]);
	float reguLamda= std::stod(argv[6]);
	float isoDScale= std::stod(argv[7]);
	float dire_thre= std::stod(argv[8]);
	float tangentThre = std::stod(argv[9]);
	int offsize= std::stoi(argv[10]);
	float considerArea= std::stoi(argv[11]);
	int comb = std::stoi(argv[12]);
	std::cerr << "select_thre " << select_thre << std::endl;
	std::cerr << "epiAreaThre " << epiAreaThre << std::endl;
	std::cerr << "reguLamda " << reguLamda << std::endl;
	std::cerr << "isoDScale " << isoDScale << std::endl;
	std::cerr << "dire_thre " << dire_thre << std::endl;
	std::cerr << "tangentThre " << tangentThre << std::endl;
	std::cerr << "offsize " << offsize << std::endl;
	std::cerr << "considerArea " << considerArea << std::endl;
	std::cerr << "comb " << comb << std::endl;

	vector<vector<int>> T = CurveOperation::CombinationChoose(comb);
	path += group_name;
	float scale = 0.3;

	/* Load Image camera parameters and pixel information*/
	vector<ProjectImage> Views(viewNum);
	for (int i = 0; i < viewNum; i++) {
		Views[i].LoadView(path, i + 1, type);
		//cout << i<<"\t"<<Views[i].groupIdxStart.size() << endl;
	}

	vector<vector<int>> takenList(Views.size());
	for (int i = 0; i < viewNum; i++) {
		vector<int> taken(Views[i].img2DPoints.size(), 0);
		takenList[i] = taken;
	}

	vector<ReconCurve> totalReconCurveList;
	vector<int> totalReconStartIdx; totalReconStartIdx.push_back(0);
	vector<int> reconNumIter;
	vector<int> partCurveNum;
	vector<int> refIdxIter(T.size());
	vector<vector<vector<int>>> totalTakenList(T.size());
	vector<vector<vector<Vector2i>>> totalMissingSE(T.size());
	ProjectImage::CalculateFirstMissingSE(Views, totalMissingSE[0]);

	ReconCurve::DeleteFileInDire(path);

	int iterNum = T.size();
	//string filecutname = path + "/cut.txt";
	//ofstream filecut(filecutname, ofstream::out);
	for (int x = 0; x < iterNum; x++) {
		/* t2 is views img name idx, t is Views idx */
		vector<int> t2 = T[x]; vector<int> t = t2;
		for (int i = 0; i < t.size(); i++) t[i]--;
		refIdxIter[x] = t[0];
		for (int i = 1; i < viewNum; i++) {
			/* Compute the Fundamental Matrix and tangent points*/
			Views[t[i]].FMtoRefView = MatchFinder::ComputeF(Views[t[0]], Views[t[i]]);
		}
		CurveOperation::FindTangentPoints_epi(Views[t[0]], Views[t[1]], considerArea, tangentThre, epiAreaThre);

		vector<ReconCurve> reconCurveList;
		for (int count = 0; count < Views[t[0]].reconNumInGroup.size(); count++) {
			vector<Vector3f> reconP;
			MatchFinder::FindIntersectForView_epi(count, Views, t, epiAreaThre, reconP);

			vector<int> startPointIdx = CurveOperation::SeperateCurve(x, count, Views[t[1]]);

			/*if (x == 0) {
				for (int p = 0; p < startPointIdx.size(); p++) {
					filecut << startPointIdx[p] + 1 << "\t";
				}
				filecut << endl;
			}*/

			vector<vector<vector<int>>> totalPairString(startPointIdx.size() - 1);
			int numCurve = CurveOperation::GenerateTrail(startPointIdx, Views[t[1]], totalPairString);
			
			float errorMax = 1.5;
			if (type == 'm') errorMax = 1000.0;
			CurveOperation::SaveTrail_epi(startPointIdx, totalPairString, Views, t, reconCurveList, path, x, count, offsize, considerArea, reguLamda, epiAreaThre, errorMax);
			//CurveOperation::SaveAllReconPoints(count, Views[t[0]], Views[t[1]], path, x);

			//if(x==0) UI::ShowProcess(count, path, t2[0], t2[1], Views[t[0]], Views[t[1]], scale); 

			Views[t[1]].refresh_eachPart();
			//getchar();
		}
		//UI::OnlyShowTangentPoints2(path, Views[t[1]], t2[1], considerArea, scale);
		Views[t[1]].refresh_eachIter();
		ReconCurve::SelectReconCurve3(Views, reconCurveList, select_thre, dire_thre, path, x);

		//UI::ShowReconCurveEachIter(x, path, Views, reconCurveList, scale);

		vector<vector<Vector2i>> missingStartEndList(viewNum);
		ReconCurve::FindMissingCurve(t[0], Views, reconCurveList, takenList, missingStartEndList, offsize);
		totalTakenList[x] = takenList;
		if (x < iterNum - 1) totalMissingSE[x + 1] = missingStartEndList;
		//UI::ShowLongMissingCurve(x, path, Views, missingStartEndList, scale);
		//UI::ShowFullMissingCurve(x, path, Views, takenList, scale);

		ReconCurve::SaveSelectCurve(reconCurveList, path, totalReconCurveList, partCurveNum);
		ReconCurve::MapToViews(Views, missingStartEndList);
		reconNumIter.push_back(reconCurveList.size());
		totalReconStartIdx.push_back(totalReconCurveList.size());
	}

	/*string endFileName = path + "/endpoints.off";
	ReconCurve::SaveEndPoints(totalReconCurveList, endFileName);*/
	//ReconCurve::ConnectNearCurves(totalReconCurveList, partCurveNum, connectDist, isoList);

	vector<vector<int>> dummyIdx;
	vector<vector<Vector3f>> dummyReconCurveList;
	string trailfile = path + "/trails.txt";
	ReconCurve::DummyCombinedCurve(totalReconCurveList, dummyReconCurveList, dummyIdx, trailfile);
	vector<int> isoList, notIsoList;
	//ReconCurve::FindIsolateCurves2(dummyReconCurveList, dummyIdx, isoD, isoList, notIsoList);
	ReconCurve::FindIsolateCurves3(dummyReconCurveList, dummyIdx, isoDScale, isoList, notIsoList);
	std::sort(isoList.begin(), isoList.end());

	/*for (int i = 0; i < isoList.size(); i++) {
		string filename = path + "/iso/" + totalReconCurveList[isoList[i]].curveName;
		ReconCurve::SavePoints(totalReconCurveList[isoList[i]].curve, filename);
	}*/

	for (int i = 0; i < notIsoList.size(); i++) {
		/*string filename = path + "/notiso/" + totalReconCurveList[notIsoList[i]].curveName;
		ReconCurve::SavePoints(totalReconCurveList[notIsoList[i]].curve, filename);*/
		string filename2 = path + "/notiso_color" + "/" + totalReconCurveList[notIsoList[i]].curveName;
		ReconCurve::SaveColorPoints(totalReconCurveList[notIsoList[i]].curve, filename2, totalReconCurveList[notIsoList[i]].color);
	}
	ReconCurve::IsoMapToInfo(totalReconCurveList, isoList, totalReconStartIdx, reconNumIter, partCurveNum);

	vector<int> combineStartIdx;

	/*This part has a bug*/
	for (int i = 0; i < iterNum; i++) {
		int ref = refIdxIter[i];
		vector<ReconCurve> totalReconIterList; totalReconIterList.assign(totalReconCurveList.begin() + totalReconStartIdx[i], totalReconCurveList.begin() + totalReconStartIdx[i + 1]);
		vector<int> gidx = ReconCurve::CombineReconClusterOneIter2(Views[ref], totalTakenList[iterNum - 1][ref], totalTakenList[i][ref], totalReconIterList, totalReconStartIdx[i], totalMissingSE[i][ref], combineStartIdx);
	}

	vector<float> C;
	vector<vector<float>> D, E;
	group_name = group_name + "/";
	ReconCurve::CalculateCDE3(combineStartIdx, totalReconCurveList, Views, partCurveNum, outpath, group_name);
	return 0;
}
