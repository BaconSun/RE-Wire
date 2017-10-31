#include "ReconCurve.h"

void ReconCurve::SelectReconCurve(vector<ProjectImage> Views, vector<ReconCurve> &reconCurveList, float select_thre) {
	/* Calculate C matrix */
	vector<Vector2f> nearestPoint(Views.size());
	vector<float> dist(Views.size());
	float sumdist, averdist;

	for (int i = reconCurveList.size() - 1; i > -1; i--) {
		sumdist = 0.0;
		for (int j = 0; j < reconCurveList[i].curve.size(); j++) {
			for (int k = 0; k < Views.size(); k++) {
				int idx;
				dist[k] = Views[k].FindClosestPointOnImgEachPoint(reconCurveList[i].curve[j], nearestPoint[k], idx);
				sumdist += dist[k];
			}
		}
		averdist = sqrt(sumdist) / (reconCurveList[i].curve.size()*Views.size());
		if (averdist > select_thre) reconCurveList.erase(reconCurveList.begin() + i);
	}
}

void ReconCurve::SelectReconCurve2(vector<ProjectImage> Views, vector<ReconCurve> &reconCurveList, float select_thre, float dire_thre, string path, int iter) {
	/* Calculate C matrix */
	vector<Vector2f> nearestPoint(Views.size());
	vector<float> dist(Views.size());
	vector<float> sumdist(Views.size());
	vector<float> averdist(Views.size());
	int removeTag;

	for (int i = reconCurveList.size() - 1; i > -1; i--) {
		removeTag = 0;
		for (int j = 0; j < Views.size(); j++) {
			sumdist[j] = 0.0; averdist[j] = 0.0;
		}

		vector<Vector2f> dire(Views.size());
		vector<Vector2f> dire2(Views.size());
		vector<float> dotP(Views.size());
		Vector4f sPoint = Vector4f::Ones(); sPoint.head(3) = reconCurveList[i].curve[0];
		Vector4f ePoint = Vector4f::Ones(); ePoint.head(3) = reconCurveList[i].curve.back();
		for (int k = 0; k < Views.size(); k++) {
			Vector3f sImgPoint = Views[k].projectionMatrix*sPoint;
			sImgPoint = sImgPoint / sImgPoint[2];

			Vector3f eImgPoint = Views[k].projectionMatrix*ePoint;
			eImgPoint = eImgPoint / eImgPoint[2];
			dire[k] = eImgPoint.head(2) - sImgPoint.head(2);
			dire[k].normalize();

			int sIdx;
			float dist = Views[k].FindClosestPointOnImgEachPoint(reconCurveList[i].curve[0], nearestPoint[k], sIdx);
			int eIdx;
			dist = Views[k].FindClosestPointOnImgEachPoint(reconCurveList[i].curve.back(), nearestPoint[k], eIdx);
			dire2[k] = Views[k].img2DPoints[eIdx] - Views[k].img2DPoints[sIdx];
			dire2[k].normalize();

			dotP[k] = dire[k].dot(dire2[k]);
		}

		for (int j = 0; j < reconCurveList[i].curve.size(); j++) {
			for (int k = 0; k < Views.size(); k++) {
				int idx;
				dist[k] = Views[k].FindClosestPointOnImgEachPoint(reconCurveList[i].curve[j], nearestPoint[k], idx);
				sumdist[k] += dist[k];
			}
		}
		for (int k = 0; k < Views.size(); k++) {
			averdist[k] = sqrt(sumdist[k]) / reconCurveList[i].curve.size();
			if (averdist[k] > select_thre) {
				removeTag = 1; break;
			}
		}

		for (int k = 0; k < Views.size(); k++) {
			if (dotP[k] < dire_thre) {
				string wrong = path + "\\wrong\\iter" + std::to_string(iter) + "_v" + std::to_string(k) + "_c" + std::to_string(i) + ".off";
				ReconCurve::SavePoints(reconCurveList[i].curve, wrong);
				removeTag = 1; break;
			}
		}
		if (removeTag == 1) {
			reconCurveList.erase(reconCurveList.begin() + i);
		}
	}
}

void ReconCurve::SelectReconCurve3(vector<ProjectImage> Views, vector<ReconCurve> &reconCurveList, float select_thre, float dire_thre, string path, int iter) {
	vector<int> removeList;

	for (int i = 0; i < reconCurveList.size(); i++) {
		vector<Vector2f> nearestPoint(Views.size());
		vector<float> dist(Views.size());
		vector<float> sumdist(Views.size(), 0.0);
		vector<float> averdist(Views.size(), 0.0);
		//float averdist = 0.0;
		vector<float> dotproducts(Views.size(), 0.0); vector<int> dotNum(Views.size(), 0);
		vector<vector<Vector2f>> projectPointOn2D(Views.size());
		int removeTag = 0;
		//std::cout << removeTag << endl;
		/* To each 3D curve, backproject it to the image */
		for (int j = 0; j < reconCurveList[i].curve.size(); j++) {
			Vector4f temp4D = Vector4f::Ones(); temp4D.head(3) = reconCurveList[i].curve[j];
			for (int k = 0; k < Views.size(); k++) {
				Vector3f temp3D = Views[k].projectionMatrix*temp4D;
				temp3D = temp3D / temp3D[2];
				projectPointOn2D[k].push_back(temp3D.head(2));
			}
		}

		int reconSize = reconCurveList[i].curve.size();
		for (int j = 0; j < reconSize; j++) {
			/* Calculate the tangent of each point on the backprojected curve */
			int normStart = max(j - 3, 0);
			int normEnd = min(j + 3, reconSize);
			vector<Vector2f> projectDires(Views.size(), Vector2f::Zero());
			for (int k = 0; k < Views.size(); k++) {
				for (int t = normStart; t < normEnd - 1; t++) {
					Vector2f dire = projectPointOn2D[k][t + 1] - projectPointOn2D[k][t];
					dire.normalize();
					projectDires[k] += dire;
				}
				projectDires[k].normalize();
			}

			for (int k = 0; k < Views.size(); k++) {
				int idx;
				dist[k] = Views[k].FindClosestPointOnImgEachPoint(reconCurveList[i].curve[j], nearestPoint[k], idx);
				sumdist[k] += dist[k];

				Vector2f dire = Views[k].img2DDires[idx];
				dotproducts[k] += fabs(dire.dot(projectDires[k])); dotNum[k]++;
			}
		}

		for (int k = 0; k < Views.size(); k++) {
			averdist[k] = sqrt(sumdist[k]) / reconCurveList[i].curve.size();
			//if (k == 2) cout << averdist[k] << endl;
			if (averdist[k] > select_thre) {
				string wrong1 = path + "\\wrong_dist\\v" + std::to_string(k) + "_" + reconCurveList[i].curveName;
				ReconCurve::SavePoints(reconCurveList[i].curve, wrong1);
				removeTag = 1;
				removeList.push_back(i); break;
			}
		}

		if (removeTag == 0) {
			for (int k = 0; k < Views.size(); k++) {
				dotproducts[k] = dotproducts[k] / dotNum[k];
				if (dotproducts[k] < dire_thre) {
					//std::cout << k << "\t" << dotproducts[k] << endl;
					string wrong2 = path + "\\wrong\\v" + std::to_string(k) + "_" + reconCurveList[i].curveName;
					//ReconCurve::SavePoints(reconCurveList[i].curve, wrong2);
					removeList.push_back(i); break;
				}
			}
		}
	}
	for (int i = removeList.size() - 1; i > -1; i--) {
		reconCurveList.erase(reconCurveList.begin() + removeList[i]);
	}
}

void ReconCurve::SelectReconCurve4(vector<ProjectImage> Views, vector<ReconCurve> &reconCurveList, float select_thre, float dire_thre, string path, int iter) {
	vector<Vector2f> nearestPoint(Views.size());
	vector<float> dist(Views.size());
	vector<float> sumdist(Views.size());
	vector<float> averdist(Views.size());
	vector<float> dotproducts(Views.size()); vector<int> dotNum(Views.size());
	vector<vector<Vector2f>> projectPointOn2D(Views.size());
	vector<vector<Vector2f>> projectPointOn2D2(Views.size());

	int removeTag;

	/*string filename1 = path + "new_points.txt";
	ifstream file1(filename1, ifstream::in);*/
	for (int i = reconCurveList.size() - 1; i > -1; i--) {
		removeTag = 0;
		for (int j = 0; j < Views.size(); j++) {
			sumdist[j] = 0.0; averdist[j] = 0.0;
		}
		for (int j = 0; j < Views.size(); j++) {
			dotproducts[j] = 0.0; dotNum[j] = 0;
		}

		/* To each 3D curve, backproject it to the image */
		for (int j = 0; j < reconCurveList[i].curve.size(); j++) {
			Vector4f temp4D = Vector4f::Ones(); temp4D.head(3) = reconCurveList[i].curve[j];
			for (int k = 0; k < Views.size(); k++) {
				Vector3f temp3D = Views[k].projectionMatrix*temp4D;
				temp3D = temp3D / temp3D[2];
				projectPointOn2D[k].push_back(temp3D.head(2));
				if (j % 5 == 2) projectPointOn2D2[k].push_back(temp3D.head(2));
			}
		}

		int reconSize = reconCurveList[i].curve.size();
		for (int j = 0; j < reconSize; j++) {
			for (int k = 0; k < Views.size(); k++) {
				int idx;
				dist[k] = Views[k].FindClosestPointOnImgEachPoint(reconCurveList[i].curve[j], nearestPoint[k], idx);
				sumdist[k] += dist[k];
			}
		}

		for (int j = 0; j < projectPointOn2D2[0].size(); j++) {
			/* Calculate the tangent of each point on the backprojected curve */
			int sampleSize = projectPointOn2D2[0].size();
			int normStart = max(j - 3, 0);
			int normEnd = min(j + 3, sampleSize);
			vector<Vector2f> projectDires(Views.size(), Vector2f::Zero());
			for (int k = 0; k < Views.size(); k++) {
				for (int t = normStart; t < normEnd - 1; t++) {
					Vector2f dire = projectPointOn2D2[k][t + 1] - projectPointOn2D2[k][t];
					dire.normalize();
					projectDires[k] += dire;
				}
				projectDires[k].normalize();
			}

			for (int k = 0; k < Views.size(); k++) {
				int idx;
				dist[k] = Views[k].FindClosestPointOnImgEachPoint(reconCurveList[i].curve[5 * j + 2], nearestPoint[k], idx);

				Vector2f dire = Views[k].img2DDires[idx];
				dotproducts[k] += fabs(dire.dot(projectDires[k])); dotNum[k]++;
			}
		}


		for (int k = 0; k < Views.size(); k++) {
			averdist[k] = sqrt(sumdist[k]) / reconCurveList[i].curve.size();
			if (averdist[k] > select_thre) {
				removeTag = 1; break;
			}
		}

		for (int k = 0; k < Views.size(); k++) {
			dotproducts[k] = dotproducts[k] / dotNum[k];
			if (dotproducts[k] < dire_thre) {
				//cout << k << "\t" << dotproducts[k] << endl;
				string wrong = path + "\\wrong\\iter" + std::to_string(iter) + "_v" + std::to_string(k) + "_c" + std::to_string(i) + ".off";
				//ReconCurve::SavePoints(reconCurveList[i].curve, wrong);
				removeTag = 1; break;
			}
		}

		if (removeTag == 1) {
			reconCurveList.erase(reconCurveList.begin() + i);
		}
	}
}

void ReconCurve::SaveSelectCurve(vector<ReconCurve> &reconCurveList, string path, vector<ReconCurve> &totalReconCurveList, vector<int> &partCurveNum) {
	int prePartIdx, preClusterIdx;
	int tempPartNum = 0;
	for (int i = 0; i < reconCurveList.size(); i++) {
		string iter = std::to_string(reconCurveList[i].iter);
		int partIdx = reconCurveList[i].partIdx;
		int clusterIdx = reconCurveList[i].clusterIdx;

		string partIdxStr;
		if (partIdx < 10) partIdxStr = "0" + std::to_string(partIdx);
		else partIdxStr = std::to_string(partIdx);

		string clusterIdxStr;
		if (clusterIdx < 10) clusterIdxStr = "00" + std::to_string(clusterIdx);
		else if (clusterIdx < 100) clusterIdxStr = "0" + std::to_string(clusterIdx);
		else clusterIdxStr = std::to_string(clusterIdx);

		string filename2 = path + "/x1/" + reconCurveList[i].curveName;
		totalReconCurveList.push_back(reconCurveList[i]);
		SavePoints(reconCurveList[i].curve, filename2);
		string filename3 = path + "/x1_color/" + reconCurveList[i].curveName;
		SaveColorPoints(reconCurveList[i].curve, filename3, reconCurveList[i].color);

		if (i == 0)  tempPartNum++;
		else {
			if (partIdx == prePartIdx && clusterIdx == preClusterIdx) tempPartNum++;
			else {
				partCurveNum.push_back(tempPartNum);
				tempPartNum = 1;
			}
		}
		prePartIdx = partIdx;
		preClusterIdx = clusterIdx;
	}
	if (tempPartNum != 0) partCurveNum.push_back(tempPartNum);
}

void ReconCurve::SavePoints(vector<Vector3f> points, string filename) {
	ofstream file(filename, ofstream::out);

	file << "OFF\n";
	file << points.size() << "\t0\t0\n";
	for (int i = 0; i < points.size(); i++) {
		file << points[i][0] << "\t" << points[i][1] << "\t" << points[i][2] << "\n";
	}
	file.close();
}

void ReconCurve::SaveColorPoints(vector<Vector3f> points, string filename, Vector3f color) {
	ofstream file(filename, ofstream::out);

	file << "COFF\n";
	file << points.size() << "\t0\t0\n";
	for (int i = 0; i < points.size(); i++) {
		file << points[i][0] << "\t" << points[i][1] << "\t" << points[i][2] << "\t" << color[0] << "\t" << color[1] << "\t" << color[2] << "\t255" << endl;
	}
	file.close();
}

void ReconCurve::SaveEndPoints(vector<ReconCurve> reconCurves, string filename) {
	ofstream file(filename, ofstream::out);

	file << "OFF\n";
	file << reconCurves.size() * 2 << "\t0\t0\n";
	for (int i = 0; i < reconCurves.size(); i++) {
		file << reconCurves[i].curve.front()[0] << "\t" << reconCurves[i].curve.front()[1] << "\t" << reconCurves[i].curve.front()[2] << "\n";
		file << reconCurves[i].curve.back()[0] << "\t" << reconCurves[i].curve.back()[1] << "\t" << reconCurves[i].curve.back()[2] << "\n";
	}
	file.close();
}

void ReconCurve::FindMissingCurve(int refIdx, vector<ProjectImage> Views, vector<ReconCurve> reconCurveList, vector<vector<int>> &takenList, vector<vector<Vector2i>> &missingStartEndList, int offsize) {
	for (int k = 0; k < Views.size(); k++) {
		if (k != refIdx) {
			for (int i = 0; i < reconCurveList.size(); i++) {
				vector<int> projIdx;
				for (int j = 0; j < reconCurveList[i].curve.size(); j++) {
					Vector2f nearest; int idx;
					float dist = Views[k].FindClosestPointOnImgEachPoint(reconCurveList[i].curve[j], nearest, idx);
					//takenList[k][idx] = 1;
					projIdx.push_back(idx);
				}
				ProjPointToTakenList(projIdx, Views[k], k, takenList);
			}
		}
		else {
			for (int i = 0; i < reconCurveList.size(); i++) {
				int start = reconCurveList[i].startRef;
				int end = reconCurveList[i].endRef;
				for (int j = start; j < end + 1; j++) {
					takenList[k][j] = 1;
				}
			}
		}
	}
	GenerateCurve(Views, takenList, missingStartEndList, offsize);
}


void ReconCurve::ProjPointToTakenList(vector<int> projIdx, ProjectImage view, int viewIdx, vector<vector<int>> &takenList) {
	int pMinOnCurve = 5;
	std::sort(projIdx.begin(), projIdx.end());
	std::unique(projIdx.begin(), projIdx.end());
	vector<int> gisINpi(view.groupIdxStart.size() + 1);
	gisINpi[0] = 0; gisINpi.back() = projIdx.size();
	for (int i = 1; i < view.groupIdxStart.size(); i++) {
		int j;
		for (j = gisINpi[i - 1]; j < projIdx.size(); j++) {
			if (projIdx[j] > view.groupIdxStart[i]) {
				gisINpi[i] = j;
				break;
			}
		}
		if (j == projIdx.size()) gisINpi[i] = projIdx.size();
	}

	for (int i = 0; i < gisINpi.size() - 1; i++) {
		int pm = min(pMinOnCurve, view.pointsNumInGroup[i] / 2);
		if (gisINpi[i + 1] - gisINpi[i] > pm) {
			for (int j = gisINpi[i]; j < gisINpi[i + 1] - 1; j++) {
				int start = projIdx[j];
				int end = projIdx[j + 1];
				if (projIdx[j + 1] - projIdx[j] < pMinOnCurve) {
					for (int k = start; k < end + 1; k++) {
						takenList[viewIdx][k] = 1;
					}
				}
				else {
					takenList[viewIdx][start] = 1;
					takenList[viewIdx][end] = 1;
				}
			}
		}
	}

}

void ReconCurve::GenerateCurve(vector<ProjectImage> Views, vector<vector<int>> takenList, vector<vector<Vector2i>> &missingStartEndList, int offsize) {
	for (int i = 0; i < Views.size(); i++) {
		vector<Vector2i> missingStartEnd;
		for (int j = 0; j < Views[i].pointsNumInGroup.size(); j++) {
			int start = Views[i].groupIdxStart[j];
			Vector2i tempMSE = Vector2i::Zero(); int pre = -1; int tag = 0;
			for (int k = 0; k < Views[i].pointsNumInGroup[j]; k++) {
				if (takenList[i][start + k] == 0) {
					if (pre == -1 && tag == 0) {
						tempMSE(0) = start + k; tag = 1;
					}
					pre = start + k;
				}
				else {
					if (tag == 1) {
						tempMSE(1) = pre;
						pre = -1; tag = 0;
						if (tempMSE(1) - tempMSE(0) > offsize) {
							missingStartEnd.push_back(tempMSE);
							tempMSE = Vector2i::Zero();
						}
					}
				}

			}
			if (tag == 1) {
				tempMSE(1) = start + Views[i].pointsNumInGroup[j] - 1;
				if (tempMSE(1) - tempMSE(0) > offsize) {
					missingStartEnd.push_back(tempMSE);
					//cout << tempMSE(0) << "\t" << tempMSE(1) << endl;
					tempMSE = Vector2i::Zero();
				}
			}
		}
		missingStartEndList[i] = missingStartEnd;
	}
}

void ReconCurve::MapToViews(vector<ProjectImage>& Views, vector<vector<Vector2i>> missingStartEndList) {
	for (int i = 0; i < Views.size(); i++) {
		Views[i].reconNumInGroup.clear();
		Views[i].reconGroupIdxStart.clear();
		for (int j = 0; j < missingStartEndList[i].size(); j++) {
			int start = missingStartEndList[i][j](0);
			int end = missingStartEndList[i][j](1);
			int tempPNIG = end - start + 1;
			Views[i].reconNumInGroup.push_back(tempPNIG);
			Views[i].reconGroupIdxStart.push_back(start);
			Vector2i se; se(0) = start; se(1) = end;
		}
	}
}

vector<int> ReconCurve::CombineReconClusterOneIter(ProjectImage view, vector<int> finalTaken, vector<int> curTaken, vector<ReconCurve> reconIter, vector<int>& combineStartIdx) {

	vector<int> gapPoints;

	for (int i = 0; i < view.pointsNumInGroup.size(); i++) {
		int start = view.groupIdxStart[i];
		for (int j = 0; j < view.pointsNumInGroup[i]; j++) {
			if (j == 0) gapPoints.push_back(start);
			else if (finalTaken[start + j] == 1 && curTaken[start + j] == 0) { gapPoints.push_back(start + j); /*cout << i << "\t";*/ }
		}
	}
	gapPoints.push_back(view.img2DPoints.size());

	/*cout << "gapPoints:\t" << gapPoints.size() << endl;

	cout << reconIter.size() << endl;*/
	vector<int> ctag(reconIter.size());
	for (int i = 0; i < reconIter.size(); i++) {
		int start = (reconIter[i].startRef + reconIter[i].endRef) / 2;
		for (int j = 0; j < gapPoints.size(); j++) {
			if (start < gapPoints[j]) {
				ctag[i] = j - 1; /*cout << ctag[i] << "\t"; */break;
			}
		}
	}
	/*cout << endl;*/
	int pretag;
	for (int i = 0; i < ctag.size(); i++) {
		if (i == 0) {
			combineStartIdx.push_back(i); pretag = ctag[i];
		}
		else {
			if (ctag[i] != pretag) {
				combineStartIdx.push_back(i); pretag = ctag[i];
			}
		}
	}
	return gapPoints;
}

vector<int> ReconCurve::CombineReconClusterOneIter2(ProjectImage view, vector<int> finalTaken, vector<int> curTaken, vector<ReconCurve> reconIter, int combineStart, vector<Vector2i> missingSE, vector<int>& combineStartIdx) {

	vector<int> gapPoints;

	for (int i = 0; i < missingSE.size(); i++) {
		int start = missingSE[i][0];
		int end = missingSE[i][1];
		for (int j = start; j < end + 1; j++) {
			if (j == start || j == end) gapPoints.push_back(j);
			else if (finalTaken[j] == 1 && curTaken[j] == 0) { gapPoints.push_back(j); /*cout << i << "\t";*/ }
		}
	}

	/*cout << "gapPoints:\t" << gapPoints.size() << endl;

	cout << reconIter.size() << endl;*/
	vector<int> ctag(reconIter.size());
	for (int i = 0; i < reconIter.size(); i++) {
		int start = (reconIter[i].startRef + reconIter[i].endRef) / 2;
		for (int j = 0; j < gapPoints.size(); j++) {
			if (start < gapPoints[j]) {
				ctag[i] = j - 1; /*cout << ctag[i] << "\t";*/ break;
			}
		}
	}

	/*cout << endl;*/
	int pretag;
	for (int i = 0; i < ctag.size(); i++) {
		if (i == 0) {
			combineStartIdx.push_back(combineStart + i); pretag = ctag[i];
		}
		else {
			if (ctag[i] != pretag) {
				combineStartIdx.push_back(combineStart + i); pretag = ctag[i];
			}
		}
	}
	return gapPoints;
}

void ReconCurve::CalculateCDE(vector<ReconCurve> subcurves, vector<ProjectImage> Views, vector<int> partCurveNum, vector<string> curvesName, string outpath, string solver, string model) {
	int curveNum = curvesName.size();
	int partNum = partCurveNum.size();
	int viewNum = Views.size();

	vector<float> C(curveNum);
	vector<vector<float>> D;
	vector<float> Drow(curveNum, 0.0);
	for (int i = 0; i < curveNum; i++) {
		D.push_back(Drow);
	}
	vector<vector<float>> E;
	vector<float> Erow(curveNum, 0.0);
	for (int i = 0; i < curveNum; i++) {
		E.push_back(Erow);
	}

	outpath = outpath + solver + model;
	string path1 = outpath + "CurveFileName.txt";
	ofstream curvefile(path1, ofstream::out);
	for (int i = 0; i < curveNum; i++) {
		curvefile << curvesName[i] << endl;
	}
	curvefile.close();

	string oppath = outpath + "outPartNum.txt";
	ofstream opfile(oppath, ofstream::out);
	int tempsum = 0;
	vector<int> partCurveStartIdx(partNum);
	for (int i = 0; i < partNum; i++) {
		opfile << partCurveNum[i] << endl;
		partCurveStartIdx[i] = tempsum;
		tempsum += partCurveNum[i];
	}
	opfile.close();


	vector<Vector2f> nearestPoint(viewNum);
	vector<float> dist(viewNum);
	vector<float> sumdist(curveNum, 0.0);

	for (int i = 0; i < curveNum; i++) {
		for (int j = 0; j < subcurves[i].curve.size(); j++) {
			for (int k = 0; k < viewNum; k++) {
				int idx;
				dist[k] = Views[k].FindClosestPointOnImgEachPoint(subcurves[i].curve[j], nearestPoint[k], idx);
				sumdist[i] += dist[k];
			}
		}
		C[i] = sqrt(sumdist[i]) / (subcurves[i].curve.size());
	}

	int endDireNum = 5;
	vector<Vector3f> endDire(curveNum * 2);
	for (int i = 0; i < curveNum; i++) {
		subcurves[i].EndDirection(endDireNum, endDire[2 * i], endDire[2 * i + 1]);
	}

	for (int i = 0; i < partNum; i++) {
		for (int j = 0; j < partCurveNum[i]; j++) {
			vector<Vector3f> curv1ends(2);
			curv1ends[0] = subcurves[partCurveStartIdx[i] + j].curve.front();
			curv1ends[1] = subcurves[partCurveStartIdx[i] + j].curve.back();
			for (int t = 0; t < partNum; t++) {
				for (int k = 0; k < partCurveNum[t]; k++) {
					if (i != t) {
						vector<Vector3f> curv2ends(2);
						vector<int> pair(2);
						curv2ends[0] = subcurves[partCurveStartIdx[t] + k].curve.front();
						curv2ends[1] = subcurves[partCurveStartIdx[t] + k].curve.back();
						float Ddist = ReconCurve::FindNearestPair(curv1ends, curv2ends, pair);

						Vector3f gapdire = curv2ends[pair[1]] - curv1ends[pair[0]];
						if (gapdire != Vector3f::Zero()) gapdire = gapdire.normalized();
						float dotproduct = (fabs(endDire[2 * (partCurveStartIdx[i] + j) + pair[0]].dot(gapdire)) + fabs(endDire[2 * (partCurveStartIdx[t] + k) + pair[1]].dot(gapdire))) / 2.0;
						E[partCurveStartIdx[i] + j][partCurveStartIdx[t] + k] = exp(-dotproduct);

						D[partCurveStartIdx[i] + j][partCurveStartIdx[t] + k] = Ddist;
					}

				}
			}

		}
	}

	string outC = outpath + "outC.txt";
	ofstream outfC(outC, ofstream::out);
	for (int i = 0; i < C.size(); i++) {
		outfC << C[i] << endl;
		/*cout << C[i] << endl;*/
	}
	outfC.close();

	string outE = outpath + "outE.txt";
	ofstream outfE(outE, ofstream::out);
	string outD = outpath + "outD.txt";
	ofstream outfD(outD, ofstream::out);
	for (int i = 0; i < C.size(); i++) {
		for (int j = 0; j < C.size(); j++) {
			outfE << E[i][j] << "\t";
			outfD << D[i][j] << "\t";
		}
		outfE << endl;
		outfD << endl;
	}
	outfE.close();
	outfD.close();
}

void ReconCurve::CalculateCDE2(vector<int> combineStartIdx, vector<ReconCurve> subcurves, vector<ProjectImage> Views, vector<int> partCurveNum, string outpath, string solver, string model, string suffix) {
	int curveNum = subcurves.size();
	int partNum = partCurveNum.size();
	int viewNum = Views.size();

	vector<float> C(curveNum);
	vector<vector<float>> D;
	vector<float> Drow(curveNum, 0.0);
	for (int i = 0; i < curveNum; i++) {
		D.push_back(Drow);
	}
	vector<vector<float>> E;
	vector<float> Erow(curveNum, 0.0);
	for (int i = 0; i < curveNum; i++) {
		E.push_back(Erow);
	}

	outpath = outpath + model;
	string path0 = outpath + "CombineStartIdx" + suffix + ".txt";
	ofstream combinefile(path0, ofstream::out);
	for (int i = 0; i < combineStartIdx.size(); i++) {
		combinefile << combineStartIdx[i] << endl;
	}
	combinefile.close();

	string path1 = outpath + "CurveFileName" + suffix + ".txt";
	ofstream curvefile(path1, ofstream::out);
	for (int i = 0; i < curveNum; i++) {
		curvefile << subcurves[i].curveName << endl;
	}
	curvefile.close();

	string oppath = outpath + "outPartNum" + suffix + ".txt";
	ofstream opfile(oppath, ofstream::out);
	int tempsum = 0;
	vector<int> partCurveStartIdx(partNum);
	for (int i = 0; i < partNum; i++) {
		opfile << partCurveNum[i] << endl;
		partCurveStartIdx[i] = tempsum;
		tempsum += partCurveNum[i];
	}
	opfile.close();


	vector<Vector2f> nearestPoint(viewNum);
	vector<float> dist(viewNum);
	vector<vector<float>> viewDists(curveNum);
	vector<float> sumdist(curveNum, 0.0);

	for (int i = 0; i < curveNum; i++) {
		vector<float> tempDists(viewNum, 0.0);
		for (int j = 0; j < subcurves[i].curve.size(); j++) {
			for (int k = 0; k < viewNum; k++) {
				int idx;
				dist[k] = Views[k].FindClosestPointOnImgEachPoint(subcurves[i].curve[j], nearestPoint[k], idx);
				tempDists[k] += dist[k];
				sumdist[i] += dist[k];
			}
		}
		C[i] = sqrt(sumdist[i] / subcurves[i].curve.size());
		viewDists[i] = tempDists;
		for (int k = 0; k < viewNum; k++) {
			viewDists[i][k] = sqrt(viewDists[i][k] / subcurves[i].curve.size());
		}
	}

	int endDireNum = 5;
	vector<Vector3f> endDire(curveNum * 2);
	for (int i = 0; i < curveNum; i++) {
		subcurves[i].EndDirection(endDireNum, endDire[2 * i], endDire[2 * i + 1]);
	}

	for (int i = 0; i < partNum; i++) {
		for (int j = 0; j < partCurveNum[i]; j++) {
			vector<Vector3f> curv1ends(2);
			curv1ends[0] = subcurves[partCurveStartIdx[i] + j].curve.front();
			curv1ends[1] = subcurves[partCurveStartIdx[i] + j].curve.back();
			for (int t = 0; t < partNum; t++) {
				for (int k = 0; k < partCurveNum[t]; k++) {
					if (i != t) {
						vector<Vector3f> curv2ends(2);
						vector<int> pair(2);
						curv2ends[0] = subcurves[partCurveStartIdx[t] + k].curve.front();
						curv2ends[1] = subcurves[partCurveStartIdx[t] + k].curve.back();
						float Ddist = ReconCurve::FindNearestPair(curv1ends, curv2ends, pair);

						Vector3f direCurve1, direCurve2;
						direCurve1 = endDire[2 * (partCurveStartIdx[i] + j) + pair[0]];
						if (pair[1] == pair[0]) {
							direCurve2 = -endDire[2 * (partCurveStartIdx[t] + k) + pair[1]];
						}
						else {
							direCurve2 = endDire[2 * (partCurveStartIdx[t] + k) + pair[1]];
						}

						float dotproduct = direCurve2.dot(direCurve1);
						E[partCurveStartIdx[i] + j][partCurveStartIdx[t] + k] = (1.0 - dotproduct) / 2.0;
						/*Vector3f gapdire = curv2ends[pair[1]] - curv1ends[pair[0]];
						if (gapdire != Vector3f::Zero()) gapdire = gapdire.normalized();
						float dotproduct = (fabs(endDire[2 * (partCurveStartIdx[i] + j) + pair[0]].dot(gapdire)) + fabs(endDire[2 * (partCurveStartIdx[t] + k) + pair[1]].dot(gapdire))) / 2.0;*/
						//E[partCurveStartIdx[i] + j][partCurveStartIdx[t] + k] = exp(-dotproduct);
						D[partCurveStartIdx[i] + j][partCurveStartIdx[t] + k] = Ddist;
					}
				}
			}

		}
	}

	string outC = outpath + "outC" + suffix + ".txt";
	ofstream outfC(outC, ofstream::out);
	for (int i = 0; i < C.size(); i++) {
		outfC << C[i] << endl;
		/*cout << C[i] << endl;*/
	}
	outfC.close();

	string outCV = outpath + "outCEachView" + suffix + ".txt";
	ofstream outfCV(outCV, ofstream::out);
	for (int i = 0; i < curveNum; i++) {
		for (int k = 0; k < viewNum; k++) {
			outfCV << viewDists[i][k] << "\t";
		}
		outfCV << endl;
	}
	outfCV.close();

	string outE = outpath + "outE" + suffix + ".txt";
	ofstream outfE(outE, ofstream::out);
	string outD = outpath + "outD" + suffix + ".txt";
	ofstream outfD(outD, ofstream::out);
	for (int i = 0; i < C.size(); i++) {
		for (int j = 0; j < C.size(); j++) {
			outfE << E[i][j] << "\t";
			outfD << D[i][j] << "\t";
		}
		outfE << endl;
		outfD << endl;
	}
	outfE.close();
	outfD.close();
}

void ReconCurve::CalculateCDE3(vector<int> combineStartIdx, vector<ReconCurve> subcurves, vector<ProjectImage> Views, vector<int> partCurveNum, string outpath, string model) {
	int curveNum = subcurves.size();
	int partNum = partCurveNum.size();
	int viewNum = Views.size();

	struct stat status;
	stat(outpath.c_str(), &status);

	if (!(status.st_mode & S_IFDIR)) {
		mkdir(outpath.c_str());
	}
	else {
		cout << "exists" << endl;
	}

	outpath = outpath + model;

	struct stat status2;
	stat(outpath.c_str(), &status2);

	if (!(status2.st_mode & S_IFDIR)) {
		mkdir(outpath.c_str());
	}
	else {
		cout << "exists2" << endl;
	}

	vector<float> C(curveNum);
	vector<vector<float>> D;
	vector<float> Drow(curveNum, 0.0);
	for (int i = 0; i < curveNum; i++) {
		D.push_back(Drow);
	}
	vector<vector<float>> E;
	vector<float> Erow(curveNum, 0.0);
	for (int i = 0; i < curveNum; i++) {
		E.push_back(Erow);
	}

	string path0 = outpath + "CombineStartIdx" + ".txt";
	ofstream combinefile(path0, ofstream::out);
	for (int i = 0; i < combineStartIdx.size(); i++) {
		combinefile << combineStartIdx[i] << endl;
	}
	combinefile.close();

	string path1 = outpath + "CurveFileName" + ".txt";
	ofstream curvefile(path1, ofstream::out);
	for (int i = 0; i < curveNum; i++) {
		curvefile << subcurves[i].curveName << endl;
	}
	curvefile.close();

	string oppath = outpath + "outPartNum" + ".txt";
	ofstream opfile(oppath, ofstream::out);
	for (int i = 0; i < partNum; i++) {
		opfile << partCurveNum[i] << endl;
	}
	opfile.close();


	vector<Vector2f> nearestPoint(viewNum);
	vector<float> dist(viewNum);
	vector<vector<float>> viewDists(curveNum);
	vector<float> sumdist(curveNum, 0.0);

	for (int i = 0; i < curveNum; i++) {
		vector<float> tempDists(viewNum, 0.0);
		for (int j = 0; j < subcurves[i].curve.size(); j++) {
			for (int k = 0; k < viewNum; k++) {
				int idx;
				dist[k] = Views[k].FindClosestPointOnImgEachPoint(subcurves[i].curve[j], nearestPoint[k], idx);
				tempDists[k] += dist[k];
				sumdist[i] += dist[k];
			}
		}
		C[i] = sqrt(sumdist[i] / subcurves[i].curve.size());
		viewDists[i] = tempDists;
		for (int k = 0; k < viewNum; k++) {
			viewDists[i][k] = sqrt(viewDists[i][k] / subcurves[i].curve.size());
		}
	}

	int endDireNum = 5;
	vector<Vector3f> endDire(curveNum * 2);
	for (int i = 0; i < curveNum; i++) {
		subcurves[i].EndDirection(endDireNum, endDire[2 * i], endDire[2 * i + 1]);
	}

	for (int i = 0; i < curveNum; i++) {
		vector<Vector3f> curv1ends(2);
		curv1ends[0] = subcurves[i].curve.front();
		curv1ends[1] = subcurves[i].curve.back();
		for (int j = 0; j < curveNum; j++) {
			if (i != j) {
				vector<Vector3f> curv2ends(2);
				vector<int> pair(2);
				curv2ends[0] = subcurves[j].curve.front();
				curv2ends[1] = subcurves[j].curve.back();
				float Ddist = ReconCurve::FindNearestPair(curv1ends, curv2ends, pair);

				Vector3f direCurve1, direCurve2;
				direCurve1 = endDire[2 * i + pair[0]];
				if (pair[1] == pair[0]) {
					direCurve2 = -endDire[2 * j + pair[1]];
				}
				else {
					direCurve2 = endDire[2 * j + pair[1]];
				}

				float dotproduct = direCurve2.dot(direCurve1);
				E[i][j] = (1.0 - dotproduct) / 2.0;
				D[i][j] = Ddist;
			}
		}
	}

	string outC = outpath + "outC" + ".txt";
	ofstream outfC(outC, ofstream::out);
	for (int i = 0; i < C.size(); i++) {
		outfC << C[i] << endl;
	}
	outfC.close();

	string outCV = outpath + "outCEachView" + ".txt";
	ofstream outfCV(outCV, ofstream::out);
	for (int i = 0; i < curveNum; i++) {
		for (int k = 0; k < viewNum; k++) {
			outfCV << viewDists[i][k] << "\t";
		}
		outfCV << endl;
	}
	outfCV.close();

	string outE = outpath + "outE" + ".txt";
	ofstream outfE(outE, ofstream::out);
	string outD = outpath + "outD" + ".txt";
	ofstream outfD(outD, ofstream::out);
	for (int i = 0; i < C.size(); i++) {
		for (int j = 0; j < C.size(); j++) {
			outfE << E[i][j] << "\t";
			outfD << D[i][j] << "\t";
		}
		outfE << endl;
		outfD << endl;
	}
	outfE.close();
	outfD.close();
}

void ReconCurve::EndDirection(int pointNum, Vector3f& end1, Vector3f& end2) {
	end1 = Vector3f::Zero();
	end2 = Vector3f::Zero();
	Vector3f tempdire;
	if (curve.size() < pointNum + 1) {
		for (int i = 0; i < curve.size() - 1; i++) {
			tempdire = curve[i + 1] - curve[i];
			if (tempdire != Vector3f::Zero()) tempdire = tempdire.normalized();
			end1 += tempdire;
		}
		if (end1 != Vector3f::Zero()) end1 = end1.normalized();
		end2 = end1;
	}
	else {
		for (int i = 0; i < pointNum; i++) {
			tempdire = curve[i + 1] - curve[i];
			if (tempdire != Vector3f::Zero()) tempdire = tempdire.normalized();
			end1 += tempdire;
		}
		if (end1 != Vector3f::Zero()) end1 = end1.normalized();

		for (int i = curve.size() - 1 - pointNum; i < curve.size() - 1; i++) {
			tempdire = curve[i + 1] - curve[i];
			if (tempdire != Vector3f::Zero()) tempdire = tempdire.normalized();
			end2 += tempdire;
		}
		if (end2 != Vector3f::Zero()) end2 = end2.normalized();
	}
}

float ReconCurve::FindNearestPair(vector<Vector3f> curv1ends, vector<Vector3f> curv2ends, vector<int>& pair) {
	float min_dist;
	for (int i = 0; i < curv1ends.size(); i++) {
		for (int j = 0; j < curv2ends.size(); j++) {
			float dist = (curv1ends[i] - curv2ends[j]).norm();
			if (i == 0 && j == 0) {
				min_dist = dist;
				pair[0] = pair[1] = 0;
			}
			else if (dist < min_dist) {
				min_dist = dist;
				pair[0] = i; pair[1] = j;
			}
		}
	}
	return min_dist;
}

void ReconCurve::FindIsolateCurves(vector<ReconCurve> &totalReconCurveList, vector<int> &totalReconStartIdx, float isoD, vector<int>& isoList) {
	ANNpoint queryPt;
	ANNidxArray nnIdx;
	ANNdistArray dists;
	double eps = 0.0;
	int sampleStep = 10;

	ANNpointArray dataPts;
	dataPts = annAllocPts(totalReconCurveList.size() * 2, 3);
	for (int i = 0; i < totalReconCurveList.size(); i++)
	{
		dataPts[i * 2][0] = (double)totalReconCurveList[i].curve.front()[0];
		dataPts[i * 2][1] = (double)totalReconCurveList[i].curve.front()[1];
		dataPts[i * 2][2] = (double)totalReconCurveList[i].curve.front()[2];
		dataPts[i * 2 + 1][0] = (double)totalReconCurveList[i].curve.back()[0];
		dataPts[i * 2 + 1][1] = (double)totalReconCurveList[i].curve.back()[1];
		dataPts[i * 2 + 1][2] = (double)totalReconCurveList[i].curve.back()[2];
	}
	ANNkd_tree* kdTree = new ANNkd_tree(dataPts, totalReconCurveList.size() * 2, 3);

	queryPt = annAllocPt(3);

	Vector2f nearestPoint; int idx;

	for (int i = 0; i < totalReconCurveList.size(); i++) {
		nnIdx = new ANNidx[3]; dists = new ANNdist[3];
		Vector3f start = totalReconCurveList[i].curve.front();
		queryPt[0] = (double)start[0]; queryPt[1] = (double)start[1]; queryPt[2] = (double)start[2];
		kdTree->annkSearch(queryPt, 3, nnIdx, dists, eps);
		float sdist;
		if (nnIdx[1] == 2 * i + 1)  sdist = dists[2];
		else sdist = dists[1];
		sdist = sqrt(sdist); //cout << sdist << endl;
		//cout << nnIdx[1] << endl;

		Vector3f end = totalReconCurveList[i].curve.back();
		queryPt[0] = (double)end[0]; queryPt[1] = (double)end[1]; queryPt[2] = (double)end[2];
		kdTree->annkSearch(queryPt, 3, nnIdx, dists, eps);
		float edist;
		if (nnIdx[1] == 2 * i)  edist = dists[2];
		else edist = dists[1];
		edist = sqrt(edist); //cout << edist << endl;

		if (sdist > isoD && edist > isoD) isoList.push_back(i);
		delete[] nnIdx;
		delete[] dists;
	}
}

void ReconCurve::FindIsolateCurves2(vector<vector<Vector3f>> &dummyReconCurveList, vector<vector<int>> dummyIdx, float isoD, vector<int>& isoList, vector<int>& notIsoList) {
	ANNpoint queryPt;
	ANNidxArray nnIdx;
	ANNdistArray dists;
	double eps = 0.0;
	int sampleStep = 10;
	vector<int> sampleNumList(dummyReconCurveList.size());
	vector<int> sampleStartIdx(dummyReconCurveList.size() + 1);
	vector<float> curveLength(dummyReconCurveList.size(), 0.0);

	int sN = 0;
	for (int i = 0; i < dummyReconCurveList.size(); i++)
	{
		sampleStartIdx[i] = sN;
		int pointNumOnCurve = dummyReconCurveList[i].size();
		sampleNumList[i] = (pointNumOnCurve - 1) / sampleStep + 1;
		if (sampleNumList[i] == 1)  sN++;
		else if ((pointNumOnCurve - 1) % sampleStep != 0) sN++;
		sN += sampleNumList[i];
	}
	sampleStartIdx.back() = sN;

	ANNpointArray dataPts;
	dataPts = annAllocPts(sN, 3);
	sN = 0;
	for (int i = 0; i < dummyReconCurveList.size(); i++)
	{
		int pointNumOnCurve = dummyReconCurveList[i].size();
		if (sampleNumList[i] == 1) {
			sampleNumList[i] = 2;
			Vector3f start = dummyReconCurveList[i].front();
			Vector3f end = dummyReconCurveList[i].back();
			dataPts[sN][0] = (double)start[0];
			dataPts[sN][1] = (double)start[1];
			dataPts[sN][2] = (double)start[2];
			sN++;
			dataPts[sN][0] = (double)end[0];
			dataPts[sN][1] = (double)end[1];
			dataPts[sN][2] = (double)end[2];
			sN++;
			curveLength[i] += (start - end).norm();
		}
		else {
			Vector3f pre, cur;
			for (int j = 0; j < sampleNumList[i]; j++) {
				cur = dummyReconCurveList[i][j*sampleStep];
				dataPts[sN][0] = (double)cur[0];
				dataPts[sN][1] = (double)cur[1];
				dataPts[sN][2] = (double)cur[2];
				sN++;
				if (j != 0) curveLength[i] += (pre - cur).norm();
				pre = cur;
			}
			if ((pointNumOnCurve - 1) % sampleStep != 0) {
				cur = dummyReconCurveList[i].back();
				sampleNumList[i]++;
				dataPts[sN][0] = (double)cur[0];
				dataPts[sN][1] = (double)cur[1];
				dataPts[sN][2] = (double)cur[2];
				sN++;
				curveLength[i] += (pre - cur).norm();
			}
		}
	}

	ANNkd_tree* kdTree = new ANNkd_tree(dataPts, sN, 3);

	queryPt = annAllocPt(3);

	for (int i = 0; i < dummyReconCurveList.size(); i++) {
		nnIdx = new ANNidx[sampleNumList[i] + 1]; dists = new ANNdist[sampleNumList[i] + 1];
		Vector3f start = dummyReconCurveList[i].front();
		queryPt[0] = (double)start[0]; queryPt[1] = (double)start[1]; queryPt[2] = (double)start[2];
		kdTree->annkSearch(queryPt, sampleNumList[i] + 1, nnIdx, dists, eps);
		float sdist;
		for (int j = 0; j < sampleNumList[i] + 1; j++) {
			if (nnIdx[j] < sampleStartIdx[i] || nnIdx[j] > sampleStartIdx[i + 1] - 1) {
				sdist = dists[j];
				break;
			}
		}
		sdist = sqrt(sdist);

		Vector3f end = dummyReconCurveList[i].back();
		queryPt[0] = (double)end[0]; queryPt[1] = (double)end[1]; queryPt[2] = (double)end[2];
		kdTree->annkSearch(queryPt, sampleNumList[i] + 1, nnIdx, dists, eps);
		float edist;
		for (int j = 0; j < sampleNumList[i] + 1; j++) {
			if (nnIdx[j] < sampleStartIdx[i] || nnIdx[j] > sampleStartIdx[i + 1] - 1) {
				edist = dists[j];
				break;
			}
		}
		edist = sqrt(edist);

		if (sdist > isoD && edist > isoD) {
			for (int j = 0; j < dummyIdx[i].size(); j++) {
				isoList.push_back(dummyIdx[i][j]);
			}
		}
		else {
			for (int j = 0; j < dummyIdx[i].size(); j++) {
				notIsoList.push_back(dummyIdx[i][j]);
			}
		}
		delete[] nnIdx;
		delete[] dists;
	}
}

void ReconCurve::FindIsolateCurves3(vector<vector<Vector3f>> &dummyReconCurveList, vector<vector<int>> dummyIdx, float isoDscale, vector<int>& isoList, vector<int>& notIsoList) {
	ANNpoint queryPt;
	ANNidxArray nnIdx;
	ANNdistArray dists;
	double eps = 0.0;
	int sampleStep = 10;
	vector<int> sampleNumList(dummyReconCurveList.size());
	vector<int> sampleStartIdx(dummyReconCurveList.size() + 1);
	vector<float> curveLength(dummyReconCurveList.size(), 0.0);

	int sN = 0;
	for (int i = 0; i < dummyReconCurveList.size(); i++)
	{
		sampleStartIdx[i] = sN;
		int pointNumOnCurve = dummyReconCurveList[i].size();
		sampleNumList[i] = (pointNumOnCurve - 1) / sampleStep + 1;
		if (sampleNumList[i] == 1)  sN++;
		else if ((pointNumOnCurve - 1) % sampleStep != 0) sN++;
		sN += sampleNumList[i];
	}
	sampleStartIdx.back() = sN;

	ANNpointArray dataPts;
	dataPts = annAllocPts(sN, 3);
	sN = 0;
	for (int i = 0; i < dummyReconCurveList.size(); i++)
	{
		int pointNumOnCurve = dummyReconCurveList[i].size();
		if (sampleNumList[i] == 1) {
			sampleNumList[i] = 2;
			Vector3f start = dummyReconCurveList[i].front();
			Vector3f end = dummyReconCurveList[i].back();
			dataPts[sN][0] = (double)start[0];
			dataPts[sN][1] = (double)start[1];
			dataPts[sN][2] = (double)start[2];
			sN++;
			dataPts[sN][0] = (double)end[0];
			dataPts[sN][1] = (double)end[1];
			dataPts[sN][2] = (double)end[2];
			sN++;
			curveLength[i] += (start - end).norm();
		}
		else {
			Vector3f pre, cur;
			for (int j = 0; j < sampleNumList[i]; j++) {
				cur = dummyReconCurveList[i][j*sampleStep];
				dataPts[sN][0] = (double)cur[0];
				dataPts[sN][1] = (double)cur[1];
				dataPts[sN][2] = (double)cur[2];
				sN++;
				if (j != 0) curveLength[i] += (pre - cur).norm();
				pre = cur;
			}
			if ((pointNumOnCurve - 1) % sampleStep != 0) {
				cur = dummyReconCurveList[i].back();
				sampleNumList[i]++;
				dataPts[sN][0] = (double)cur[0];
				dataPts[sN][1] = (double)cur[1];
				dataPts[sN][2] = (double)cur[2];
				sN++;
				curveLength[i] += (pre - cur).norm();
			}
		}
	}

	ANNkd_tree* kdTree = new ANNkd_tree(dataPts, sN, 3);

	queryPt = annAllocPt(3);

	for (int i = 0; i < dummyReconCurveList.size(); i++) {
		nnIdx = new ANNidx[sampleNumList[i] + 1]; dists = new ANNdist[sampleNumList[i] + 1];
		Vector3f start = dummyReconCurveList[i].front();
		queryPt[0] = (double)start[0]; queryPt[1] = (double)start[1]; queryPt[2] = (double)start[2];
		kdTree->annkSearch(queryPt, sampleNumList[i] + 1, nnIdx, dists, eps);
		float sdist;
		for (int j = 0; j < sampleNumList[i] + 1; j++) {
			if (nnIdx[j] < sampleStartIdx[i] || nnIdx[j] > sampleStartIdx[i + 1] - 1) {
				sdist = dists[j];
				break;
			}
		}
		sdist = sqrt(sdist);

		Vector3f end = dummyReconCurveList[i].back();
		queryPt[0] = (double)end[0]; queryPt[1] = (double)end[1]; queryPt[2] = (double)end[2];
		kdTree->annkSearch(queryPt, sampleNumList[i] + 1, nnIdx, dists, eps);
		float edist;
		for (int j = 0; j < sampleNumList[i] + 1; j++) {
			if (nnIdx[j] < sampleStartIdx[i] || nnIdx[j] > sampleStartIdx[i + 1] - 1) {
				edist = dists[j];
				break;
			}
		}
		edist = sqrt(edist);

		if (sdist > isoDscale*curveLength[i] && edist > isoDscale*curveLength[i]) {
			for (int j = 0; j < dummyIdx[i].size(); j++) {
				isoList.push_back(dummyIdx[i][j]);
			}
		}
		else {
			for (int j = 0; j < dummyIdx[i].size(); j++) {
				notIsoList.push_back(dummyIdx[i][j]);
			}
		}
		delete[] nnIdx;
		delete[] dists;
	}
}

void ReconCurve::DummyCombinedCurve(vector<ReconCurve> totalReconCurveList, vector<vector<Vector3f>> &dummyReconCurveList, vector<vector<int>> &dummyIdx, string filename) {
	float closeThre = 6.0; float farscale = 2.0;
	ANNpoint queryPt;
	ANNidxArray nnIdx;
	ANNdistArray dists;
	double eps = 0.0;

	ANNpointArray dataPts;
	dataPts = annAllocPts(totalReconCurveList.size() * 2, 3);

	for (int i = 0; i < totalReconCurveList.size(); i++) {
		Vector3f start = totalReconCurveList[i].curve.front();
		Vector3f end = totalReconCurveList[i].curve.back();
		dataPts[i * 2][0] = (double)start[0];
		dataPts[i * 2][1] = (double)start[1];
		dataPts[i * 2][2] = (double)start[2];
		dataPts[i * 2 + 1][0] = (double)end[0];
		dataPts[i * 2 + 1][1] = (double)end[1];
		dataPts[i * 2 + 1][2] = (double)end[2];
	}
	ANNkd_tree* kdTree = new ANNkd_tree(dataPts, totalReconCurveList.size() * 2, 3);

	queryPt = annAllocPt(3);

	vector<vector<int>> connectIdx(totalReconCurveList.size());
	for (int i = 0; i < totalReconCurveList.size(); i++) {
		nnIdx = new ANNidx[5]; dists = new ANNdist[5];
		Vector3f start = totalReconCurveList[i].curve.front();
		queryPt[0] = (double)start[0]; queryPt[1] = (double)start[1]; queryPt[2] = (double)start[2];
		kdTree->annkSearch(queryPt, 5, nnIdx, dists, eps);
		vector<int> firstTwoIdx(2); vector<float> firstTwoDist(2); int nn = 0, firstIdx = -1;
		for (int j = 1; j < 5; j++) {
			nnIdx[j] = nnIdx[j] / 2;
			if (nnIdx[j] != i && nnIdx[j] != firstIdx) {
				firstTwoIdx[nn] = nnIdx[j];
				firstTwoDist[nn] = sqrt(dists[j]);
				nn++;
				if (nn == 2) break;
				firstIdx = nnIdx[j];
			}
		}
		int pre;
		if (firstTwoDist[0] < closeThre && firstTwoDist[1]>firstTwoDist[0] * farscale) {
			connectIdx[i].push_back(firstTwoIdx[0]);
			pre = firstTwoIdx[0];
		}

		Vector3f end = totalReconCurveList[i].curve.back();
		queryPt[0] = (double)end[0]; queryPt[1] = (double)end[1]; queryPt[2] = (double)end[2];
		kdTree->annkSearch(queryPt, 5, nnIdx, dists, eps);
		nn = 0, firstIdx = -1;
		for (int j = 1; j < 5; j++) {
			nnIdx[j] = nnIdx[j] / 2;
			if (nnIdx[j] != i && nnIdx[j] != firstIdx) {
				firstTwoIdx[nn] = nnIdx[j];
				firstTwoDist[nn] = sqrt(dists[j]);
				nn++;
				if (nn == 2) break;
				firstIdx = nnIdx[j];
			}
		}
		if (firstTwoIdx[0] != pre && firstTwoDist[0] < closeThre && firstTwoDist[1]>firstTwoDist[0] * farscale) {
			connectIdx[i].push_back(firstTwoIdx[0]);
		}
	}

	vector<int> tag(totalReconCurveList.size(), 0);
	vector<vector<int>> trails;
	for (int i = 0; i < totalReconCurveList.size(); i++) {
		if (tag[i] == 0) {
			if (connectIdx[i].size() > 0) {
				for (int k = 0; k < connectIdx[i].size(); k++) {
					vector<int> trail; int pre;
					tag[i] = 1;
					pre = i;
					trail.push_back(i);
					while (1) {
						int found = 0;
						for (int j = 0; j < connectIdx[pre].size(); j++) {
							if (tag[connectIdx[pre][j]] == 0) {
								trail.push_back(connectIdx[pre][j]);
								tag[connectIdx[pre][j]] = 1;
								pre = connectIdx[pre][j];
								found = 1;
								break;
							}
						}
						if (found == 0) break;
					}
					trails.push_back(trail);
				}
			}
			else {
				vector<int> trail;
				trail.push_back(i);
				trails.push_back(trail);
				tag[i] = 1;
			}
		}
	}

	vector<int> taken(trails.size(), 0);
	for (int i = 0; i < trails.size(); i++) {
		if (taken[i] == 0) {
			int conti = 0;
			for (int j = i + 1; j < trails.size(); j++) {
				if (trails[i][0] == trails[j][0]) {
					vector<int> temptrail = trails[i];
					std::reverse(temptrail.begin(), temptrail.end());
					temptrail.insert(temptrail.end(), trails[j].begin() + 1, trails[j].end());
					dummyIdx.push_back(temptrail);
					conti = 1;
					taken[j] = 1;
					break;
				}
			}
			if (conti == 0) dummyIdx.push_back(trails[i]);
			taken[i] = 1;
		}
	}

	//ofstream file(filename, ofstream::out);
	//for (int i = 0; i < dummyIdx.size(); i++) {
	//	//cout << i << ":\t";
	//	for (int j = 0; j < dummyIdx[i].size(); j++) {
	//		file << dummyIdx[i][j] + 1 << "\t";
	//	}
	//	file << endl;
	//}

	for (int i = 0; i < dummyIdx.size(); i++) {
		if (dummyIdx[i].size() == 1) dummyReconCurveList.push_back(totalReconCurveList[dummyIdx[i][0]].curve);
		else {
			vector<Vector3f> tempCC;
			vector<Vector3f> pre;
			for (int j = 1; j < dummyIdx[i].size(); j++) {
				vector<Vector3f> cur = totalReconCurveList[dummyIdx[i][j]].curve;
				if (j == 1) {
					pre = totalReconCurveList[dummyIdx[i][0]].curve;

					vector<Vector3f> preEnds(2);
					preEnds[0] = pre.front(); preEnds[1] = pre.back();
					vector<Vector3f> curEnds(2);
					curEnds[0] = cur.front(); curEnds[1] = cur.back();
					Vector2i nearestPair; float nearestDist = 1000.0;
					for (int p = 0; p < 2; p++) {
						for (int q = 0; q < 2; q++) {
							float dist = (preEnds[p] - curEnds[q]).norm();
							if (dist < nearestDist) { nearestPair(0) = p; nearestPair(1) = q; nearestDist = dist; }
						}
					}
					if (nearestPair(0) != 1) std::reverse(pre.begin(), pre.end());
					if (nearestPair(1) != 0) std::reverse(cur.begin(), cur.end());
					tempCC.insert(tempCC.end(), pre.begin(), pre.end());
					tempCC.insert(tempCC.end(), cur.begin(), cur.end());
					pre = cur;
				}
				else {
					Vector3f preEnd = pre.back();
					vector<Vector3f> curEnds(2);
					curEnds[0] = cur.front(); curEnds[1] = cur.back();
					vector<float> dists(2);
					for (int q = 0; q < 2; q++) {
						dists[q] = (curEnds[q] - preEnd).norm();
					}
					if (dists[0] > dists[1]) std::reverse(cur.begin(), cur.end());
					tempCC.insert(tempCC.end(), cur.begin(), cur.end());
					pre = cur;
				}
			}
			dummyReconCurveList.push_back(tempCC);
		}
	}
}

void ReconCurve::IsoMapToInfo(vector<ReconCurve> &totalReconCurveList, vector<int> isoList, vector<int>& totalReconStartIdx, vector<int>& reconNumIter, vector<int>& partCurveNum) {
	for (int i = 0; i < isoList.size(); i++) {
		int iso = isoList[i];
		for (int j = 1; j < totalReconStartIdx.size(); j++) {
			if (iso < totalReconStartIdx[j]) {
				reconNumIter[j - 1]--; break;
			}
		}
	}

	for (int i = isoList.size() - 1; i > -1; i--) {
		totalReconCurveList.erase(totalReconCurveList.begin() + isoList[i]);
	}

	totalReconStartIdx.clear();
	totalReconStartIdx.push_back(0);
	int numSum = 0;
	for (int i = 0; i < reconNumIter.size(); i++) {
		numSum += reconNumIter[i];
		totalReconStartIdx.push_back(numSum);
	}

	vector<int> partNumStartIdx; partNumStartIdx.push_back(0);
	numSum = 0;
	for (int i = 0; i < partCurveNum.size(); i++) {
		numSum += partCurveNum[i];
		partNumStartIdx.push_back(numSum);
	}

	for (int i = 0; i < isoList.size(); i++) {
		int iso = isoList[i];
		for (int j = 1; j < partNumStartIdx.size(); j++) {
			if (iso < partNumStartIdx[j]) {
				partCurveNum[j - 1]--;
				break;
			}
		}
	}
}

void ReconCurve::DeleteFileInDire(string path) {
	string strPath = path + "\\notiso_color";


	struct stat status;
	stat(strPath.c_str(), &status);

	if (status.st_mode & S_IFDIR) {
		string snotiso2 = "del " + strPath + "\\*.off /a /s /f /q > log.txt";
		system(snotiso2.c_str());
	}
	else {
		mkdir(strPath.c_str());
	}
}

vector<vector<Vector3f>> ReconCurve::ReconEpiArea(int start, int startPointIdx, int endPointIdx, vector<ProjectImage> &views, vector<int> t, float considerArea, float epiAreaThre, string path, int iter, int count, int which, int num, float lamda, int fittingIterNum, Vector3f color, float errorMax) {
	vector<vector<Vector3f>> reconEpiCurves;

	float eps = 10e-6;

	Vector3f tempStart = Vector3f::Ones();
	tempStart.head(2) = views[t[0]].img2DPoints[start + startPointIdx];
	Vector3f paraStart = views[t[1]].FMtoRefView*tempStart;
	paraStart = paraStart.normalized();

	Vector3f tempEnd = Vector3f::Ones();
	tempEnd.head(2) = views[t[0]].img2DPoints[start + endPointIdx];
	Vector3f paraEnd = views[t[1]].FMtoRefView*tempEnd;
	paraEnd = paraEnd.normalized();

	vector<int> allEpiAreaTag(views[t[1]].tangentPoints.size(), 0);
	vector<int> allEpiAreaIdx;
	for (int i = startPointIdx; i < endPointIdx + 1; i++) {
		for (int j = 0; j < views[t[1]].epiAreaIdx[i].size(); j++) {
			allEpiAreaTag[views[t[1]].epiAreaIdx[i][j]]++;
		}
	}

	for (int i = 0; i < allEpiAreaTag.size(); i++) {
		if (allEpiAreaTag[i] > 5) {
			Vector3f tangentPara = views[t[1]].tangentPara[i];
			Vector2f startTestPoint;
			startTestPoint[0] = views[t[1]].tangentPoints[i][0];
			startTestPoint[1] = -(paraStart[0] * startTestPoint[0] + paraStart[2]) / paraStart[1];
			Vector2f endTestPoint;
			endTestPoint[0] = views[t[1]].tangentPoints[i][0];
			endTestPoint[1] = -(paraEnd[0] * endTestPoint[0] + paraEnd[2]) / paraEnd[1];
			float signStart = tangentPara[0] * startTestPoint[0] + tangentPara[1] * startTestPoint[1] + tangentPara[2];
			float signEnd = tangentPara[0] * endTestPoint[0] + tangentPara[1] * endTestPoint[1] + tangentPara[2];
			if (signStart*views[t[1]].tangentSign[i] > eps && signEnd*views[t[1]].tangentSign[i] > eps) {
				allEpiAreaIdx.push_back(i);
			}
			else {
				allEpiAreaTag[i] = 0;
			}
		}
	}
	vector<Vector2f> intersectStart;
	int internumStart = MatchFinder::FindIntersectForOnePoint_epi(paraStart, views[t[1]], intersectStart);
	vector<int> intersectEpiAreaTagStart(intersectStart.size(), 0);
	vector<Vector2i> intersectForEpiAreaStart(allEpiAreaIdx.size());


	vector<Vector2f> intersectEnd;
	int internumEnd = MatchFinder::FindIntersectForOnePoint_epi(paraEnd, views[t[1]], intersectEnd);
	//vector<int> intersectEpiAreaTagEnd(intersectEnd.size(), 0);
	vector<Vector2i> intersectForEpiAreaEnd(allEpiAreaIdx.size());
	for (int i = 0; i < allEpiAreaIdx.size(); i++) {
		int left = views[t[1]].tangentLeftRight[allEpiAreaIdx[i]][0];
		Vector2f leftPoint = views[t[1]].img2DPoints[left];
		int right = views[t[1]].tangentLeftRight[allEpiAreaIdx[i]][1];
		Vector2f rightPoint = views[t[1]].img2DPoints[right];
		int minIdxLeft, minIdxRight;
		float minDistLeft, minDistRight;
		for (int j = 0; j < intersectStart.size(); j++) {
			float distLeft = (intersectStart[j] - leftPoint).norm();
			if (distLeft < minDistLeft || j == 0) {
				minIdxLeft = j; minDistLeft = distLeft;
			}
			float distRight = (intersectStart[j] - rightPoint).norm();
			if (distRight < minDistRight || j == 0) {
				minIdxRight = j; minDistRight = distRight;
			}
		}
		intersectForEpiAreaStart[i][0] = minIdxLeft;
		intersectForEpiAreaStart[i][1] = minIdxRight;
		intersectEpiAreaTagStart[minIdxLeft] = 1;
		intersectEpiAreaTagStart[minIdxRight] = 1;

		int minIdxLeft2, minDistLeft2, minIdxRight2, minDistRight2;
		for (int j = 0; j < intersectEnd.size(); j++) {
			float distLeft = (intersectEnd[j] - leftPoint).norm();
			if (distLeft < minDistLeft2 || j == 0) {
				minIdxLeft2 = j; minDistLeft2 = distLeft;
			}
			float distRight = (intersectEnd[j] - rightPoint).norm();
			if (distRight < minDistRight2 || j == 0) {
				minIdxRight2 = j; minDistRight2 = distRight;
			}
		}
		intersectForEpiAreaEnd[i][0] = minIdxLeft2;
		intersectForEpiAreaEnd[i][1] = minIdxRight2;
		//intersectEpiAreaTagEnd[minIdxLeft2] = 1;
		//intersectEpiAreaTagEnd[minIdxRight2] = 1;
	}

	int nn = 0;
	for (int i = 0; i < allEpiAreaIdx.size(); i++) {
		Vector3f start3D = MatchFinder::computeClosest3DPoint(views[t[0]].img2DPoints[start + startPointIdx], intersectStart[intersectForEpiAreaStart[i][0]], views[t[0]].projectionMatrix, views[t[1]].projectionMatrix);
		Vector3f end3D = MatchFinder::computeClosest3DPoint(views[t[0]].img2DPoints[start + endPointIdx], intersectEnd[intersectForEpiAreaEnd[i][1]], views[t[0]].projectionMatrix, views[t[1]].projectionMatrix);
		vector<Vector3f> fittingPoints;
		bool RV = ReconCurve::FittingEpiArea32(fittingPoints, path, iter, count, which, nn, start3D, end3D, allEpiAreaIdx[i], paraStart, views, t, num, lamda, fittingIterNum, color, errorMax);
		//bool RV = ReconCurve::FittingEpiArea5(fittingPoints, path, iter, count, which, nn, start + startPointIdx, start + endPointIdx, start3D, end3D, paraStart, allEpiAreaIdx[i], views, t, num, lamda, fittingIterNum, color);
		//bool RV = ReconCurve::FittingEpiArea7(fittingPoints, path, iter, count, which, nn, start + startPointIdx, start + endPointIdx, start3D, end3D, paraStart, allEpiAreaIdx[i], views, t, num, lamda, fittingIterNum, color);
		if (RV == true) {
			reconEpiCurves.push_back(fittingPoints);
			nn++;
		}

		Vector3f start3D2 = MatchFinder::computeClosest3DPoint(views[t[0]].img2DPoints[start + startPointIdx], intersectStart[intersectForEpiAreaStart[i][1]], views[t[0]].projectionMatrix, views[t[1]].projectionMatrix);
		Vector3f end3D2 = MatchFinder::computeClosest3DPoint(views[t[0]].img2DPoints[start + endPointIdx], intersectEnd[intersectForEpiAreaEnd[i][0]], views[t[0]].projectionMatrix, views[t[1]].projectionMatrix);
		vector<Vector3f> fittingPoints2;
		bool RV2 = ReconCurve::FittingEpiArea32(fittingPoints2, path, iter, count, which, nn, start3D2, end3D2, allEpiAreaIdx[i], paraStart, views, t, num, lamda, fittingIterNum, color, errorMax);
		//bool RV2 = ReconCurve::FittingEpiArea5(fittingPoints2, path, iter, count, which, nn, start + startPointIdx, start + endPointIdx, start3D2, end3D2, paraStart, allEpiAreaIdx[i], views, t, num, lamda, fittingIterNum, color);
		//bool RV2 = ReconCurve::FittingEpiArea7(fittingPoints2, path, iter, count, which, nn, start + startPointIdx, start + endPointIdx, start3D2, end3D2, paraStart, allEpiAreaIdx[i], views, t, num, lamda, fittingIterNum, color);
		if (RV2 == true) {
			reconEpiCurves.push_back(fittingPoints2);
			nn++;
		}
	}


	vector<vector<Vector2f>> intersectEpiArea(endPointIdx - startPointIdx + 1);
	vector<Vector2f> tempIntersectStart;
	for (int i = 0; i < intersectEpiAreaTagStart.size(); i++) {
		if (intersectEpiAreaTagStart[i] == 0) tempIntersectStart.push_back(intersectStart[i]);
	}
	intersectEpiArea[0] = tempIntersectStart;

	/*vector<Vector2f> tempIntersectEnd;
	for (int i = 0; i < intersectEpiAreaTagEnd.size(); i++) {
		if (intersectEpiAreaTagEnd[i] == 0) tempIntersectEnd.push_back(intersectEnd[i]);
	}
	intersectEpiArea.back() = tempIntersectEnd;*/

	vector<vector<int>> trails(intersectEpiArea[0].size());
	for (int i = 0; i < intersectEpiArea[0].size(); i++) {
		vector<int> tempTrail(endPointIdx - startPointIdx + 1, -1);
		tempTrail[0] = i;
		trails[i] = tempTrail;
	}

	for (int i = 1; i < endPointIdx - startPointIdx + 1; i++) {
		/* Calculate its epipolar line */
		Vector3f temp = Vector3f::Ones();
		temp.head(2) = views[t[0]].img2DPoints[start + startPointIdx + i];
		Vector3f para = views[t[1]].FMtoRefView*temp;
		//cout << "para1111\t" << para[0] << "\t" << para[1] << "\t" << para[2] << endl;

		int internum = MatchFinder::FindIntersectForOnePoint_epi(para, views[t[1]], intersectEpiArea[i]);

		vector<int> tag(intersectEpiArea[i].size(), 0);
		for (int j = 0; j < trails.size(); j++) {
			int minIdx; float minDist = 10000.0;
			if (trails[j][i - 1] != -1) {
				for (int k = 0; k < intersectEpiArea[i].size(); k++) {
					if (tag[k] == 0) {
						int dist = (intersectEpiArea[i - 1][trails[j][i - 1]] - intersectEpiArea[i][k]).norm();
						if (dist < minDist) {
							minIdx = k;
							minDist = dist;
						}
					}
				}
				if (minDist < 2.0) {
					tag[minIdx] = 1;
					trails[j][i] = minIdx;
				}
			}
		}

	}

	for (int i = 0; i < trails.size(); i++) {
		vector<Vector3f> tempRecon;
		for (int j = 0; j < endPointIdx - startPointIdx + 1; j++) {
			if (trails[i][j] != -1) {
				Vector3f temp = MatchFinder::computeClosest3DPoint(views[t[0]].img2DPoints[start + startPointIdx + j], intersectEpiArea[j][trails[i][j]], views[t[0]].projectionMatrix, views[t[1]].projectionMatrix);
				tempRecon.push_back(temp);
			}
		}
		reconEpiCurves.push_back(tempRecon);
	}

	return reconEpiCurves;
}

vector<Vector3f> ReconCurve::FittingEpiArea(string path, int iter, int which, int nn, Vector3f startPoint, Vector3f endPoint, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color) {
	vector<Vector3f> samplePointOnCurve;

	float step = 1.0 / num;
	for (int i = 0; i < num - 1; i++) {
		Vector3f tempP = (1 - step*i)*startPoint + (step*i)*endPoint;
		samplePointOnCurve.push_back(tempP);
	}
	samplePointOnCurve.push_back(endPoint);

	string beforename = path + "\\before\\iter" + std::to_string(iter) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, beforename, color);
	int count = 0;

	while (count < fittingIterNum) {
		vector<Vector2f> nearestPointList1(num), nearestPointList2(num);
		//VectorXf controlPoints(3 * samplePointOnCurve.size());
		MatrixXf para_B1(2 * (num - 2), 3 * num);
		VectorXf para_b1(2 * (num - 2));
		MatrixXf para_B2(2 * (num - 2), 3 * num);
		VectorXf para_b2(2 * (num - 2));

		/*Find the closest points of the sample points on the image*/
		for (int i = 1; i < num - 1; i++) {
			int idx;
			views[t[0]].FindClosestPointOnImgEachPoint2(samplePointOnCurve[i], nearestPointList1[i], idx);
			views[t[1]].FindClosestPointOnImgEachPoint2(samplePointOnCurve[i], nearestPointList2[i], idx);
		}

		/*Calculate the parameters of the objective function*/
		views[t[0]].RcCalculateB(num, nearestPointList1, para_B1, para_b1);
		views[t[1]].RcCalculateB(num, nearestPointList2, para_B2, para_b2);

		/*MatrixXf reguInb(3 * reconstructPointList.size(), 1);
		for (int i = 0; i < reconstructPointList.size(); i++) {
		reguInb(i * 3) = reconstructPointList[i][0];
		reguInb(i * 3 + 1) = reconstructPointList[i][1];
		reguInb(i * 3 + 2) = reconstructPointList[i][2];
		}*/

		MatrixXf paraBAll(7 * (num - 2) + 6, 3 * num);
		VectorXf parabAll(7 * (num - 2) + 6);

		paraBAll.block(0, 0, 2 * (num - 2), 3 * num) = para_B1;
		paraBAll.block(2 * (num - 2), 0, 2 * (num - 2), 3 * num) = para_B2;

		parabAll.block(0, 0, 2 * (num - 2), 1) = para_b1;
		parabAll.block(2 * (num - 2), 0, 2 * (num - 2), 1) = para_b2;


		/*For Regularization*/

		MatrixXf reguB = MatrixXf::Zero(3 * (num - 2), 3 * num);
		for (int i = 0; i < 3 * (num - 2); i++) {
			reguB(i, i) = lamda * 1; reguB(i, i + 3) = lamda*(-2); reguB(i, i + 6) = lamda * 1;
		}
		VectorXf regub = VectorXf::Zero(3 * (num - 2));

		paraBAll.block(4 * (num - 2), 0, 3 * (num - 2), 3 * num) = reguB;
		parabAll.block(4 * (num - 2), 0, 3 * (num - 2), 1) = regub;

		/*For X1 and Xn*/
		MatrixXf X1B = MatrixXf::Zero(3, 3 * num);
		X1B.block(0, 0, 3, 3) = Matrix3f::Identity();
		Vector3f X1b = samplePointOnCurve[0];

		MatrixXf XnB = MatrixXf::Zero(3, 3 * num);
		XnB.block(0, 3 * (num - 1), 3, 3) = Matrix3f::Identity();
		Vector3f Xnb = samplePointOnCurve[num - 1];

		paraBAll.block(7 * (num - 2), 0, 3, 3 * num) = X1B;
		paraBAll.block(7 * (num - 2) + 3, 0, 3, 3 * num) = XnB;

		parabAll.block(7 * (num - 2), 0, 3, 1) = X1b;
		parabAll.block(7 * (num - 2) + 3, 0, 3, 1) = Xnb;

		VectorXf before(3 * num);
		for (int i = 0; i < num; i++) {
			before[i * 3] = samplePointOnCurve[i][0];
			before[i * 3 + 1] = samplePointOnCurve[i][1];
			before[i * 3 + 2] = samplePointOnCurve[i][2];
		}
		VectorXf errorvec = paraBAll*before - parabAll;
		float error = errorvec.norm();
		//cout << "The error before optimization:\t" << error << endl;

		/*Minimize the objective function, update the location of reconstructed points*/
		VectorXf controlPoints = MatchFinder::SolveForCP(paraBAll, parabAll, samplePointOnCurve);

		VectorXf after(3 * num);
		for (int i = 0; i < num; i++) {
			after[i * 3] = samplePointOnCurve[i][0];
			after[i * 3 + 1] = samplePointOnCurve[i][1];
			after[i * 3 + 2] = samplePointOnCurve[i][2];
		}

		errorvec = paraBAll*controlPoints - parabAll;
		error = errorvec.norm();
		//cout << "The error after optimization:\t" << error << endl;
		count++;
	}

	string aftername = path + "\\after\\iter" + std::to_string(iter) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, aftername, color);
	return samplePointOnCurve;
}

vector<Vector3f> ReconCurve::FittingEpiArea2(string path, int iter, int which, int nn, Vector3f startPoint, Vector3f endPoint, vector<ProjectImage> views, vector<int> t, int numori, float lamda, int fittingIterNum, Vector3f color) {
	vector<Vector3f> samplePointOnCurve;
	int num; int sampleStep = 5;
	if (numori > sampleStep * 5) num = numori / sampleStep;
	else num = numori;
	float step = 1.0 / num;
	for (int i = 0; i < num - 1; i++) {
		Vector3f tempP = (1 - step*i)*startPoint + (step*i)*endPoint;
		samplePointOnCurve.push_back(tempP);
	}
	samplePointOnCurve.push_back(endPoint);

	string beforename = path + "\\before\\iter" + std::to_string(iter) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, beforename, color);
	int count = 0;

	while (count < fittingIterNum) {
		vector<Vector2f> nearestPointList1(num), nearestPointList2(num);
		//VectorXf controlPoints(3 * samplePointOnCurve.size());
		MatrixXf para_B1(2 * (num - 2), 3 * num);
		VectorXf para_b1(2 * (num - 2));
		MatrixXf para_B2(2 * (num - 2), 3 * num);
		VectorXf para_b2(2 * (num - 2));

		/*Find the closest points of the sample points on the image*/
		for (int i = 1; i < num - 1; i++) {
			int idx;
			views[t[0]].FindClosestPointOnImgEachPoint(samplePointOnCurve[i], nearestPointList1[i], idx);
			views[t[1]].FindClosestPointOnImgEachPoint(samplePointOnCurve[i], nearestPointList2[i], idx);
		}

		/*Calculate the parameters of the objective function*/
		views[t[0]].RcCalculateB(num, nearestPointList1, para_B1, para_b1);
		views[t[1]].RcCalculateB(num, nearestPointList2, para_B2, para_b2);

		/*MatrixXf reguInb(3 * reconstructPointList.size(), 1);
		for (int i = 0; i < reconstructPointList.size(); i++) {
		reguInb(i * 3) = reconstructPointList[i][0];
		reguInb(i * 3 + 1) = reconstructPointList[i][1];
		reguInb(i * 3 + 2) = reconstructPointList[i][2];
		}*/

		MatrixXf paraBAll(7 * (num - 2) + 6, 3 * num);
		VectorXf parabAll(7 * (num - 2) + 6);

		paraBAll.block(0, 0, 2 * (num - 2), 3 * num) = para_B1;
		paraBAll.block(2 * (num - 2), 0, 2 * (num - 2), 3 * num) = para_B2;

		parabAll.block(0, 0, 2 * (num - 2), 1) = para_b1;
		parabAll.block(2 * (num - 2), 0, 2 * (num - 2), 1) = para_b2;


		/*For Regularization*/

		MatrixXf reguB = MatrixXf::Zero(3 * (num - 2), 3 * num);
		for (int i = 0; i < 3 * (num - 2); i++) {
			reguB(i, i) = lamda * 1; reguB(i, i + 3) = lamda*(-2); reguB(i, i + 6) = lamda * 1;
		}
		VectorXf regub = VectorXf::Zero(3 * (num - 2));

		paraBAll.block(4 * (num - 2), 0, 3 * (num - 2), 3 * num) = reguB;
		parabAll.block(4 * (num - 2), 0, 3 * (num - 2), 1) = regub;

		/*For X1 and Xn*/
		MatrixXf X1B = MatrixXf::Zero(3, 3 * num);
		X1B.block(0, 0, 3, 3) = Matrix3f::Identity();
		Vector3f X1b = samplePointOnCurve[0];

		MatrixXf XnB = MatrixXf::Zero(3, 3 * num);
		XnB.block(0, 3 * (num - 1), 3, 3) = Matrix3f::Identity();
		Vector3f Xnb = samplePointOnCurve[num - 1];

		paraBAll.block(7 * (num - 2), 0, 3, 3 * num) = X1B;
		paraBAll.block(7 * (num - 2) + 3, 0, 3, 3 * num) = XnB;

		parabAll.block(7 * (num - 2), 0, 3, 1) = X1b;
		parabAll.block(7 * (num - 2) + 3, 0, 3, 1) = Xnb;

		VectorXf before(3 * num);
		for (int i = 0; i < num; i++) {
			before[i * 3] = samplePointOnCurve[i][0];
			before[i * 3 + 1] = samplePointOnCurve[i][1];
			before[i * 3 + 2] = samplePointOnCurve[i][2];
		}
		VectorXf errorvec = paraBAll*before - parabAll;
		float error = errorvec.norm();
		//cout << "The error before optimization:\t" << error << endl;

		/*Minimize the objective function, update the location of reconstructed points*/
		VectorXf controlPoints = MatchFinder::SolveForCP(paraBAll, parabAll, samplePointOnCurve);

		VectorXf after(3 * num);
		for (int i = 0; i < num; i++) {
			after[i * 3] = samplePointOnCurve[i][0];
			after[i * 3 + 1] = samplePointOnCurve[i][1];
			after[i * 3 + 2] = samplePointOnCurve[i][2];
		}

		errorvec = paraBAll*controlPoints - parabAll;
		error = errorvec.norm();
		//cout << "The error after optimization:\t" << error << endl;
		count++;
	}

	vector<Vector3f> newSamplePointOnCurve;
	Vector3f tempSample;
	float sstep = 1.0 / sampleStep;
	if (numori > sampleStep * 5) {
		for (int i = 0; i < samplePointOnCurve.size() - 1; i++) {
			for (int j = 0; j < sampleStep; j++) {
				tempSample = (1 - j*sstep)*samplePointOnCurve[i] + (j*sstep)*samplePointOnCurve[i + 1];
				newSamplePointOnCurve.push_back(tempSample);
			}
		}
		newSamplePointOnCurve.push_back(samplePointOnCurve.back());
	}
	else {
		newSamplePointOnCurve = samplePointOnCurve;
	}
	string aftername = path + "\\after\\iter" + std::to_string(iter) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(newSamplePointOnCurve, aftername, color);
	return newSamplePointOnCurve;
}

bool ReconCurve::FittingEpiArea3(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, Vector3f startPoint, Vector3f endPoint, int tangentIdx, Vector3f para, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color) {
	float errorMax = 1.2;
	bool returnValue = true;
	float eps = 10e-6;
	int numNearest = 100;
	float step = 1.0 / num;
	Vector2f tangentPoint = views[t[1]].tangentPoints[tangentIdx];
	int sign = tangentPoint[0] * para[0] + tangentPoint[1] * para[1] + para[2];

	for (int i = 0; i < num - 1; i++) {
		Vector3f tempP = (1 - step*i)*startPoint + (step*i)*endPoint;
		samplePointOnCurve.push_back(tempP);
	}
	samplePointOnCurve.push_back(endPoint);

	string beforename = path + "\\before\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, beforename, color);
	int iternn = 0;

	while (iternn < fittingIterNum) {
		vector<Vector2f> nearestPointList1(num), nearestPointList2(num);
		//VectorXf controlPoints(3 * samplePointOnCurve.size());
		MatrixXf para_B1(2 * (num - 2), 3 * num);
		VectorXf para_b1(2 * (num - 2));
		MatrixXf para_B2(2 * (num - 2), 3 * num);
		VectorXf para_b2(2 * (num - 2));

		/*Find the closest points of the sample points on the image*/
		for (int i = 1; i < num - 1; i++) {
			int idx1;
			views[t[0]].FindClosestPointOnImgEachPoint2(samplePointOnCurve[i], nearestPointList1[i], idx1);

			int idx2;
			views[t[1]].FindClosestPointOnImgEachPoint2(samplePointOnCurve[i], nearestPointList2[i], idx2);

			/*vector<int> idx2(numNearest); vector<Vector2f> tempNearestPoints2(numNearest);
			views[t[1]].FindN_NearestPointOnImgEachPoint(samplePointOnCurve[i], tempNearestPoints2, idx2, numNearest);
			nearestPointList2[i] = tempNearestPoints2[0];*/
			/*int ftag = 0;
			for (int j = 0; j < numNearest; j++) {
				float tempSign = tempNearestPoints2[j][0] * para[0] + tempNearestPoints2[j][1] * para[1] + para[2];
				if (tempSign*sign > eps) {
					nearestPointList2[i] = tempNearestPoints2[j];
					ftag = 1;
					break;
				}
			}
			if (ftag == 0) {
				cout << "In the fitting function, cannot find the valid nearest point!" << endl;
				returnValue = false;
			}*/
		}

		/*Calculate the parameters of the objective function*/
		views[t[0]].RcCalculateB(num, nearestPointList1, para_B1, para_b1);
		views[t[1]].RcCalculateB(num, nearestPointList2, para_B2, para_b2);

		/*MatrixXf reguInb(3 * reconstructPointList.size(), 1);
		for (int i = 0; i < reconstructPointList.size(); i++) {
		reguInb(i * 3) = reconstructPointList[i][0];
		reguInb(i * 3 + 1) = reconstructPointList[i][1];
		reguInb(i * 3 + 2) = reconstructPointList[i][2];
		}*/

		MatrixXf paraBAll(7 * (num - 2) + 6, 3 * num);
		VectorXf parabAll(7 * (num - 2) + 6);

		paraBAll.block(0, 0, 2 * (num - 2), 3 * num) = para_B1;
		paraBAll.block(2 * (num - 2), 0, 2 * (num - 2), 3 * num) = para_B2;

		parabAll.block(0, 0, 2 * (num - 2), 1) = para_b1;
		parabAll.block(2 * (num - 2), 0, 2 * (num - 2), 1) = para_b2;


		/*For Regularization*/

		MatrixXf reguB = MatrixXf::Zero(3 * (num - 2), 3 * num);
		for (int i = 0; i < 3 * (num - 2); i++) {
			reguB(i, i) = lamda * 1; reguB(i, i + 3) = lamda*(-2); reguB(i, i + 6) = lamda * 1;
		}
		VectorXf regub = VectorXf::Zero(3 * (num - 2));

		paraBAll.block(4 * (num - 2), 0, 3 * (num - 2), 3 * num) = reguB;
		parabAll.block(4 * (num - 2), 0, 3 * (num - 2), 1) = regub;

		/*For X1 and Xn*/
		MatrixXf X1B = MatrixXf::Zero(3, 3 * num);
		X1B.block(0, 0, 3, 3) = Matrix3f::Identity();
		Vector3f X1b = samplePointOnCurve[0];

		MatrixXf XnB = MatrixXf::Zero(3, 3 * num);
		XnB.block(0, 3 * (num - 1), 3, 3) = Matrix3f::Identity();
		Vector3f Xnb = samplePointOnCurve[num - 1];

		paraBAll.block(7 * (num - 2), 0, 3, 3 * num) = X1B;
		paraBAll.block(7 * (num - 2) + 3, 0, 3, 3 * num) = XnB;

		parabAll.block(7 * (num - 2), 0, 3, 1) = X1b;
		parabAll.block(7 * (num - 2) + 3, 0, 3, 1) = Xnb;

		VectorXf before(3 * num);
		for (int i = 0; i < num; i++) {
			before[i * 3] = samplePointOnCurve[i][0];
			before[i * 3 + 1] = samplePointOnCurve[i][1];
			before[i * 3 + 2] = samplePointOnCurve[i][2];
		}
		VectorXf errorvec = (paraBAll*before - parabAll) / 1000.0;
		float error = errorvec.norm();
		//cout << "The error before optimization:\t" << error << endl;

		if (error > errorMax) {
			returnValue = false;
			break;
		}

		/*Minimize the objective function, update the location of reconstructed points*/
		VectorXf controlPoints = MatchFinder::SolveForCP(paraBAll, parabAll, samplePointOnCurve);

		VectorXf after(3 * num);
		for (int i = 0; i < num; i++) {
			after[i * 3] = samplePointOnCurve[i][0];
			after[i * 3 + 1] = samplePointOnCurve[i][1];
			after[i * 3 + 2] = samplePointOnCurve[i][2];
		}

		errorvec = (paraBAll*controlPoints - parabAll) / 1000.0;
		error = errorvec.norm();
		//cout << "The error after optimization:\t" << error << endl;
		iternn++;
	}

	string aftername = path + "\\after\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, aftername, color);
	return returnValue;
}

bool ReconCurve::FittingEpiArea32(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, Vector3f startPoint, Vector3f endPoint, int tangentIdx, Vector3f para, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color, float errorMax) {
	bool returnValue = true;
	float eps = 10e-6;
	int numNearest = 100;
	float step = 1.0 / num;
	Vector2f tangentPoint = views[t[1]].tangentPoints[tangentIdx];
	int sign = tangentPoint[0] * para[0] + tangentPoint[1] * para[1] + para[2];

	for (int i = 0; i < num - 1; i++) {
		Vector3f tempP = (1 - step*i)*startPoint + (step*i)*endPoint;
		samplePointOnCurve.push_back(tempP);
	}
	samplePointOnCurve.push_back(endPoint);

	string beforename = path + "\\before\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	//ReconCurve::SaveColorPoints(samplePointOnCurve, beforename, color);
	int iternn = 0;

	while (iternn < fittingIterNum) {
		vector<Vector2f> nearestPointList1(num), nearestPointList2(num);
		//VectorXf controlPoints(3 * samplePointOnCurve.size());
		MatrixXf para_B1(2 * (num - 2), 3 * num);
		VectorXf para_b1(2 * (num - 2));
		MatrixXf para_B2(2 * (num - 2), 3 * num);
		VectorXf para_b2(2 * (num - 2));

		vector<Vector2f> sampleOnV1(num), sampleOnV2(num);

		for (int p = 0; p < samplePointOnCurve.size(); p++) {
			Vector4f tempPoint3D; Vector3f imgPoint;
			tempPoint3D.head(3) = samplePointOnCurve[p]; tempPoint3D[3] = 1.0;
			imgPoint = views[t[0]].projectionMatrix*tempPoint3D;
			imgPoint = imgPoint / imgPoint[2];
			sampleOnV1[p] = imgPoint.head(2);

			imgPoint = views[t[1]].projectionMatrix*tempPoint3D;
			imgPoint = imgPoint / imgPoint[2];
			sampleOnV2[p] = imgPoint.head(2);
		}

		/*Find the closest points of the sample points on the image*/
		for (int i = 1; i < num - 1; i++) {
			vector<int> idx1(numNearest); vector<Vector2f> tempNearestPoints1(numNearest);
			views[t[0]].FindN_NearestPointOnImgEachPoint(samplePointOnCurve[i], tempNearestPoints1, idx1, numNearest);

			float min_arclength1; int min_arcidx1;
			for (int j = 0; j < numNearest; j++) {
				float min_dist; int min_idx;
				for (int k = 0; k < num; k++) {
					float dist = (tempNearestPoints1[j] - sampleOnV1[k]).norm();
					if (k == 0 || dist < min_dist) {
						min_dist = dist;
						min_idx = k;
					}
				}
				float arclength = 0.0;
				if (min_idx < i) {
					for (int p = min_idx; p < i; p++) {
						arclength += (sampleOnV1[p] - sampleOnV1[p + 1]).norm();
					}
				}
				else {
					for (int p = i; p < min_idx; p++) {
						arclength += (sampleOnV1[p] - sampleOnV1[p + 1]).norm();
					}
				}
				if (j == 0 || arclength < min_arclength1) {
					min_arclength1 = arclength;
					min_arcidx1 = j;
				}
			}
			nearestPointList1[i] = tempNearestPoints1[min_arcidx1];

			vector<int> idx2(numNearest); vector<Vector2f> tempNearestPoints2(numNearest);
			views[t[1]].FindN_NearestPointOnImgEachPoint(samplePointOnCurve[i], tempNearestPoints2, idx2, numNearest);

			float min_arclength2; int min_arcidx2;
			for (int j = 0; j < numNearest; j++) {
				float min_dist; int min_idx;
				for (int k = 0; k < num; k++) {
					float dist = (tempNearestPoints2[j] - sampleOnV2[k]).norm();
					if (k == 0 || dist < min_dist) {
						min_dist = dist;
						min_idx = k;
					}
				}
				float arclength = 0.0;
				if (min_idx < i) {
					for (int p = min_idx; p < i; p++) {
						arclength += (sampleOnV2[p] - sampleOnV2[p + 1]).norm();
					}
				}
				else {
					for (int p = i; p < min_idx; p++) {
						arclength += (sampleOnV2[p] - sampleOnV2[p + 1]).norm();
					}
				}
				if (j == 0 || arclength < min_arclength2) {
					min_arclength2 = arclength;
					min_arcidx2 = j;
				}
			}
			nearestPointList2[i] = tempNearestPoints2[min_arcidx2];
		}

		/*Calculate the parameters of the objective function*/
		views[t[0]].RcCalculateB(num, nearestPointList1, para_B1, para_b1);
		views[t[1]].RcCalculateB(num, nearestPointList2, para_B2, para_b2);

		/*MatrixXf reguInb(3 * reconstructPointList.size(), 1);
		for (int i = 0; i < reconstructPointList.size(); i++) {
		reguInb(i * 3) = reconstructPointList[i][0];
		reguInb(i * 3 + 1) = reconstructPointList[i][1];
		reguInb(i * 3 + 2) = reconstructPointList[i][2];
		}*/

		MatrixXf paraBAll(7 * (num - 2) + 6, 3 * num);
		VectorXf parabAll(7 * (num - 2) + 6);

		paraBAll.block(0, 0, 2 * (num - 2), 3 * num) = para_B1;
		paraBAll.block(2 * (num - 2), 0, 2 * (num - 2), 3 * num) = para_B2;

		parabAll.block(0, 0, 2 * (num - 2), 1) = para_b1;
		parabAll.block(2 * (num - 2), 0, 2 * (num - 2), 1) = para_b2;


		/*For Regularization*/

		MatrixXf reguB = MatrixXf::Zero(3 * (num - 2), 3 * num);
		for (int i = 0; i < 3 * (num - 2); i++) {
			reguB(i, i) = lamda * 1; reguB(i, i + 3) = lamda*(-2); reguB(i, i + 6) = lamda * 1;
		}
		VectorXf regub = VectorXf::Zero(3 * (num - 2));

		paraBAll.block(4 * (num - 2), 0, 3 * (num - 2), 3 * num) = reguB;
		parabAll.block(4 * (num - 2), 0, 3 * (num - 2), 1) = regub;

		/*For X1 and Xn*/
		MatrixXf X1B = MatrixXf::Zero(3, 3 * num);
		X1B.block(0, 0, 3, 3) = Matrix3f::Identity();
		Vector3f X1b = samplePointOnCurve[0];

		MatrixXf XnB = MatrixXf::Zero(3, 3 * num);
		XnB.block(0, 3 * (num - 1), 3, 3) = Matrix3f::Identity();
		Vector3f Xnb = samplePointOnCurve[num - 1];

		paraBAll.block(7 * (num - 2), 0, 3, 3 * num) = X1B;
		paraBAll.block(7 * (num - 2) + 3, 0, 3, 3 * num) = XnB;

		parabAll.block(7 * (num - 2), 0, 3, 1) = X1b;
		parabAll.block(7 * (num - 2) + 3, 0, 3, 1) = Xnb;

		VectorXf before(3 * num);
		for (int i = 0; i < num; i++) {
			before[i * 3] = samplePointOnCurve[i][0];
			before[i * 3 + 1] = samplePointOnCurve[i][1];
			before[i * 3 + 2] = samplePointOnCurve[i][2];
		}
		VectorXf errorvec = (paraBAll*before - parabAll) / 1000.0;
		float error = errorvec.norm();
		//cout << "The error before optimization:\t" << error << endl;

		if (error > errorMax) {
			returnValue = false;
			break;
		}

		/*Minimize the objective function, update the location of reconstructed points*/
		VectorXf controlPoints = MatchFinder::SolveForCP(paraBAll, parabAll, samplePointOnCurve);

		VectorXf after(3 * num);
		for (int i = 0; i < num; i++) {
			after[i * 3] = samplePointOnCurve[i][0];
			after[i * 3 + 1] = samplePointOnCurve[i][1];
			after[i * 3 + 2] = samplePointOnCurve[i][2];
		}

		errorvec = (paraBAll*controlPoints - parabAll) / 1000.0;
		error = errorvec.norm();
		//cout << "The error after optimization:\t" << error << endl;
		iternn++;
	}

	string aftername = path + "\\after\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	//ReconCurve::SaveColorPoints(samplePointOnCurve, aftername, color);
	return returnValue;
}

/* Wrong */
bool ReconCurve::FittingEpiArea5(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, int startV1Idx, int endV1Idx, Vector3f startPoint, Vector3f endPoint, Vector3f para, int epiAreaIdx, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color) {
	float errorMax = 2000.0;
	bool returnValue = true;
	vector<Vector2f> sampleV1Points, sampleV2Points;

	float eps = 10e-6;
	int numNearest = 100;
	float step = 1.0 / num;
	/* Sample points on the line segment and calculate the backprojection on the images */
	for (int i = 0; i < num - 1; i++) {
		Vector3f tempP = (1 - step*i)*startPoint + (step*i)*endPoint;
		samplePointOnCurve.push_back(tempP);

		Vector4f tempPoint3D = Vector4f::Ones(); tempPoint3D.head(3) = tempP;
		Vector3f V1Point = views[t[0]].projectionMatrix*tempPoint3D;
		V1Point = V1Point / V1Point[2];
		sampleV1Points.push_back(V1Point.head(2));

		Vector3f V2Point = views[t[1]].projectionMatrix*tempPoint3D;
		V2Point = V2Point / V2Point[2];
		sampleV2Points.push_back(V2Point.head(2));
	}
	samplePointOnCurve.push_back(endPoint);
	Vector4f tempEnd3D = Vector4f::Ones(); tempEnd3D.head(3) = endPoint;
	Vector3f V1EndPoint = views[t[0]].projectionMatrix*tempEnd3D;
	V1EndPoint = V1EndPoint / V1EndPoint[2];
	sampleV1Points.push_back(V1EndPoint.head(2));

	Vector3f V2EndPoint = views[t[1]].projectionMatrix*tempEnd3D;
	V2EndPoint = V2EndPoint / V2EndPoint[2];
	sampleV2Points.push_back(V2EndPoint.head(2));

	/*vector<Vector2f> V1ImgMapPoints;
	for (int i = 0; i < endV1Idx - startV1Idx + 1; i++) {
		V1ImgMapPoints.push_back(views[t[0]].img2DPoints[i]);
	}*/



	string beforename = path + "\\before\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, beforename, color);
	int iternn = 0;

	while (iternn < fittingIterNum) {
		vector<Vector2f> nearestPointList1(num), nearestPointList2(num);
		//VectorXf controlPoints(3 * samplePointOnCurve.size());
		MatrixXf para_B1(2 * (num - 2), 3 * num);
		VectorXf para_b1(2 * (num - 2));
		MatrixXf para_B2(2 * (num - 2), 3 * num);
		VectorXf para_b2(2 * (num - 2));

		/*Find the closest points of the sample points on the image*/
		for (int i = 1; i < num - 1; i++) {
			int idx1; float minDist1;
			for (int j = startV1Idx; j < endV1Idx + 1; j++) {
				float dist = (views[t[0]].img2DPoints[j] - sampleV1Points[i]).norm();
				if (dist > minDist1 || j == startV1Idx) {
					minDist1 = dist;
					idx1 = j;
				}
			}
			nearestPointList1[i] = views[t[0]].img2DPoints[idx1];


			vector<int> idx2(numNearest); vector<Vector2f> tempNearestPoints2(numNearest);
			int sign = views[t[1]].tangentSign[epiAreaIdx];
			views[t[1]].FindN_NearestPointOnImgEachPoint(samplePointOnCurve[i], tempNearestPoints2, idx2, numNearest);
			for (int j = 0; j < numNearest; j++) {
				float tempSign = tempNearestPoints2[j][0] * para[0] + tempNearestPoints2[j][1] * para[1] + para[2];
				if (tempSign*sign < -eps) {
					nearestPointList2[i] = tempNearestPoints2[j]; break;
				}
			}
			/*int idx2; float minDist2;
			for (int j = views[t[1]].tangentLeftRight[epiAreaIdx][0]; j < views[t[1]].tangentLeftRight[epiAreaIdx][1]+1; j++) {
				float dist = (views[t[1]].img2DPoints[j] - sampleV2Points[i]).norm();
				if (dist > minDist2 || j == views[t[1]].tangentLeftRight[epiAreaIdx][0]) {
					minDist2 = dist;
					idx2 = j;
				}
			}
			nearestPointList2[i] = views[t[1]].img2DPoints[idx2];*/
		}

		/*Calculate the parameters of the objective function*/
		views[t[0]].RcCalculateB(num, nearestPointList1, para_B1, para_b1);
		views[t[1]].RcCalculateB(num, nearestPointList2, para_B2, para_b2);

		/*MatrixXf reguInb(3 * reconstructPointList.size(), 1);
		for (int i = 0; i < reconstructPointList.size(); i++) {
		reguInb(i * 3) = reconstructPointList[i][0];
		reguInb(i * 3 + 1) = reconstructPointList[i][1];
		reguInb(i * 3 + 2) = reconstructPointList[i][2];
		}*/

		MatrixXf paraBAll(7 * (num - 2) + 6, 3 * num);
		VectorXf parabAll(7 * (num - 2) + 6);

		paraBAll.block(0, 0, 2 * (num - 2), 3 * num) = para_B1;
		paraBAll.block(2 * (num - 2), 0, 2 * (num - 2), 3 * num) = para_B2;

		parabAll.block(0, 0, 2 * (num - 2), 1) = para_b1;
		parabAll.block(2 * (num - 2), 0, 2 * (num - 2), 1) = para_b2;


		/*For Regularization*/

		MatrixXf reguB = MatrixXf::Zero(3 * (num - 2), 3 * num);
		for (int i = 0; i < 3 * (num - 2); i++) {
			reguB(i, i) = lamda * 1; reguB(i, i + 3) = lamda*(-2); reguB(i, i + 6) = lamda * 1;
		}
		VectorXf regub = VectorXf::Zero(3 * (num - 2));

		paraBAll.block(4 * (num - 2), 0, 3 * (num - 2), 3 * num) = reguB;
		parabAll.block(4 * (num - 2), 0, 3 * (num - 2), 1) = regub;

		/*For X1 and Xn*/
		MatrixXf X1B = MatrixXf::Zero(3, 3 * num);
		X1B.block(0, 0, 3, 3) = Matrix3f::Identity();
		Vector3f X1b = samplePointOnCurve[0];

		MatrixXf XnB = MatrixXf::Zero(3, 3 * num);
		XnB.block(0, 3 * (num - 1), 3, 3) = Matrix3f::Identity();
		Vector3f Xnb = samplePointOnCurve[num - 1];

		paraBAll.block(7 * (num - 2), 0, 3, 3 * num) = X1B;
		paraBAll.block(7 * (num - 2) + 3, 0, 3, 3 * num) = XnB;

		parabAll.block(7 * (num - 2), 0, 3, 1) = X1b;
		parabAll.block(7 * (num - 2) + 3, 0, 3, 1) = Xnb;

		VectorXf before(3 * num);
		for (int i = 0; i < num; i++) {
			before[i * 3] = samplePointOnCurve[i][0];
			before[i * 3 + 1] = samplePointOnCurve[i][1];
			before[i * 3 + 2] = samplePointOnCurve[i][2];
		}
		VectorXf errorvec = paraBAll*before - parabAll;
		float error = errorvec.norm();
		//cout << "The error before optimization:\t" << error << endl;

		if (error > errorMax) {
			returnValue = false;
			break;
		}
		/*Minimize the objective function, update the location of reconstructed points*/
		VectorXf controlPoints = MatchFinder::SolveForCP(paraBAll, parabAll, samplePointOnCurve);

		VectorXf after(3 * num);
		for (int i = 0; i < num; i++) {
			after[i * 3] = samplePointOnCurve[i][0];
			after[i * 3 + 1] = samplePointOnCurve[i][1];
			after[i * 3 + 2] = samplePointOnCurve[i][2];
		}

		errorvec = paraBAll*controlPoints - parabAll;
		error = errorvec.norm();
		//cout << "The error after optimization:\t" << error << endl;
		iternn++;
	}

	string aftername = path + "\\after\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, aftername, color);
	return returnValue;
}

/* Wrong */
bool ReconCurve::FittingEpiArea6(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, int startV1Idx, int endV1Idx, Vector3f startPoint, Vector3f endPoint, Vector3f para, int epiAreaIdx, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color) {
	float errorMax = 2000.0;
	bool returnValue = true;
	vector<Vector2f> sampleV1Points, sampleV2Points;

	float eps = 10e-6;
	int numNearest = 50;
	float step = 1.0 / num;
	/* Sample points on the line segment and calculate the backprojection on the images */
	for (int i = 0; i < num - 1; i++) {
		Vector3f tempP = (1 - step*i)*startPoint + (step*i)*endPoint;
		samplePointOnCurve.push_back(tempP);

		Vector4f tempPoint3D = Vector4f::Ones(); tempPoint3D.head(3) = tempP;
		Vector3f V1Point = views[t[0]].projectionMatrix*tempPoint3D;
		V1Point = V1Point / V1Point[2];
		sampleV1Points.push_back(V1Point.head(2));

		Vector3f V2Point = views[t[1]].projectionMatrix*tempPoint3D;
		V2Point = V2Point / V2Point[2];
		sampleV2Points.push_back(V2Point.head(2));
	}
	samplePointOnCurve.push_back(endPoint);
	Vector4f tempEnd3D = Vector4f::Ones(); tempEnd3D.head(3) = endPoint;
	Vector3f V1EndPoint = views[t[0]].projectionMatrix*tempEnd3D;
	V1EndPoint = V1EndPoint / V1EndPoint[2];
	sampleV1Points.push_back(V1EndPoint.head(2));

	Vector3f V2EndPoint = views[t[1]].projectionMatrix*tempEnd3D;
	V2EndPoint = V2EndPoint / V2EndPoint[2];
	sampleV2Points.push_back(V2EndPoint.head(2));

	/*vector<Vector2f> V1ImgMapPoints;
	for (int i = 0; i < endV1Idx - startV1Idx + 1; i++) {
	V1ImgMapPoints.push_back(views[t[0]].img2DPoints[i]);
	}*/



	string beforename = path + "\\before\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, beforename, color);
	int iternn = 0;

	while (iternn < fittingIterNum) {
		vector<Vector2f> nearestPointList1(num), nearestPointList2(num);
		//VectorXf controlPoints(3 * samplePointOnCurve.size());
		MatrixXf para_B1(2 * (num - 2), 3 * num);
		VectorXf para_b1(2 * (num - 2));
		MatrixXf para_B2(2 * (num - 2), 3 * num);
		VectorXf para_b2(2 * (num - 2));

		/*Find the closest points of the sample points on the image*/
		for (int i = 1; i < num - 1; i++) {
			int idx1; float minDist1;
			for (int j = startV1Idx; j < endV1Idx + 1; j++) {
				float dist = (views[t[0]].img2DPoints[j] - sampleV1Points[i]).norm();
				if (dist > minDist1 || j == startV1Idx) {
					minDist1 = dist;
					idx1 = j;
				}
			}
			nearestPointList1[i] = views[t[0]].img2DPoints[idx1];



			vector<int> idx2(numNearest); vector<Vector2f> tempNearestPoints2(numNearest);
			int sign = views[t[1]].tangentSign[epiAreaIdx];
			views[t[1]].FindN_NearestPointOnImgEachPoint(samplePointOnCurve[i], tempNearestPoints2, idx2, numNearest);
			for (int j = 0; j < numNearest; j++) {
				int minBackIdx;  float minBackDist;
				for (int k = 0; k < sampleV2Points.size(); k++) {
					float backDist = (sampleV2Points[k] - tempNearestPoints2[j]).norm();
					if (backDist < minBackDist || k == 0) {
						minBackIdx = k;
						minBackDist = backDist;
					}
				}
				if (abs(i - minBackIdx) < 3) {
					nearestPointList2[i] = tempNearestPoints2[j]; break;
				}
			}

		}

		/*Calculate the parameters of the objective function*/
		views[t[0]].RcCalculateB(num, nearestPointList1, para_B1, para_b1);
		views[t[1]].RcCalculateB(num, nearestPointList2, para_B2, para_b2);

		/*MatrixXf reguInb(3 * reconstructPointList.size(), 1);
		for (int i = 0; i < reconstructPointList.size(); i++) {
		reguInb(i * 3) = reconstructPointList[i][0];
		reguInb(i * 3 + 1) = reconstructPointList[i][1];
		reguInb(i * 3 + 2) = reconstructPointList[i][2];
		}*/

		MatrixXf paraBAll(7 * (num - 2) + 6, 3 * num);
		VectorXf parabAll(7 * (num - 2) + 6);

		paraBAll.block(0, 0, 2 * (num - 2), 3 * num) = para_B1;
		paraBAll.block(2 * (num - 2), 0, 2 * (num - 2), 3 * num) = para_B2;

		parabAll.block(0, 0, 2 * (num - 2), 1) = para_b1;
		parabAll.block(2 * (num - 2), 0, 2 * (num - 2), 1) = para_b2;


		/*For Regularization*/

		MatrixXf reguB = MatrixXf::Zero(3 * (num - 2), 3 * num);
		for (int i = 0; i < 3 * (num - 2); i++) {
			reguB(i, i) = lamda * 1; reguB(i, i + 3) = lamda*(-2); reguB(i, i + 6) = lamda * 1;
		}
		VectorXf regub = VectorXf::Zero(3 * (num - 2));

		paraBAll.block(4 * (num - 2), 0, 3 * (num - 2), 3 * num) = reguB;
		parabAll.block(4 * (num - 2), 0, 3 * (num - 2), 1) = regub;

		/*For X1 and Xn*/
		MatrixXf X1B = MatrixXf::Zero(3, 3 * num);
		X1B.block(0, 0, 3, 3) = Matrix3f::Identity();
		Vector3f X1b = samplePointOnCurve[0];

		MatrixXf XnB = MatrixXf::Zero(3, 3 * num);
		XnB.block(0, 3 * (num - 1), 3, 3) = Matrix3f::Identity();
		Vector3f Xnb = samplePointOnCurve[num - 1];

		paraBAll.block(7 * (num - 2), 0, 3, 3 * num) = X1B;
		paraBAll.block(7 * (num - 2) + 3, 0, 3, 3 * num) = XnB;

		parabAll.block(7 * (num - 2), 0, 3, 1) = X1b;
		parabAll.block(7 * (num - 2) + 3, 0, 3, 1) = Xnb;

		VectorXf before(3 * num);
		for (int i = 0; i < num; i++) {
			before[i * 3] = samplePointOnCurve[i][0];
			before[i * 3 + 1] = samplePointOnCurve[i][1];
			before[i * 3 + 2] = samplePointOnCurve[i][2];
		}
		VectorXf errorvec = paraBAll*before - parabAll;
		float error = errorvec.norm();
		//cout << "The error before optimization:\t" << error << endl;

		if (error > errorMax) {
			returnValue = false;
			break;
		}
		/*Minimize the objective function, update the location of reconstructed points*/
		VectorXf controlPoints = MatchFinder::SolveForCP(paraBAll, parabAll, samplePointOnCurve);

		VectorXf after(3 * num);
		for (int i = 0; i < num; i++) {
			after[i * 3] = samplePointOnCurve[i][0];
			after[i * 3 + 1] = samplePointOnCurve[i][1];
			after[i * 3 + 2] = samplePointOnCurve[i][2];
		}

		errorvec = paraBAll*controlPoints - parabAll;
		error = errorvec.norm();
		//cout << "The error after optimization:\t" << error << endl;
		iternn++;
	}

	string aftername = path + "\\after\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, aftername, color);
	return returnValue;
}

bool ReconCurve::FittingEpiArea7(vector<Vector3f>& samplePointOnCurve, string path, int iter, int count, int which, int nn, int startV1Idx, int endV1Idx, Vector3f startPoint, Vector3f endPoint, Vector3f para, int epiAreaIdx, vector<ProjectImage> views, vector<int> t, int num, float lamda, int fittingIterNum, Vector3f color) {
	bool returnValue = true;
	int numNearest = 100;
	int sign = views[t[1]].tangentSign[epiAreaIdx];
	float eps = 10e-6;
	float step = 1.0 / num;

	vector<Vector3f> beforeSamplePoint;
	for (int i = 0; i < num - 1; i++) {
		Vector3f tempP = (1 - step*i)*startPoint + (step*i)*endPoint;
		beforeSamplePoint.push_back(tempP);
	}
	beforeSamplePoint.push_back(endPoint);

	string beforename = path + "\\before\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(beforeSamplePoint, beforename, color);

	Vector3f midPoint = (startPoint + endPoint) / 2.0;

	Vector4f tempPoint3D = Vector4f::Ones(); tempPoint3D.head(3) = midPoint;
	Vector3f midPointProjection = views[t[0]].projectionMatrix*tempPoint3D;
	midPointProjection = midPointProjection / midPointProjection[2];

	Vector2f nearestPointV1, nearestPointV2;
	int idx1; float minDist1;
	for (int j = startV1Idx; j < endV1Idx + 1; j++) {
		float dist = (views[t[0]].img2DPoints[j] - midPointProjection.head(2)).norm();
		if (dist > minDist1 || j == startV1Idx) {
			minDist1 = dist;
			idx1 = j;
		}
	}
	nearestPointV1 = views[t[0]].img2DPoints[idx1];

	vector<int> idx2(numNearest); vector<Vector2f> tempNearestPoints2(numNearest);
	views[t[1]].FindN_NearestPointOnImgEachPoint(midPoint, tempNearestPoints2, idx2, numNearest);
	for (int j = 0; j < numNearest; j++) {
		float tempSign = tempNearestPoints2[j][0] * para[0] + tempNearestPoints2[j][1] * para[1] + para[2];
		if (tempSign*sign < -eps) {
			nearestPointV2 = tempNearestPoints2[j]; break;
		}
	}

	Vector3f newMidPoint = MatchFinder::computeClosest3DPoint(nearestPointV1, nearestPointV2, views[t[0]].projectionMatrix, views[t[1]].projectionMatrix);

	vector<Vector3f> midPointIter2(2);
	vector<Vector2f> nearestPointV1Iter2(2), nearestPointV2Iter2(2);
	midPointIter2[0] = (startPoint + newMidPoint) / 2.0;
	midPointIter2[1] = (newMidPoint + endPoint) / 2.0;

	vector<Vector3f> newMidPointIter2(2);
	for (int i = 0; i < 2; i++) {
		Vector4f temp3D = Vector4f::Ones(); temp3D.head(3) = midPointIter2[i];
		Vector3f midPointIter2Projection = views[t[0]].projectionMatrix*temp3D;
		midPointIter2Projection = midPointIter2Projection / midPointIter2Projection[2];
		int idx1; float minDist1;
		for (int j = startV1Idx; j < endV1Idx + 1; j++) {
			float dist = (views[t[0]].img2DPoints[j] - midPointIter2Projection.head(2)).norm();
			if (dist > minDist1 || j == startV1Idx) {
				minDist1 = dist;
				idx1 = j;
			}
		}
		nearestPointV1Iter2[i] = views[t[0]].img2DPoints[idx1];

		vector<int> idx2(numNearest); vector<Vector2f> tempNearestPoints2(numNearest);
		views[t[1]].FindN_NearestPointOnImgEachPoint(midPointIter2[i], tempNearestPoints2, idx2, numNearest);
		for (int j = 0; j < numNearest; j++) {
			float tempSign = tempNearestPoints2[j][0] * para[0] + tempNearestPoints2[j][1] * para[1] + para[2];
			if (tempSign*sign < -eps) {
				nearestPointV2Iter2[i] = tempNearestPoints2[j]; break;
			}
		}
		newMidPointIter2[i] = MatchFinder::computeClosest3DPoint(nearestPointV1Iter2[i], nearestPointV2Iter2[i], views[t[0]].projectionMatrix, views[t[1]].projectionMatrix);
	}

	step = step * 4;
	for (int i = 0; i < num / 4; i++) {
		Vector3f tempP = (1 - step*i)*startPoint + (step*i)*newMidPointIter2[0];
		samplePointOnCurve.push_back(tempP);
	}
	for (int i = 0; i < num / 4; i++) {
		Vector3f tempP = (1 - step*i)*newMidPointIter2[0] + (step*i)*newMidPoint;
		samplePointOnCurve.push_back(tempP);
	}
	for (int i = 0; i < num / 4; i++) {
		Vector3f tempP = (1 - step*i)*newMidPoint + (step*i)*newMidPointIter2[1];
		samplePointOnCurve.push_back(tempP);
	}
	for (int i = 0; i < num / 4; i++) {
		Vector3f tempP = (1 - step*i)*newMidPointIter2[1] + (step*i)*endPoint;
		samplePointOnCurve.push_back(tempP);
	}
	samplePointOnCurve.push_back(endPoint);

	string aftername = path + "\\after\\iter" + std::to_string(iter) + "_part" + std::to_string(count) + "_sub" + std::to_string(which) + "_t" + std::to_string(nn) + ".off";
	ReconCurve::SaveColorPoints(samplePointOnCurve, aftername, color);
	return returnValue;
}