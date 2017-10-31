#pragma once
#include "LibParty.h"
#include "ReconCurve.h"

class CurveOperation {
public:
	static vector<vector<int>> CombinationChoose(int comb) {
		vector<vector<int>> T;
		switch(comb) {
		case 1:
			 T = { { 2,1,3 },{ 3,1,2 },{ 3,2,1 },{ 2,3,1 },{ 1,2,3 },{ 1,3,2 } };//1
			 break;
		case 2:
			 T = { { 3,1,2 },{ 3,2,1 },{ 2,3,1 },{ 1,2,3 },{ 1,3,2 },{ 2,1,3 } };//2
			 break;
		case 3:
			 T = { { 3,2,1 },{ 2,3,1 },{ 1,2,3 },{ 1,3,2 },{ 2,1,3 },{ 3,1,2 } };//3
			 break;
		case 4:
			 T = { { 2,3,1 },{ 1,2,3 },{ 1,3,2 },{ 2,1,3 },{ 3,1,2 },{ 3,2,1 } };//4
			 break;
		case 5:
			 T = { { 1,2,3 },{ 1,3,2 },{ 2,1,3 },{ 3,1,2 },{ 3,2,1 },{ 2,3,1 } };//5
			 break;
		case 6:
			 T = { { 1,3,2 },{ 2,1,3 },{ 3,1,2 },{ 3,2,1 },{ 2,3,1 },{ 1,2,3 } };//6
			 break;
		};
		return T;
	}

	static vector<int> SeperateCurve(int x, int count, ProjectImage view) {
		int fixnum = view.intersectNum[0];
		vector<int> startPointIdx;
		startPointIdx.push_back(0);
		for (int i = 1; i < view.intersectNum.size(); i++) {
			if (view.intersectNum[i] != fixnum) { startPointIdx.push_back(i); fixnum = view.intersectNum[i]; }
		}
		startPointIdx.push_back(view.intersectNum.size());
		return startPointIdx;
	}

	static int GenerateTrail(vector<int> startPointIdx, ProjectImage view, vector<vector<vector<int>>> &totalPairString) {
		int numCurve = 0;
		/* Calculate the distance between the intersection points of two neighboring points in view1 */
		for (int t = 0; t < startPointIdx.size() - 1; t++) {
			/* Record the string of intersection points on a subcurve */
			int stringSpace = startPointIdx[t + 1] - startPointIdx[t];
			if (stringSpace > 0 && view.intersectNum[startPointIdx[t]] != 0) {
				vector<vector<int>> pairString; vector<int> tempString(stringSpace, -1);
				for (int i = 0; i < view.intersectNum[startPointIdx[t]]; i++) {
					tempString[0] = view.intersectIdxStart[startPointIdx[t]] + i;
					pairString.push_back(tempString);
				}
				for (int i = startPointIdx[t]; i < startPointIdx[t + 1] - 1; i++) {
					/* For the intersection points of two neighbor points i and i+1 */
					vector<int> curtag(view.intersectNum[i + 1], 0);
					for (int p = 0; p < pairString.size(); p++) {
						if (pairString[p][i - startPointIdx[t]] != -1) {
							Vector2f prepoint = view.intersect[pairString[p][i - startPointIdx[t]]];
							Vector2f presumdire = Vector2f::Zero();

							int whichq = -1;
							float mindist = 5.0;
							for (int q = 0; q < view.intersectNum[i + 1]; q++) {
								Vector2f curpoint = view.intersect[view.intersectIdxStart[i + 1] + q];
								if (curtag[q] == 0) {
									float curdist = (curpoint - prepoint).norm();
									if (curdist < mindist) {
										whichq = q;
										mindist = curdist;
									}
								}
							}
							if (whichq != -1) {
								pairString[p][i + 1 - startPointIdx[t]] = view.intersectIdxStart[i + 1] + whichq;
								curtag[whichq] = 1;
							}
						}
					}
				}
				totalPairString[t] = pairString;
				numCurve += view.intersectNum[startPointIdx[t]];
			}
		}
		return numCurve;
	}

	//static int GenerateTrail_inEpiArea(vector<int> startPointIdx, ProjectImage view, vector<vector<vector<int>>> &totalPairString) {
	//	int numCurve = 0;
	//	/* Calculate the distance between the intersection points of two neighboring points in view1 */
	//	for (int t = 0; t < startPointIdx.size() - 1; t++) {
	//		/* Record the string of intersection points on a subcurve */
	//		int stringSpace = startPointIdx[t + 1] - startPointIdx[t];
	//		if (stringSpace > 0 && view.intersectNum[startPointIdx[t]] != 0) {
	//			vector<vector<int>> pairString; vector<int> tempString(stringSpace, -1);
	//			for (int i = 0; i < view.intersectNum[startPointIdx[t]]; i++) {
	//				tempString[0] = view.intersectIdxStart[startPointIdx[t]] + i;
	//				pairString.push_back(tempString);
	//			}
	//			for (int i = startPointIdx[t]; i < startPointIdx[t + 1] - 1; i++) {
	//				/* For the intersection points of two neighbor points i and i+1 */
	//				vector<int> curtag(view.intersectNum[i + 1], 0);
	//				for (int p = 0; p < pairString.size(); p++) {
	//					if (pairString[p][i - startPointIdx[t]] != -1) {
	//						Vector2f prepoint = view.intersect[pairString[p][i - startPointIdx[t]]];
	//						Vector2f presumdire = Vector2f::Zero();

	//						int whichq = -1;
	//						float mindist = 1000.0; float maxdot = -1;
	//						for (int q = 0; q < view.intersectNum[i + 1]; q++) {
	//							Vector2f curpoint = view.intersect[view.intersectIdxStart[i + 1] + q];
	//							if (curtag[q] == 0) {
	//								float curdist = (curpoint - prepoint).norm();
	//								if (curdist < mindist) {
	//									whichq = q;
	//									mindist = curdist;
	//								}
	//							}
	//						}
	//						if (whichq != -1) {
	//							pairString[p][i + 1 - startPointIdx[t]] = view.intersectIdxStart[i + 1] + whichq;
	//							curtag[whichq] = 1;
	//						}
	//					}
	//				}
	//			}
	//			totalPairString[t] = pairString;
	//			numCurve += view.intersectNum[startPointIdx[t]];
	//		}
	//	}
	//	return numCurve;
	//}

	static int GenerateLongTrail(vector<int> startPointIdx, vector<ProjectImage> &views, ProjectImage &refView, vector<vector<Vector3f>> reconList, float onecurve_dist, vector<vector<Vector3f>> &newReconList) {
		for (int i = 0; i < views[0].intersectNum[0]; i++) {
			newReconList.push_back(reconList[i]);
		}

		for (int i = 0; i < startPointIdx.size() - 1; i++) {
			vector<vector<Vector3f>> tempNewReconList;
			for (int j = 0; j < newReconList.size(); j++) {
				Vector3f firstEnd = newReconList[j].back();
				Vector4f temp3D = Vector4f::Ones();
				temp3D.head(3) = firstEnd;
				for (int k = 0; k < views[0].intersectNum[startPointIdx[i + 1]]; k++) {
					vector<Vector3f> secondTrail = reconList[views[0].intersectIdxStart[i + 1] + k];
					Vector3f secondStart = secondTrail.front();
					Vector4f temp3D2 = Vector4f::Ones();
					temp3D2.head(3) = secondStart;
					float max_dist;
					for (int x = 0; x < views.size(); x++) {
						Vector3f temp2D = views[x].projectionMatrix*temp3D;
						temp2D = temp2D / temp2D[2];
						Vector2f firstEndPixel = temp2D.head(2);
						Vector3f temp2D2 = views[x].projectionMatrix*temp3D2;
						temp2D2 = temp2D2 / temp2D2[2];
						Vector2f secondStartPixel = temp2D2.head(2);
						float dist = (firstEndPixel - secondStartPixel).norm();
						if (x == 0) max_dist = dist;
						else if (dist > max_dist) max_dist = dist;
					}

					if (max_dist < onecurve_dist) {
						vector<Vector3f> tempNewRecon = newReconList[j];
						tempNewRecon.insert(tempNewRecon.end(), secondTrail.begin(), secondTrail.end());
						tempNewReconList.push_back(tempNewRecon);
					}
				}

			}
			//cout << "aa" << i << endl;
			//cout << tempNewReconList.size() << endl;
			newReconList.insert(newReconList.end(), tempNewReconList.begin(), tempNewReconList.end());
			//cout << newReconList.size() << endl;
			tempNewReconList.clear();
			//cout << tempNewReconList.size() << endl;
		}
		return newReconList.size();
	}

	static void SaveTrail(vector<int> startPointIdx, vector<vector<vector<int>>> totalPairString, vector<ProjectImage> views, vector<int> t, vector<ReconCurve> &reconCurveList, string path, int iter, int count, int offsize) {
		int start = views[t[0]].reconGroupIdxStart[count];
		vector<Vector3f> colors = { Vector3f(255, 0, 0), Vector3f(0, 255, 0), Vector3f(0, 0, 255), Vector3f(255, 255, 0), Vector3f(255, 0, 255), Vector3f(0, 255, 255) };
		for (int i = 0; i < startPointIdx.size() - 1; i++) {
			int stringSpace = startPointIdx[i + 1] - startPointIdx[i];
			Vector3f color = colors[i % 6];
			/*color[0] = rand() % 256; color[1] = rand() % 256; color[2] = rand() % 256;*/
			if (stringSpace > 0 && views[t[1]].intersectNum[startPointIdx[i]] != 0) {
				for (int j = 0; j < views[t[1]].intersectNum[startPointIdx[i]]; j++) {
					string full_count;
					if (count < 10) full_count = "0" + std::to_string(count);
					else full_count = std::to_string(count);

					string full_i;
					if (i < 10) full_i = "00" + std::to_string(i);
					else if (i < 100) full_i = "0" + std::to_string(i);
					else full_i = std::to_string(i);

					string filename = "iter" + std::to_string(iter) + "_part" + full_count + "_sub" + full_i + "_t" + std::to_string(j) + ".off";
					string recon = path + "/recon/" + filename;
					string recon_color = path + "/recon_color/" + filename;

					ReconCurve tempCurve;
					vector<Vector3f> reconstructPointList;
					for (int k = 0; k < stringSpace; k++) {
						if (totalPairString[i][j][k] != -1) {
							Vector3f temprecon = MatchFinder::computeClosest3DPoint(views[t[0]].img2DPoints[start + startPointIdx[i] + k], views[t[1]].intersect[totalPairString[i][j][k]], views[t[0]].projectionMatrix, views[t[1]].projectionMatrix);
							reconstructPointList.push_back(temprecon);
						}
					}
					int startRef = start + startPointIdx[i];
					int endRef = start + startPointIdx[i + 1] - 1;

					if (reconstructPointList.size() > offsize) {
						ReconCurve::SavePoints(reconstructPointList, recon);
						ReconCurve::SaveColorPoints(reconstructPointList, recon_color, color);
						tempCurve.curve = reconstructPointList;
						tempCurve.curveName = filename;
						tempCurve.iter = iter;
						tempCurve.partIdx = count;
						tempCurve.clusterIdx = i;
						tempCurve.idxInCluster = j;
						tempCurve.refIdx = t[0];
						tempCurve.startRef = startRef;
						tempCurve.endRef = endRef;
						tempCurve.color = color;
						reconCurveList.push_back(tempCurve);
					}
				}
			}
		}
	}

	static void SaveTrail2(vector<vector<Vector3f>> trails, string path, int count) {
		string full_count;
		if (count < 10) full_count = "0" + std::to_string(count);
		else full_count = std::to_string(count);

		for (int j = 0; j < trails.size(); j++) {
			string full_j;
			if (j < 10) full_j = "00" + std::to_string(j);
			else if (j < 100) full_j = "0" + std::to_string(j);
			else full_j = std::to_string(j);

			string recon = path + "/recon_long/string" + full_count + "_trail" + full_j + ".off";
			ofstream file(recon, ofstream::out);
			file << "OFF\n";
			file << trails[j].size() << "\t0\t0\n";
			for (int i = 0; i < trails[j].size(); i++) {
				file << trails[j][i][0] << "\t" << trails[j][i][1] << "\t" << trails[j][i][2] << "\n";
			}
			file.close();
		}
	}

	static void SaveTrail_epi(vector<int> startPointIdx, vector<vector<vector<int>>> totalPairString, vector<ProjectImage> &views, vector<int> t, vector<ReconCurve> &reconCurveList, string path, int iter, int count, int offsize, float considerArea, float lamda, float epiAreaThre, float errorMax) {
		//vector<Vector3f> colors = { Vector3f(255, 0, 0), Vector3f(0, 255, 0), Vector3f(0, 0, 255), Vector3f(255, 255, 0), Vector3f(255, 0, 255), Vector3f(0, 255, 255) };
		int fittingIterNum = 2;
		int start = views[t[0]].reconGroupIdxStart[count];
		for (int i = 0; i < startPointIdx.size() - 1; i++) {
			int stringSpace = startPointIdx[i + 1] - startPointIdx[i];
			//Vector3f color = colors[i % 6];
			Vector3f color;
			color[0] = rand() % 256; color[1] = rand() % 256; color[2] = rand() % 256;
			if (stringSpace > offsize && views[t[1]].epiAreaTag[startPointIdx[i]] == 0) {
				for (int j = 0; j < views[t[1]].intersectNum[startPointIdx[i]]; j++) {
					string full_count;
					if (count < 10) full_count = "0" + std::to_string(count);
					else full_count = std::to_string(count);

					string full_i;
					if (i < 10) full_i = "00" + std::to_string(i);
					else if (i < 100) full_i = "0" + std::to_string(i);
					else full_i = std::to_string(i);

					string filename = "iter" + std::to_string(iter) + "_part" + full_count + "_sub" + full_i + "_t" + std::to_string(j) + ".off";
					string recon = path + "/recon/" + filename;
					string recon_color = path + "/recon_color/" + filename;

					ReconCurve tempCurve;
					vector<Vector3f> reconstructPointList;
					for (int k = 0; k < stringSpace; k++) {
						if (totalPairString[i][j][k] != -1) {
							Vector3f temprecon = MatchFinder::computeClosest3DPoint(views[t[0]].img2DPoints[start + startPointIdx[i] + k], views[t[1]].intersect[totalPairString[i][j][k]], views[t[0]].projectionMatrix, views[t[1]].projectionMatrix);
							reconstructPointList.push_back(temprecon);
						}
					}
					int startRef = start + startPointIdx[i];
					int endRef = start + startPointIdx[i] + stringSpace - 1;

					if (reconstructPointList.size() > offsize) {
						/*ReconCurve::SavePoints(reconstructPointList, recon);
						ReconCurve::SaveColorPoints(reconstructPointList, recon_color, color);*/
						tempCurve.curve = reconstructPointList;
						tempCurve.curveName = filename;
						tempCurve.curve = reconstructPointList;
						tempCurve.iter = iter;
						tempCurve.partIdx = count;
						tempCurve.clusterIdx = i;
						tempCurve.idxInCluster = j;
						tempCurve.refIdx = t[0];
						tempCurve.startRef = startRef;
						tempCurve.endRef = endRef;
						tempCurve.color = color;
						reconCurveList.push_back(tempCurve);
					}
				}
				//pairString[i].clear();
			}
			else if (stringSpace > offsize && views[t[1]].epiAreaTag[startPointIdx[i]] == 1) {
				string full_count;
				if (count < 10) full_count = "0" + std::to_string(count);
				else full_count = std::to_string(count);

				string full_i;
				if (i < 10) full_i = "00" + std::to_string(i);
				else if (i < 100) full_i = "0" + std::to_string(i);
				else full_i = std::to_string(i);

				int sampleNum = startPointIdx[i + 1] - startPointIdx[i];

				vector<vector<Vector3f>> reconEpiAreaCurve = ReconCurve::ReconEpiArea(start, startPointIdx[i], startPointIdx[i + 1] - 1, views, t, considerArea, epiAreaThre, path, iter, count, i, sampleNum, lamda, fittingIterNum, color, errorMax);

				for (int s = 0; s < reconEpiAreaCurve.size(); s++) {
					ReconCurve tempCurve;
					if (reconEpiAreaCurve[s].size() > offsize) {
						string filename = "iter" + std::to_string(iter) + "_part" + full_count + "_sub" + full_i + "_t" + std::to_string(s) + ".off";
						string recon = path + "/recon/" + filename;
						string recon_color = path + "/recon_color/" + filename;
						string fitting = path + "/fitting/" + filename;
						/*ReconCurve::SavePoints(reconEpiAreaCurve[s], recon);
						ReconCurve::SaveColorPoints(reconEpiAreaCurve[s], recon_color, color);
						ReconCurve::SaveColorPoints(reconEpiAreaCurve[s], fitting, color);*/
						tempCurve.curve = reconEpiAreaCurve[s];
						tempCurve.curveName = filename;
						tempCurve.iter = iter;
						tempCurve.partIdx = count;
						tempCurve.clusterIdx = i;
						tempCurve.idxInCluster = s;
						tempCurve.refIdx = t[0];
						tempCurve.startRef = start + startPointIdx[i];
						tempCurve.endRef = start + startPointIdx[i + 1] - 1;
						tempCurve.color = color;
						reconCurveList.push_back(tempCurve);
					}
				}

			}
		}
	}

	static vector<Vector2f> FindTangentPoints(ProjectImage &refView, ProjectImage &neighborView, int offstep) {
		CalculateEpipole(refView, neighborView);
		vector<Vector2f> tangent_points;
		float eps = 10e-6;
		Vector2f epipole = neighborView.epipole;
		for (int i = 0; i < neighborView.groupIdxStart.size(); i++) {
			int startIdx = neighborView.groupIdxStart[i];
			int endIdx = startIdx + neighborView.pointsNumInGroup[i];
			for (int j = 0; j < neighborView.pointsNumInGroup[i]; j++) {
				Vector3f paraline;
				//cout << "epi:\n"<<epipole << endl;
				paraline[0] = epipole[1] - neighborView.img2DPoints[startIdx + j][1];
				paraline[1] = neighborView.img2DPoints[startIdx + j][0] - epipole[0];
				paraline[2] = neighborView.img2DPoints[startIdx + j][1] * epipole[0] - epipole[1] * neighborView.img2DPoints[startIdx + j][0];
				//cout << "para before:\n" << paraline << endl;
				paraline = paraline.normalized();
				//cout << "para:\n"<<paraline << endl;
				int sign = 0; int tag = 1;
				int start = max(j - offstep, 0);
				int end = min(j + offstep, neighborView.pointsNumInGroup[i]);

				//cout << "begin:" << j << endl;
				for (int k = start; k < end; k++) {
					float tempsign = paraline[0] * neighborView.img2DPoints[startIdx + k][0] + paraline[1] * neighborView.img2DPoints[startIdx + k][1] + paraline[2];
					//cout << tempsign << endl;
					if (sign == 0) {
						if (tempsign < -eps) sign = -1;
						else if (tempsign > eps) sign = 1;
					}
					else {
						if (tempsign*sign < -eps) {
							tag = 0; break;
						}
					}
				}
				if (tag == 1) {
					Vector2f para;
					para = paraline.head(2).normalized();
					float ang = fabs(para[0] * neighborView.img2DDires[startIdx + j][0] + para[1] * neighborView.img2DDires[startIdx + j][1]);
					if (ang < 0.02) {
						//cout << "start\t" << i << "\t" << j << "\t" << ang << "\t" << neighborView.img2DDires[startIdx + j][0] << "\t" << neighborView.img2DDires[startIdx + j][1] << endl;
						tangent_points.push_back(neighborView.img2DPoints[startIdx + j]);
					}
				}
			}
		}
		return tangent_points;
	}

	static void FindTangentPoints_epi(ProjectImage &refView, ProjectImage &neighborView, int offstep, float angMatchThre, float epiAreaThre) {
		CalculateEpipole(refView, neighborView);
		vector<Vector2f> tangent_points;
		double eps = 10e-8;
		Vector2f epipole = neighborView.epipole;
		for (int i = 0; i < neighborView.groupIdxStart.size(); i++) {
			int startIdx = neighborView.groupIdxStart[i];
			for (int j = offstep; j < neighborView.pointsNumInGroup[i] - offstep + 1; j++) {
				Vector3f paraline;
				//cout << "epi:\n"<<epipole << endl;
				paraline[0] = epipole[1] - neighborView.img2DPoints[startIdx + j][1];
				paraline[1] = neighborView.img2DPoints[startIdx + j][0] - epipole[0];
				paraline[2] = neighborView.img2DPoints[startIdx + j][1] * epipole[0] - epipole[1] * neighborView.img2DPoints[startIdx + j][0];
				paraline = paraline.normalized();
				int sign = 0; int tag = 1; int istangent = 0;
				int start = j - offstep;
				int end = j + offstep;

				//cout << "begin:" << j << endl;
				for (int k = start; k < end; k++) {
					if (k != j) {
						float tempsign = paraline[0] * neighborView.img2DPoints[startIdx + k][0] + paraline[1] * neighborView.img2DPoints[startIdx + k][1] + paraline[2];
						if (sign == 0) {
							if (tempsign < -eps) sign = -1;
							else if (tempsign > eps) sign = 1;
						}
						else {
							if (tempsign*sign < -eps) {
								tag = 0; break;
							}
						}
					}
				}

				if (tag == 1) {
					Vector2f para;
					para = paraline.head(2).normalized();
					float ang = fabs(para[0] * neighborView.img2DDires[startIdx + j][0] + para[1] * neighborView.img2DDires[startIdx + j][1]);
					if (ang < angMatchThre) {
						//cout << neighborView.tangentPoints.size() << "\t" << paraline[0] << "\t" << paraline[1] << "\t" << paraline[2] << endl;
						//cout << "start\t" << i<<"\t"<<j << "\t" << ang <<"\t"<< neighborView.img2DNorms[startIdx + j][0]<< "\t" << neighborView.img2DNorms[startIdx + j][1] << endl;
						Vector2f tangentPoint = neighborView.img2DPoints[startIdx + j];
						neighborView.tangentPoints.push_back(neighborView.img2DPoints[startIdx + j]);
						neighborView.tangentIdx.push_back(startIdx + j);
						neighborView.tangentPara.push_back(paraline);
						neighborView.tangentSign.push_back(sign);

						Vector2i tangentLR;
						for (int k = start; k < j; k++) {
							Vector3f tempPara;
							//cout << "epi:\n"<<epipole << endl;
							tempPara[0] = epipole[1] - neighborView.img2DPoints[startIdx + k][1];
							tempPara[1] = neighborView.img2DPoints[startIdx + k][0] - epipole[0];
							tempPara[2] = neighborView.img2DPoints[startIdx + k][1] * epipole[0] - epipole[1] * neighborView.img2DPoints[startIdx + k][0];

							double numer = tempPara[0] * tangentPoint[0] + tempPara[1] * tangentPoint[1] + tempPara[2];
							double distp2l = fabs(numer) / sqrt(tempPara[0] * tempPara[0] + tempPara[1] * tempPara[1]);

							if (distp2l < epiAreaThre) {
								tangentLR[0] = startIdx + k;
								break;
							}
						}
						for (int k = end-1; k > j; k--) {
							Vector3f tempPara;
							//cout << "epi:\n"<<epipole << endl;
							tempPara[0] = epipole[1] - neighborView.img2DPoints[startIdx + k][1];
							tempPara[1] = neighborView.img2DPoints[startIdx + k][0] - epipole[0];
							tempPara[2] = neighborView.img2DPoints[startIdx + k][1] * epipole[0] - epipole[1] * neighborView.img2DPoints[startIdx + k][0];

							double numer = tempPara[0] * tangentPoint[0] + tempPara[1] * tangentPoint[1] + tempPara[2];
							double distp2l = fabs(numer) / sqrt(tempPara[0] * tempPara[0] + tempPara[1] * tempPara[1]);

							if (distp2l < epiAreaThre) {
								tangentLR[1] = startIdx + k;
								break;
							}
						}
						neighborView.tangentLeftRight.push_back(tangentLR);
					}
				}
			}
		}
	}

	static void FindTangentPoints_epi2(ProjectImage &refView, ProjectImage &neighborView, int offstep, float angMatchThre) {
		neighborView.tangentPoints.clear();
		neighborView.tangentIdx.clear();
		CalculateEpipole(refView, neighborView);
		float eps = 10e-6;
		Vector2f epipole = neighborView.epipole;
		for (int i = 0; i < neighborView.groupIdxStart.size(); i++) {
			int startIdx = neighborView.groupIdxStart[i];
			for (int j = offstep; j < neighborView.pointsNumInGroup[i] - offstep + 1; j++) {
				Vector3f paraline;
				//cout << "epi:\n"<<epipole << endl;
				paraline[0] = epipole[1] - neighborView.img2DPoints[startIdx + j][1];
				paraline[1] = neighborView.img2DPoints[startIdx + j][0] - epipole[0];
				paraline[2] = neighborView.img2DPoints[startIdx + j][1] * epipole[0] - epipole[1] * neighborView.img2DPoints[startIdx + j][0];
				//cout << "para before:\n" << paraline << endl;
				paraline = paraline.normalized();
				//cout << "para:\n"<<paraline << endl;
				int sign = 0; int tag = 1;
				int start = j - offstep;
				int end = j + offstep;

				//cout << "begin:" << j << endl;
				for (int k = start; k < end; k++) {
					float tempsign = paraline[0] * neighborView.img2DPoints[startIdx + k][0] + paraline[1] * neighborView.img2DPoints[startIdx + k][1] + paraline[2];
					//cout << tempsign << endl;
					if (sign == 0) {
						if (tempsign < -eps) sign = -1;
						else if (tempsign > eps) sign = 1;
					}
					else {
						if (tempsign*sign < -eps) {
							tag = 0; break;
						}
					}
				}
				if (tag == 1) {
					Vector2f para;
					para = paraline.head(2).normalized();
					float ang = fabs(para[0] * neighborView.img2DDires[startIdx + j][0] + para[1] * neighborView.img2DDires[startIdx + j][1]);
					if (ang < angMatchThre) {
						//cout << "start\t" << i<<"\t"<<j << "\t" << ang <<"\t"<< neighborView.img2DNorms[startIdx + j][0]<< "\t" << neighborView.img2DNorms[startIdx + j][1] << endl;
						neighborView.tangentPoints.push_back(neighborView.img2DPoints[startIdx + j]);
						neighborView.tangentIdx.push_back(startIdx + j);
					}
				}
			}
		}
	}

	static void CalculateEpipole(ProjectImage &refView, ProjectImage &neighborView) {
		Matrix3f rot;
		rot = neighborView.rot*refView.rot.transpose();
		Vector3f trans;
		trans = neighborView.trans - rot*refView.trans;

		Vector3f e_v2 = neighborView.K*trans;
		e_v2 = e_v2 / e_v2[2];
		neighborView.epipole = e_v2.head(2);
	}

	static void SaveAllReconPoints(int count, ProjectImage refView, ProjectImage neighborView, string path, int x) {
		Vector3f color;
		color[0] = rand() % 256; color[1] = rand() % 256; color[2] = rand() % 256;
		vector<Vector3f> reconstructPointList_all(neighborView.intersect.size());
		int start = refView.reconGroupIdxStart[count];
		int pointsNum = refView.reconNumInGroup[count];
		for (int i = 0; i < pointsNum; i++) {
			for (int j = 0; j < neighborView.intersectNum[i]; j++) {
				reconstructPointList_all[neighborView.intersectIdxStart[i] + j] = MatchFinder::computeClosest3DPoint(refView.img2DPoints[start + i], neighborView.intersect[neighborView.intersectIdxStart[i] + j], refView.projectionMatrix, neighborView.projectionMatrix);
			}
		}
		string path1 = path + "/recon_all/iter" + std::to_string(x) + "_part" + std::to_string(count) + ".off";
		ReconCurve::SavePoints(reconstructPointList_all, path);
		string path2 = path + "/recon_all_color/iter" + std::to_string(x) + "_part" + std::to_string(count) + ".off";
		ReconCurve::SaveColorPoints(reconstructPointList_all, path2, color);
	}
};

