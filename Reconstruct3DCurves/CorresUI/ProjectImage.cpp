#include "ProjectImage.h"

void ProjectImage::LoadView(string path, int Idxnum, char type) {
	string Idx = std::to_string(Idxnum);
	string cam_name = path + "/camera" + Idx + ".txt";
	LoadCamera(cam_name, type);
	string img_name = path + "/curvepoints" + Idx + ".txt";
	LoadImagePoints3(img_name);
}

void ProjectImage::LoadCamera(string filename, char type) {
	if (type == 'r') {
		ifstream file(filename, ifstream::in);
		int lineSize = 1024;
		char *buffer = new char[lineSize];
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\n", &focalLength);
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\n", &cameraCenter[0], &cameraCenter[1]);
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &trans[0], &trans[1], &trans[2]);
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &rot(0, 0), &rot(0, 1), &rot(0, 2));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &rot(1, 0), &rot(1, 1), &rot(1, 2));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &rot(2, 0), &rot(2, 1), &rot(2, 2));

		K << focalLength, 0.0, cameraCenter[0],
			0.0, focalLength, cameraCenter[1],
			0.0, 0.0, 1.0;

		MatrixXf RT(3, 4);
		RT.block(0, 0, 3, 3) = rot;
		RT.block(0, 3, 3, 1) = trans;

		projectionMatrix = K*RT;
		/*cout << "PM:\n" << projectionMatrix << endl;
		cout << "K:\n" << K << endl;
		cout << "rot:\n" << rot << endl;
		cout << "trans:\n" << trans << endl;*/
	}
	else if (type == 'm') {
		ifstream file(filename, ifstream::in);
		int lineSize = 1024;
		char *buffer = new char[lineSize];
		MatrixXf PM(3, 4);
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\t%f\n", &PM(0, 0), &PM(0, 1), &PM(0, 2), &PM(0, 3));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\t%f\n", &PM(1, 0), &PM(1, 1), &PM(1, 2), &PM(1, 3));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\t%f\n", &PM(2, 0), &PM(2, 1), &PM(2, 2), &PM(2, 3));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &K(0, 0), &K(0, 1), &K(0, 2));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &K(1, 0), &K(1, 1), &K(1, 2));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &K(2, 0), &K(2, 1), &K(2, 2));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &rot(0, 0), &rot(0, 1), &rot(0, 2));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &rot(1, 0), &rot(1, 1), &rot(1, 2));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &rot(2, 0), &rot(2, 1), &rot(2, 2));
		file.getline(buffer, lineSize);
		sscanf(buffer, "%f\t%f\t%f\n", &trans[0], &trans[1], &trans[2]);
		projectionMatrix = PM;
		/*cout << "PM:\n" << PM << endl;
		cout << "K:\n" << K << endl;
		cout << "rot:\n" << rot << endl;
		cout << "trans:\n" << trans << endl;*/
	}
	else { cout << "please input a right type!" << endl; }
}

void ProjectImage::LoadImagePoints(string filename) {

	ifstream file(filename, ifstream::in);
	int lineSize = 1024; int tempgroup;
	char *buffer = new char[lineSize];
	Vector2i tempPointOnCurve; Vector2f tempPointOnCurvefloat;

	file.getline(buffer, lineSize);
	int groupNum;
	sscanf(buffer, "%d\n", &groupNum);

	//cout << groupNum << endl;
	file.getline(buffer, lineSize);
	sscanf(buffer, "%d\t%d\t%d\n", &tempgroup, &tempPointOnCurve[1], &tempPointOnCurve[0]);

	int count = 0;
	int pointsIdx = 0;
	while (count < groupNum) {
		int tempGroupNum = 0;
		do {
			tempPointOnCurvefloat[0] = (float)tempPointOnCurve[0]; tempPointOnCurvefloat[1] = (float)tempPointOnCurve[1];
			img2DPoints.push_back(tempPointOnCurvefloat);
			tempGroupNum++;

			file.getline(buffer, lineSize);
			sscanf(buffer, "%d\t%d\t%d\n", &tempgroup, &tempPointOnCurve[1], &tempPointOnCurve[0]);
			//cout << tempPointOnCurve[0] << "\t" << tempPointOnCurve[1] << endl;

		} while ((tempgroup == count + 1) && !file.eof());
		pointsNumInGroup.push_back(tempGroupNum);
		groupIdxStart.push_back(pointsIdx);
		pointsIdx += tempGroupNum;
		count++;
	}
	file.close();
	BuildKdTree(img2DPoints);
}

/* No smooth*/
void ProjectImage::LoadImagePoints2(string filename) {

	ifstream file(filename, ifstream::in);
	int lineSize = 1024; int tempgroup;
	char *buffer = new char[lineSize];
	Vector2i tempPointOnCurve; Vector2f tempPointOnCurvefloat;
	Vector2f tempNormal;

	file.getline(buffer, lineSize);
	int groupNum;
	sscanf(buffer, "%d\n", &groupNum);

	//cout << groupNum << endl;
	file.getline(buffer, lineSize);
	sscanf(buffer, "%d\t%d\t%d\t%f\t%f\n", &tempgroup, &tempPointOnCurve[1], &tempPointOnCurve[0], &tempNormal[1], &tempNormal[0]);

	int count = 0;
	int pointsIdx = 0;
	while (count < groupNum) {
		int tempGroupNum = 0;
		do {
			tempPointOnCurvefloat[0] = (float)tempPointOnCurve[0]; tempPointOnCurvefloat[1] = (float)tempPointOnCurve[1];
			img2DPoints.push_back(tempPointOnCurvefloat);
			tempGroupNum++;

			file.getline(buffer, lineSize);
			sscanf(buffer, "%d\t%d\t%d\t%f\t%f\n", &tempgroup, &tempPointOnCurve[1], &tempPointOnCurve[0], &tempNormal[1], &tempNormal[0]);
			//cout << tempPointOnCurve[0] << "\t" << tempPointOnCurve[1] << endl;

		} while ((tempgroup == count + 1) && !file.eof());
		pointsNumInGroup.push_back(tempGroupNum);
		groupIdxStart.push_back(pointsIdx);
		pointsIdx += tempGroupNum;
		//cout << tempGroupNum << endl;
		count++;
	}
	file.close();
	reconNumInGroup = pointsNumInGroup;
	reconGroupIdxStart = groupIdxStart;
	BuildKdTree(img2DPoints);
}

/* Smooth the curve and calculate the tangent for each point*/
void ProjectImage::LoadImagePoints3(string filename) {

	ifstream file(filename, ifstream::in);
	int lineSize = 1024; int tempgroup;
	char *buffer = new char[lineSize];
	Vector2i tempPointOnCurve; 
	Vector2f tempNormal;
	float eps = 1e-6;

	file.getline(buffer, lineSize);
	int groupNum;
	sscanf(buffer, "%d\n", &groupNum);

	vector<Vector2i> temp2DPoints;
	vector<int> tempPointsNumInGroup;
	vector<int> tempGroupIdxStart;

	/* Import point data */
	file.getline(buffer, lineSize);
	sscanf(buffer, "%d\t%d\t%d\t%f\t%f\n", &tempgroup, &tempPointOnCurve[1], &tempPointOnCurve[0], &tempNormal[1], &tempNormal[0]);

	int count = 0;
	int pointsIdx = 0;
	while (count < groupNum) {
		int tempGroupNum = 0;
		do {
			//tempPointOnCurvefloat[0] = (float)tempPointOnCurve[0]; tempPointOnCurvefloat[1] = (float)tempPointOnCurve[1];
			temp2DPoints.push_back(tempPointOnCurve);
			tempGroupNum++;
			file.getline(buffer, lineSize);
			sscanf(buffer, "%d\t%d\t%d\t%f\t%f\n", &tempgroup, &tempPointOnCurve[1], &tempPointOnCurve[0], &tempNormal[1], &tempNormal[0]);
		} while ((tempgroup == count + 1) && !file.eof());
		tempPointsNumInGroup.push_back(tempGroupNum);
		tempGroupIdxStart.push_back(pointsIdx);
		pointsIdx += tempGroupNum;
		count++;
	}
	file.close();

	//cout << "aa:\t"<<temp2DPoints.size() << endl;

	/* Find Turning Points */
	Vector2f tempFloat;
	vector<Vector2f> temp2DPoints2;
	vector<int> tempPointsNumInGroup2;
	vector<int> tempGroupIdxStart2;
	float preDire, nextDire;
	pointsIdx = 0;
	vector<vector<int>> turnPoints2(groupNum);
	for (int i = 0; i < groupNum; i++) {
		int start = tempGroupIdxStart[i];
		int tempGroupNum2 = 0;
		for (int j = 0; j < tempPointsNumInGroup[i]; j++) {
			if (j == 0) {
				tempFloat(0) = (float)temp2DPoints[start](0); tempFloat(1) = (float)temp2DPoints[start](1);
				temp2DPoints2.push_back(tempFloat);
				tempGroupNum2++;
				turnPoints2[i].push_back(j);
				if (temp2DPoints[start + j + 1](0) == temp2DPoints[start + j](0)) nextDire = 10000.0;
				else nextDire = (float)(temp2DPoints[start + j + 1](1) - temp2DPoints[start + j](1)) / (float)(temp2DPoints[start + j + 1](0) - temp2DPoints[start + j](0));
			}
			else if (j == tempPointsNumInGroup[i] - 1) {
				tempFloat(0) = (float)temp2DPoints[start + j](0); tempFloat(1) = (float)temp2DPoints[start + j](1);
				temp2DPoints2.push_back(tempFloat);
				tempGroupNum2++;
				turnPoints2[i].push_back(j);
			}
			else {
				preDire = nextDire;
				if (temp2DPoints[start + j + 1](0) == temp2DPoints[start + j](0)) nextDire = 10000.0;
				else nextDire = (float)(temp2DPoints[start + j + 1](1) - temp2DPoints[start + j](1)) / (float)(temp2DPoints[start + j + 1](0) - temp2DPoints[start + j](0));
				if (fabs(preDire - nextDire) < eps) {
					tempFloat(0) = (float)temp2DPoints[start + j](0); tempFloat(1) = (float)temp2DPoints[start + j](1);
					temp2DPoints2.push_back(tempFloat);
					tempGroupNum2++;
					turnPoints2[i].push_back(j);
				}
			}
		}
		tempPointsNumInGroup2.push_back(tempGroupNum2);
		tempGroupIdxStart2.push_back(pointsIdx);
		pointsIdx += tempGroupNum2;
	}

	//cout << "bb:\t" << temp2DPoints2.size() << endl;

	/* Calculate the mid points on each line segements */
	pointsIdx = 0;
	Vector2f tempPoint3;
	vector<vector<int>> turnPoints3(groupNum);
	for (int i = 0; i < groupNum; i++) {
		int start = tempGroupIdxStart2[i];
		int tempGroupNum3 = 0;
		for (int j = 0; j < tempPointsNumInGroup2[i]; j++) {
			if (j == 0 || j == tempPointsNumInGroup2[i] - 1) {
				neiImg2DPoints.push_back(temp2DPoints2[start + j]);
				tempGroupNum3++;
				turnPoints3[i].push_back(turnPoints2[i][j]);
			}
			else if (j < tempPointsNumInGroup2[i] - 2) {
				tempPoint3 = (temp2DPoints2[start + j] + temp2DPoints2[start + j + 1]) / 2.0;
				neiImg2DPoints.push_back(tempPoint3);
				tempGroupNum3++;
				float tempTurnIdx = (turnPoints2[i][j] + turnPoints2[i][j + 1]) / 2.0;
				turnPoints3[i].push_back(tempTurnIdx);
			}
		}
		neiPointsNumInGroup.push_back(tempGroupNum3);
		neiGroupIdxStart.push_back(pointsIdx);
		pointsIdx += tempGroupNum3;
		//cout << tempPointsNumInGroup[i] - 1 << endl;
	}

	//cout << "cc:\t" << neiImg2DPoints.size() << endl;

	/* Map the original points to the line segments */
	pointsIdx = 0;
	for (int i = 0; i < groupNum; i++) {
		int start = neiGroupIdxStart[i];
		int tempGroupNum4 = 0;
		for (int j = 0; j < neiPointsNumInGroup[i]-1; j++) {
			int num = round(turnPoints3[i][j + 1]) - round(turnPoints3[i][j]);
			float step = 1.0 / num;
			for (int t = 0; t < num; t++) {
				Vector2f inter = (1.0 - step*t)*neiImg2DPoints[start+j] + (step*t)*neiImg2DPoints[start+j + 1];
				img2DPoints.push_back(inter);
				tempGroupNum4++;
			}
		}
		img2DPoints.push_back(neiImg2DPoints[start + neiPointsNumInGroup[i] - 1]);
		tempGroupNum4++;
		pointsNumInGroup.push_back(tempGroupNum4);
		groupIdxStart.push_back(pointsIdx);
		pointsIdx += tempGroupNum4;
	}
	//cout << "dd:\t" << img2DPoints.size() << endl;

	/* Calculate tangents for each point */
	for (int i = 0; i < groupNum; i++) {
		int start = groupIdxStart[i];
		int tempGroupNum5 = 0;
		for (int j = 0; j < pointsNumInGroup[i]; j++) {
			int normStart = max(j - 3, 0);
			int normEnd = min(j + 3, pointsNumInGroup[i]);
			Vector2f sumDire = Vector2f::Zero();
			for (int t = normStart; t < normEnd - 1; t++) {
				Vector2f dire = img2DPoints[start + t+1] - img2DPoints[start + t];
				dire.normalize();
				sumDire += dire;
			}
			sumDire.normalize();
			img2DDires.push_back(sumDire);
		}
	}
	reconNumInGroup = pointsNumInGroup;
	reconGroupIdxStart = groupIdxStart;
	
	//cout << "ee" << endl;
	BuildKdTree(img2DPoints);
}

void ProjectImage::LoadImagePoints4(string filename) {
	ifstream file(filename, ifstream::in);
	int lineSize = 1024; int tempgroup;
	char *buffer = new char[lineSize];
	Vector2i tempPointOnCurve; Vector2f tempPointOnCurvefloat;
	Vector2f tempNormal;
	float eps = 1e-6;

	file.getline(buffer, lineSize);
	int groupNum;
	sscanf(buffer, "%d\n", &groupNum);

	vector<Vector2f> temp2DPoints;
	vector<Vector2f> temp2DNorms;
	vector<int> tempPointsNumInGroup;
	vector<int> tempGroupIdxStart;

	file.getline(buffer, lineSize);
	sscanf(buffer, "%d\t%d\t%d\t%f\t%f\n", &tempgroup, &tempPointOnCurve[1], &tempPointOnCurve[0], &tempNormal[1], &tempNormal[0]);

	int count = 0;
	int pointsIdx = 0;
	while (count < groupNum) {
		int tempGroupNum = 0;
		do {
			tempPointOnCurvefloat[0] = (float)tempPointOnCurve[0]; tempPointOnCurvefloat[1] = (float)tempPointOnCurve[1];
			temp2DPoints.push_back(tempPointOnCurvefloat);
			temp2DNorms.push_back(tempNormal);
			tempGroupNum++;
			file.getline(buffer, lineSize);
			sscanf(buffer, "%d\t%d\t%d\t%f\t%f\n", &tempgroup, &tempPointOnCurve[1], &tempPointOnCurve[0], &tempNormal[1], &tempNormal[0]);
			//cout << tempPointOnCurve[0] << "\t" << tempPointOnCurve[1] << endl;

		} while ((tempgroup == count + 1) && !file.eof());
		tempPointsNumInGroup.push_back(tempGroupNum);
		tempGroupIdxStart.push_back(pointsIdx);
		pointsIdx += tempGroupNum;
		count++;
	}
	file.close();



	//cout << "aa:\t" << temp2DPoints.size() << endl;

	vector<Vector2f> temp2DPoints2;
	vector<Vector2f> temp2DNorms2;
	vector<int> tempPointsNumInGroup2;
	vector<int> tempGroupIdxStart2;
	float preDire, nextDire;
	pointsIdx = 0;
	vector<vector<int>> turnPoints(groupNum);
	for (int i = 0; i < groupNum; i++) {
		int start = tempGroupIdxStart[i];
		int tempGroupNum2 = 0;
		for (int j = 0; j < tempPointsNumInGroup[i]; j++) {
			if (j == 0) {
				temp2DPoints2.push_back(temp2DPoints[start]);
				temp2DNorms2.push_back(temp2DNorms[start]);
				tempGroupNum2++;
				turnPoints[i].push_back(j);
				if (temp2DPoints[start + j](0) == temp2DPoints[start + j - 1](0)) preDire = 10000.0;
				else preDire = (temp2DPoints[start + j](1) - temp2DPoints[start + j - 1](1)) / (temp2DPoints[start + j](0) - temp2DPoints[start + j - 1](0));
			}
			else if (j == tempPointsNumInGroup[i] - 1) {
				temp2DPoints2.push_back(temp2DPoints[start + j]);
				temp2DNorms2.push_back(temp2DNorms[start + j]);
				tempGroupNum2++;
				turnPoints[i].push_back(j);
			}
			else {
				preDire = nextDire;
				if (temp2DPoints[start + j + 1](0) == temp2DPoints[start + j](0)) nextDire = 10000.0;
				else nextDire = (temp2DPoints[start + j + 1](1) - temp2DPoints[start + j](1)) / (temp2DPoints[start + j + 1](0) - temp2DPoints[start + j](0));
				if (fabs(preDire - nextDire) < eps) {
					temp2DPoints2.push_back(temp2DPoints[start + j]);
					temp2DNorms2.push_back(temp2DNorms[start + j]);
					tempGroupNum2++;
					turnPoints[i].push_back(j);
				}
			}
		}
		tempPointsNumInGroup2.push_back(tempGroupNum2);
		tempGroupIdxStart2.push_back(pointsIdx);
		pointsIdx += tempGroupNum2;
	}

	//cout << "bb:\t" << temp2DPoints2.size() << endl;

	pointsIdx = 0;
	Vector2f tempPoint2;
	vector<vector<int>> turnPoints2(groupNum);
	for (int i = 0; i < groupNum; i++) {
		int start = tempGroupIdxStart2[i];
		int tempGroupNum3 = 0;
		for (int j = 0; j < tempPointsNumInGroup2[i]; j++) {
			if (j == 0 || j == tempPointsNumInGroup2[i] - 1) {
				neiImg2DPoints.push_back(temp2DPoints2[start + j]);
				tempGroupNum3++;
				turnPoints2[i].push_back(turnPoints[i][j]);
			}
			else if (j < tempPointsNumInGroup2[i] - 2) {
				tempPoint2 = (temp2DPoints2[start + j] + temp2DPoints2[start + j + 1]) / 2.0;
				neiImg2DPoints.push_back(tempPoint2);
				tempGroupNum3++;
				float tempTurnIdx = (turnPoints[i][j] + turnPoints[i][j + 1]) / 2.0;
				turnPoints2[i].push_back(tempTurnIdx);
			}
		}
		neiPointsNumInGroup.push_back(tempGroupNum3);
		neiGroupIdxStart.push_back(pointsIdx);
		pointsIdx += tempGroupNum3;
		//cout << tempPointsNumInGroup[i] - 1 << endl;
	}

	//cout << "cc:\t" << neiImg2DPoints.size() << endl;

	pointsIdx = 0;
	for (int i = 0; i < groupNum; i++) {
		int start = neiGroupIdxStart[i];
		int tempGroupNum4 = 0;
		for (int j = 0; j < neiPointsNumInGroup[i] - 1; j++) {
			int num = round(turnPoints2[i][j + 1]) - round(turnPoints2[i][j]);
			float step = 1.0 / num;
			for (int t = 0; t < num; t++) {
				Vector2f inter = (1.0 - step*t)*neiImg2DPoints[start + j] + (step*t)*neiImg2DPoints[start + j + 1];
				img2DPoints.push_back(inter);
				tempGroupNum4++;
			}
		}
		img2DPoints.push_back(neiImg2DPoints[start + neiPointsNumInGroup[i] - 1]);
		tempGroupNum4++;
		pointsNumInGroup.push_back(tempGroupNum4);
		groupIdxStart.push_back(pointsIdx);
		pointsIdx += tempGroupNum4;
	}
	//cout << "dd:\t" << img2DPoints.size() << endl;

	reconNumInGroup = pointsNumInGroup;
	reconGroupIdxStart = groupIdxStart;

	BuildKdTree(img2DPoints);
}

//void ProjectImage::DrawProjectImage(vector<Vector3f> points3D, string name) {
//	Mat img(height, width, CV_8UC3, Scalar::all(UCHAR_MAX));
//	Vector3f pixel; Vector4f temp;
//
//	for (int i = 0; i < points3D.size(); i++) {
//		temp.head(3) = points3D[i]; temp[3] = 1.0;
//		pixel = projectionMatrix*temp;
//		pixel = pixel / pixel[2];
//		for (int off = -4; off < 4; off++) {
//			Vec3b& bgr = img.at<Vec3b>((int)(pixel[1] + off), (int)(pixel[0] + off));
//			bgr[0] = 0; bgr[1] = 0; bgr[2] = 0;
//		}
//	}
//	Mat img_s;
//	resize(img, img_s, Size(), scale, scale);
//	imwrite(name, img_s);
//}

int ProjectImage::removeOneIntersect(int part, int k) {
	int start = intersectIdxStart[part];
	intersect.erase(intersect.begin() + start + k);
	//cout << start + k << endl;
	intersectNum[part]--;
	for (int i = part + 1; i < intersectIdxStart.size(); i++) {
		intersectIdxStart[i]--;
	}
	return 0;
}

void ProjectImage::BuildKdTree(vector<Vector2f> imgPoints) {
	ANNpointArray dataPts;
	dataPts = annAllocPts(imgPoints.size(), 2);
	for (int j = 0; j < imgPoints.size(); j++)
	{
		dataPts[j][0] = (double)imgPoints[j][0];
		dataPts[j][1] = (double)imgPoints[j][1];
	}

	kdTree = std::make_shared<ANNkd_tree>(dataPts, imgPoints.size(), 2);
}

float ProjectImage::FindClosestPointOnImgEachPoint(Vector3f point3D, Vector2f& nearestPoint, int &idx) {
	ANNpoint queryPt;
	ANNidxArray nnIdx;
	ANNdistArray dists;
	Vector3f imgPoint;
	double eps = 0.0;

	queryPt = annAllocPt(2);
	nnIdx = new ANNidx[1]; dists = new ANNdist[1];

	Vector4f tempPoint3D;
	tempPoint3D[0] = point3D[0]; tempPoint3D[1] = point3D[1]; tempPoint3D[2] = point3D[2]; tempPoint3D[3] = 1.0;
	imgPoint = projectionMatrix*tempPoint3D;
	imgPoint = imgPoint / imgPoint[2];
	
	queryPt[0] = (double)imgPoint[0];
	queryPt[1] = (double)imgPoint[1];
	/*cout << "yes" << endl;
	cout << imgPoint[0] << "\t" << imgPoint[1] << endl;*/

	kdTree->annkSearch(queryPt, 1, nnIdx, dists, eps);
	nearestPoint[0] = img2DPoints[nnIdx[0]][0];
	nearestPoint[1] = img2DPoints[nnIdx[0]][1];
	idx = nnIdx[0];

	delete nnIdx; // clean things up
	delete dists;
	//delete kdTree;
	//annClose(); // done with ANN

	return dists[0];
}

float ProjectImage::FindClosestPointOnImgEachPoint2(Vector3f point3D, Vector2f& nearestPoint, int &idx) {
	ANNpoint queryPt;
	ANNidxArray nnIdx;
	ANNdistArray dists;
	Vector3f imgPoint;
	double eps = 0.0;

	queryPt = annAllocPt(2);
	nnIdx = new ANNidx[1]; dists = new ANNdist[1];

	Vector4f tempPoint3D;
	tempPoint3D[0] = point3D[0]; tempPoint3D[1] = point3D[1]; tempPoint3D[2] = point3D[2]; tempPoint3D[3] = 1.0;
	imgPoint = projectionMatrix*tempPoint3D;
	imgPoint = imgPoint / imgPoint[2];

	queryPt[0] = (double)imgPoint[0];
	queryPt[1] = (double)imgPoint[1];
	/*cout << "yes" << endl;
	cout << imgPoint[0] << "\t" << imgPoint[1] << endl;*/

	kdTree->annkSearch(queryPt, 1, nnIdx, dists, eps);
	nearestPoint[0] = img2DPoints[nnIdx[0]][0];
	nearestPoint[1] = img2DPoints[nnIdx[0]][1];
	idx = nnIdx[0];

	delete nnIdx; // clean things up
	delete dists;
	//delete kdTree;
	//annClose(); // done with ANN

	return dists[0];
}

void ProjectImage::FindN_NearestPointOnImgEachPoint(Vector3f point3D, vector<Vector2f>& nearestPoint, vector<int> &idx, int numNearest) {
	ANNpoint queryPt;
	ANNidxArray nnIdx;
	ANNdistArray dists;
	Vector3f imgPoint;
	double eps = 0.0;

	queryPt = annAllocPt(2);
	nnIdx = new ANNidx[numNearest]; dists = new ANNdist[numNearest];

	Vector4f tempPoint3D;
	tempPoint3D[0] = point3D[0]; tempPoint3D[1] = point3D[1]; tempPoint3D[2] = point3D[2]; tempPoint3D[3] = 1.0;
	imgPoint = projectionMatrix*tempPoint3D;
	imgPoint = imgPoint / imgPoint[2];

	queryPt[0] = (double)imgPoint[0];
	queryPt[1] = (double)imgPoint[1];
	/*cout << "yes" << endl;
	cout << imgPoint[0] << "\t" << imgPoint[1] << endl;*/

	kdTree->annkSearch(queryPt, numNearest, nnIdx, dists, eps);

	for (int i = 0; i < numNearest; i++) {
		nearestPoint[i][0] = img2DPoints[nnIdx[i]][0];
		nearestPoint[i][1] = img2DPoints[nnIdx[i]][1];
		idx[i] = nnIdx[i];
	}
	

	delete nnIdx; // clean things up
	delete dists;
	//delete kdTree;
	//annClose(); // done with ANN
}


void ProjectImage::CalculateFirstMissingSE(vector<ProjectImage> views, vector<vector<Vector2i>>& firstMissingSE) {
	for (int i = 0; i < views.size(); i++) {
		vector<Vector2i> mse;
		for (int j = 0; j < views[i].pointsNumInGroup.size(); j++) {
			Vector2i se;
			se(0) = views[i].groupIdxStart[j];
			se(1) = se(0) + views[i].pointsNumInGroup[j] - 1;

			mse.push_back(se);
		}
		firstMissingSE.push_back(mse);
	}
}

void ProjectImage::RcCalculateB(int num, vector<Vector2f> nearestPoint, MatrixXf& para_B, VectorXf& para_b) {
	MatrixXf tempB(2, 3 * num); Vector2f tempb;

	for (int i = 1; i < num - 1; i++) {
		RcCalculateEachPointB(nearestPoint[i], i, num, tempB, tempb);
		para_B.block(2 * (i - 1), 0, 2, 3 * num) = tempB;
		para_b.block(2 * (i - 1), 0, 2, 1) = tempb;
	}
}

void ProjectImage::RcCalculateEachPointB(Vector2f nP, int i, int num, MatrixXf& para_B, Vector2f& para_b) {
	MatrixXf A(2, 3);
	A(0, 0) = projectionMatrix(0, 0) - nP[0] * projectionMatrix(2, 0);
	A(1, 0) = projectionMatrix(1, 0) - nP[1] * projectionMatrix(2, 0);
	A(0, 1) = projectionMatrix(0, 1) - nP[0] * projectionMatrix(2, 1);
	A(1, 1) = projectionMatrix(1, 1) - nP[1] * projectionMatrix(2, 1);
	A(0, 2) = projectionMatrix(0, 2) - nP[0] * projectionMatrix(2, 2);
	A(1, 2) = projectionMatrix(1, 2) - nP[1] * projectionMatrix(2, 2);

	MatrixXf T = MatrixXf::Zero(3, 3 * num);
	T(0, i * 3) = T(1, i * 3 + 1) = T(2, i * 3 + 2) = 1.0;

	para_B = A*T;

	//cout << "para_B:\n" << para_B << endl;
	para_b[0] = nP[0] * projectionMatrix(2, 3) - projectionMatrix(0, 3);
	para_b[1] = nP[1] * projectionMatrix(2, 3) - projectionMatrix(1, 3);
}

void ProjectImage::refresh_eachPart() {
	intersect.clear();
	intersectNum.clear();
	intersectIdxStart.clear();
	epiAreaTag.clear();
	epiAreaIdx.clear();
}

void ProjectImage::refresh_eachIter() {
	tangentPoints.clear();
	tangentSign.clear();
	tangentIdx.clear();
	tangentLeftRight.clear();
	tangentPara.clear();
}
