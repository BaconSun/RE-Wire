//#include "stdafx.h"
#include "MatchFinder.h"


Vector3f MatchFinder::computeClosest3DPoint(Vector2f refPixel, Vector2f neighborPixel, MatrixXf refCam, MatrixXf neighborCam)
{
	//compute the 3D point that fits best the two image projections
	//minimize: ||P*X - px||^2 where P is the projection matrix, X is the 3D point, px is the image projection

	MatrixXf A(4, 3);
	Vector4f b;
	Vector3f point3D;

	//first column
	A(0, 0) = refCam(0, 0) - refCam(2, 0)*refPixel[0];
	A(1, 0) = refCam(1, 0) - refCam(2, 0)*refPixel[1];
	A(2, 0) = neighborCam(0, 0) - neighborCam(2, 0)*neighborPixel[0];
	A(3, 0) = neighborCam(1, 0) - neighborCam(2, 0)*neighborPixel[1];

	//second column
	A(0, 1) = refCam(0, 1) - refCam(2, 1)*refPixel[0];
	A(1, 1) = refCam(1, 1) - refCam(2, 1)*refPixel[1];
	A(2, 1) = neighborCam(0, 1) - neighborCam(2, 1)*neighborPixel[0];
	A(3, 1) = neighborCam(1, 1) - neighborCam(2, 1)*neighborPixel[1];

	//third column
	A(0, 2) = refCam(0, 2) - refCam(2, 2)*refPixel[0];
	A(1, 2) = refCam(1, 2) - refCam(2, 2)*refPixel[1];
	A(2, 2) = neighborCam(0, 2) - neighborCam(2, 2)*neighborPixel[0];
	A(3, 2) = neighborCam(1, 2) - neighborCam(2, 2)*neighborPixel[1];


	//b
	b[0] = refCam(2, 3)*refPixel[0] - refCam(0, 3);
	b[1] = refCam(2, 3)*refPixel[1] - refCam(1, 3);
	b[2] = neighborCam(2, 3)*neighborPixel[0] - neighborCam(0, 3);
	b[3] = neighborCam(2, 3)*neighborPixel[1] - neighborCam(1, 3);

	point3D = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);

	return point3D;
}

Matrix3f MatchFinder::ComputeF(ProjectImage &refView, ProjectImage &neighborView) {
	Matrix3f rot;
	rot = neighborView.rot*refView.rot.transpose();
	Vector3f trans;
	trans = neighborView.trans - rot*refView.trans;

	Matrix3f tx;
	Matrix3f essentialMatrix;
	tx << 0.0, -trans[2], trans[1],
		trans[2], 0.0, -trans[0],
		-trans[1], trans[0], 0.0;
	essentialMatrix = tx*rot;
	Matrix3f fundamentalMatrix;
	fundamentalMatrix = ((neighborView.K.inverse()).transpose())*essentialMatrix*(refView.K.inverse());
	return fundamentalMatrix;
}

int MatchFinder::FindIntersectPoints(Vector3f para, vector<Vector2f> &pointsOnCurve2, vector<Vector2f> &intersect, int intersectSum, float dist_thre, float points_thre) {
	vector<float> distList(pointsOnCurve2.size());
	/* The number of intersection points (View2) for each point (View1) */
	int internum = 0;
	/* Calculate each point in view2 to the epipolar line*/
	for (int i = 0; i < pointsOnCurve2.size(); i++) {
		float numer = fabs(para[0] * pointsOnCurve2[i][0] + para[1] * pointsOnCurve2[i][1] + para[2]);
		float denomi = sqrt(para[0] * para[0] + para[1] * para[1]);
		float dist = numer / denomi;
		distList[i] = dist;
	}
	vector<int> idxList(pointsOnCurve2.size());
	for (int i = 0; i < idxList.size(); i++) {
		idxList[i] = i;
	}
	Sorting::QuickSort(distList, idxList, 0, distList.size() - 1);
	for (int i = 0; i < idxList.size(); i++) {
		if (!internum) { intersect.push_back(pointsOnCurve2[idxList[i]]); internum++; }
		else if (distList[idxList[i]] > dist_thre) break;
		else {
			int flag = 1;
			for (int j = 0; j < internum; j++) {
				if ((pointsOnCurve2[idxList[i]] - intersect[intersectSum + j]).norm() < points_thre) flag = 0;
			}
			if (flag) { intersect.push_back(pointsOnCurve2[idxList[i]]); internum++; }
		}
	}
	return internum;
}

int MatchFinder::FindIntersectPoints2(Vector3f para, ProjectImage& view) {
	/* The number of intersection points (View2) for each point (View1) */
	int internum = 0;
	double eps = 1e-10;
	float prenumer;
	/* Calculate each point in view2 to the epipolar line*/
	for (int j = 0; j < view.pointsNumInGroup.size(); j++) {
		for (int i = 0; i < view.pointsNumInGroup[j]; i++) {
			double numer = para[0] * view.img2DPoints[view.groupIdxStart[j] + i][0] + para[1] * view.img2DPoints[view.groupIdxStart[j] + i][1] + para[2];
			double dist = fabs(numer) / sqrt(para[0] * para[0] + para[1] * para[1]);
			if (dist < eps) {
				view.intersect.push_back(view.img2DPoints[view.groupIdxStart[j] + i]); internum++; //cout << startIdx[j] << "\t" << i << "\t x" << endl; 
			}
			if (i != 0) {
				Vector2f pre = view.img2DPoints[view.groupIdxStart[j] + i - 1];
				Vector2f cur = view.img2DPoints[view.groupIdxStart[j] + i];

				if (prenumer*numer < 0) {
					if (fabs(cur[0] - pre[0]) < eps) {
						Vector2f inter;
						inter[0] = cur[0];
						inter[1] = -(para[0] * inter[0] + para[2]) / para[1];
						view.intersect.push_back(inter); internum++;
					}
					else if (fabs(cur[1] - pre[1]) < eps) {
						Vector2f inter;
						inter[1] = cur[1];
						inter[0] = -(para[1] * inter[1] + para[2]) / para[0];

						view.intersect.push_back(inter); internum++;
					}
					else {
						float k = (cur[1] - pre[1]) / (cur[0] - pre[0]);
						float q = (cur[0] * pre[1] - pre[0] * cur[1]) / (cur[0] - pre[0]);
						Vector2f inter;
						inter[0] = -(para[1] * q + para[2]) / (para[0] + para[1] * k);
						inter[1] = k*inter[0] + q;
						view.intersect.push_back(inter); internum++;
					}
				}

			}
			prenumer = numer;
		}
	}
	return internum;
}

//int MatchFinder::FindIntersectPoints3(Vector3f para, ProjectImage& view, float p2ldist) {
//	/* The number of intersection points (View2) for each point (View1) */
//	int internum = 0;
//	double eps = 1e-10;
//	float prenumer;
//	/* Calculate each point in view2 to the epipolar line*/
//	for (int j = 0; j < view.pointsNumInGroup.size(); j++) {
//		for (int i = 0; i < view.pointsNumInGroup[j]; i++) {
//			double numer = para[0] * view.img2DPoints[view.groupIdxStart[j] + i][0] + para[1] * view.img2DPoints[view.groupIdxStart[j] + i][1] + para[2];
//			double dist = fabs(numer) / sqrt(para[0] * para[0] + para[1] * para[1]);
//			if (dist < eps) {
//				if (CheckIntersection(para, view.img2DPoints[view.groupIdxStart[j] + i], view.img2DNorms[view.groupIdxStart[j] + i], p2ldist, view.tangentPoints) == true) {
//					view.intersect.push_back(view.img2DPoints[view.groupIdxStart[j] + i]); internum++; //cout << startIdx[j] << "\t" << i << "\t x" << endl; 
//				}
//			}
//			if (i != 0) {
//				Vector2f pre = view.img2DPoints[view.groupIdxStart[j] + i - 1];
//				Vector2f cur = view.img2DPoints[view.groupIdxStart[j] + i];
//
//				if (prenumer*numer < 0) {
//					if (fabs(cur[0] - pre[0]) < eps) {
//						Vector2f inter;
//						inter[0] = cur[0];
//						inter[1] = -(para[0] * inter[0] + para[2]) / para[1];
//						if (CheckIntersection(para, view.img2DPoints[view.groupIdxStart[j] + i], view.img2DNorms[view.groupIdxStart[j] + i], p2ldist, view.tangentPoints) == true) {
//							view.intersect.push_back(inter); internum++;
//						}
//					}
//					else if (fabs(cur[1] - pre[1]) < eps) {
//						Vector2f inter;
//						inter[1] = cur[1];
//						inter[0] = -(para[1] * inter[1] + para[2]) / para[0];
//						if (CheckIntersection(para, view.img2DPoints[view.groupIdxStart[j] + i], view.img2DNorms[view.groupIdxStart[j] + i], p2ldist, view.tangentPoints) == true) {
//							view.intersect.push_back(inter); internum++;
//						}
//					}
//					else {
//						float k = (cur[1] - pre[1]) / (cur[0] - pre[0]);
//						float q = (cur[0] * pre[1] - pre[0] * cur[1]) / (cur[0] - pre[0]);
//						Vector2f inter;
//						inter[0] = -(para[1] * q + para[2]) / (para[0] + para[1] * k);
//						inter[1] = k*inter[0] + q;
//						if (CheckIntersection(para, view.img2DPoints[view.groupIdxStart[j] + i], view.img2DNorms[view.groupIdxStart[j] + i], p2ldist, view.tangentPoints) == true) {
//							view.intersect.push_back(inter); internum++;
//						}
//					}
//				}
//
//			}
//			prenumer = numer;
//		}
//	}
//	return internum;
//}
//
//bool MatchFinder::CheckIntersection(Vector3f para, Vector2f intersect, Vector2f norm, float p2ldist, vector<Vector2f> tangentPoints) {
//	Vector2f lineNorm = para.head(2).normalized();
//	float eps = 1e-6;
//	for (int i = 0; i < tangentPoints.size(); i++) {
//		double numer = para[0] * tangentPoints[i][0] + para[1] * tangentPoints[i][1] + para[2];
//		double distp2l = fabs(numer) / sqrt(para[0] * para[0] + para[1] * para[1]);
//		if (distp2l < p2ldist) {
//			double distp2p = (tangentPoints[i] - intersect).norm();
//			if (distp2p < p2ldist) {
//				float ang = fabs(lineNorm[0] * norm[0] + lineNorm[1] * norm[1]);
//				if (ang > 0.8) {
//					if ((tangentPoints[i] - intersect).norm() < eps) return true;
//					else {
//						//cout << "false:\t" << dist << "\t" << ang << endl;
//						return false;
//					}
//				}
//			}
//
//		}
//	}
//	return true;
//}


void MatchFinder::FindIntersectForView(int count, ProjectImage &refView, ProjectImage &neighborView, vector<ProjectImage> views, vector<int> t, float project_mindist) {
	int start = refView.reconGroupIdxStart[count];
	int pointsNum = refView.reconNumInGroup[count];
	int intersectSum = 0;
	int internum;
	/* For each point on a subcurve in view1 */
	for (int j = 0; j < pointsNum; j++) {
		/* Calculate its epipolar line */
		Vector3f temp = Vector3f::Ones();
		temp.head(2) = refView.img2DPoints[start + j];
		Vector3f para = neighborView.FMtoRefView*temp;
		
		//internum[i] = MatchFinder::FindIntersectPoints(para, Views[i].img2DPoints, Views[i].intersect, intersectSum[i], dist_thre, points_thre);
		internum = MatchFinder::FindIntersectPoints2(para, neighborView);
		
		//float p2ldist = 10.0;
		//internum[i] = MatchFinder::FindIntersectPoints3(para, Views[i], p2ldist);
		neighborView.intersectNum.push_back(internum);
		neighborView.intersectIdxStart.push_back(intersectSum);
		intersectSum += internum;
	}

	//for (int i = pointsNum-1; i > -1; i--) {
	//	int startIdx = neighborView.intersectIdxStart[i];
	//	for (int j = neighborView.intersectNum[i] - 1; j >-1; j--) {
	//		Vector3f temprecon = MatchFinder::computeClosest3DPoint(refView.img2DPoints[start + i], neighborView.intersect[startIdx + j], refView.projectionMatrix, neighborView.projectionMatrix);
	//		float max_dist;
	//		for (int k = 0; k < views.size(); k++) {
	//			Vector2f nearest; int idx;
	//			float dist = views[t[k]].FindClosestPointOnImgEachPoint(temprecon, nearest, idx);
	//			if (k == 0) max_dist = dist;
	//			else if (dist > max_dist) max_dist = dist;
	//		}
	//		if (sqrt(max_dist) > project_mindist) {
	//			/*cout << i << "\t" << j << endl;*/ neighborView.removeOneIntersect(i, j);
	//		}
	//	}
	//}
}

void MatchFinder::FindIntersectForView_epi(int count, vector<ProjectImage>& views, vector<int> t, float epiAreaThre, vector<Vector3f> &reconP) {
	int start = views[t[0]].reconGroupIdxStart[count];
	int pointsNum = views[t[0]].reconNumInGroup[count];
	int intersectSum = 0;
	int internum;
	vector<Vector2f> tangentPoints = views[t[1]].tangentPoints;
	/* For each point on a subcurve in view1 */
	for (int j = 0; j < pointsNum; j++) {
		/* Calculate its epipolar line */
		Vector3f temp = Vector3f::Ones();
		temp.head(2) = views[t[0]].img2DPoints[start + j];
		Vector3f para = views[t[1]].FMtoRefView*temp;
		para = para.normalized();

		int tempTag = 0;
		vector<int> tempEpiAreaIdx;
		for (int k = 0; k < tangentPoints.size(); k++) {
			double numer = para[0] * tangentPoints[k][0] + para[1] * tangentPoints[k][1] + para[2];
			double distp2l = fabs(numer) / sqrt(para[0] * para[0] + para[1] * para[1]);
			if (distp2l < epiAreaThre) {
				tempTag = 1; 
				tempEpiAreaIdx.push_back(k);
			}
		}
		views[t[1]].epiAreaTag.push_back(tempTag);
		views[t[1]].epiAreaIdx.push_back(tempEpiAreaIdx);

		if (tempTag == 0) {
			//internum[i] = MatchFinder::FindIntersectPoints(para, Views[i].img2DPoints, Views[i].intersect, intersectSum[i], dist_thre, points_thre);
			internum = MatchFinder::FindIntersectPoints2(para, views[t[1]]);

			views[t[1]].intersectNum.push_back(internum);
			views[t[1]].intersectIdxStart.push_back(intersectSum);
			intersectSum += internum;
		}
		else {
			internum = 0;
			views[t[1]].intersectNum.push_back(internum);
			views[t[1]].intersectIdxStart.push_back(intersectSum);
			intersectSum += internum;
		}
	}
}

int MatchFinder::FindIntersectForOnePoint_epi(Vector3f para, ProjectImage view, vector<Vector2f>& intersect) {
	/* The number of intersection points (View2) for each point (View1) */
	int internum = 0;
	double eps = 1e-6;
	float prenumer;
	/* Calculate each point in view2 to the epipolar line*/
	for (int j = 0; j < view.pointsNumInGroup.size(); j++) {
		for (int i = 0; i < view.pointsNumInGroup[j]; i++) {
			double numer = para[0] * view.img2DPoints[view.groupIdxStart[j] + i][0] + para[1] * view.img2DPoints[view.groupIdxStart[j] + i][1] + para[2];
			double dist = fabs(numer) / sqrt(para[0] * para[0] + para[1] * para[1]);
			if (dist < eps) {
				intersect.push_back(view.img2DPoints[view.groupIdxStart[j] + i]); internum++; //cout << startIdx[j] << "\t" << i << "\t x" << endl; 
			}
			if (i != 0) {
				Vector2f pre = view.img2DPoints[view.groupIdxStart[j] + i - 1];
				Vector2f cur = view.img2DPoints[view.groupIdxStart[j] + i];

				if (prenumer*numer < 0) {
					if (fabs(cur[0] - pre[0]) < eps) {
						Vector2f inter;
						inter[0] = cur[0];
						inter[1] = -(para[0] * inter[0] + para[2]) / para[1];
						intersect.push_back(inter); internum++;
					}
					else if (fabs(cur[1] - pre[1]) < eps) {
						Vector2f inter;
						inter[1] = cur[1];
						inter[0] = -(para[1] * inter[1] + para[2]) / para[0];

						intersect.push_back(inter); internum++;
					}
					else {
						float k = (cur[1] - pre[1]) / (cur[0] - pre[0]);
						float q = (cur[0] * pre[1] - pre[0] * cur[1]) / (cur[0] - pre[0]);
						Vector2f inter;
						inter[0] = -(para[1] * q + para[2]) / (para[0] + para[1] * k);
						inter[1] = k*inter[0] + q;
						intersect.push_back(inter); internum++;
					}
				}

			}
			prenumer = numer;
		}
	}
	return internum;
}

VectorXf MatchFinder::SolveForCP(MatrixXf paraBAll, VectorXf parabAll, vector<Vector3f>& updateCP) {
	VectorXf controlPoints(3 * updateCP.size());

	controlPoints = paraBAll.jacobiSvd(ComputeThinU | ComputeThinV).solve(parabAll);

	for (int i = 0; i < updateCP.size(); i++) {
		updateCP[i][0] = controlPoints[3 * i];
		updateCP[i][1] = controlPoints[3 * i + 1];
		updateCP[i][2] = controlPoints[3 * i + 2];
	}

	return controlPoints;
}

void MatchFinder::ShowTangentCut(string path, ProjectImage view, int vidx, float scale, Vector3f paraStart, Vector3f paraEnd, int vRefIdx, Vector3f refStart, Vector3f refEnd) {
	float scaleinv = 1.0 / scale;
	string imgname = path + "/view" + std::to_string(vidx) + "_bd1" + ".png";
	Mat img;
	Mat imS;

	img = cv::imread(imgname).clone();
	cv::resize(img, imS, Size(), scale, scale);
	for (int j = 0; j < view.tangentPoints.size(); j++) {
		Vector2i tLR = view.tangentLeftRight[j];
		cv::circle(imS, Point((int)(view.tangentPoints[j][0] / scaleinv), (int)(view.tangentPoints[j][1] / scaleinv)), 2, Scalar(0, 0, 255), -1, 8);
		cv::circle(imS, Point((int)(view.img2DPoints[tLR[0]][0] / scaleinv), (int)(view.img2DPoints[tLR[0]][1] / scaleinv)), 2, Scalar(255, 0, 0), -1, 8);
		cv::circle(imS, Point((int)(view.img2DPoints[tLR[1]][0] / scaleinv), (int)(view.img2DPoints[tLR[1]][1] / scaleinv)), 2, Scalar(255, 0, 0), -1, 8);
	}
	int colnum = imS.cols;
	Point start; start.x = 0; start.y = -(paraStart[2]/scaleinv) / paraStart[1];
	Point end; end.x = colnum; end.y = -(paraStart[0] * end.x + (paraStart[2] / scaleinv)) / paraStart[1];
	line(imS, start, end, Scalar(0, 255, 0), 1, 8, 0);

	start.x = 0; start.y = -(paraEnd[2]/scaleinv) / paraEnd[1];
	end.x = colnum; end.y = -(paraEnd[0] * end.x + (paraEnd[2] / scaleinv)) / paraEnd[1];
	line(imS, start, end, Scalar(255, 0, 0), 1, 8, 0);
	string title = "Tangent Cut: imgNei" + std::to_string(vidx);
	imshow(title, imS);

	string imgnameb = path + "/view" + std::to_string(vRefIdx) + "_bd1" + ".png";
	Mat imgb;
	Mat imSb;

	imgb = cv::imread(imgnameb).clone();
	cv::resize(imgb, imSb, Size(), scale, scale);
	cv::circle(imSb, Point((int)(refStart[0] / scaleinv), (int)(refStart[1] / scaleinv)), 2, Scalar(0, 255, 0), -1, 8);
	cv::circle(imSb, Point((int)(refEnd[0] / scaleinv), (int)(refEnd[1] / scaleinv)), 2, Scalar(255, 0, 0), -1, 8);
	string titleb = "Tangent Cut: imgRef" + std::to_string(vRefIdx);
	imshow(titleb, imSb);
	waitKey(0);
}
