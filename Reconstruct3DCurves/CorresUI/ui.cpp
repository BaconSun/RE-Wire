#include "ui.h"
Mat UI::img_a_;
Mat UI::img_b_;
vector<cv::Point2i>UI::handles_a_;
vector<cv::Point2i>UI::handles_b_;
float UI::scale_;
vector<Vector2f> UI::corres_a_;
vector<Vector2f> UI::corres_b_;
Matrix3f UI::em_a_;
vector<Point2i> UI::points_a_;
vector<Point2i> UI::points_b_;
int idx_a;
float linethick = 1;
vector<int> UI::corres_num_;
int interp = 0;
vector<int> startIdx;
int pstep = 10;
string imgIdx = "1";

void UI::SetImages(Mat img_a, Mat img_b, float scale, Matrix3f& em_a, vector<Vector2f> &pointsOnCurve, vector<Vector2f> &intersection, vector<int> &intersectStart, vector<int> interNum)
{
	corres_a_.clear();
	corres_b_.clear();

	points_a_.clear();
	points_b_.clear();

	img_a_ = img_a;
	img_b_ = img_b;
	handles_b_.clear();
	handles_a_.clear();
	scale_ = scale;
	em_a_ = em_a;
	corres_a_ = pointsOnCurve;
	corres_b_ = intersection;
	idx_a = 0;
	startIdx = intersectStart;
	corres_num_ = interNum;
}

void UI::SetImages(Mat img_a, Mat img_b, float scale, Matrix3f& em_a, vector<Vector2f>& corres_a, vector<Vector2f>& corres_b)
{
	img_a_ = img_a;
	img_b_ = img_b;
	handles_b_.clear();
	handles_a_.clear();
	scale_ = scale;
	em_a_ = em_a;
	corres_a_ = corres_a;
	corres_b_ = corres_b;
}

void UI::MouseEventb(int event, int x, int y, int flags, void* param)
{
	switch (event)
	{
	case CV_EVENT_LBUTTONDOWN:
	{
		if (idx_a < corres_a_.size()) {
			circle(img_a_, Point((int)(corres_a_[idx_a][0] / scale_), (int)(corres_a_[idx_a][1] / scale_)), linethick, Scalar(0, 0, 255), -1, 8);
			//points_a_.push_back(Point((int)(corres_a_[pidx_][0] / scale_), (int)(corres_a_[pidx_][1] / scale_)));
			imshow("imga", img_a_);

			vector<Vector2f> temp2points;
			EndPointsOnLine(corres_a_[idx_a], temp2points, 'a');
			Vector2f start = temp2points[0], end = temp2points[1];
			Mat imgbtemp = img_b_.clone();
			line(imgbtemp, Point(start[0], start[1]), Point(end[0], end[1]), Scalar(0, 255, 0), linethick, 8, 0);
			imshow("imgb", imgbtemp);
			for (int i = 0; i < corres_num_[idx_a]; i++) {
				circle(imgbtemp, Point((int)(corres_b_[startIdx[idx_a] + i][0] / scale_), (int)(corres_b_[startIdx[idx_a] + i][1] / scale_)), linethick, Scalar(0, 0, 255), -1, 8);
				imshow("imgb", imgbtemp);
			}
			idx_a = idx_a + pstep;
		}
	}
	}
}

void UI::Show()
{
	cv::imshow("imga", img_a_);
	cv::imshow("imgb", img_b_);

	cv::setMouseCallback("imgb", MouseEventb);
	cv::waitKey(0);
}

void UI::OnlyShowEndPoints()
{
	/*cv::imshow("imga", img_a_);
	cv::imshow("imgb", img_b_);*/
	//cout << corres_a_.size() << endl;

	circle(img_a_, Point((int)(corres_a_.front()[0] / scale_), int(corres_a_.front()[1] / scale_)), 4, Scalar(0, 0, 255), -1, 8);
	vector<Vector2f> temp2points;
	EndPointsOnLine(corres_a_.front(), temp2points, 'a');
	Vector2f start = temp2points[0], end = temp2points[1];
	//Mat imgbtemp = img_b_.clone();
	line(img_b_, Point(start[0], start[1]), Point(end[0], end[1]), Scalar(0, 255, 0), linethick, 8, 0);

	circle(img_a_, Point((int)(corres_a_.back()[0] / scale_), int(corres_a_.back()[1] / scale_)), 4, Scalar(0, 0, 255), -1, 8);
	vector<Vector2f> temp2points2;
	EndPointsOnLine(corres_a_.back(), temp2points2, 'a');
	Vector2f start2 = temp2points2[0], end2 = temp2points2[1];
	line(img_b_, Point(start2[0], start2[1]), Point(end2[0], end2[1]), Scalar(0, 255, 0), linethick, 8, 0);
	imshow("imga", img_a_);
	imshow("imgb", img_b_);

	//cv::setMouseCallback("imgb", MouseEventb);


	cv::waitKey(0);
}

void UI::OnlyShowTangentPoints()
{
	/*cv::imshow("imga", img_a_);
	cv::imshow("imgb", img_b_);*/
	//cout << corres_a_.size() << endl;

	circle(img_a_, Point((int)(corres_a_.front()[0] / scale_), int(corres_a_.front()[1] / scale_)), 4, Scalar(0, 0, 255), -1, 8);
	vector<Vector2f> temp2points;
	EndPointsOnLine(corres_a_.front(), temp2points, 'a');
	Vector2f start = temp2points[0], end = temp2points[1];
	//Mat imgbtemp = img_b_.clone();
	line(img_b_, Point(start[0], start[1]), Point(end[0], end[1]), Scalar(0, 255, 0), linethick, 8, 0);

	circle(img_a_, Point((int)(corres_a_.back()[0] / scale_), int(corres_a_.back()[1] / scale_)), 4, Scalar(0, 0, 255), -1, 8);
	vector<Vector2f> temp2points2;
	EndPointsOnLine(corres_a_.back(), temp2points2, 'a');
	Vector2f start2 = temp2points2[0], end2 = temp2points2[1];
	line(img_b_, Point(start2[0], start2[1]), Point(end2[0], end2[1]), Scalar(0, 255, 0), linethick, 8, 0);
	imshow("imga", img_a_);
	imshow("imgb", img_b_);

	//cv::setMouseCallback("imgb", MouseEventb);


	cv::waitKey(0);
}

void UI::OnlyShow()
{
	cv::imshow("imga_again", img_a_);
	cv::imshow("imgb_again", img_b_);

	for (int i = 0; i < corres_a_.size(); i++) {
		circle(img_a_, Point((int)(corres_a_[i][0] / scale_), int(corres_a_[i][1] / scale_)), linethick, Scalar(0, 0, 255), -1, 8);
		imshow("imga_again", img_a_);
	}

	for (int i = 0; i < corres_b_.size(); i++) {
		circle(img_b_, Point((int)(corres_b_[i][0] / scale_), int(corres_b_[i][1] / scale_)), linethick, Scalar(0, 0, 255), -1, 8);
		imshow("imgb_again", img_b_);
	}

	cv::waitKey(0);
}

void UI::EndPointsOnLine(Vector2f& refpoint, vector<Vector2f>& endpoints, char which) {
	Matrix3f em;
	float max_y, max_x;
	if (which == 'a') {
		em = em_a_;
		max_y = img_a_.rows;
		max_x = img_a_.cols;
	}
	else {
		cout << "wrong type!" << endl;
	}
	Vector3f temp = Vector3f::Ones();
	temp.head(2) = refpoint;
	Vector3f para = em*temp;
	para[2] = para[2] / scale_;


	int count = 0;
	float value;
	Vector2f tempEP;
	value = -para[2] / para[1];
	if (!(value > max_y || value < 0.0 || count>1)) {
		tempEP << 0.0, value;
		endpoints.push_back(tempEP); count++;
	}
	value = -(para[2] + para[0] * max_x) / para[1];
	if (!(value > max_y || value < 0.0 || count>1)) {
		tempEP << max_x, value;
		endpoints.push_back(tempEP); count++;
	}
	value = -para[2] / para[0];
	if (!(value > max_x || value < 0.0 || count>1)) {
		tempEP << value, 0.0;
		endpoints.push_back(tempEP); count++;
	}
	value = -(para[2] + para[1] * max_y) / para[0];
	if (!(value > max_x || value < 0.0 || count>1)) {
		tempEP << value, max_y;
		endpoints.push_back(tempEP); count++;
	}
}

void UI::ShowProcess(int count, string path, int refIdxnum, int neiIdxnum, ProjectImage refView, ProjectImage neighborView, float scale) {
	float scaleinv = 1.0 / scale;
	int start = refView.reconGroupIdxStart[count];
	int pointsNum = refView.reconNumInGroup[count];
	cv::Mat imS_a, imS_b;

	string refIdx = std::to_string(refIdxnum);
	string neiIdx = std::to_string(neiIdxnum);
	string imga_name = path + "/view" + refIdx + "_bd" + imgIdx + ".png";
	string imgb_name = path + "/view" + neiIdx + "_bd" + imgIdx + ".png";

	/*string imga_name = path + "/view" + refIdx + ".jpg";
	string imgb_name = path + "/view" + neiIdx + ".jpg";*/

	Mat img_a = cv::imread(imga_name).clone();
	Mat img_b = cv::imread(imgb_name).clone();

	cv::resize(img_a, imS_a, Size(), scale, scale);
	cv::resize(img_b, imS_b, Size(), scale, scale);
	for (int i = 0; i < neighborView.tangentPoints.size(); i++) {
		cv::circle(imS_b, Point((int)(neighborView.tangentPoints[i][0] / scaleinv), (int)(neighborView.tangentPoints[i][1] / scaleinv)), 3, Scalar(255, 0, 255), -1, 8);
	}
	UI x;
	vector<Vector2f> tempRefP; tempRefP.assign(refView.img2DPoints.begin() + start, refView.img2DPoints.begin() + start + pointsNum);
	x.SetImages(imS_a, imS_b, scaleinv, neighborView.FMtoRefView, tempRefP, neighborView.intersect, neighborView.intersectIdxStart, neighborView.intersectNum);
	x.Show();
}

void UI::ShowLongMissingCurve(int iter, string path, vector<ProjectImage> views, vector<vector<Vector2i>> missingStartEndList, float scale) {
	float scaleinv = 1.0 / scale;
	vector<Scalar> color3 = { Scalar(0, 0, 255), Scalar(255, 0, 0), Scalar(0, 255, 0), Scalar(255, 255, 0) };
	vector<string> imgnames(views.size());
	vector<Mat> imgs(views.size());
	vector<Mat> imSs(views.size());
	for (int i = 0; i < views.size(); i++) {
		imgnames[i] = path + "/view" + std::to_string(i + 1) + "_bd" + imgIdx + ".png";
		imgs[i] = cv::imread(imgnames[i]).clone();
		cv::resize(imgs[i], imSs[i], Size(), scale, scale);
	}
	for (int i = 0; i < views.size(); i++) {
		for (int j = 0; j < missingStartEndList[i].size(); j++) {
			int start = missingStartEndList[i][j](0);
			int end = missingStartEndList[i][j](1);
			for (int k = start; k < end + 1; k++) {
				cv::circle(imSs[i], Point((int)(views[i].img2DPoints[k][0] / scaleinv), (int)(views[i].img2DPoints[k][1] / scaleinv)), 1, color3[j % 4], -1, 8);
			}
		}
		string imgname = "Long Missing Curves: iter" + std::to_string(iter + 1) + "_img" + std::to_string(i + 1);
		imshow(imgname, imSs[i]);
	}
	waitKey(0);
}

void UI::ShowFullMissingCurve(int iter, string path, vector<ProjectImage> views, vector<vector<int>> takenList, float scale) {
	float scaleinv = 1.0 / scale;

	vector<string> imgnames(views.size());
	vector<Mat> imgs(views.size());
	vector<Mat> imSs(views.size());
	for (int i = 0; i < views.size(); i++) {
		imgnames[i] = path + "/view" + std::to_string(i + 1) + "_bd" + imgIdx + ".png";
		imgs[i] = cv::imread(imgnames[i]).clone();
		cv::resize(imgs[i], imSs[i], Size(), scale, scale);
	}
	for (int i = 0; i < views.size(); i++) {
		for (int j = 0; j < views[i].img2DPoints.size(); j++) {
			if (takenList[i][j] == 0) cv::circle(imSs[i], Point((int)(views[i].img2DPoints[j][0] / scaleinv), (int)(views[i].img2DPoints[j][1] / scaleinv)), 1, Scalar(0, 0, 255), -1, 8);
		}
		string imgname = "Full Missing Curves: iter" + std::to_string(iter + 1) + "_img" + std::to_string(i + 1);
		imshow(imgname, imSs[i]);
	}
	waitKey(0);
}

void UI::ShowReconCurveEachIter(int iter, string path, vector<ProjectImage> views, vector<ReconCurve> reconCurveList, float scale) {
	float scaleinv = 1.0 / scale;
	vector<Scalar> color3 = { Scalar(0, 0, 255), Scalar(255, 0, 0), Scalar(0, 255, 0), Scalar(255, 255, 0) };
	vector<string> imgnames(views.size());
	vector<Mat> imgs(views.size());
	vector<Mat> imSs(views.size());
	for (int i = 0; i < views.size(); i++) {
		imgnames[i] = path + "/view" + std::to_string(i + 1) + "_bd" + imgIdx + ".png";
		imgs[i] = cv::imread(imgnames[i]).clone();
		cv::resize(imgs[i], imSs[i], Size(), scale, scale);
	}

	vector<vector<int>> takenList(views.size());
	for (int i = 0; i < views.size(); i++) {
		vector<int> taken(views[i].img2DPoints.size(), 0);
		takenList[i] = taken;
	}

	for (int k = 0; k < views.size(); k++) {
		for (int i = 0; i < reconCurveList.size(); i++) {
			vector<int> projIdx;

			for (int j = 0; j < reconCurveList[i].curve.size(); j++) {
				Vector2f nearest; int idx;
				float dist = views[k].FindClosestPointOnImgEachPoint(reconCurveList[i].curve[j], nearest, idx);
				takenList[k][idx] = 1;

			}

		}
	}

	for (int i = 0; i < views.size(); i++) {
		for (int j = 0; j < views[i].img2DPoints.size(); j++) {
			if (takenList[i][j] == 1) cv::circle(imSs[i], Point((int)(views[i].img2DPoints[j][0] / scaleinv), (int)(views[i].img2DPoints[j][1] / scaleinv)), 1, Scalar(255, 0, 0), -1, 8);
		}
		string imgname = "Full Recon Curves: iter" + std::to_string(iter + 1) + "_img" + std::to_string(i + 1);
		imshow(imgname, imSs[i]);
	}
	waitKey(0);
}

void UI::ShowTakenList(string path, int vidx, ProjectImage view, vector<int> gidx, vector<int> curTakenList, float scale) {
	float scaleinv = 1.0 / scale;

	string imgname = path + "/view" + std::to_string(vidx + 1) + "_bd" + imgIdx + ".png";
	Mat img = cv::imread(imgname).clone();
	Mat imS;
	cv::resize(img, imS, Size(), scale, scale);


	for (int j = 0; j < view.img2DPoints.size(); j++) {
		if (curTakenList[j] == 1) cv::circle(imS, Point((int)(view.img2DPoints[j][0] / scaleinv), (int)(view.img2DPoints[j][1] / scaleinv)), 1, Scalar(255, 0, 0), -1, 8);
	}
	for (int j = 0; j < gidx.size(); j++) {
		cv::circle(imS, Point((int)(view.img2DPoints[gidx[j]][0] / scaleinv), (int)(view.img2DPoints[gidx[j]][1] / scaleinv)), 1, Scalar(0, 0, 255), -1, 8);
	}
	string title = "Full Missing Curves: img" + std::to_string(vidx + 1);
	imshow(title, imS);

	waitKey(0);
}

//void UI::MouseEvent2(int event, int x, int y, int flags, void* param)
//{
//	switch (event)
//	{
//	case CV_EVENT_LBUTTONDOWN:
//	{
//		cv::resize(imgs[i], imSs[i], Size(), scale, scale);
//	}
//	}
//
//}

void UI::ZoomImage(string path, vector<ProjectImage> views, float scale) {
	float scaleinv = 1.0 / scale;
	vector<Scalar> color3 = { Scalar(0, 0, 255), Scalar(255, 0, 0), Scalar(0, 255, 0), Scalar(255, 255, 0) };
	vector<string> imgnames(views.size());
	vector<Mat> imgs(views.size());
	vector<Mat> imSs(views.size());
	for (int i = 0; i < views.size(); i++) {
		imgnames[i] = path + "/view" + std::to_string(i + 1) + "_bd" + imgIdx + ".png";
		imgs[i] = cv::imread(imgnames[i]).clone();
		cv::resize(imgs[i], imSs[i], Size(), scale, scale);
	}
	/*setMouseCallback("My Window", MouseEvent2, NULL);
	waitKey(0);*/
	Mat src, dst, tmp, tempdst;
	src = imSs[0];
	tmp = src;
	dst = tmp;
	tempdst = dst;
	while (true)
	{
		int c;
		c = waitKey(10);

		if ((char)c == 27)
		{
			break;
		}
		if ((char)c == 'u')
		{
			cv::resize(tmp, tempdst, Size(), 2.0, 2.0);
			/*cout << tmp.size().width << tmp.size().height;
			cout << tempdst.size().width << tempdst.size().height;*/
			Rect roi = Rect(tmp.size().width / 4, tmp.size().height / 4, tmp.size().width / 4 + tmp.size().width, tmp.size().height / 4 + tmp.size().height);
			dst = tempdst(roi);
			printf("** Zoom In: Image x 2 \n");
		}
		else if ((char)c == 'd')
		{
			cv::resize(tmp, dst, Size(), 0.5, 0.5);
			//pyrDown(tmp, dst, Size(tmp.cols / 2, tmp.rows / 2));
			printf("** Zoom Out: Image / 2 \n");
		}

		imshow("Zoom", dst);
		tmp = dst;
	}
	cv::waitKey(0);

	/*for (int i = 0; i < views.size(); i++) {
		for (int j = 0; j < views[i].img2DPoints.size(); j++) {
			if (takenList[i][j] == 1) cv::circle(imSs[i], Point((int)(views[i].img2DPoints[j][0] / scaleinv), (int)(views[i].img2DPoints[j][1] / scaleinv)), 1, Scalar(255, 0, 0), -1, 8);
		}
		string imgname = "Full Recon Curves: iter" + std::to_string(iter + 1) + "_img" + std::to_string(i + 1);
		imshow(imgname, imSs[i]);
	}
	waitKey(0);*/
}

void UI::ShowTangentOnPoints(string path, vector<ProjectImage> views, float scale) {
	float scaleinv = 1.0 / scale;
	vector<Scalar> color3 = { Scalar(0, 0, 255), Scalar(255, 0, 0), Scalar(0, 255, 0), Scalar(255, 255, 0) };
	vector<string> imgnames(views.size());
	vector<Mat> imgs(views.size());
	vector<Mat> imSs(views.size());
	for (int i = 0; i < views.size(); i++) {
		imgnames[i] = path + "/view" + std::to_string(i + 1) + "_bd" + imgIdx + ".png";
		imgs[i] = cv::imread(imgnames[i]).clone();
		cv::resize(imgs[i], imSs[i], Size(), scale, scale);
	}

	for (int i = 0; i < views.size(); i++) {
		for (int j = 0; j < views[i].img2DPoints.size(); j++) {
			Vector2f start = views[i].img2DPoints[j];
			Vector2f end = start + views[i].img2DDires[j] * 10;
			//circle(imSs[i], Point((int)(start[0] / scaleinv), (int)(start[1] / scaleinv)), 1, Scalar(0, 0, 255), -1, 8);
			line(imSs[i], Point(start[0] / scaleinv, start[1] / scaleinv), Point(end[0] / scaleinv, end[1] / scaleinv), Scalar(0, 255, 0), 1, 8, 0);
		}
		string imgname = "Img" + std::to_string(i + 1);
		imshow(imgname, imSs[i]);
	}
	waitKey(0);
}

void UI::OnlyShowTangentPoints(string path, vector<ProjectImage> views, float scale) {
	float scaleinv = 1.0 / scale;
	vector<Scalar> color3 = { Scalar(0, 0, 255), Scalar(255, 0, 0), Scalar(0, 255, 0), Scalar(255, 255, 0) };
	vector<string> imgnames(views.size());
	vector<Mat> imgs(views.size());
	vector<Mat> imSs(views.size());
	for (int i = 0; i < views.size(); i++) {
		imgnames[i] = path + "/view" + std::to_string(i + 1) + "_bd" + imgIdx + ".png";
		imgs[i] = cv::imread(imgnames[i]).clone();
		cv::resize(imgs[i], imSs[i], Size(), scale, scale);
		for (int j = 0; j < views[i].tangentPoints.size(); j++) {
			cv::circle(imSs[i], Point((int)(views[i].tangentPoints[j][0] / scaleinv), (int)(views[i].tangentPoints[j][1] / scaleinv)), 3, Scalar(255, 0, 255), -1, 8);
		}
		string imgname = "Tangent Points: img" + std::to_string(i + 1);
		imshow(imgname, imSs[i]);
	}
	waitKey(0);
}

void UI::OnlyShowTangentPoints2(string path, ProjectImage view, int vidx, int offstep, float scale) {
	float scaleinv = 1.0 / scale;
	string imgname = path + "/view" + std::to_string(vidx) + "_bd" + imgIdx + ".png";
	Mat img;
	Mat imS;

	vector<Scalar> color7 = { Scalar(0, 0, 255), Scalar(0, 255, 0), Scalar(255, 0, 0), Scalar(0, 255, 255), Scalar(255, 255, 0), Scalar(255, 0, 255), Scalar(0, 0, 0) };

	img = cv::imread(imgname).clone();
	cv::resize(img, imS, Size(), scale, scale);
	for (int j = 0; j < view.tangentPoints.size(); j++) {
		cv::circle(imS, Point((int)(view.tangentPoints[j][0] / scaleinv), (int)(view.tangentPoints[j][1] / scaleinv)), 2, color7[j%3], -1, 8);
	}
	//for (int j = 0; j < view.tangentPoints.size(); j++) {
	//	cv::circle(imS, Point((int)(view.img2DPoints[view.tangentLeftRight[j][0]][0] / scaleinv), (int)(view.img2DPoints[view.tangentLeftRight[j][0]][1] / scaleinv)), 2, Scalar(255, 0, 0), -1, 8);
	//	cv::circle(imS, Point((int)(view.img2DPoints[view.tangentLeftRight[j][1]][0] / scaleinv), (int)(view.img2DPoints[view.tangentLeftRight[j][1]][1] / scaleinv)), 2, Scalar(255, 0, 0), -1, 8);

	//	/*cv::circle(imS, Point((int)(view.img2DPoints[view.tangentIdx[j]-offstep][0] / scaleinv), (int)(view.img2DPoints[view.tangentIdx[j] - offstep][1] / scaleinv)), 3, Scalar(0, 0, 255), -1, 8);
	//	cv::circle(imS, Point((int)(view.img2DPoints[view.tangentIdx[j] + offstep][0] / scaleinv), (int)(view.img2DPoints[view.tangentIdx[j] + offstep][1] / scaleinv)), 3, Scalar(0, 0, 255), -1, 8);*/
	//}
	string title = "Tangent Points: img" + std::to_string(vidx);
	imshow(title, imS);

	waitKey(0);
}

