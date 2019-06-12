#include"union_smooth.h"
const double pi = 3.1415926;
const double e = 2.718282;

using namespace std;
using namespace cv;
using namespace Bundled;

union_smooth::union_smooth(string open_org, string save_dest, string save_cut_dst, set_of_path org_paths, std::vector<std::vector<int>>index, int video_h, int video_w, double qw, double qh, double a, double b, double c, double d, double e)
{
	this->index = index;
	this->h = video_h;
	this->w = video_w;
	this->qh = qh;
	this->qw = qw;
	this->org_paths = org_paths;
	this->a = a;
	this->b = b;
	this->c = c;
	this->d = d;
	this->e = e;
	this->open_org_video_path = open_org;
	this->save_dest_video_path = save_dest;
	this->save_cut_dst_video_path = save_cut_dst;

	source = new Bundled::Mesh(h, w, qw, qh);
	destin = new Bundled::Mesh(h, w, qw, qh);
	mesh_height = source->GetMeshHeight();
	mesh_width = source->GetMeshWidth();
}




void union_smooth::TransformByHomo(cv::Mat & mat, double x, double y, double & dx, double & dy)
{
	double * Mi0 = mat.ptr<double>(0);
	double * Mi1 = mat.ptr<double>(1);
	double * Mi2 = mat.ptr<double>(2);


	dx = Mi0[0] * x + Mi0[1] * y + Mi0[2];
	dy = Mi1[0] * x + Mi1[1] * y + Mi1[2];
	double dz = Mi2[0] * x + Mi2[1] * y + Mi2[2];

	dx = dx / dz;
	dy = dy / dz;
}
void union_smooth::set_control_points(const vector<std::vector<v2d>> frame_points)
{

	this->frame_points = frame_points;
	frame_num = frame_points.size();
	num_of_feature_points = 0;
	for (int i = 0; i < frame_points.size(); i++)
	{
		num_of_feature_points += frame_points[i].size();
	}


	int frame_num = frame_points.size();
	vector<double> V00;
	vector<double> V01;
	vector<double> V10;
	vector<double> V11;
	TensorInt           histogram;
	TensorDou           histogram2;
	int                 maxHistogram;
	histogram.set(mesh_width - 1, mesh_height - 1);//一个Int类型的矩阵
	histogram2.set(mesh_width, mesh_height);
	histogram.fill(0);
	histogram2.fill(0);

	//为每一帧设置控制点
	for (int i = 0; i < frame_num; i++)
	{

		int len = frame_points[i].size();
		for (int j = 0; j < len; j++)
		{
			double x, y;
			x = frame_points[i][j][0];
			y = frame_points[i][j][1];
			v2d point = _v2d_(x, y);

			int mesh_i = floor(point[1] / qh);
			int mesh_j = floor(point[0] / qw);



			if (mesh_i >= mesh_height - 1)
			{
				mesh_i -= 1;
			}

			if (mesh_j >= mesh_width - 1)
			{
				mesh_j -= 1;
			}

			Bundled::Quad qd = source->getQuad(mesh_i, mesh_j);



			double coe[4];

			qd.GetBilinearCoordinates(point, coe);//通过双线性插值获得四个顶点的坐标

			V00.push_back(coe[0]);
			V01.push_back(coe[1]);
			V10.push_back(coe[2]);
			V11.push_back(coe[3]);

			histogram.at(mesh_i, mesh_j) = histogram.at(mesh_i, mesh_j) + 1;//这里是什么意思，对这个矩阵的每一个坐标值加一
		}

		dataterm_element_V00.push_back(V00);
		dataterm_element_V01.push_back(V01);
		dataterm_element_V10.push_back(V10);
		dataterm_element_V11.push_back(V11);


		V00.clear();
		V01.clear();
		V10.clear();
		V11.clear();


		for (int i = 0; i < mesh_height; i++)
		{
			for (int j = 0; j < mesh_width; j++)
			{
				double fT = 0;
				double fD = 0;

				if (i - 1 >= 0 && j < mesh_width - 1) //fT是一个2*2矩阵的和
				{
					fT += histogram.at(i - 1, j);
					fD++;
				}

				if (j - 1 >= 0 && i < mesh_height - 1)
				{
					fT += histogram.at(i, j - 1);
					fD++;
				}

				if (i - 1 >= 0 && j - 1 >= 0)
				{
					fT += histogram.at(i - 1, j - 1);
					fD++;
				}

				if (i < mesh_height - 1 && j < mesh_width - 1)
				{
					fT += histogram.at(i, j);
					fD++;
				}

				histogram2.at(i, j) = fT / fD;//histogram2是相邻四个坐标值的平均数
			}

		}

		maxHistogram = histogram.at(0, 0);

		for (int i = 0; i < mesh_height - 1; i++)
		{
			for (int j = 0; j < mesh_width - 1; j++)
			{
				if (histogram.at(i, j) > maxHistogram)
					maxHistogram = histogram.at(i, j);
			}
		}

		each_frame_histogram.push_back(histogram);
		each_frame_histogram2.push_back(histogram2);
		each_frame_maxHistogram.push_back(maxHistogram);


	}

	histogram.erase();
	histogram2.erase();

	V00.clear();
	V01.clear();
	V10.clear();
	V11.clear();
	vector<double>(V00).swap(V00);
	vector<double>(V01).swap(V01);
	vector<double>(V10).swap(V10);
	vector<double>(V11).swap(V11);
}

double union_smooth::GetWeight(int frame, int i, int j, double c)
{
	double nn = each_frame_histogram[frame].at(i, j);
	double sigma = each_frame_maxHistogram[frame] * 0.5;
	double alpha = c*exp(-0.5*nn*nn / (sigma*sigma));

	return alpha;
}
void union_smooth::GetSmoothWeight(v2d V1, v2d V2, v2d V3, double uv[2])
{
	v2d V21 = _v2d_(V1[0] - V2[0], V1[1] - V2[1]);
	v2d V23 = _v2d_(V3[0] - V2[0], V3[1] - V2[1]);

	double d1 = sqrt(V21[0] * V21[0] + V21[1] * V21[1]);//
	double d3 = sqrt(V23[0] * V23[0] + V23[1] * V23[1]);

	double cosin = V21[0] * V23[0] + V21[1] * V23[1];
	cosin = cosin / (d1*d3);

	double u_dis = cosin * d1;
	uv[0] = u_dis / d3;

	double v_dis = sqrt(d1 * d1 - u_dis * u_dis);
	uv[1] = v_dis / d3;
}

void union_smooth::solve_by_Eigen(bool ratation)
{
	//Firstly,establishing the optimization equation.
	//There are five items in total.
	//items                path_data           path_smooth            mesh_data                 mesh_smooth                                       mesh_time_smooth
	//parameter            All feature points  All feature points     All feature points        Number of mesh vertices*number of frames          Number of mesh vertices*number of frames

	num_of_path_data_term = num_of_feature_points * 2;
	num_of_path_smooth_term_2 = (num_of_feature_points - org_paths.size()) * 2;

	num_of_path_smooth_term = num_of_feature_points * 2 - org_paths.size() * 2 * 2;
	num_of_mesh_data_term = num_of_feature_points * 2;
	num_of_mesh_smooth_term = (mesh_height - 1)*(mesh_width - 1)*frame_num * 2 * 2;
	num_of_mesh_path_smooth_term = (frame_num - 2)*mesh_height*mesh_width * 2;
	//num_of_mesh_time_smooth_term = mesh_height*mesh_width*frame_num*2;

	//A*x=B
	Bundled::ArrayDou B;
	B.set(num_of_path_data_term + num_of_path_smooth_term_2 + num_of_path_smooth_term + num_of_mesh_data_term + num_of_mesh_smooth_term + num_of_mesh_path_smooth_term);
	B.fill(0);
	row_count = 0;


	x_index.set(num_of_feature_points + mesh_width*mesh_height*frame_num);
	y_index.set(num_of_feature_points + mesh_width*mesh_height*frame_num);

	x_index.fill(0);
	y_index.fill(0);

	for (int i = 0; i < num_of_feature_points; i++)
	{
		x_index[i] = i;
		y_index[i] = i + num_of_feature_points;
	}
	int j = 0;
	for (int i = num_of_feature_points; i <num_of_feature_points + mesh_height*mesh_width*frame_num; i++)
	{
		x_index[i] = num_of_feature_points * 2 + j;
		y_index[i] = num_of_feature_points * 2 + mesh_height*mesh_width*frame_num + j;
		j++;
	}


	//path_data_term 
	Creat_path_data_term(B);
	cout << row_count << endl;
	cout << num_of_mesh_data_term << endl << endl;

	//path_smooth_term a

	Creat_path_smooth_term(B);
	cout << row_count << endl;
	cout << num_of_mesh_data_term + num_of_path_smooth_term << endl << endl;

	//path_smooth_term2 b

	Creat_path_smooth_term2(B, ratation);
	cout << row_count << endl;
	cout << num_of_mesh_data_term + num_of_path_smooth_term + num_of_path_smooth_term_2 << endl << endl;

	//mesh_data_term c

	Creat_mash_data_term(B);
	cout << row_count << endl;
	cout << num_of_mesh_data_term + num_of_path_smooth_term + num_of_path_smooth_term_2 + num_of_mesh_data_term << endl << endl;

	//mesh_shap_smooth_term d

	Creat_mesh_smooth_term(B);
	cout << row_count << endl;
	cout << num_of_mesh_data_term + num_of_path_smooth_term + num_of_path_smooth_term_2 + num_of_mesh_data_term + num_of_mesh_smooth_term << endl << endl;

	//mesh_path_smooth_term e
	Creat_mesh_path_smooth_term(B);
	cout << row_count << endl;
	cout << num_of_mesh_data_term + num_of_path_smooth_term + num_of_path_smooth_term_2 + num_of_mesh_data_term + num_of_mesh_smooth_term + num_of_mesh_path_smooth_term << endl << endl;






	//jacobi 
	int j_rows;
	int j_cols;

	//number of term
	j_rows = num_of_path_data_term +
		num_of_path_smooth_term_2 +
		num_of_path_smooth_term +
		num_of_mesh_data_term +
		num_of_mesh_smooth_term +
		num_of_mesh_path_smooth_term;
	//number of variable
	j_cols = num_of_feature_points * 2 + frame_num*mesh_width*mesh_height * 2;


	typedef Eigen::Triplet<double> T;
	vector<T> tripletlist;

	Eigen::SparseMatrix<double> matj_sparse(j_rows, j_cols);
	//Eigen::MatrixXd matJ(j_rows,j_cols);
	int num_of_term_check = 0;

	Eigen::VectorXd bJ(j_rows);

	int NN = (int)Path_data_term.size();
	num_of_term_check += NN;
	for (int i = 0; i < NN; i++)
	{
		int row = Path_data_term[i][0];
		int col = Path_data_term[i][1];
		double v = Path_data_term[i][2];

		tripletlist.push_back(T(row, col, v));
	}

	NN = (int)Path_smooth_term.size();
	num_of_term_check += NN;
	for (int i = 0; i<NN; i++)
	{
		int row = Path_smooth_term[i][0];
		int col = Path_smooth_term[i][1];
		double v = Path_smooth_term[i][2];

		tripletlist.push_back(T(row, col, v));
	}


	NN = (int)Path_smooth_term2.size();
	num_of_term_check += NN;
	for (int i = 0; i < NN; i++)
	{
		int row = Path_smooth_term2[i][0];
		int col = Path_smooth_term2[i][1];
		double v = Path_smooth_term2[i][2];

		tripletlist.push_back(T(row, col, v));
	}


	NN = (int)Mesh_data_term.size();
	num_of_term_check += NN;
	for (int i = 0; i < NN; i++)
	{
		int row = Mesh_data_term[i][0];
		int col = Mesh_data_term[i][1];
		double v = Mesh_data_term[i][2];


		tripletlist.push_back(T(row, col, v));
	}

	NN = (int)Mesh_smooth_term.size();
	num_of_term_check += NN;
	for (int i = 0; i < NN; i++)
	{
		int row = Mesh_smooth_term[i][0];
		int col = Mesh_smooth_term[i][1];
		double v = Mesh_smooth_term[i][2];


		tripletlist.push_back(T(row, col, v));
	}

	NN = (int)Mesh_path_smooth_term.size();
	num_of_term_check += NN;
	for (int i = 0; i < NN; i++)
	{
		int row = Mesh_path_smooth_term[i][0];
		int col = Mesh_path_smooth_term[i][1];
		double v = Mesh_path_smooth_term[i][2];

		tripletlist.push_back(T(row, col, v));
	}

	for (int i = 0; i<j_rows; i++)
	{
		bJ(i) = B[i];

	}

	//sparse mat

	//Eigen::SparseMatrix<double> matj_sparse;
	matj_sparse.setFromTriplets(tripletlist.begin(), tripletlist.end());

	Eigen::VectorXd x_sparse(j_cols);

	Eigen::SparseMatrix<double> mat_t = matj_sparse.transpose()*matj_sparse;
	Eigen::VectorXd b_t = matj_sparse.transpose()*bJ;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	solver.compute(mat_t);
	x_sparse = solver.solve(b_t);





	set_of_path dst_path;
	dst_path = org_paths;

	int count = 0;
	for (int i = 0; i < dst_path.size(); i++)
	{
		for (int j = 0; j < dst_path[i].PointPath.size(); j++)
		{
			dst_path[i].PointPath[j][0] = x_sparse[count];
			dst_path[i].PointPath[j][1] = x_sparse[count + num_of_feature_points];
			count++;
		}
	}



	tripletlist.clear();
	dst_path.clear();
	vector<T>(tripletlist).swap(tripletlist);
	set_of_path(dst_path).swap(dst_path);


	cout << endl;

	cv::Mat Point_imge = cv::Mat(1000, 1000, CV_8UC3, cv::Scalar(0, 0, 0));

	//test code
	/*int index;
	for (int i = 0; i < org_paths.size(); i++)
	{
	if (org_paths[i].start_frame == 1&&org_paths[i].end_frame==112)
	index = i;
	}
	int index_path=0;
	for (int i = 0; i < index; i++)
	{
	index_path += org_paths[i].PointPath.size();
	}

	for (int i = index_path; i < index_path+org_paths[index].PointPath.size(); i++)
	{

	cv::circle(Point_imge, Point2d(x_sparse[i], x_sparse[i + num_of_feature_points]), 1, cv::Scalar(255, 0, 0));
	}


	for (int i = 0; i < org_paths[index].PointPath.size(); i++)
	{
	cout << org_paths[i].PointPath[index][0] << " " << org_paths[i].PointPath[index][1] << endl;
	cv::circle(Point_imge, Point2d(org_paths[index].PointPath[i][0]+500, org_paths[index].PointPath[i][1]), 1, cv::Scalar(0, 255, 0));
	}

	namedWindow("Point");
	imshow("Point", Point_imge);
	waitKey(0);
	*/

	/*save mesh*/
	//system("cls");
	double xx, yy;


	for (int t = 0; t < frame_num; t++)
	{
		for (int i = 0; i < mesh_height; i++)
		{
			for (int j = 0; j < mesh_width; j++)
			{
				xx = x_sparse[num_of_feature_points * 2 + t*mesh_height*mesh_width + i*mesh_width + j];
				yy = x_sparse[num_of_feature_points * 2 + frame_num*mesh_height*mesh_width + t*mesh_height*mesh_width + i*mesh_width + j];

				destin->SetVertex(i, j, _v2d_(xx, yy));


			}
		}

		each_frame_destin.push_back(*destin);

	}

	//Release memory
	matj_sparse.resize(0, 0);
	mat_t.resize(0, 0);
	b_t.resize(0);
	x_sparse.resize(0);
	bJ.resize(0);
	x_index.erase();
	y_index.erase();
	B.erase();

	//show the mesh
	//for (int t = 0; t < frame_num; t++)
	//{
	//	cv::Mat big_image = cv::Mat(1500, 700, CV_8UC3, cv::Scalar(0, 0, 0));
	//	for (int row = 0; row < each_frame_destin[t].GetMeshHeight(); row++)
	//	{
	//		for (int col = 0; col < each_frame_destin[t].GetMeshWidth(); col++)
	//		{
	//			cv::circle(big_image, Point2d(each_frame_destin[t].GetVertex(row, col)[0], each_frame_destin[t].GetVertex(row, col)[1]), 2, cv::Scalar(255, 0, 0));
	//		}
	//	}

	//	namedWindow("text");
	//	imshow("text", big_image);
	//	waitKey(20);
	//}


}

void union_smooth::Creat_path_data_term(ArrayDou& B)
{
	int len = num_of_feature_points;
	int term_count = 0;
	int col_count = 0;

	Path_data_term.erase();
	Path_data_term.set(num_of_feature_points * 2);
	Path_data_term.fill(_v3d_(0, 0, 0));

	for (int i = 0; i < org_paths.size(); i++)
	{
		for (int j = 0; j < org_paths[i].PointPath.size(); j++)
		{

			Path_data_term[term_count][0] = row_count;
			Path_data_term[term_count][1] = x_index(col_count);
			Path_data_term[term_count][2] = 1;
			term_count++;
			B[row_count++] = org_paths[i].PointPath[j][0];

			Path_data_term[term_count][0] = row_count;
			Path_data_term[term_count][1] = y_index(col_count);
			Path_data_term[term_count][2] = 1;
			term_count++;
			B[row_count++] = org_paths[i].PointPath[j][1];

			col_count++;
		}
	}
}

void union_smooth::Creat_path_smooth_term2(ArrayDou& B,bool rotation)
{
	int term_count = 0;
	int col_count = 0;

	Path_smooth_term2.erase();
	Path_smooth_term2.set((num_of_feature_points - org_paths.size()) * 2 * 2);
	Path_smooth_term2.fill(_v3d_(0, 0, 0));

	for (int i = 0; i < org_paths.size(); i++)
	{
		for (int j = 0; j < org_paths[i].PointPath.size(); j++)
		{
			if (j == 0)
			{
				col_count++;
				continue;
			}

			int x_balance_parameter = 1;
			int y_balance_parameter = 1;
			//computer the x and y distance betwen two neighbor feature points
			double x_distance = abs(org_paths[i].PointPath[j][0] - org_paths[i].PointPath[j-1][0]);
			double y_distance = abs(org_paths[i].PointPath[j][1] - org_paths[i].PointPath[j-1][1]);
			//when the distance of neirghbor feature points more than 5 pixel, it must be adjust.

			if (rotation)
			{
				x_distance > 15 ? x_balance_parameter = b : x_balance_parameter = 1;
				y_distance > 15 ? y_balance_parameter = b : y_balance_parameter = 1;
			}
			else
			{
				x_balance_parameter = 1;
				y_balance_parameter = 1;
			}

			//x
			Path_smooth_term2[term_count][0] = row_count;
			Path_smooth_term2[term_count][1] = x_index(col_count);
			Path_smooth_term2[term_count][2] = 1 * b/ x_balance_parameter;
			term_count++;

			Path_smooth_term2[term_count][0] = row_count;
			Path_smooth_term2[term_count][1] = x_index(col_count - 1);
			Path_smooth_term2[term_count][2] = -1 * b / x_balance_parameter;
			term_count++;
			row_count++;

			//y
			Path_smooth_term2[term_count][0] = row_count;
			Path_smooth_term2[term_count][1] = y_index(col_count);
			Path_smooth_term2[term_count][2] = 1 * b / y_balance_parameter;
			term_count++;

			Path_smooth_term2[term_count][0] = row_count;
			Path_smooth_term2[term_count][1] = y_index(col_count - 1);
			Path_smooth_term2[term_count][2] = -1 * b / y_balance_parameter;
			term_count++;
			row_count++;

			col_count++;

		}

	}

}

void union_smooth::Creat_path_smooth_term(ArrayDou& B)
{
	int term_count = 0;
	int col_count = 0;
	Path_smooth_term.erase();
	Path_smooth_term.set(num_of_feature_points * 2 * 3 - org_paths.size() * 2 * 2 * 3);
	Path_smooth_term.fill(_v3d_(0, 0, 0));


	for (int i = 0; i < org_paths.size(); i++)
	{
		for (int j = 0; j < org_paths[i].PointPath.size(); j++)
		{
			if (j == 0)
			{
				col_count++;

				continue;
			}

			if (j == org_paths[i].PointPath.size() - 1)
			{
				col_count++;

				continue;
			}



			else
			{

				//x
				Path_smooth_term[term_count][0] = row_count;
				Path_smooth_term[term_count][1] = x_index(col_count - 1);
				Path_smooth_term[term_count][2] = -1 * a;
				term_count++;


				Path_smooth_term[term_count][0] = row_count;
				Path_smooth_term[term_count][1] = x_index(col_count);
				Path_smooth_term[term_count][2] = 2 * a;
				term_count++;

				Path_smooth_term[term_count][0] = row_count;
				Path_smooth_term[term_count][1] = x_index(col_count + 1);
				Path_smooth_term[term_count][2] = -1 * a;
				term_count++;

				row_count++;

				//y
				Path_smooth_term[term_count][0] = row_count;
				Path_smooth_term[term_count][1] = y_index(col_count - 1);
				Path_smooth_term[term_count][2] = -1 * a;
				term_count++;

				Path_smooth_term[term_count][0] = row_count;
				Path_smooth_term[term_count][1] = y_index(col_count);
				Path_smooth_term[term_count][2] = 2 * a;
				term_count++;

				Path_smooth_term[term_count][0] = row_count;
				Path_smooth_term[term_count][1] = y_index(col_count + 1);
				Path_smooth_term[term_count][2] = -1 * a;
				term_count++;




				row_count++;
				col_count++;



			}

		}
	}
}

void union_smooth::Creat_mash_data_term(ArrayDou& B)
{
	Mesh_data_term.erase();
	Mesh_data_term.set(num_of_feature_points * 2 * 5);
	Mesh_data_term.fill(_v3d_(0, 0, 0));
	int term_count = 0;



	for (int i = 0; i < frame_points.size(); i++)
	{
		for (int j = 0; j < frame_points[i].size(); j++)
		{
			double w00 = dataterm_element_V00[i][j];
			double w01 = dataterm_element_V01[i][j];
			double w10 = dataterm_element_V10[i][j];
			double w11 = dataterm_element_V11[i][j];

			int mesh_i, mesh_j;

			mesh_i = floor(frame_points[i][j][1] / qh);
			mesh_j = floor(frame_points[i][j][0] / qw);

			if (mesh_i >= mesh_height - 1)
			{
				mesh_i -= 1;
			}

			if (mesh_j >= mesh_width - 1)
			{
				mesh_j -= 1;
			}



			//x
			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = x_index(num_of_feature_points + i*mesh_height*mesh_width + mesh_i*mesh_width + mesh_j);
			Mesh_data_term[term_count][2] = w00*c;
			term_count++;


			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = x_index(num_of_feature_points + i*mesh_height*mesh_width + mesh_i*mesh_width + mesh_j + 1);
			Mesh_data_term[term_count][2] = w01*c;
			term_count++;


			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = x_index(num_of_feature_points + i*mesh_height*mesh_width + (mesh_i + 1)*mesh_width + mesh_j);
			Mesh_data_term[term_count][2] = w10*c;
			term_count++;


			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = x_index(num_of_feature_points + i*mesh_height*mesh_width + (mesh_i + 1)*mesh_width + (mesh_j + 1));
			Mesh_data_term[term_count][2] = w11*c;
			term_count++;


			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = x_index(index[i][j]);
			Mesh_data_term[term_count][2] = -1 * c;
			term_count++;

			row_count++;

			//y
			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = y_index(num_of_feature_points + i*mesh_height*mesh_width + mesh_i*mesh_width + mesh_j);
			Mesh_data_term[term_count][2] = w00*c;
			term_count++;

			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = y_index(num_of_feature_points + i*mesh_height*mesh_width + mesh_i*mesh_width + mesh_j + 1);
			Mesh_data_term[term_count][2] = w01*c;
			term_count++;


			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = y_index(num_of_feature_points + i*mesh_height*mesh_width + (mesh_i + 1)*mesh_width + mesh_j);
			Mesh_data_term[term_count][2] = w10*c;
			term_count++;


			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = y_index(num_of_feature_points + i*mesh_height*mesh_width + (mesh_i + 1)*mesh_width + mesh_j + 1);
			Mesh_data_term[term_count][2] = w11*c;
			term_count++;


			Mesh_data_term[term_count][0] = row_count;
			Mesh_data_term[term_count][1] = y_index(index[i][j]);
			Mesh_data_term[term_count][2] = -1 * c;
			term_count++;

			row_count++;

		}

	}
}

void union_smooth::Creat_mesh_smooth_term(ArrayDou& b)
{
	Mesh_smooth_term.erase();
	Mesh_smooth_term.set(frame_num*(mesh_width - 1)*(mesh_height - 1) * 2 * 2 * 5);
	Mesh_smooth_term.fill(_v3d_(0, 0, 0));
	int term_count = 0;

	for (int t = 0; t < frame_num; t++)
	{
		for (int i = 0; i < mesh_height - 1; i++)
		{
			for (int j = 0; j < mesh_width - 1; j++)
			{
				double aveg_c = GetWeight(t, i, j, d);

				//if (weight_of_all_vertexs_similarity[t][i][j] < 0.5)
				//	aveg_c *= 0.5;
				//else
				//	aveg_c = aveg_c*weight_of_all_vertexs_similarity[t][i][j];

				/*cout << weight_of_all_vertexs_similarity[t][i][j] << endl;*/

				SetUpTriangle(t, i, j, aveg_c, term_count);

				SetDownTriangle(t, i, j, aveg_c, term_count);


			}
		}
	}


}

void union_smooth::SetUpTriangle(int frame, int i, int j, double aveg_c, int& term_count)
{
	//v1-----v2
	//	     |
	//	     |
	//	    v3

	v2d V1 = source->GetVertex(i, j);
	v2d V2 = source->GetVertex(i, j + 1);
	v2d V3 = source->GetVertex(i + 1, j + 1);

	//int weight_v1 = aveg_c*weight_of_all_vertexs[frame][i][j];
	//int weight_v2 = aveg_c*weight_of_all_vertexs[frame][i][j + 1];
	//int weight_v3 = aveg_c*weight_of_all_vertexs[frame][i + 1][j + 1];


	double uv[2];

	GetSmoothWeight(V1, V2, V3, uv);
	double u = uv[0];
	double v = uv[1];

	//cout << 1 - u << " " << u << " " << v << " " << endl;
	//x
	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + i*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = (1.0 - u)*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = u*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + i*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = v*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = (-1.0 * v)*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + i*mesh_width + j];
	Mesh_smooth_term[term_count][2] = (-1.0)*aveg_c;
	term_count++;

	row_count++;

	//y
	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + i*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = (1.0 - u)*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = u*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = v*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + i*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = (-1.0 * v)*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + i*mesh_width + j];
	Mesh_smooth_term[term_count][2] = (-1.0)*aveg_c;
	term_count++;

	row_count++;
}

void union_smooth::SetDownTriangle(int frame, int i, int j, double aveg_c, int& term_count)
{
	/**
	*  V1
	*  |
	*  |
	*  V2____V3
	*/

	v2d V1 = source->GetVertex(i, j);
	v2d V2 = source->GetVertex(i + 1, j);
	v2d V3 = source->GetVertex(i + 1, j + 1);

	//int weight_v1 = aveg_c*weight_of_all_vertexs[frame][i][j];
	//int weight_v2 = aveg_c*weight_of_all_vertexs[frame][i + 1][j];
	//int weight_v3 = aveg_c*weight_of_all_vertexs[frame][i + 1][j + 1];

	double uv[2];
	GetSmoothWeight(V1, V2, V3, uv);
	double u = uv[0];
	double v = uv[1];

	//cout << 1 - u << " " << u << " " << v << " " << endl;
	//x
	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j];
	Mesh_smooth_term[term_count][2] = (1.0 - u)*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = u*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = v*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j];
	Mesh_smooth_term[term_count][2] = (-1.0*v)*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + i*mesh_width + j];
	Mesh_smooth_term[term_count][2] = (-1.0)*aveg_c;
	term_count++;

	row_count++;

	//y

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j];
	Mesh_smooth_term[term_count][2] = (1.0 - u)*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = u*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j];
	Mesh_smooth_term[term_count][2] = v*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = x_index[num_of_feature_points + frame*mesh_height*mesh_width + (i + 1)*mesh_width + j + 1];
	Mesh_smooth_term[term_count][2] = (-1.0*v)*aveg_c;
	term_count++;

	Mesh_smooth_term[term_count][0] = row_count;
	Mesh_smooth_term[term_count][1] = y_index[num_of_feature_points + frame*mesh_height*mesh_width + i*mesh_width + j];
	Mesh_smooth_term[term_count][2] = (-1.0)*aveg_c;
	term_count++;

	row_count++;


}

void union_smooth::Creat_mesh_path_smooth_term(ArrayDou& b)
{
	int term_count = 0;
	Mesh_path_smooth_term.erase();
	Mesh_path_smooth_term.set((frame_num - 2)*(mesh_height)*(mesh_width) * 2 * 3);
	Mesh_path_smooth_term.fill(_v3d_(0, 0, 0));

	for (int t = 1; t < frame_num - 1; t++)
	{
		for (int i = 0; i < mesh_height; i++)
		{
			for (int j = 0; j < mesh_width; j++)
			{
				double smooth_e;
				if (weight_of_all_vertexs_smooth[t][i][j] < 1)
					smooth_e = 1 * e;
				else
				{
					smooth_e = weight_of_all_vertexs_smooth[t][i][j] * e;
				}

				smooth_e = weight_of_all_vertexs_smooth[t][i][j] * e;
				smooth_e = e;
				//x
				Mesh_path_smooth_term[term_count][0] = row_count;
				Mesh_path_smooth_term[term_count][1] = x_index(num_of_feature_points + (t - 1)*mesh_height*mesh_width + i*mesh_width + j);
				Mesh_path_smooth_term[term_count][2] = -1 * smooth_e;
				term_count++;

				Mesh_path_smooth_term[term_count][0] = row_count;
				Mesh_path_smooth_term[term_count][1] = x_index(num_of_feature_points + (t)*mesh_height*mesh_width + i*mesh_width + j);
				Mesh_path_smooth_term[term_count][2] = 2 * smooth_e;
				term_count++;

				Mesh_path_smooth_term[term_count][0] = row_count;
				Mesh_path_smooth_term[term_count][1] = x_index(num_of_feature_points + (t + 1)*mesh_height*mesh_width + i*mesh_width + j);
				Mesh_path_smooth_term[term_count][2] = -1 * smooth_e;
				term_count++;

				row_count++;

				//y
				Mesh_path_smooth_term[term_count][0] = row_count;
				Mesh_path_smooth_term[term_count][1] = y_index(num_of_feature_points + (t - 1)*mesh_height*mesh_width + i*mesh_width + j);
				Mesh_path_smooth_term[term_count][2] = -1 * smooth_e;
				term_count++;

				Mesh_path_smooth_term[term_count][0] = row_count;
				Mesh_path_smooth_term[term_count][1] = y_index(num_of_feature_points + (t)*mesh_height*mesh_width + i*mesh_width + j);
				Mesh_path_smooth_term[term_count][2] = 2 * smooth_e;
				term_count++;

				Mesh_path_smooth_term[term_count][0] = row_count;
				Mesh_path_smooth_term[term_count][1] = y_index(num_of_feature_points + (t + 1)*mesh_height*mesh_width + i*mesh_width + j);
				Mesh_path_smooth_term[term_count][2] = -1 * smooth_e;
				term_count++;

				row_count++;

			}
		}

	}
}

void union_smooth::Creat_mesh_time_smooth_term(ArrayDou& b)
{

	int term_count = 0;
	Mesh_time_smooth_term.erase();
	Mesh_time_smooth_term.set((frame_num)*mesh_height*mesh_width * 2 * 2 - mesh_height*mesh_width * 2);
	Mesh_time_smooth_term.fill(_v3d_(0, 0, 0));

	for (int t = 0; t < frame_num; t++)
	{
		for (int i = 0; i < mesh_height; i++)
		{
			for (int j = 0; j < mesh_width; j++)
			{
				if (t == 0)
				{
					//x
					Mesh_time_smooth_term[term_count][0] = row_count;
					Mesh_time_smooth_term[term_count][1] = x_index[num_of_feature_points + t*mesh_width*mesh_height + i*mesh_width + j];
					Mesh_time_smooth_term[term_count][2] = 1;
					term_count++;

					b[row_count] = source->GetVertex(i, j)[0];
					row_count++;

					//y
					Mesh_time_smooth_term[term_count][0] = row_count;
					Mesh_time_smooth_term[term_count][1] = y_index[num_of_feature_points + t*mesh_width*mesh_height + i*mesh_width + j];
					Mesh_time_smooth_term[term_count][2] = 1;
					term_count++;

					b[row_count] = source->GetVertex(i, j)[1];
					row_count++;

					continue;
				}

				//x
				Mesh_time_smooth_term[term_count][0] = row_count;
				Mesh_time_smooth_term[term_count][1] = x_index[num_of_feature_points + t*mesh_width*mesh_height + i*mesh_width + j];
				Mesh_time_smooth_term[term_count][2] = 1;
				term_count++;

				Mesh_time_smooth_term[term_count][0] = row_count;
				Mesh_time_smooth_term[term_count][1] = x_index[num_of_feature_points + (t - 1)*mesh_height*mesh_width + i*mesh_width + j];
				Mesh_time_smooth_term[term_count][2] = -1;
				term_count++;
				row_count++;

				//y
				Mesh_time_smooth_term[term_count][0] = row_count;
				Mesh_time_smooth_term[term_count][1] = y_index[num_of_feature_points + t*mesh_width*mesh_height + i*mesh_width + j];
				Mesh_time_smooth_term[term_count][2] = 1;
				term_count++;

				Mesh_time_smooth_term[term_count][0] = row_count;
				Mesh_time_smooth_term[term_count][1] = y_index[num_of_feature_points + (t - 1)*mesh_width*mesh_height + i*mesh_width + j];
				Mesh_time_smooth_term[term_count][2] = -1;
				term_count++;
				row_count++;
			}
		}
	}

}

void union_smooth::retern_mesh(vector<Bundled::Mesh>& mesh)
{
	mesh = each_frame_destin;
}

void union_smooth::change_mesh()
{
	for (int t = 0; t < frame_num; t++)
	{
		for (int i = 0; i < mesh_height - 1; i++)
		{
			for (int j = 0; j < mesh_width - 1; j++)
			{

				v2d p0 = each_frame_destin[t].GetVertex(i, j);
				v2d p1 = each_frame_destin[t].GetVertex(i, j + 1);
				v2d p2 = each_frame_destin[t].GetVertex(i + 1, j + 1);
				v2d p3 = each_frame_destin[t].GetVertex(i + 1, j);

				if (p0[0] > p1[0])
				{
					double temp = p0[0];
					p0[0] = p1[0];
					p1[0] = temp;

				}

				if (p3[0] > p2[0])
				{
					double temp = p3[0];
					p3[0] = p2[0];
					p2[0] = temp;

				}

				if (p0[1] > p3[1])
				{
					double temp = p0[1];
					p0[1] = p3[1];
					p3[1] = temp;

				}

				if (p1[1] > p2[1])
				{
					double temp = p1[1];
					p1[1] = p2[1];
					p2[1] = temp;

				}

				if (abs(p0[0] - p1[0]) < 2)
				{
					p1[0] += 2;
				}

				if (abs(p3[0] - p2[0]) < 2)
				{
					p2[0] += 2;
				}



				each_frame_destin[t].SetVertex(i, j, p0);
				each_frame_destin[t].SetVertex(i, j + 1, p1);
				each_frame_destin[t].SetVertex(i + 1, j + 1, p2);
				each_frame_destin[t].SetVertex(i + 1, j, p3);
			

			}
		}
	}
}

void union_smooth::find_homograph_from_warped_to_org()
{
	Bundled::TensorMatrix3d homos_avg;
	homos_avg.set(mesh_width, mesh_height);

	

	for (int t = 0; t < frame_num; t++)
	{
		for (int i = 0; i < mesh_height - 1; i++)
		{
			for (int j = 0; j < mesh_width - 1; j++)
			{
				double minx = j * qw;
				double miny = i * qh;
				double maxx = (j + 1) * qw;
				double maxy = (i + 1) * qh;

				cv::Mat s_(4, 2, CV_64F);
				cv::Mat d_(4, 2, CV_64F);

				v2d p0 = each_frame_destin[t].GetVertex(i, j);
				v2d p1 = each_frame_destin[t].GetVertex(i, j + 1);
				v2d p2 = each_frame_destin[t].GetVertex(i + 1, j + 1);
				v2d p3 = each_frame_destin[t].GetVertex(i + 1, j);

				
				

				s_.at<double>(0, 0) = minx;
				s_.at<double>(0, 1) = miny;

				s_.at<double>(1, 0) = maxx;
				s_.at<double>(1, 1) = miny;

				s_.at<double>(2, 0) = maxx;
				s_.at<double>(2, 1) = maxy;

				s_.at<double>(3, 0) = minx;
				s_.at<double>(3, 1) = maxy;

				d_.at<double>(0, 0) = p0[0];
				d_.at<double>(0, 1) = p0[1];

				d_.at<double>(1, 0) = p1[0];
				d_.at<double>(1, 1) = p1[1];

				d_.at<double>(2, 0) = p2[0];
				d_.at<double>(2, 1) = p2[1];

				d_.at<double>(3, 0) = p3[0];
				d_.at<double>(3, 1) = p3[1];
				
			
				homos_avg.at(i, j) = cv::findHomography(d_, s_);

				if (homos_avg.at(i, j).empty())
				{
					Mat a(3, 3, CV_32FC1);
					homos_avg.at(i, j) = a;
				}


			}
		}
		homograph_from_warped_to_org.push_back(homos_avg);
	}
}

void union_smooth::warping_video(double fps, bool border)
{

	if (border)
	{
		cv::VideoWriter vw;
		vw.open(save_dest_video_path, CV_FOURCC('X', 'V', 'I', 'D'),
			fps,
			cv::Size(w + 240, h + 240),
			true);

		//save the result video with mesh
		cv::VideoWriter vw_mesh;
		vw_mesh.open(save_cut_dst_video_path, CV_FOURCC('X', 'V', 'I', 'D'),
			fps,
			cv::Size(w + 240, h + 240),
			true);




		//warping every mesh

		/*cv::Mat map_x(h, w, CV_32FC1, cv::Scalar(-1));*/
		//cv::Mat map_y(h, w, CV_32FC1, cv::Scalar(-1));
		cv::VideoCapture vc1(open_org_video_path);

		//test

		Mat map_x_single(h, w, CV_32FC1, cv::Scalar(-1));
		Mat map_y_single(h, w, CV_32FC1, cv::Scalar(-1));
		Mat org_single;
		Mat dst_single;

		for (int t = 0; t < frame_num - 1; t++)
		{

			cv::Mat map_x(h + 240, w + 240, CV_32FC1, cv::Scalar(-1));
			cv::Mat map_y(h + 240, w + 240, CV_32FC1, cv::Scalar(-1));

			Mat org;
			Mat biger_org(h + 240, w + 240, CV_8UC3, cv::Scalar(0, 0, 0));
			Mat dst(h + 240, w + 240, CV_8SC3, cv::Scalar(0, 0, 0));
			vc1 >> org;
			org.copyTo(biger_org(Range(0, h), Range(0, w)));

			for (int i = 0; i < mesh_height - 1; i++)
			{
				for (int j = 0; j < mesh_width - 1; j++)
				{


					double ox1 = j * qw;
					double oy1 = i * qh;
					double ox2 = ox1 + qw;
					double oy2 = oy1 + qh;

					v2d point1 = each_frame_destin[t].GetVertex(i, j);
					v2d point2 = each_frame_destin[t].GetVertex(i, j + 1);
					v2d point3 = each_frame_destin[t].GetVertex(i + 1, j + 1);
					v2d point4 = each_frame_destin[t].GetVertex(i + 1, j);

					double x1 = point1[0];
					double x2 = point2[0];
					double x3 = point3[0];
					double x4 = point4[0];

					double y1 = point1[1];
					double y2 = point2[1];
					double y3 = point3[1];
					double y4 = point4[1];

					


					

					



					/*double x1 = avg_warped_meshes[t].at(i, j)[0];
					double x2 = avg_warped_meshes[t].at(i, j + 1)[0];
					double x3 = avg_warped_meshes[t].at(i + 1, j + 1)[0];
					double x4 = avg_warped_meshes[t].at(i + 1, j)[0];

					double y1 = avg_warped_meshes[t].at(i, j)[1];
					double y2 = avg_warped_meshes[t].at(i, j + 1)[1];
					double y3 = avg_warped_meshes[t].at(i + 1, j + 1)[1];
					double y4 = avg_warped_meshes[t].at(i + 1, j)[1];*/

					//finding min x
					double minx = x1 < x2 ? x1 : x2; minx = minx < x3 ? minx : x3; minx = minx < x4 ? minx : x4;
					double maxx = x1 > x2 ? x1 : x2; maxx = maxx > x3 ? maxx : x3; maxx = maxx > x4 ? maxx : x4;
					//finding min y
					double miny = y1 < y2 ? y1 : y2; miny = miny < y3 ? miny : y3; miny = miny < y4 ? miny : y4;
					double maxy = y1 > y2 ? y1 : y2; maxy = maxy > y3 ? maxy : y3; maxy = maxy > y4 ? maxy : y4;



					for (int ii = floor(miny); ii <= ceil(maxy); ii++)
					{
						for (int jj = floor(minx); jj <= ceil(maxx); jj++)
						{

							double ux, uy;
							TransformByHomo(homograph_from_warped_to_org[t].at(i, j), jj, ii, ux, uy);
							//TransformByHomo(home[t].at(i, j), jj, ii, ux, uy);


							if (ux >= ox1 && ux < ox2 && uy >= oy1 && uy < oy2)
							{
								map_x.at<float>(ii + 120, jj + 120) = ux;
								map_y.at<float>(ii + 120, jj + 120) = uy;
							}

						}
					}

				}
			}


			cv::remap(biger_org, dst, map_x, map_y, cv::INTER_LINEAR, cv::BORDER_CONSTANT, Scalar(0, 0, 0));


			vw << dst;



			for (int h = 0; h < each_frame_destin[t].GetMeshHeight(); h++)
				for (int w = 0; w < each_frame_destin[t].GetMeshWidth(); w++)
				{
					for (int x = -1; x <= 1; x++)
					{

						if ((x + w) < 0 || (x + w) >= each_frame_destin[t].GetMeshWidth())
						{
							continue;
						}

						line(dst, Point(each_frame_destin[t].GetVertex(h, w)[0] + 120, each_frame_destin[t].GetVertex(h, w)[1] + 120),
							Point(each_frame_destin[t].GetVertex(h, w + x)[0] + 120, each_frame_destin[t].GetVertex(h, w + x)[1] + 120), Scalar(255, 255, 0), 2);
					}
					for (int y = -1; y <= 1; y++)
					{
						if ((y + h) < 0 || (y + h) >= each_frame_destin[t].GetMeshHeight())
						{
							continue;
						}
						line(dst, Point(each_frame_destin[t].GetVertex(h, w)[0] + 120, each_frame_destin[t].GetVertex(h, w)[1] + 120),
							Point(each_frame_destin[t].GetVertex(h + y, w)[0] + 120, each_frame_destin[t].GetVertex(h + y, w)[1] + 120), Scalar(255, 255, 0), 2);
					}
				}

			vw_mesh << dst;

			namedWindow("test");
			imshow("test", dst);


		}

	}
	else
	{
		cv::VideoWriter vw;
		vw.open(save_dest_video_path, CV_FOURCC('X', 'V', 'I', 'D'),
			fps,
			cv::Size(w, h),
			true);

		//save mesh
		cv::VideoWriter vw_mesh;
		vw_mesh.open(save_cut_dst_video_path, CV_FOURCC('X', 'V', 'I', 'D'),
			fps,
			cv::Size(w, h),
			true);




		//warping every mesh

		/*cv::Mat map_x(h, w, CV_32FC1, cv::Scalar(-1));*/
		//cv::Mat map_y(h, w, CV_32FC1, cv::Scalar(-1));
		cv::VideoCapture vc1(open_org_video_path);

		//test

		for (int t = 0; t < frame_num - 1; t++)
		{

			cv::Mat map_x(h, w, CV_32FC1, cv::Scalar(-1));
			cv::Mat map_y(h, w, CV_32FC1, cv::Scalar(-1));

			Mat org;
			Mat dst(h, w, CV_8SC3, cv::Scalar(0, 0, 0));
			vc1 >> org;

			for (int i = 0; i < mesh_height - 1; i++)
			{
				for (int j = 0; j < mesh_width - 1; j++)
				{


					double ox1 = j * qw;
					double oy1 = i * qh;
					double ox2 = ox1 + qw;
					double oy2 = oy1 + qh;

					v2d point1 = each_frame_destin[t].GetVertex(i, j);
					v2d point2 = each_frame_destin[t].GetVertex(i, j + 1);
					v2d point3 = each_frame_destin[t].GetVertex(i + 1, j + 1);
					v2d point4 = each_frame_destin[t].GetVertex(i + 1, j);

					double x1 = point1[0];
					double x2 = point2[0];
					double x3 = point3[0];
					double x4 = point4[0];

					double y1 = point1[1];
					double y2 = point2[1];
					double y3 = point3[1];
					double y4 = point4[1];

					

					//finding min x
					double minx = x1 < x2 ? x1 : x2; minx = minx < x3 ? minx : x3; minx = minx < x4 ? minx : x4;
					double maxx = x1 > x2 ? x1 : x2; maxx = maxx > x3 ? maxx : x3; maxx = maxx > x4 ? maxx : x4;
					//finding min y
					double miny = y1 < y2 ? y1 : y2; miny = miny < y3 ? miny : y3; miny = miny < y4 ? miny : y4;
					double maxy = y1 > y2 ? y1 : y2; maxy = maxy > y3 ? maxy : y3; maxy = maxy > y4 ? maxy : y4;



					for (int ii = floor(miny); ii <= ceil(maxy); ii++)
					{
						for (int jj = floor(minx); jj <= ceil(maxx); jj++)
						{
							if (ii >= 0 && ii < h &&
								jj >= 0 && jj < w)
							{
								double ux, uy;
								TransformByHomo(homograph_from_warped_to_org[t].at(i, j), jj, ii, ux, uy);
								//TransformByHomo(home[t].at(i, j), jj, ii, ux, uy);


								if (ux >= ox1 && ux < ox2 && uy >= oy1 && uy < oy2)
								{
									map_x.at<float>(ii, jj) = ux;
									map_y.at<float>(ii, jj) = uy;
								}
							}
						}
					}

				}
			}


			cv::remap(org, dst, map_x, map_y, cv::INTER_LINEAR, cv::BORDER_CONSTANT, Scalar(255, 255, 255));


			vw << dst;


			for (int h = 0; h < each_frame_destin[t].GetMeshHeight(); h++)
				for (int w = 0; w < each_frame_destin[t].GetMeshWidth(); w++)
				{
					for (int x = -1; x <= 1; x++)
					{

						if ((x + w) < 0 || (x + w) >= each_frame_destin[t].GetMeshWidth())
						{
							continue;
						}

						line(dst, Point(each_frame_destin[t].GetVertex(h, w)[0], each_frame_destin[t].GetVertex(h, w)[1]),
							Point(each_frame_destin[t].GetVertex(h, w + x)[0], each_frame_destin[t].GetVertex(h, w + x)[1]), Scalar(255, 255, 0), 2);
					}
					for (int y = -1; y <= 1; y++)
					{
						if ((y + h) < 0 || (y + h) >= each_frame_destin[t].GetMeshHeight())
						{
							continue;
						}
						line(dst, Point(each_frame_destin[t].GetVertex(h, w)[0], each_frame_destin[t].GetVertex(h, w)[1]),
							Point(each_frame_destin[t].GetVertex(h + y, w)[0], each_frame_destin[t].GetVertex(h + y, w)[1]), Scalar(255, 255, 0), 2);
					}
				}

			vw_mesh << dst;

			namedWindow("test");
			imshow("test", dst);
			waitKey(10);
		
		}





	}
}





inline double computer_sigma(double value)
{
	double sigma = sqrt(-(value*value) / (2 * log(0.5)));
	return sigma;
}

inline double Gaussian(double sigma, double value)
{
	/*double a = 1 / (sigma*sqrt(2 * pi));*/
	double a = 3;
	double result = a*exp(-(value*value) / (2 * sigma*sigma));
	return result;
}

void union_smooth::ComputeMeshPointsAndWeight(vector<std::vector<v2d>> frame_points)
{



	//Traversing through all the features of a frame
	for (int t = 0; t < frame_points.size(); t++)
	{
		Points_num_in_mesh.clear();
		Points_num_in_mesh.resize(mesh_height + 1);

		for (int i = 0; i < mesh_height + 1; i++)
			Points_num_in_mesh[i].resize(mesh_width + 1);

		for (int i = 0; i < frame_points[t].size(); i++)
		{
			int mesh_x = ceil(frame_points[t][i][0] / qw);
			int mesh_y = ceil(frame_points[t][i][1] / qh);
			mesh_x >(mesh_width - 1) ? mesh_x = (mesh_width - 1) : mesh_x = mesh_x;
			mesh_y >(mesh_height - 1) ? mesh_y = (mesh_height - 1) : mesh_y = mesh_y;

			Points_num_in_mesh[mesh_y][mesh_x]++;

		}

		Points_num_in_all_mesh.push_back(Points_num_in_mesh);
	}


	cout << "compute point number over" << endl;
	get_Median();
	if (weight_mesh_similarty < 1)
		weight_mesh_similarty = 1;
	if (weight_mesh_smooth < 1)
		weight_mesh_smooth = 1;
	//caught
	cout << weight_mesh_smooth << " " << weight_mesh_similarty << endl;


	//weight of mesh smooth
	for (int t = 0; t < frame_num; t++)
	{
		vector<vector<double>> weight_of_vertex_smooth;
		weight_of_vertex_smooth.resize(mesh_height);

		vector<vector<double>> weight_of_vertex_similarity;
		weight_of_vertex_similarity.resize(mesh_height - 1);

		for (int h = 0; h<mesh_height; h++)
			for (int w = 0; w < mesh_width; w++)
			{


				double Point_num_smooth = Points_num_in_all_mesh[t][h + 1][w + 1] +
					Points_num_in_all_mesh[t][h][w + 1] +
					Points_num_in_all_mesh[t][h][w] +
					Points_num_in_all_mesh[t][h + 1][w];
				double weight_smooth = Gaussian(weight_mesh_smooth, Point_num_smooth);
				weight_of_vertex_smooth[h].push_back(weight_smooth);
			}

		weight_of_all_vertexs_smooth.push_back(weight_of_vertex_smooth);

		//weight of mesh similarity
		for (int h = 0; h < mesh_height - 1; h++)
		{
			for (int w = 0; w < mesh_width - 1; w++)
			{


				double Point_num_similarity = Points_num_in_all_mesh[t][h + 1][w + 1] +
					Points_num_in_all_mesh[t][h + 1][w + 2] +
					Points_num_in_all_mesh[t][h + 2][w + 2] +
					Points_num_in_all_mesh[t][h + 2][w + 1] +
					Points_num_in_all_mesh[t][h + 2][w] +
					Points_num_in_all_mesh[t][h + 1][w] +
					Points_num_in_all_mesh[t][h][w] +
					Points_num_in_all_mesh[t][h][w + 1] +
					Points_num_in_all_mesh[t][h][w + 2];

				double weight_similarity = Gaussian(weight_mesh_similarty, Point_num_similarity);
				weight_of_vertex_similarity[h].push_back(weight_similarity);




			}
		}

		weight_of_all_vertexs_similarity.push_back(weight_of_vertex_similarity);

		weight_of_vertex_smooth.clear();
		weight_of_vertex_smooth.clear();

		vector<vector<double>>(weight_of_vertex_smooth).swap(weight_of_vertex_smooth);
		vector<vector<double>>(weight_of_vertex_similarity).swap(weight_of_vertex_similarity);
	}



	cout << "ComputeMeshPointsAndWeight over" << endl;
}

void union_smooth::get_Median()
{
	//smooth
	int* feature = new int[frame_num*(mesh_height)*(mesh_width)];

	for (int t = 0; t < frame_num; t++)
	{
		for (int i = 0; i<mesh_height; i++)
			for (int j = 0; j < mesh_width; j++)
			{
				double Point_num_smooth = Points_num_in_all_mesh[t][i + 1][j + 1] +
					Points_num_in_all_mesh[t][i][j + 1] +
					Points_num_in_all_mesh[t][i][j] +
					Points_num_in_all_mesh[t][i + 1][j];

				feature[t*mesh_height*mesh_width + i*(mesh_width)+j] = Point_num_smooth;

				Point_num_smooth = 0;
			}
	}


	std::sort(feature, feature + frame_num*(mesh_height)*(mesh_width));
	int median_smooth = (frame_num*(mesh_height)*(mesh_width) / 4);
	weight_mesh_smooth = feature[median_smooth];




	delete[] feature;

	feature = new int[frame_num*(mesh_height - 1)*(mesh_width - 1)];

	//similarity
	for (int t = 0; t < frame_num; t++)
	{
		for (int i = 0; i<mesh_height - 1; i++)
			for (int j = 0; j < mesh_width - 1; j++)
			{
				double Point_num_similarity = Points_num_in_all_mesh[t][i + 1][j + 1] +
					Points_num_in_all_mesh[t][i + 1][j + 2] +
					Points_num_in_all_mesh[t][i + 2][j + 2] +
					Points_num_in_all_mesh[t][i + 2][j + 1] +
					Points_num_in_all_mesh[t][i + 2][j] +
					Points_num_in_all_mesh[t][i + 1][j] +
					Points_num_in_all_mesh[t][i][j] +
					Points_num_in_all_mesh[t][i][j + 1] +
					Points_num_in_all_mesh[t][i][j + 2];

				feature[t*(mesh_height - 1)*(mesh_width - 1) + (i)*(mesh_width - 1) + j] = Point_num_similarity;

				Point_num_similarity = 0;
			}
	}

	std::sort(feature, feature + frame_num*(mesh_height - 1)*(mesh_width - 1));
	int median_similarity = (frame_num*(mesh_height - 1)*(mesh_width - 1) / 3);
	weight_mesh_similarty = feature[median_similarity];


	delete[] feature;
	feature = nullptr;
}


