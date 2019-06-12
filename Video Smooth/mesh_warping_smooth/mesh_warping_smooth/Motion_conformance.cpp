#include"Motion_conformance.h"

Motion_conformantion::Motion_conformantion(set_of_path org_paths, vector<vector<v2d>> frame_points, double alpha)
{
	this->org_path = org_paths;
	this->frame_points = frame_points;
	this->alpha = alpha;

	Point_num = 0;
	frame_num = frame_points.size();

	for (int i = 0; i < frame_points.size(); i++)
	{
		Point_num += frame_points[i].size();
	}
	
	x_index.resize(Point_num, 0);
	y_index.resize(Point_num, 0);

	for (int i = 0; i < Point_num; i++)
	{
		x_index[i] = i;
		y_index[i] = i + Point_num;
	}

}

void Motion_conformantion::sreach_neighbor(int path_index, int point_index, vector<int>& neighbor_paths,vector<int>& neighbor_points)
{
	
	int t;

	for (t = path_index + 1; t < org_path.size(); t++)
	{

		if (org_path[t].start_frame <= point_index + org_path[path_index].start_frame <= org_path[t].end_frame)
		{
			neighbor_paths.push_back(t);
			neighbor_points.push_back(point_index + org_path[path_index].start_frame - org_path[t].start_frame);
		}
	}

		
		for (t = path_index - 1; t >= 0; t--)
		{
			if (org_path[t].start_frame <= point_index + org_path[path_index].start_frame <= org_path[t].end_frame)
			{
				neighbor_paths.push_back(t);
				neighbor_points.push_back(point_index + org_path[path_index].start_frame - org_path[t].start_frame);
			}
		}


}

int Motion_conformantion::get_index(int i, int j)
{
	int index=0;
	for (int x = 0; x < i;  x++)
	{
		index += org_path[x].PointPath.size();
	}

	for (int x = 0; x < j; x++)
	{
		index++;
	}

	return index;
}

void Motion_conformantion::save_conformantion(Eigen::VectorXd path)
{
	for (int i = 0; i < org_path.size(); i++)
	{
		for (int j = 0; j < org_path[i].PointPath.size(); j++)
		{
			org_path[i].PointPath[j][0] = path[get_index(i, j)];
			org_path[i].PointPath[j][1] = path[Point_num + get_index(i, j)];
		}
	 }
}

set_of_path Motion_conformantion::return_org_path()
{
	return org_path;
}

void Motion_conformantion::solve_by_eigen()
{

	vector<double> B;

	matJ_cols = Point_num*2;
	matJ_rows = 0;

 

	//data_term
    set_data_term(B);


	//conformantion_term
	set_conformantion_term(B);

	//for (int i = Point_num; i < term.size(); i++)
	//{
	//	cout << term[i].row() << " " << term[i].col() << " " << term[i].value() << endl;
	//}


	
	
	matJ_rows = term.size();

	Eigen::SparseMatrix<double> matj_sparse(matJ_rows, matJ_cols);

	Eigen::VectorXd bJ(matJ_rows);

	for (int i = 0; i < B.size(); i++)
	{
		bJ[i] = B[i];
	}
	for (int i = B.size(); i < matJ_rows; i++)
	{
		bJ[i] = 0;
	}

	matj_sparse.setFromTriplets(term.begin(), term.end());

	Eigen::VectorXd x_sparse(matJ_cols);


	Eigen::SparseMatrix<double> mat_t = matj_sparse.transpose()*matj_sparse;
	Eigen::VectorXd b_t = matj_sparse.transpose()*bJ;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	solver.compute(mat_t);
	x_sparse = solver.solve(b_t);


	save_conformantion(x_sparse);


}

void Motion_conformantion::set_data_term(vector<double>& B)
{
	
	int row;
	int col;
	double value;

	int count = 0;

	for (int i = 0; i < org_path.size(); i++)
	{
		for (int j = 0; j < org_path[i].PointPath.size(); j++)
		{

			//x
			row = matJ_rows;
			col = x_index[count];
			value = 1;
			term.push_back(Eigen::Triplet<double>(row,col,value));
			B.push_back(org_path[i].PointPath[j][0]);
			matJ_rows++;

			//y
			row = matJ_rows;
			col = y_index[count];
			value = 1;
			term.push_back(Eigen::Triplet<double>(row, col, value));
			B.push_back(org_path[i].PointPath[j][1]);
			
			
			count++;
			matJ_rows++;
		}
	}

}

void Motion_conformantion::set_conformantion_term(vector<double>& B)
{

	int count = 0;
	vector<int> neighbor_path;
	vector<int> neighbor_point;
	int row;
	int col;
	double value;



	for (int i = 0; i < org_path.size(); i++)
	{
		for (int j = 0; j < org_path[i].PointPath.size(); j++)
		{

			sreach_neighbor(i, j, neighbor_path, neighbor_point);

			cout << neighbor_path.size() << endl;
			
			if (j == 0)
			{
				count++;
				continue;
			}

			//x
			row = matJ_rows;
			col = x_index[count];
			value = neighbor_path.size() * alpha;
			term.push_back(Eigen::Triplet<double>(row, col, value));

			row = matJ_rows;
			col = x_index[count - 1];
			value = -1*neighbor_path.size() * alpha;
			term.push_back(Eigen::Triplet<double>(row, col, value));



			for (int t = 0; t < neighbor_path.size(); t++)
			{
				//x
				
				row = matJ_rows;
				col = x_index[get_index(neighbor_path[t], neighbor_point[t])];
				value = -1 * alpha;
				term.push_back(Eigen::Triplet<double>(row, col, value));

				row = matJ_rows;
				col = x_index[get_index(neighbor_path[t], neighbor_point[t] - 1)];
				value = 1 * alpha;
				term.push_back(Eigen::Triplet<double>(row, col, value));


			}
			matJ_rows++;

			//y
			row = matJ_rows;
			col = y_index[count];
			value = neighbor_path.size() * alpha;
			term.push_back(Eigen::Triplet<double>(row, col, value));

			row = matJ_rows;
			col = y_index[count - 1];
			value = -1*neighbor_path.size() * alpha;
			term.push_back(Eigen::Triplet<double>(row, col, value));

			for (int t = 0; t < neighbor_path.size(); t++)
			{
				//y

				row = matJ_rows;
				col = y_index[get_index(neighbor_path[t], neighbor_point[t])];
				value = -1 * alpha;
				term.push_back(Eigen::Triplet<double>(row, col, value));

				row = matJ_rows;
				col = y_index[get_index(neighbor_path[t], neighbor_point[t] - 1)];
				value = 1 * alpha;
				term.push_back(Eigen::Triplet<double>(row, col, value));

				
				
			}
			matJ_rows++;
			count++;

			neighbor_path.clear();
			neighbor_point.clear();
		}
	}
}