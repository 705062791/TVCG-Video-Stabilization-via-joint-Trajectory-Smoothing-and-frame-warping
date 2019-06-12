#pragma once
#include"common.h"
#include <Eigen/Dense>
#include<Eigen/Sparse>
#include<iostream>
using namespace std;

	class Motion_conformantion
	{
	public:

		Motion_conformantion(set_of_path org_paths, vector<vector<v2d>> frame_points,double alpha);
		//~Motion_conformantion();
		void solve_by_eigen();
		set_of_path return_org_path();
	private:
		int Point_num;
		int matJ_cols;
		int matJ_rows;
		int frame_num;
		int alpha;

		vector<int> x_index;
		vector<int> y_index;
		vector<vector<v2d>> frame_points;
		set_of_path org_path;
		set_of_path motion_conformantion_path;

		vector<Eigen::Triplet<double>> term;

		void set_data_term(vector<double>& B);
		void set_conformantion_term(vector<double>& B);
		void sreach_neighbor(int path_index, int point_index,vector<int>& neighbor_path, vector<int>& neighbor_point);
		void save_conformantion(Eigen::VectorXd path);
		int get_index(int i, int j);
	};

