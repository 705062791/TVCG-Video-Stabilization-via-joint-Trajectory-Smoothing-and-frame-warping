#pragma once
#include"common.h"
#include"commonio.h"
#include"Mesh.h"
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<cmath>
#include<memory>
#include<algorithm>
namespace Bundled
{
	class union_smooth
	{
	public:


		union_smooth(string open_org,string save_dest ,string save_cut_dst,set_of_path org_paths, std::vector<std::vector<int>>index, int video_h, int video_w, double qw, double qh, double a, double b, double c, double d,double e);

		void change_mesh();

		void set_control_points(vector<std::vector<v2d>> frame_points);

		void solve_by_Eigen(bool rotaion);

		void retern_mesh(vector<Bundled::Mesh>& mesh);

		void find_homograph_from_warped_to_org();

		void warping_video(double fps,bool border);

		void TransformByHomo(cv::Mat & mat, double x, double y, double & dx, double & dy);

		void CutSave();

		void ComputeMeshPointsAndWeight(vector<std::vector<v2d>> frame_points);

		~union_smooth()
		{
			Path_data_term.erase();
			Path_smooth_term2.erase();
			Path_smooth_term.erase();
			Mesh_data_term.erase();
			Mesh_smooth_term.erase();
			Mesh_time_smooth_term.erase();
			Mesh_path_smooth_term.erase();

			x_index, y_index;
			each_frame_destin.clear();
			frame_points.clear();
			homograph_from_warped_to_org.clear();
			each_frame_histogram.clear();
			each_frame_histogram2.clear();
			each_frame_maxHistogram.clear();
			Points_num_in_mesh.clear();
			Points_num_in_all_mesh.clear();
			dataterm_element_V00.clear();
			dataterm_element_V01.clear();
			dataterm_element_V10.clear();
			dataterm_element_V11.clear();
			weight_of_all_vertexs_similarity.clear();
			weight_of_all_vertexs_smooth.clear();
			index.clear();

			set_of_path(org_paths).swap(org_paths);
			vector<Bundled::Mesh>(each_frame_destin).swap(each_frame_destin);
			vector<std::vector<v2d>>(frame_points).swap(frame_points);
			vector<Bundled::TensorMatrix3d>(homograph_from_warped_to_org).swap(homograph_from_warped_to_org);
			vector<TensorInt>(each_frame_histogram).swap(each_frame_histogram);
			vector<TensorDou>(each_frame_histogram2).swap(each_frame_histogram2);
			vector<int>(each_frame_maxHistogram).swap(each_frame_maxHistogram);
			vector<vector<int>>(Points_num_in_mesh).swap(Points_num_in_mesh);
			vector<vector<vector<int>>>(Points_num_in_all_mesh).swap(Points_num_in_all_mesh);

			vector<vector<double>>(dataterm_element_V00).swap(dataterm_element_V00);
			vector<vector<double>>(dataterm_element_V01).swap(dataterm_element_V01);
			vector<vector<double>>(dataterm_element_V10).swap(dataterm_element_V10);
			vector<vector<double>>(dataterm_element_V11).swap(dataterm_element_V11);
			vector<vector<vector<double>>>(weight_of_all_vertexs_similarity).swap(weight_of_all_vertexs_similarity);
			vector<vector<vector<double>>>(weight_of_all_vertexs_smooth).swap(weight_of_all_vertexs_smooth);
			std::vector<std::vector<int>>(index).swap(index);
		}
	private:

		Bundled::Mesh       *source ,*destin;
		vector<Bundled::Mesh> each_frame_destin;
		Mat all_path;//all feature in each frame were saved in path
		set_of_path org_paths;//all paths 
		vector<std::vector<v2d>> frame_points;
		vector<Bundled::TensorMatrix3d> homograph_from_warped_to_org;

		int mesh_height;
		int mesh_width;
		int  h;
		int  w;
		int num_of_feature_points;
		int frame_num;
		int row_count;

		int num_of_path_data_term;
		int num_of_path_smooth_term_2;
		int num_of_path_smooth_term;
		int num_of_mesh_data_term;
		int num_of_mesh_smooth_term;
		int num_of_mesh_time_smooth_term;
		int num_of_mesh_path_smooth_term;

		double qw;
		double qh;
		double a, b, c, d, e;  //weight for each term
		double weight_mesh_smooth;
		double weight_mesh_similarty;
		

		string open_org_video_path;
		string save_dest_video_path;
		string save_cut_dst_video_path;

		vector<TensorInt>           each_frame_histogram;
		vector<TensorDou>           each_frame_histogram2;
		vector<int>                 each_frame_maxHistogram;
		vector<vector<int>>         Points_num_in_mesh;
		vector<vector<vector<int>>> Points_num_in_all_mesh;

		vector<vector<double>>            dataterm_element_V00; //¿ØÖÆµã
		vector<vector<double>>            dataterm_element_V01;
		vector<vector<double>>            dataterm_element_V10;
		vector<vector<double>>            dataterm_element_V11;
		vector<vector<vector<double>>>    weight_of_all_vertexs_similarity;
		vector<vector<vector<double>>>    weight_of_all_vertexs_smooth;

		ArrayV3d            Path_data_term;
		ArrayV3d            Path_smooth_term2;
		ArrayV3d            Path_smooth_term;
		ArrayV3d            Mesh_data_term;
		ArrayV3d            Mesh_smooth_term;
		ArrayV3d            Mesh_time_smooth_term;
		ArrayV3d            Mesh_path_smooth_term;
	 
		ArrayInt            x_index, y_index;

		Rect cut;

		std::vector<std::vector<int>>index;
		double GetWeight(int frame,int i, int h, double c);

		void Creat_path_data_term(ArrayDou& b);

		void Creat_path_smooth_term2(ArrayDou& b,bool);

		void Creat_path_smooth_term(ArrayDou& b);

		void Creat_mash_data_term(ArrayDou& b);

		void Creat_mesh_smooth_term(ArrayDou& b);

		void Creat_mesh_path_smooth_term(ArrayDou& b);

		void SetUpTriangle(int frame, int i, int j, double aveg_c,int& term_count);
		void SetDownTriangle(int frame, int i, int j, double aveg_c,int& term_count);
		void GetSmoothWeight(v2d v1, v2d v2, v2d v3, double uv[2]);

		void Creat_mesh_time_smooth_term(ArrayDou& b);

		void get_Median();


	};
}
