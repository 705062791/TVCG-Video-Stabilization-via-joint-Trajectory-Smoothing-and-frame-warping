#include "common.h"

void read_computer_path(string path, set_of_path& PathSet)
{
	point_path temp_path_path;
	string temp;
	ifstream read(path);
	if (!read)
	{
		cout << "open file false" << endl;
	}
	double start=0, end=0;
	double x, y;
	while (getline(read, temp) && temp != " ")
	{

		istringstream(temp) >> start >> end;;
		int lengh = end - start+1;
		temp_path_path.start_frame = start;
		temp_path_path.end_frame = end;

		for (int i = 0; i < lengh; i++)
		{
			getline(read, temp);
			istringstream(temp) >> x >> y;
			v2d point_;
			point_[0] = x;
			point_[1] = y;
			temp_path_path.PointPath.push_back(point_);
		}

		PathSet.push_back(temp_path_path);
		temp_path_path.PointPath.clear();
	}


	temp_path_path.PointPath.clear();
	vector<LibIV::Math::v2d>(temp_path_path.PointPath).swap(temp_path_path.PointPath);

	cout << "we have" << PathSet.size() << "paths in totall" << endl;;


}

void writr_computered_path(string path, set_of_path PathSet)
{
	string temp, temp2;
	ofstream write(path);
	if (!write)
	{
		cout << "open write file false" << endl;
	}

	stringstream writer_in_txt;

	for (int i = 0; i < PathSet.size(); i++)
	{
		writer_in_txt.clear();
		writer_in_txt << PathSet[i].start_frame << " " << PathSet[i].end_frame;
		writer_in_txt >> temp;
		write << temp;
		write << " ";
		writer_in_txt >> temp;
		write << temp;
		write << endl;
		for (int k = 0; k < PathSet[i].PointPath.size(); k++)
		{
			writer_in_txt.clear();
			//cout << (double)PathSet[i].PointPath[k][0] << " " << (double)PathSet[i].PointPath[k][1] << endl;
			writer_in_txt << (double)PathSet[i].PointPath[k][0] << " " << (double)PathSet[i].PointPath[k][1];
			writer_in_txt >> temp2;
			write << temp2;
			write << " ";
			writer_in_txt >> temp2;
			write << temp2;
			write << endl;
		}
		//cout << endl;
	}
	write.close();
}

void ComputePointsOfFrame(int frame_num, std::vector<std::vector<v2d>> & points_of_frames, vector<std::vector<int>>& index, const set_of_path & paths)
{

	points_of_frames.clear();
	points_of_frames.resize(frame_num);

	index.resize(frame_num);

	int point_index=0;
	for (size_t i = 0; i < paths.size(); i++)
	{

		for (size_t j = 0; j < paths[i].PointPath.size(); j++)
		{
			points_of_frames[paths[i].start_frame + j - 1].push_back(paths[i].PointPath[j]);
			
	        
			index[paths[i].start_frame + j - 1].push_back(point_index + j);
			
		}

		point_index += paths[i].PointPath.size();
	}

	
}

cv::Mat read_path_and_save_in_mat(string path, set_of_path& org_paths, int frame_num)
{
	
	read_computer_path(path, org_paths);
	int num_of_path = org_paths.size();


	Mat all_path = Mat::zeros(frame_num * 2,num_of_path,  CV_64F);
	cout << num_of_path << " " << frame_num * 2 << endl;

	v2d temp;
	int row_count = 0;

	for (size_t i = 0; i < org_paths.size(); i++)
	{
		row_count++;
		int start = org_paths[i].start_frame;
		for (size_t j = 0; j < org_paths[i].PointPath.size(); j++)
		{

			temp = org_paths[i].PointPath[j];

			all_path.at<double>((start - 1 + j) * 2, row_count) = temp[0];
			all_path.at<double>((start - 1 + j) * 2 + 1, row_count) = temp[1];
		}
	}
	return all_path;
}