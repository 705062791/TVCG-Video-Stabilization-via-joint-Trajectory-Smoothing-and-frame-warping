#pragma once
#include"common.h"

class Median_Filter
{
public:

	Median_Filter(set_of_path org_path,int frame_num,int window);
	void median_filter();
	set_of_path return_org_path();
private:
	int window;  
	int path_num;
	int frame_num;

	vector<vector<Point2d>> motion;
	vector<vector<Point2d>> frame_motion;
	vector<vector<Point2i>> index_frame_motion;

	set_of_path org_path;

	
	void get_motion();
	void find_frame_motion();
	void exchang(vector<Point2d>& input, int i, int j);
	void sort(vector<Point2d> input);
};

