#pragma once
#include "common.h"

void read_computer_path(string path, set_of_path& PathSet);

void writr_computered_path(string path, set_of_path PathSet);

void ComputePointsOfFrame(int frame_num, std::vector<std::vector<v2d>> & points_of_frames, vector<std::vector<int>>& index, const set_of_path & paths);

cv::Mat read_path_and_save_in_mat(string path, set_of_path& org_paths, int frame_num);