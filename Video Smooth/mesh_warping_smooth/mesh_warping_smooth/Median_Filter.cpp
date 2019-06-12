#include"Median_Filter.h"

Median_Filter::Median_Filter(set_of_path org_path,int frame_num, int window)
{
	this->frame_num = frame_num;
	this->org_path = org_path;
	this->window = window;
	path_num = org_path.size();
}

set_of_path Median_Filter::return_org_path()
{
	return org_path;
}

void Median_Filter::get_motion()
{
	motion.resize(org_path.size()-1);
	for (int i = 0; i < org_path.size()-1; i++)
	{
		for (int j = 1; j < org_path[i].PointPath.size(); j++)
		{
			double v_x = org_path[i].PointPath[j][0] - org_path[i].PointPath[j - 1][0];
			double v_y = org_path[i].PointPath[j][1] - org_path[i].PointPath[j - 1][1];

			motion[i].push_back(Point2d(v_x, v_y));
		}
	}
}

void Median_Filter::find_frame_motion()
{
	frame_motion.resize(frame_num-1);
	index_frame_motion.resize(frame_num-1);

	for (int i = 0; i < motion.size(); i++)
	{
		for (int j = 0; j < motion[i].size(); j++)
		{
			frame_motion[org_path[i].start_frame - 1 + j].push_back(motion[i][j]);
			index_frame_motion[org_path[i].start_frame - 1 + j].push_back(Point2i(j, i));
		}
	}
}

void inline Median_Filter::exchang(vector<Point2d>& input,int i,int j)
{
	Point2d temp;
	input[i] = temp;
	input[i] = input[j];
	input[j] = temp;
}

void Median_Filter::sort(vector<Point2d> input)
{
	for (int i = 0; i < input.size()/2; i++)
	{
		for (int j = i+1; j < input.size(); j++)
		{
			if ((input[i].x*input[i].x + input[i].y*input[i].y) < (input[j].x*input[j].x + input[j].y*input[j].y))
			{
				exchang(input, i, j);
			}
		}
	}
}

void Median_Filter::median_filter()
{
	vector<Point2d> filter_window;

	get_motion();

	find_frame_motion();
	
	for (int iteration = 0; iteration < 1; iteration++)
	{
		for (int i = 0; i < frame_motion.size(); i++)
		{
			for (int j = 0; j < frame_motion[i].size(); j++)
			{

				filter_window.clear();

			
				if (j < window / 2)
				{
					for (int m = 0; m < window / 2 - j; m++)
					{
						filter_window.push_back(frame_motion[i][frame_motion[i].size() - 1 - m]);
					}

					for (int m = 0; m < (window / 2) + 1 + j; m++)
					{
						filter_window.push_back(frame_motion[i][m]);
					}
				}

				else if (j > frame_motion[i].size() - window / 2)
				{

					for (int m = 0; m < window / 2 - (frame_motion[i].size() - 1 - j); m++)
					{

						filter_window.push_back(frame_motion[i][m]);
					}

					for (int m = 0; m < window / 2 + (frame_motion[i].size() - 1 - j); m++)
					{
						filter_window.push_back(frame_motion[i][j - window / 2 + m]);
					}

				}

		

				else
				{
					for (int m = 0; m < window; m++)
					{
						filter_window.push_back(frame_motion[i][j - window / 2]);
					}
				}

				sort(filter_window);

				
				Point2d median = filter_window[filter_window.size() / 2];



				int path = index_frame_motion[i][j].y;
				int point = index_frame_motion[i][j].x;

				org_path[path].PointPath[point + 1][0] = org_path[path].PointPath[point][0] + median.x;
				org_path[path].PointPath[point + 1][1] = org_path[path].PointPath[point][1] + median.y;

			

				frame_motion[i][j] = median;






			}
		}
	}

}