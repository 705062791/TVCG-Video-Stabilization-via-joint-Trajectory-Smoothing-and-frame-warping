# TVCG-Video-Stabilization-via-joint-Trajectory-Smoothing-and-frame-warping
Source code for "Effective Video stabilization via joint trajectory smoothing and Frame warping" on TVCG
  
--------
  
## 1.Environment
 
  - `window10` 
  - `vs2015`  
  
## 2.You need

  - `opencv 3.4.2`   
  - `matlab R2015b`  
  - `Eigen`   
  - `libivcore`  
  - `c++`  
  
  - We have provided `Eigen` and `libivcore` in `\lib`
  
## 3.How to start  

  This program was divided into two parts, the `\Trajectory tracking` and `\Video Smooth`. Each part can run independently
  
  ### Trajectory tracking
  
  `input`: Source Video (source_video.mp4)   
  `output`: trajectorys （videoXX_path,txt）  
  
1. open `\Trajectory tracking\my_tracker.m`  
  
2. Modify some main variables  
    - `caseFile` :The folder where the source video is located   
    - `saveFile` :A folder for storing trajectory date 
    - `casename`:Name of your source video
       
3. running  
   
### Video smooth  
   
   `input`:trajectorys （videoXX_path,txt)  
   `output`:stable video (stable_video_XX.mp4)
   
1. Open the project in `/Video smooth/` by visual studio, and Adding dependencies to projects
   
2. Modify some main variables
   
   - `video_name`: Name of your source video
       
   -  `video_file`: The folder where the source video is located
       
   -  `path_file`: A folder for storing trajectory data from the previous step
       
   -  `save_path`: A folder for storing stable video results
       
   -  {`alpha_1`,`alpha_2`,`Beta_1`,`Beta_2`,`gamma`}: Five parameters of energy equation (detailed reference paper)
       
3. running
### Ultimately, you'll get stable videos in the specified folder（save_path）
   
----

### Please visit [http://nieyongwei.net/](http://nieyongwei.net/) for more information
### Please cite paper:  
    T. Ma, Y. Nie, Q. Zhang, Z. Zhang, H. Sun and G. Li, "Effective Video Stabilization via Joint Trajectory Smoothing and Frame Warping," in IEEE Transactions on Visualization and Computer Graphics.
    doi: 10.1109/TVCG.2019.2923196
    keywords: {Trajectory;Smoothing methods;Cameras;Two dimensional displays;Three-dimensional displays;Feature extraction;Streaming media;Video stabilization;trajectory smoothing;mesh warping;optimization},
    URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8737754&isnumber=4359476
