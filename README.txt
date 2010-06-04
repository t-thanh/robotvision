For an installation instruction, please look at "INSTALL.txt".

If RobotVision is already install, you can run the demo application
"rss2010_demo" in order to recreate the figures of the paper 

>H. Strasdat, J.M.M. Montiel, A.J. Davision: "Scale Drift-Aware
 Large Scale Monocular SLAM", Proc. of Robotics: Science and
 Systems (RSS), 2010.<

"./rss2010_demo 2" creates Figure 2. Here, the Keble College trajectory before
loop closure is first read from a test file (black curve), then 6 DoF
optimisation (green curve) and 7 DoF optimisation is performed (red curve).

"./rss2010_demo 3" creates Figure 3. Attention, this will take one to two
hours! However, the whole process is visualised.

"./rss2010_demo 4" creates Figure 4. This is a simulation of an aircraft 
flying over a sphere. It performs monocular visual odometry using a 
downward looking camera. In the end, several loop closure constraints 
are identified and 7 DoF optimisation is performed.


Version 1.0 of RobotVision offers "Bundle Adjustment" (bundle_adjuster.h) 
and "Pose-graph Optimisation" (graph_optimizer.h) among other things. 

The core code of BundleAdjuster is in the method "calcFull". This method is 
rather complex because of the exploitation of its underlying sparseness.
For a more compact implementation of Levenberg-Marquardt, please refer to
"calcStructOnly" and "calcMotionOnly".

BundleAduster and GraphOptimiser are highly generic. They can be used for
various different kind of problems. For instance BundleAdjuster generalises
over different prediction functions using different kinds of transformations, 
different landmark types and different observation dimensions.
An example prediction SE3XYZ which depends of 3D rigid transformation SE3,
3D Euclidean points and 2D observation is provided in "transformations.h".
This could be seen as the default prediction for monocular bundle adjustment.
Writing a different prediction class (as a  sub-class of AbstractPrediction)
-- such as for stereo-camera SLAM or 2D range-bearing SLAM - should be 
straight-forward.

