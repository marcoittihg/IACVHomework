# IACVHomework
Homework for the Image Analysis and Computer Vision (IACV) course

For a more in-depth explanation read  [Report.pdf](Report.pdf)

## Input image
<img src="First.jpg" width="70%" class="center">

### Features:
##### 1) Edges detection.
<img src="Images/CannyImp1.jpg" width="45%"> <img src="Images/CannyImp2.jpg" width="45%"> <img src="Images/CannyImp3.jpg" width="45%"> 
##### 2) Corners detection.
<img src="Images/HarrisImp1.png" width="45%"> <img src="Images/HarrisImp2.png" width="45%"> <img src="Images/HarrisImp3.png" width="45%"> 
##### 3) Lines detection.
<img src="Images/LineDec1.png" width="45%"> <img src="Images/LineDec2.png" width="45%"> <img src="Images/LineDec3.png" width="45%"> 

### Geometry:
##### 4.1) Find two horizontal vanishing points and one vertical vanishing point.
<img src="Images/VanishingPoints.png" width="50%" class="center">

##### 4.2) Look at the image below and consider the horizontal section of facades 1, 2 and 6 , depicted in yellow: metrically reconstruct this horizontal section*, so as to determine the relative coordinates of features points indicated in blue. These point features are placed (i) at the intersections between facades 1 and 2 and between facades 6 and 7, and (ii) in correspondence of borders of the windows.
<img src="Images/MetricReconstruction.png" width="50%" class="center">

##### 4.3) Estimate the calibration matrix K of the camera. Assume the camera is zero-skew, but not natural.
<img src="Images/KValues.png" width="50%">

##### 4.4) Use the knowledge of K to rectify also a vertical facade, as, e.g., facade 1 or 4 or 2+6.
<img src="Images/VerticalRectification.png" width="35%" class="center">

##### 4.5) Fix a suitable reference frame attached to the building, and localize the camera relative to the fixed reference.
<img src="Images/PositionAndRotation.png" width="50%">