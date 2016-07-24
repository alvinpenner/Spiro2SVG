# Spiro2SVG  
Convert spirographs and roulettes to SVG format using Bezier curves.  
<center><a href="https://github.com/alvinpenner/Spiro2SVG/tree/master/svg">
<img src="https://github.com/alvinpenner/Spiro2SVG/tree/master/svg/01_Multiple.svg"
 width="25%" height="25%" title="01_Multiple.svg" />
<img src="https://github.com/alvinpenner/Spiro2SVG/tree/master/svg/decor3.svg"
 width="25%" height="25%" title="decor3.svg" />
<img src="https://github.com/alvinpenner/Spiro2SVG/tree/master/svg/meduza.svg"
 width="25%" height="25%" title="meduza.svg" /></a></center>

The Bézier curves are calculated by matching both the slope and the curvature of the spirograph at the endpoints, which leads to a very accurate fit using very few curves.  
The input spirograph data files are from the program Spirograph v-1.0.2.1, available at:  
http://mathiversity.com/online-spirograph,  
or from the program SpiroJ, available at:  
http://sourceforge.net/projects/spiroj.  
(also see the .spiro and .xml files in the Samples folder).  

to run from the gui: double-click on dist/Spiro2SVG.jar  
from the command line, use: java -jar Spiro2SVG.jar -?  
