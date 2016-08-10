# Spiro2SVG  
Convert spirographs and roulettes to SVG format using Bezier curves.  
<center><a href="https://github.com/alvinpenner/Spiro2SVG/tree/master/svg">
<img src="http://vaxxine.com/apenner/images/01_Multiple222.svg"
 width="24%" height="24%" title="01_Multiple.svg" />
<img src="http://vaxxine.com/apenner/images/decor3222.svg"
 width="24%" height="24%" title="decor3.svg" />
<img src="http://vaxxine.com/apenner/images/meduza222.svg"
 width="24%" height="24%" title="meduza.svg" />
<img src="http://vaxxine.com/apenner/images/Farris_1_7_-17_222.svg"
 width="24%" height="24%" title="Farris_1_7_-17.svg" /></a></center>

The Bézier curves are calculated by matching both the slope and the curvature of the spirograph at the endpoints, which leads to an accurate fit using very few points.  
The input spirograph data files can come from four sources:  
<ul>
<li>program Spirograph v-1.0.2.1, available at:
<a target="_blank" href="http://mathiversity.com/downloads">http://mathiversity.com</a>.<br>
A gallery of these images is at:
<a target="_blank" href="http://mathiversity.com/online-spirograph/gallery">http://mathiversity.com/.../gallery</a>.</li>
<li>program SpiroJ, available at:
<a target="_blank" href="http://sourceforge.net/projects/spiroj/">http://sourceforge.net/projects/spiroj</a>.</li>
<li>download file ShapesDemo.zip, at:
<a target="_blank" href="http://www.codeproject.com/Articles/76878/Spirograph-Shapes-WPF-Bezier-Shapes-from-Math-Form">Spirograph-Shapes-WPF-Bezier-Shapes</a>.<br>
see also Farris Wheels at: <a target="_blank" href="http://scholarcommons.scu.edu/cgi/viewcontent.cgi?article=1004&context=math_compsci">Wheels on Wheels on Wheels - Surprising Symmetry</a>.</li>
<li>see also the .spiro and .xml files in the Spiro2SVG/Samples folder.</li>
</ul>
to run from the gui: double-click on dist/Spiro2SVG.jar  
from the command line, use: java -jar Spiro2SVG.jar -?  
