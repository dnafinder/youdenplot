# youdenplot
Plot a Youden plot<br/>
The Youden plot is a graphical method to analyse inter-laboratory data,
where all laboratories have analysed 2 samples. The plot visualises
within-laboratory variability as well as between-laboratory variability.  
For the original Youden plot (Youden WJ (1959) Graphical diagnosis of
interlaboratory test results. Industrial Quality Control, 15, 24-28.) the
two  samples must be similar and reasonably close in the magnitude of the
property evaluated.  
The axes in this plot are drawn on the same scale: one unit on the x-axis
has the same length as one unit on the y-axis. 
Each point in the plot corresponds to the results of one laboratory and
is defined by a first response variable on the horizontal axis (i.e. run
1 or product 1 response value) and a second response variable 2 (i.e.,
run 2 or product 2 response value) on the vertical axis.   
A horizontal median line is drawn parallel to the x-axis so that there
are as many points above the line as there are below it. A second median
line is drawn parallel to the y-axis so that there are as many points on
the left as there are on the right of this line. Outliers are not used in
determining the position of the median lines. The intersection of the two
median lines is called the Manhattan median.     
A circle is drawn that should include a percent (usually 95%) of the
laboratories if individual constant errors could be eliminated. 
A 45-degree reference line is drawn through the Manhattan median.
Moreover, there are two tangents to the circle and parallel to the 45°
line.

Interpretation
Points that lie into the circle: only random error
Points that lie outside the circle but inside the tangents: systematic error
Points that lie near the 45-degree reference: very precise results 
Points that lie near the 45-degree reference but outside the circle: very precise results but systematic error
Points that lie outside the tangents: gross errors

Syntax: Out=youdenplot(data,group,alpha,verbose)

Input: X - This is a Nx2 data matrix. This input is mandatory<br/>
       Group - specifies one or more grouping variables G, producing a
       separate scatter plot for each set of X values sharing the same G
       value or values. Grouping variables must have one column per element
       of X. Specify a single grouping variable in G by using a vector.
       By default, G=1:1:N;<br/>
       ALPHA - significance level (default 0.05)<br/>
       VERBOSE - if you want to see report (0-no; 1-yes by default);<br/>

Output: if verbose = 0
        the Youden plots
        if verbose = 1
        A table with Total, Random, Systematic errors, and if the point
        is within the circle, within the tangents or outside.
        
        if Out is declared, you will have a struct:
        Out.m=manhattan media;
        Out.s=standard deviations;
        OUT.confidence=(1-alpha)*100 confidence 
        Out.r=circle radius;
        Out.stats=[Total_error Random_Error Systematic_error In_circle In_tangent Out_of_all];


          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2014) Youden's Plot: compute the Youden's plot for laboratories variability
http://www.mathworks.com/matlabcentral/fileexchange/48039
