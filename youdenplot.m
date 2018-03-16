function varargout=youdenplot(data,varargin)
% Youden Plot - Plot a Youden plot
% The Youden plot is a graphical method to analyse inter-laboratory data,
% where all laboratories have analysed 2 samples. The plot visualises
% within-laboratory variability as well as between-laboratory variability.  
% For the original Youden plot (Youden WJ (1959) Graphical diagnosis of
% interlaboratory test results. Industrial Quality Control, 15, 24-28.) the
% two  samples must be similar and reasonably close in the magnitude of the
% property evaluated.  
% The axes in this plot are drawn on the same scale: one unit on the x-axis
% has the same length as one unit on the y-axis. 
% Each point in the plot corresponds to the results of one laboratory and
% is defined by a first response variable on the horizontal axis (i.e. run
% 1 or product 1 response value) and a second response variable 2 (i.e.,
% run 2 or product 2 response value) on the vertical axis.   
% A horizontal median line is drawn parallel to the x-axis so that there
% are as many points above the line as there are below it. A second median
% line is drawn parallel to the y-axis so that there are as many points on
% the left as there are on the right of this line. Outliers are not used in
% determining the position of the median lines. The intersection of the two
% median lines is called the Manhattan median.     
% A circle is drawn that should include a percent (usually 95%) of the
% laboratories if individual constant errors could be eliminated. 
% A 45-degree reference line is drawn through the Manhattan median.
% Moreover, there are two tangents to the circle and parallel to the 45°
% line.
% 
% Interpretation
% Points that lie into the circle: only random error
% Points that lie outside the circle but inside the tangents: systematic error
% Points that lie near the 45-degree reference: very precise results 
% Points that lie near the 45-degree reference but outside the circle: very precise results but systematic error
% Points that lie outside the tangents: gross errors
% 
% Syntax: Out=youdenplot(data,group,alpha,verbose)
% 
% Input: X - This is a Nx2 data matrix. This input is mandatory
%        Group - specifies one or more grouping variables G, producing a
%        separate scatter plot for each set of X values sharing the same G
%        value or values. Grouping variables must have one column per element
%        of X. Specify a single grouping variable in G by using a vector.
%        By default, G=1:1:N;
%        ALPHA - significance level (default 0.05)
%        VERBOSE - if you want to see report (0-no; 1-yes by default);
% 
% Output: if verbose = 0
%         the Youden plots
%         if verbose = 1
%         A table with Total, Random, Systematic errors, and if the point
%         is within the circle, within the tangents or outside.
%         
%         if Out is declared, you will have a struct:
%         Out.m=manhattan media;
%         Out.s=standard deviations;
%         OUT.confidence=(1-alpha)*100 confidence 
%         Out.r=circle radius;
%         Out.stats=[Total_error Random_Error Systematic_error In_circle In_tangent Out_of_all];
% 
% 
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
% 
% To cite this file, this would be an appropriate format:
% Cardillo G. (2014) Youden's Plot: compute the Youden's plot for laboratories variability
% http://www.mathworks.com/matlabcentral/fileexchange/48039

%Input error handling
p = inputParser;
addRequired(p,'data',@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','nonempty','ncols',2}));
addOptional(p,'g',[],@(x) isempty(x) || (~isnan(x) && isnumeric(x) && isreal(x) && isfinite(x) && isrow(x)));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'verbose',1, @(x) ~isnan(x) && isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==1 || x==0));
parse(p,data,varargin{:});
data=p.Results.data; g=p.Results.g; alpha=p.Results.alpha; verbose=p.Results.verbose; 
if isempty(g)
    g=1:1:length(data);
else
    assert(size(data,1)==length(g),'Warning: data and g must have the same length')
end
clear p

g=g(:); m=median(data); dm=diff(m); n=length(data); mm=repmat(m,n,1);
colors=distinguishable_colors(max(g));
%Pythagora's theorem
ERR=realsqrt(sum((data-mm).^2,2)); %total magnitude of error (distance of the data point from the manhattan median)
int45=(repmat(sum(data,2),1,2)+[mm(:,1)-mm(:,2) +mm(:,2)-mm(:,1)])./2; %intercept coordinates of perpendicular from datapoints to 45° reference line
ran=realsqrt(sum((data-int45).^2,2)); %Random Distances from Points to Intercepts
sys=realsqrt(sum((mm-int45).^2,2)); %Systematic Distances from Origin to Intercepts 
clear mm int45
%The components of systematic and random error must be fitted on to the
%total error
c=ERR./(sys+ran);
SYS=sys.*c; SYS(sys==0)=0;
RAN=ran.*c; RAN(ran==0)=0;
clear c sys ran 
%The radius of the circle is based on the standard deviation of the lab
%random errors.  This standard deviation of the random errors is then
%multiplied by a confidence factor that should make the circle contain
%(1-alpha)% of all points IF the systematic errors were eliminated.   
%The radius of the circle is T times Student's t for a two tailed (1-alpha)
%confidence interval with n-1 degrees of freedom.
df=n-1; %degrees of freedom
s=realsqrt(sum(RAN.^2)/df); %Standard deviation of random error
r=s*tinv(1-alpha/2,df); %Radius of the circle
clear df
%Youden's Plot
hFig=figure; %Create a maximized figure window
set(hFig,'units','normalized','outerposition',[0 0 1 1]);
set(hFig, 'Color', 'white'); % sets the color to white
axis equal
hold on
%plot the circle
%x=x0+r*cos(theta)   y=y0+r*sin(theta)
t=linspace(0,2*pi,500); ct=cos(t); st=sin(t);
C0=plot(m(1)+r.*ct,m(2)+r.*st,'g-','Linewidth',2);
%plot the points
gscatter(data(:,1),data(:,2),g,colors,[],[],'off');
clear colors t st ct
%Find the correct limits for the square plot
Ax=get(gca,'XLim'); Ay=get(gca,'YLim');
%Now plot the 45% reference lines
f=@(x,m,k)x+dm+k; %equation of the 45° sheaf of lines
%the line that pass through the centre of the circle has k=0
C1=plot(Ax,feval(f,Ax,m,0),'Color',[255 153 0]./255,'Linewidth',2);
%the upper tangent parallel to 45° reference line has k=r*realsqrt(2)
%the lower tangent parallel to 45° reference line has k=-r*realsqrt(2)
k=r*realsqrt(2);
plot(Ax,feval(f,Ax,dm,k),'Color','r','Linewidth',2)
plot(Ax,feval(f,Ax,dm,-k),'Color','r','Linewidth',2)
set(gca,'Xlim',Ax,'Ylim',Ay)
%plot the axes passing through the Manhattan median
C2=plot([m(1) m(1)],Ay,'k-','Linewidth',2);
plot(Ax,[m(2) m(2)],'k-','Linewidth',2)
hold off
clear Ax Ay f 
legend([C0 C1 C2],[num2str((1-alpha)*100) '% circle of random error'],...
     'Precision line',['Manhattan median: ' num2str(m)],...
     'Location','NorthEastOutside');
set(gca,'FontSize',14)
title('Youden''s plot','FontSize',16);
xlabel('First measure','FontSize',16);
ylabel('Second measure','FontSize',16);
clear C0 C1 C2
Y.m=m; Y.s=s; Y.confidence=(1-alpha)*100; Y.r=r; 
Inc=ERR<=r; %points within the circle ERR<radius
Int=zeros(n,1); %points between the tangents
%find the parallels to 45° reference line passing through the
%datapoints outside the circle =>y=x+delta=x+yp-xp
delta=data(Inc==0,2)-data(Inc==0,1);
%the points between the tangents
Int(Inc==0)=delta<=dm+k & delta>=dm-k;
Out=~Inc & ~Int; %points outside tangents
Y.stats=[ERR RAN SYS Inc Int Out];
clear dm k m s alpha r n delta
if nargout
    varargout=Y;
end
if verbose==1
    table(g,data,ERR,RAN,SYS,Inc,logical(Int),Out,'VariableNames',{'G','Data','Total_error','Random_error','Systematic_error','Within_circle','Inside_tangents','Outside_tangents'})
end