function varargout = youdenplot(data, varargin)
%YOUDENPLOT Youden's Plot for inter-laboratory variability.
%
%   Y = YOUDENPLOT(DATA) produces a Youden Plot for an Nx2 matrix DATA,
%   where each row corresponds to one laboratory and the two columns are
%   the responses on sample 1 and sample 2.
%
%   Y = YOUDENPLOT(DATA, G) also specifies a grouping vector G (one value
%   per row of DATA). Each distinct value of G corresponds to a laboratory
%   and is used only for coloring / labeling in the plot. By default,
%   G = 1:N.
%
%   Y = YOUDENPLOT(DATA, G, ALPHA) specifies the significance level ALPHA
%   for the confidence circle (default ALPHA = 0.05). The circle radius is
%   based on the standard deviation of the random error component, scaled
%   by the appropriate Student's t quantile for a two-tailed confidence
%   interval with N−1 degrees of freedom.
%
%   Y = YOUDENPLOT(DATA, G, ALPHA, VERBOSE) controls the table output:
%       VERBOSE = 0 → no table is displayed.
%       VERBOSE = 1 → a classification table is printed (default).
%
%   ------------------------------------------------------------------
%   Syntax:
%       Y = youdenplot(DATA)
%       Y = youdenplot(DATA, G)
%       Y = youdenplot(DATA, G, ALPHA)
%       Y = youdenplot(DATA, G, ALPHA, VERBOSE)
%
%   Inputs:
%       DATA    - Nx2 data matrix (real, finite, non-NaN).
%                 Column 1 = first measurement (sample/run 1)
%                 Column 2 = second measurement (sample/run 2)
%
%       G       - (optional) grouping vector, length N.
%                 Used only for colouring/identification of laboratories.
%                 Default: G = 1:N.
%
%       ALPHA   - (optional) significance level for the confidence circle:
%                 (1 − ALPHA)*100% of random-error-only points should lie
%                 inside the circle.
%                 Default: 0.05.
%
%       VERBOSE - (optional) flag:
%                 0 → do not display the summary table
%                 1 → display summary table (default).
%
%   Output:
%       Y       - structure with the following fields:
%                   Y.m          = [m1 m2] Manhattan median of DATA
%                   Y.s          = std. dev. of random error
%                   Y.confidence = (1 − ALPHA)*100
%                   Y.r          = circle radius
%                   Y.stats      = [Total Random Systematic InCircle InTang OutTang]
%                                 where:
%                                 Total     = total error (distance from median)
%                                 Random    = random component
%                                 Systematic= systematic component
%                                 InCircle  = logical flag (inside circle)
%                                 InTang    = logical flag (between tangents)
%                                 OutTang   = logical flag (outside tangents)
%
%   ------------------------------------------------------------------
%   Strategy:
%   ------------------------------------------------------------------
%   1) Compute the Manhattan median M = (m1, m2) of the two columns of DATA.
%
%   2) For each laboratory i:
%         - compute the total error as the Euclidean distance from M;
%         - decompose the error into:
%               * random component  (distance from point to its 45°-line
%                                   intercept)
%               * systematic component (distance from median M to that
%                                   intercept)
%
%   3) Estimate the standard deviation S of the random error component and
%      derive the circle radius:
%               r = S * tinv(1 − ALPHA/2, N−1)
%
%   4) Draw:
%         - the confidence circle centered at M with radius r,
%         - the 45° precision line through M,
%         - two tangents to the circle parallel to the 45° line,
%         - the vertical and horizontal lines through M.
%
%      Classify each laboratory as:
%         - inside the circle      → mostly random error
%         - outside circle but
%           between tangents       → systematic error
%         - outside tangents       → gross or out-of-control error
%
%   ------------------------------------------------------------------
%   Author and citation:
%   ------------------------------------------------------------------
%   Created by:  Giuseppe Cardillo
%   E-mail:      giuseppe.cardillo.75@gmail.com
%
%   To cite this file:
%   Cardillo G. (2014) Youden''s Plot: compute the Youden''s plot for
%   laboratories variability.
%
%   GitHub: https://github.com/dnafinder/youdenplot
%   GitHub namespace: https://github.com/dnafinder/
%
%   ------------------------------------------------------------------

%% Input parsing
p = inputParser;
p.FunctionName = 'youdenplot';

addRequired(p,'data',@(x) validateattributes( ...
    x,{'numeric'},{'real','finite','nonnan','nonempty','ncols',2}));

addOptional(p,'g',[],@(x) isempty(x) || ...
    (isnumeric(x) && isvector(x) && all(isfinite(x(:)))));

addOptional(p,'alpha',0.05,@(x) validateattributes( ...
    x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));

addOptional(p,'verbose',1,@(x) isnumeric(x) && isscalar(x) && ismember(x,[0 1]));

parse(p,data,varargin{:});
data    = p.Results.data;
g       = p.Results.g;
alpha   = p.Results.alpha;
verbose = p.Results.verbose;

%% Group handling
n = size(data,1);
if isempty(g)
    g = (1:n).';
else
    g = g(:);
    assert(numel(g)==n,'DATA and G must have the same number of rows.');
end

%% Manhattan median
m  = median(data,1);
dm = diff(m);
mm = repmat(m,n,1);

%% Required dependency
assert(exist('distinguishable_colors','file')==2, ...
    'youdenplot:MissingDependency', ...
    'distinguishable_colors.m must be on the MATLAB path.');

colors = distinguishable_colors(max(g));

%% Error components
% Total error: distance from Manhattan median
ERR = realsqrt(sum((data-mm).^2,2));

% Intercepts on the 45° reference line
sumd  = sum(data,2);
int45 = (repmat(sumd,1,2) + ...
        [mm(:,1)-mm(:,2), mm(:,2)-mm(:,1)]) ./ 2;

% Random and systematic distances
ran = realsqrt(sum((data-int45).^2,2));  % point → intercept
sys = realsqrt(sum((mm-int45).^2,2));    % median → intercept

% Fit random & systematic to total error
c   = ERR ./ (sys + ran);
SYS = sys .* c;  SYS(sys==0) = 0;
RAN = ran .* c;  RAN(ran==0) = 0;

%% Circle radius
df = n - 1;
s  = realsqrt(sum(RAN.^2) / df);
r  = s * tinv(1 - alpha/2, df);

%% Plot
hFig = figure;
set(hFig,'Units','normalized','OuterPosition',[0 0 1 1],'Color','white');
axis equal; hold on;

% Circle
t  = linspace(0,2*pi,500);
ct = cos(t); st = sin(t);
C0 = plot(m(1)+r.*ct, m(2)+r.*st,'g-','LineWidth',2);

% Points
gscatter(data(:,1), data(:,2), g, colors, [], [], 'off');

% Axis limits
Ax = get(gca,'XLim');
Ay = get(gca,'YLim');

% 45° sheaf: y = x + dm + k
f = @(x,dm,k) x + dm + k;

% Precision line (k = 0)
C1 = plot(Ax, f(Ax, dm, 0), 'Color', [255 153 0]/255, 'LineWidth',2);

% Tangents (k = ± r*sqrt(2))
k = r * realsqrt(2);
plot(Ax, f(Ax, dm,  k),'r','LineWidth',2);
plot(Ax, f(Ax, dm, -k),'r','LineWidth',2);

% Manhattan median axes
set(gca,'XLim',Ax,'YLim',Ay);
C2 = plot([m(1) m(1)],Ay,'k-','LineWidth',2);
plot(Ax,[m(2) m(2)],'k-','LineWidth',2);

hold off;

legend([C0 C1 C2], ...
    [num2str((1-alpha)*100) '% circle'], ...
    'Precision line', ...
    ['Manhattan median: ' num2str(m)], ...
    'Location','NorthEastOutside');

set(gca,'FontSize',14)
title('Youden''s plot','FontSize',16);
xlabel('First measure','FontSize',16);
ylabel('Second measure','FontSize',16);

%% Classification
InCircle = ERR <= r;

delta  = data(~InCircle,2) - data(~InCircle,1);
InTang = false(n,1);
InTang(~InCircle) = (delta <= dm + k) & (delta >= dm - k);

OutTang = ~InCircle & ~InTang;

%% Output struct
Y.m          = m;
Y.s          = s;
Y.confidence = (1 - alpha) * 100;
Y.r          = r;
Y.stats      = [ERR RAN SYS InCircle InTang OutTang];

if nargout
    varargout{1} = Y;
end

%% Verbose table
if verbose
    table(g, data, ERR, RAN, SYS, InCircle, logical(InTang), OutTang, ...
        'VariableNames', {'G','Data','Total_error','Random_error', ...
        'Systematic_error','Within_circle','Inside_tangents','Outside_tangents'})
end

end
