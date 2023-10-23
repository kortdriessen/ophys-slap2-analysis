function hLine= shadedErrorPlot(x,y,e, colors, hAx)
if nargin<5
    hAx = axes;
end

assert(isvector(x));
x = x(:);
npoints = length(x);
szy = size(y);
if szy(2)==npoints
    y= y';
    e = e';
end
if size(e,2)==1 && size(y,2)>1
    e = repmat(e, 1, size(y,2));
end
if size(colors,1)==1
    colors = repmat(colors, size(y,2), 1);
end
assert(all(size(e)==size(y)));

for lnum = 1:size(y,2)
    cc = colors(lnum,:);
    yy = y(:,lnum);
    ee = e(:,lnum);
    patch([x ; flipud(x)], [yy+ee ; flipud(yy-ee)], cc, 'edgecolor', 'none', 'facealpha', 0.5, 'parent', hAx)
    hold on
    plot(hAx, x,y, 'color', cc, 'linewidth', 2)
end