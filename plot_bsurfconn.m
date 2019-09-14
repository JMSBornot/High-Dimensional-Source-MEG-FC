function plot_bsurfconn(surf, anatcol, Conn, Opt)
% Example 1:
% addpath('E:\MATLAB\fieldtrip-20160521\external\freesurfer')
% Opt.nodesize = 8;
% Opt.nodestyle = 'ok';
% Opt.nodefacecolor = 'k';
% Opt.edgesize = 2;
% Opt.edgestyle = '-';
% Opt.edgecolor = 'b';
% dirname = 'E:\DATOS\Template\SPM\';
% lhcurv = read_curv([dirname 'surf\lh.curv']);
% rhcurv = read_curv([dirname 'surf\rh.curv']);
% anatcol = sigmf(-[lhcurv;rhcurv],[10 0])*ones(1,3);
% [LH.pnt, LH.tri] = read_surf([dirname 'surf\lh.white']);
% [RH.pnt, RH.tri] = read_surf([dirname 'surf\rh.white']);
% LH.tri = LH.tri + 1;
% RH.tri = RH.tri + 1;
% nl = size(LH.pnt,1);
% LH.pntcol = anatcol(1:nl,:);
% RH.pntcol = anatcol(nl+1:end,:);
% np = length(LH.pnt) + length(RH.pnt);
% Conn = sparse([100000;250000], [200000;300000], ones(2,1), np, np);
% LH.Conn = Conn(1:nl, 1:nl);
% RH.Conn = Conn(nl+1:end, nl+1:end);
% figure; plot_bsurfconn(LH, LH.pntcol, LH.Conn, Opt);
% figure; plot_bsurfconn(RH, RH.pntcol, RH.Conn, Opt);
% surf.pnt = [LH.pnt; RH.pnt];
% surf.tri = [LH.tri; RH.tri + nl];
% figure; plot_bsurfconn(surf, anatcol, Conn, Opt);

if isfield(surf, 'face')
    face = surf.face;
elseif isfield(surf, 'tri')
    face = surf.tri;
end
if isfield(surf, 'pnt')
    pnt = surf.pnt;
elseif isfield(surf, 'pos')
    pnt = surf.pos;
elseif isfield(surf, 'vert')
    pnt = surf.vert;
end
Nvert = length(pnt);
if isempty(anatcol) || isempty(Conn)
    error('You should provide anatomical and connectivity data');
end
% check options
if isfield(Opt, 'alpha') 
    alpha = Opt.alpha;
else
    alpha = 0.5;
end
if isempty(alpha)
    alpha = 0.5;
end
if issymmetric(Conn)
    Conn = triu(Conn,1);
end
ncateg = max(max(Conn));
if (length(Opt.edgesize) == 1) % size
    indcateg_edgesize = ones(ncateg,1);
elseif (length(Opt.edgesize) < ncateg)
    error('This value should be at least equals to the number of categories');
else
    indcateg_edgesize = 1:ncateg;
end
if (length(Opt.edgestyle) == 1) % style
    indcateg_edgestyle = ones(ncateg,1);
    Opt.edgestyle = cellstr(Opt.edgestyle);
elseif (length(Opt.edgestyle) < ncateg)
    error('This value should be at least equals to the number of categories');
else
    indcateg_edgestyle = 1:ncateg;
end
if (length(Opt.edgecolor) == 1) % color
    indcateg_edgecolor = ones(ncateg,1);
    Opt.edgecolor = cellstr(Opt.edgecolor);
elseif (length(Opt.edgecolor) < ncateg)
    error('This value should be at least equals to the number of categories');
else
    indcateg_edgecolor = 1:ncateg;
end
% create connection lines
nedges = nnz(Conn);
xline = NaN(2,nedges);
yline = NaN(2,nedges);
zline = NaN(2,nedges);
[ii,jj,val] = find(Conn);
indnodes = unique([ii;jj]);
xyzcent = pnt(indnodes,:);
for it = 1:nedges
    ind = [find(indnodes == ii(it)); find(indnodes == jj(it))];
    xline(:,it) = xyzcent(ind,1);
    yline(:,it) = xyzcent(ind,2);
    zline(:,it) = xyzcent(ind,3);
end
% plot anatomical surface
if (isvector(anatcol) && length(anatcol)==3) % RGB color
    patch('faces', face, 'vertices', pnt, 'FaceAlpha', alpha, ...
        'EdgeColor', 'none', 'FaceColor', anatcol);
elseif any(size(pnt) ~= size(anatcol)) % one RGB color per point (vertice)
    error('Provide one RBG code for anatomy surface coloring or one RGB code (rows) per vertex');
else
    patch('faces', face, 'vertices', pnt, 'FaceVertexCData', anatcol, ...
        'FaceAlpha', alpha, 'EdgeColor', 'none', 'FaceColor', 'interp');
end
hold on
for it = 1:ncateg % draw edges
    flag = (val == it);
    plot3(xline(:,flag), yline(:,flag), zline(:,flag), ...
        'LineStyle', Opt.edgestyle{indcateg_edgestyle(it)}, ...
        'LineWidth', Opt.edgesize(indcateg_edgesize(it)), ...
        'Color', Opt.edgecolor{indcateg_edgecolor(it)});
end
plot3(xyzcent(:,1), xyzcent(:,2), xyzcent(:,3), Opt.nodestyle, ... % draw nodes
    'MarkerFaceColor', Opt.nodefacecolor, 'MarkerSize', Opt.nodesize);
hold off
axis image
view([0 0  1]); [AZ,EL]=view; lightangle(AZ,EL); lighting GOURAUD;
view([0 0 -1]); [AZ,EL]=view; lightangle(AZ,EL); lighting GOURAUD;
view([1 0 0]);