function hax = plot_bsurfconn_spm(anatfunc, code, Conn, Opt)
% See also: plot_bsurfconn_template

% Default settings
if ~exist('Opt','var')
    Opt = [];
end
if ~isfield(Opt,'nodesize'), Opt.nodesize = 8; end
if ~isfield(Opt,'nodestyle'), Opt.nodestyle = 'ok'; end
if ~isfield(Opt,'nodefacecolor'), Opt.nodefacecolor = 'k'; end
if ~isfield(Opt,'edgesize'), Opt.edgesize = 2; end
if ~isfield(Opt,'edgestyle'), Opt.edgestyle = '-'; end
if ~isfield(Opt,'edgecolor'), Opt.edgecolor = 'b'; end
if ~isfield(Opt,'alpha'), Opt.alpha = 0.5; end

% --- Check connectivity matrix
if size(Conn,1) ~= size(Conn,2)
    error('The connectivity matrix should be square.');
end
% only consider the upper triangular part if the matrix is symmetric
if issymmetric(Conn)
    Conn = triu(Conn,1);
end
% re-write Conn matrix (if only specified for Opt.indsel indices) as sparse
% for the whole mesh space as determined by Opt.Nvert points of the mesh.
if isfield(Opt, 'indsel')
    if (length(Opt.indsel) ~= length(Conn))
        error('The number of selected indices has to be equal the Conn matrix length.');
    end
    if isfield(Opt,'Nvert')
        Nvert = Opt.Nvert;
    else
        error(['The length for the connectivity matrix should be specified ' ...
            'corresponding to the mesh that contains the source elements.']);
    end
    % --- Create the sparse Conn matrix in the highres mesh
    [ii,jj,kk] = find(Conn);
    ii = Opt.indsel(ii);
    jj = Opt.indsel(jj);
    Conn = sparse(ii,jj,kk,Nvert,Nvert);
else
    Nvert = length(Conn);
end
if (Nvert == 5124)
    mesh = spm_eeg_inv_mesh([], 1);
elseif (Nvert == 8196)
    mesh = spm_eeg_inv_mesh([], 2);
elseif (Nvert == 20484)
    mesh = spm_eeg_inv_mesh([], 3);
else
    error('The number of vertices must correspond to those is one of the three SPM templates.');
end
surf = export(gifti(mesh.tess_ctx), 'ft');
% setup left and right hemispheres struct
nl = Nvert/2;
flag = all(surf.tri <= nl, 2);
LH.tri = surf.tri(flag,:);
ind = unique(LH.tri(:));
LH.pnt = surf.pnt(ind,:);
RH.tri = surf.tri(~flag,:);
ind = unique(RH.tri(:));
RH.pnt = surf.pnt(ind,:);
RH.tri = RH.tri - nl;
LH.Conn = Conn(1:nl, 1:nl);
RH.Conn = Conn(nl+1:end, nl+1:end);
% set anatomical space mesh and read curvature for coloring it.
switch anatfunc
    case 'brain'
        % ignore code value
        anatcol = [238 203 193]/255;
    case 'uniform'
        % the RGB color is in code
        if ~(isvector(code) && length(code)==3)
            error('Code should be an RGB code value.');
        end
        anatcol = code;
    otherwise
        anatcol = [238 203 193]/255;
end
if isempty(anatcol)
    error('unexpected');
elseif isvector(anatcol)
    LH.anatcol = anatcol;
    RH.anatcol = anatcol;
else %  ismatrix(anatcol) 
    LH.anatcol = anatcol(1:nl,:);
    RH.anatcol = anatcol(nl+1:end,:);
end
% template surface has RAS orientation
hsep = 30;
vsep = 20;
range = max(max(LH.pnt),max(RH.pnt)) - min(min(LH.pnt),min(RH.pnt));
horz = [hsep range(1) hsep range(1) hsep range(2) hsep range(2) hsep];
horz = horz/sum(horz);
vert = [4*vsep range(3) vsep range(3) vsep];
vert = vert/sum(vert);
hpos = cumsum(horz);
vpos = cumsum(vert);
lag = 0.4*min(horz(1), vert(1));
% front view
hax = zeros(1,8);
hax(1) = axes('position', [hpos(1)-lag vpos(3)-lag horz(2)+2*lag vert(2)+2*lag]);
plot_bsurfconn(surf, anatcol, Conn, Opt); view([ 0  1  0]); axis off;
% posterior view
hax(2) = axes('position', [hpos(1)-lag vpos(1)-lag horz(2)+2*lag vert(2)+2*lag]);
plot_bsurfconn(surf, anatcol, Conn, Opt); view([ 0 -1  0]); axis off;
% superior view
hax(3) = axes('position', [hpos(3) vpos(3) horz(2) vert(2)]);
plot_bsurfconn(surf, anatcol, Conn, Opt); view([ 0  0  1]); axis off;
% inferior view
hax(4) = axes('position', [hpos(3) vpos(1) horz(2) vert(2)]);
plot_bsurfconn(surf, anatcol, Conn, Opt); view([ 0  0 -1]); axis off;
% LH - left view
hax(5) = axes('position', [hpos(5)-lag vpos(3)-lag horz(6)+2*lag vert(2)+2*lag]);
plot_bsurfconn(LH, LH.anatcol, LH.Conn, Opt); view([-1 0 0]); axis off;
% LH - right view
hax(6) = axes('position', [hpos(5)-lag vpos(1)-lag horz(6)+2*lag vert(2)+2*lag]);
plot_bsurfconn(LH, LH.anatcol, LH.Conn, Opt); view([ 1 0 0]); axis off;
% RH - right view
hax(7) = axes('position', [hpos(7)-lag vpos(3)-lag horz(6)+2*lag vert(2)+2*lag]);
plot_bsurfconn(RH, RH.anatcol, RH.Conn, Opt); view([ 1 0 0]); axis off;
% RH - left view
hax(8) = axes('position', [hpos(7)-lag vpos(1)-lag horz(6)+2*lag vert(2)+2*lag]);
plot_bsurfconn(RH, RH.anatcol, RH.Conn, Opt); view([-1 0 0]); axis off;