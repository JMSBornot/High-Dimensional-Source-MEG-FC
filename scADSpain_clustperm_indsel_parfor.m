% scADSpain_clustperm_indsel_parfor
warning off
dirname = '\\sceis_cl1fs\shared\PhD Students\JoseSanchez';
cd(dirname);

%% setting path to main tools
addpath('C:\WORK\MATLAB\Utiles');
addpath(genpath('C:\WORK\MATLAB\spm12'));
addpath('C:\WORK\MATLAB\fieldtrip-20160521\external\freesurfer');
addpath('C:\WORK\DATOS\Template\SPM\surf\');
addpath('C:\WORK\DATOS\Template\SPM\label');

%% control panel

% fcompute_block = true;
fcompute_block = false;

% itype = 1; % ranksum
% itype = 2; % Spearman's correlation (MMSE)
% itype = 3; % Spearman's correlation (IRM)
itype = 4; % Spearman's correlation (DRM)
% itype = 5; % correlation (MMSE)
% itype = 6; % correlation (IRM)
% itype = 7; % correlation (DRM)

% fcompute_stat = true;
fcompute_stat = false;

fvisualize_stat = true;

strmeth = {'ranksum' 'SpcorrMMSE' 'SpcorrIRM' 'SpcorrDRM' 'corrMMSE' 'corrIRM' 'corrDRM'};

%% Loading dataset header info and setting main variables
[ndata, fieldnames] = xlsread(fullfile(dirname, 'Resting_eyes_closed_30vs30.xlsx'));
id_subj = ndata(:, strcmp(fieldnames, 'ID_meg'));
Group = ndata(:,strcmp(fieldnames,'Group'));
GroupLabels = {'HC'; 'MCI'}; % corresponding to Group's codes 1 or 2
GroupLabels = GroupLabels(Group);
MMSE = ndata(:,strcmp(fieldnames,'MMSE'));
flagMMSE = ~isnan(MMSE);
IRM = ndata(:,strcmp(fieldnames,'Immediate recall'));
flagIRM = ~isnan(IRM);
DRM = ndata(:,strcmp(fieldnames,'Delayed recall'));
flagDRM = ~isnan(DRM);
bandlab = {'delta'; 'thetaI'; 'thetaII'; 'alphaI'; 'alphaII'; 'betaI'; 'betaII'; 'gamma'};
nband = length(bandlab);
isize = [500*ones(1,15) 696];
Ndip = sum(isize);
nFC = nchoosek(Ndip,2);
nblock = length(isize);
Niter = nchoosek(nblock+1,2);
nsubj = length(id_subj);
% setting matrix block subindices
itrow = zeros(Niter+1,1); %#ok<*UNRCH>
itcol = [nblock; zeros(Niter,1)];
for cont = 1:Niter
    if (itcol(cont) == nblock)
        itrow(cont+1) = itrow(cont) + 1;
        itcol(cont+1) = itrow(cont+1);
    else
        itrow(cont+1) = itrow(cont);
        itcol(cont+1) = itcol(cont) + 1;
    end
end
itrow(1) = [];
itcol(1) = [];

%% Compute indices for nonparametric cluster permutation statistics per blocks
if fcompute_block
    Nr = 1001;
    nX = nnz(strcmp(GroupLabels,'HC'));
    nY = nnz(strcmp(GroupLabels,'MCI'));
    M = nX + nY;
    rng('default'); % reset random number generator to guarantee replicability
    indperm = zeros(M, Nr);
    indperm(:,1) = 1:M;
    for it = 2:Nr
        indperm(:,it) = randperm(M);
    end
    if (itype == 2) || (itype == 5)
        y = MMSE(flagMMSE);
        indperm(indperm > nnz(flagMMSE)) = -50; % tag elements to remove with any negative value
    elseif (itype == 3) || (itype == 6)
        y = IRM(flagIRM);
        indperm(indperm > nnz(flagIRM)) = -50;
    elseif (itype == 4) || (itype == 7)
        y = DRM(flagDRM);
        indperm(indperm > nnz(flagDRM)) = -50;
    end
    indperm(indperm < 0) = []; % remove tagged elements and then reshape
    indperm = reshape(indperm, [], Nr);
    if (itype == 1)
        leftThresh = [579 603 629];
        rightThresh = [1251 1227 1201];
        % leftThresh = [710 705 700];
        % rightThresh = [1120 1110 1100];
    else
        leftThresh = [-0.62 -0.58 -0.53];
        rightThresh = [0.62  0.58  0.53];
    end
    disp('Computing nonparametric statistics...');
    parfor it = 1:Niter
        format compact
        disp([it Niter]);
        for itb = 1:nband
            data = cell(nsubj,1);
            for k = 1:nsubj
                tag = num2str(id_subj(k)); %#ok<*PFBNS>
                fname = fullfile(dirname, 'RESULTS', tag, sprintf('EIC_rn%dcn%d_%s.mat',itrow(it),itcol(it),bandlab{itb}));
                block = load(fname);
                data{k} = block.block;
            end
            data = double(cell2mat(data));
            ncol = size(data,2);
            if (itype == 1)
                X = data(strcmp(GroupLabels,'HC'),:);
                Y = data(strcmp(GroupLabels,'MCI'),:);
                % leftThresh = [710 705 700];
                % rightThresh = [1120 1110 1100];
                % tic; cellind = compute_indices_ranksum(X, Y, indperm, leftThresh, rightThresh); toc
                % Z = [X;Y]; t=2;
                % for ii=1:10
                %     [pval,H,STATS] = ranksum_mc(Z(indperm(1:30,ii),:),Z(indperm(31:end,ii),:));
                %     ind = find(STATS.ranksum>rightThresh(t));
                %     [max(abs(uint64(ind)-cellind{ii,2,t}-1)) length(ind)]
                % end
                cellind = compute_indices_ranksum(X, Y, indperm, leftThresh, rightThresh);
                fname = fullfile(dirname, 'RESULTS', 'STATS', sprintf('block_clustperm_indsel_rn%dcn%d_%s.mat',itrow(it),itcol(it),bandlab{itb}));
            else
                if (itype == 2) || (itype == 5)
                    data = data(flagMMSE,:);
                    fname = fullfile(dirname, 'RESULTS', 'STATS', sprintf('block_clustperm_%s_indsel_rn%dcn%d_%s.mat',strmeth{itype},itrow(it),itcol(it),bandlab{itb}));
                elseif (itype == 3) || (itype == 6)
                    data = data(flagIRM,:);
                    fname = fullfile(dirname, 'RESULTS', 'STATS', sprintf('block_clustperm_%s_indsel_rn%dcn%d_%s.mat',strmeth{itype},itrow(it),itcol(it),bandlab{itb}));
                elseif (itype == 4) || (itype == 7)
                    data = data(flagDRM,:);
                    fname = fullfile(dirname, 'RESULTS', 'STATS', sprintf('block_clustperm_%s_indsel_rn%dcn%d_%s.mat',strmeth{itype},itrow(it),itcol(it),bandlab{itb}));
                end
% leftThresh = [-0.50 -0.45 -0.40];
% rightThresh = [0.50  0.45  0.40];
% tic; cellind = compute_indices_correlation(data, y, indperm, 'Spearman', leftThresh, rightThresh); toc
% t=2;
% for ii=1:10
%      coeff = corr(data, y(indperm(:,ii)), 'Type', 'Spearman');
%      ind = find(coeff>rightThresh(t));
%      [max(abs(uint64(ind')-cellind{ii,2,t}-1)) length(ind)]
%      % ind = find(coeff<leftThresh(t));
%      % [max(abs(uint64(ind')-cellind{ii,1,t}-1)) length(ind)]
% end
                if ((itype >= 2) && (itype <= 4))
                    cellind = compute_indices_correlation(data, y, indperm, 'Spearman', leftThresh, rightThresh);
                else
                    cellind = compute_indices_correlation(data, y, indperm, [], leftThresh, rightThresh);
                end
            end
            save_cellind_ncol_Thresh(fname, cellind, ncol, leftThresh, rightThresh);
        end
    end
    
    keyboard;
end
clear cellind ncol;

%% creating auxiliar variables
% indices for mapping from FC indices to Ndip x Ndip connectivity matrix
IND = mat2cell(triu(ones(Ndip),1), isize, isize);
cont = 0;
for itr = 1:nblock
    for itc = itr:nblock
        cont = cont+1;
        IND{itr,itc} = cont*IND{itr,itc};
    end
end
IND = cell2mat(IND);
indtriu = cell(Niter, 1);
for it = 1:Niter
    indtriu{it} = find(IND == it);
end
indtriu = cell2mat(indtriu);
% auxiliar data structure
mesh = spm_eeg_inv_mesh([], 2);
surf = export(gifti(mesh.tess_ctx), 'ft');
TR = triangulation(surf.tri, surf.pnt);
E = edges(TR);
E = sortrows([E; fliplr(E)], [1 2]);
[~,ivertNeighBeg] = ismember(1:Ndip, E(:,1));
ivert2irow = repmat((1:Ndip)',1,Ndip);
ivert2irow = ivert2irow(indtriu);
ivert2icol = repmat(1:Ndip,Ndip,1);
ivert2icol = ivert2icol(indtriu);


%% Compute nonparametric cluster permutation statistics
if fcompute_stat
    accumind = zeros(1,Niter);
    for itb = 1:nband
        fprintf('Analyzing frequency band %d of %d ...\n', itb, nband);
        disp('Concatenating the indices of significant FC links along blocks...');
        Nindsel = 0; % counting the number of significant FC links for each MC resampling and threshold
        nbase = 0;
        for it = 1:Niter
            if (itype == 1)
                fname = fullfile(dirname, 'RESULTS', 'STATS', sprintf('block_clustperm_indsel_rn%dcn%d_%s.mat',itrow(it),itcol(it),bandlab{itb}));
            elseif (itype == 2) || (itype == 5)
                fname = fullfile(dirname, 'RESULTS', 'STATS', sprintf('block_clustperm_%s_indsel_rn%dcn%d_%s.mat',strmeth{itype},itrow(it),itcol(it),bandlab{itb}));
            elseif (itype == 3) || (itype == 6)
                fname = fullfile(dirname, 'RESULTS', 'STATS', sprintf('block_clustperm_%s_indsel_rn%dcn%d_%s.mat',strmeth{itype},itrow(it),itcol(it),bandlab{itb}));
            elseif (itype == 4) || (itype == 7)
                fname = fullfile(dirname, 'RESULTS', 'STATS', sprintf('block_clustperm_%s_indsel_rn%dcn%d_%s.mat',strmeth{itype},itrow(it),itcol(it),bandlab{itb}));
            else
                error('Unexpected');
            end
            load(fname);
            if (itb == 1)
                accumind(it) = ncol; % collecting number of FC links per block
                Nindsel = Nindsel + cellfun(@length, cellind);
            end
            if (itb == 1) && (it == 1)
                % initialization
                [Nr,~,Nth] = size(Nindsel);
                maxClusterExt = zeros(Nr,2,Nth,nband);
                clust_cell = cell(2,Nth,nband);
                clustext_cell = cell(2,Nth,nband);
                IndselBands_cell = cell(2,Nth,nband);
            end
            if (it == 1)
                tmp = num2cell(repmat(Niter, [Nr 2 Nth]));
                Indsel_cell = cellfun(@(N) cell(1,N), tmp, 'UniformOutput', false);
            end
            for it1 = 1:Nr
                for it2 = 1:2
                    for it3 = 1:Nth
                        % indsel = Indsel_cell{it1,it2,it3} + 1; % add 1 so indices start at 1 following Matlab convention
                        Indsel_cell{it1,it2,it3}{it} = nbase + cellind{it1,it2,it3} + 1;
                    end
                end
            end
            nbase = nbase + accumind(it);
        end
        Indsel_cell = cellfun(@cell2mat, Indsel_cell, 'UniformOutput', false); % concatenate the internal cell arrays
        for it2 = 1:2
            for it3 = 1:Nth
                IndselBands_cell{it2,it3,itb} = Indsel_cell{1,it2,it3};
            end
        end
        disp('Now, computing the clusters and the maximum cluster statistics...');
        for it1 = 1:Nr
            disp([it1 Nr]);
            for it2 = 1:2
                for it3 = 1:Nth
                    indsel = Indsel_cell{it1,it2,it3};
                    if ~isempty(indsel)
                        [clust, clust_ext] = compute_cluster(E, ivert2irow, ivert2icol, ivertNeighBeg, indsel);
                        if (it1 == 1)
                            clust_cell{it2,it3,itb} = clust;
                            clustext_cell{it2,it3,itb} = clust_ext;
                        end
                        maxClusterExt(it1,it2,it3,itb) = max(clust_ext);
                    end
                end
            end
        end
    end
    if (itype == 1)
        save clustperm_statistics_ranksum IndselBands_cell clust_cell clustext_cell maxClusterExt accumind Nth -v7.3
    else
        save(sprintf('clustperm_statistics_%s.mat',strmeth{itype}), 'IndselBands_cell', 'clust_cell', 'clustext_cell', 'maxClusterExt', 'accumind', 'Nth', '-v7.3');
    end
end

%% Summarize and plot results
if (itype == 1)
    load('clustperm_statistics_ranksum.mat');
else
    load(sprintf('clustperm_statistics_%s.mat', strmeth{itype}));
end
if (itype == 1)
    leftThresh = [579 603 629];
    rightThresh = [1251 1227 1201];
else
    leftThresh = [-0.62 -0.58 -0.53];
    rightThresh = [0.62  0.58  0.53];
end
ThMat = [leftThresh; rightThresh];
% summarize cluster stat across all bands and considering together the
% clusters separately obtained for the lower and upper tails
maxClusterExt_stat = squeeze(max(max(maxClusterExt,[],4),[],2));
% Threshold estimation
ThrClusterSize = zeros(1,Nth);
for it = 1:Nth
    ThrClusterSize(it) = prctile(maxClusterExt_stat(2:end,it),95);
end
disp('Threshold level estimated from the 95% percentile of the maximum cluster-permutation stat');
disp(ThrClusterSize);
% maximum cluster for our statistical analysis
stat = squeeze(maxClusterExt(1,:,:,:));
disp('Size of the maximum cluster statistics in our analysis, separately per band:')
for itb = 1:nband
    disp(bandlab{itb});
    disp(stat(:,:,itb));
end
% significant clusters
disp('Maximum cluster size of clusters that survive thresholding by cluster permutation statistics:')
for itb = 1:nband
    disp(bandlab{itb});
    disp(stat(:,:,itb).*(stat(:,:,itb) > ThrClusterSize));
end
% putting together all significant FC links for all frequency bands for a
% frequency band of interest (FBOI)

% itb = 1;
itb = 3;

Indsel_FClinks = cell(2,Nth);
for it1 = 1:2
    for it2 = 1:Nth
        indsel = [];
        ind = find(clustext_cell{it1,it2,itb} > ThrClusterSize(it2));
        for ic = ind
            indsel = union(indsel, IndselBands_cell{it1,it2,itb}(clust_cell{it1,it2,itb}==ic));
        end
        Indsel_FClinks{it1,it2} = indsel;
    end
end
% visualize the cluster(s) of FC links that survived the corresponding
% cluster-permutation statistical analysis for this FBOI
Optconn = [];
Optconn.alpha = 0.5;
Optconn.edgecolor = 'b';
strtype = {'Hypersynchronization', 'Hyposynchronization'};
for it1 = 1:2
    for it2 = 1:Nth
        indsel = Indsel_FClinks{it1,it2};
        if ~isempty(indsel)
            [ii,jj] = ind2sub([Ndip Ndip], indtriu(indsel));
            Conn = sparse(ii,jj,1,Ndip,Ndip);
            figure; plot_bsurfconn_spm('brain', [], Conn, Optconn);
            annotation('textbox', [0.2 0.05 0.6 0.05], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'normalized', ...
                'String', sprintf('%s, %s, %s: %d connections (%.2f, cluster stat)', strtype{it1}, bandlab{itb}, strmeth{itype}, length(indsel), ThMat(it1,it2)), 'FontUnits', 'normalized', 'FitBoxToText', 'on');
        end
    end
end
% visualize all FC links below the p-value threshold
for it1 = 1:2
    for it2 = 1:min(Nth,2)
        [ii,jj] = ind2sub([Ndip Ndip], indtriu(IndselBands_cell{it1,it2,itb}));
        Conn = sparse(ii,jj,1,Ndip,Ndip);
        figure; plot_bsurfconn_spm('brain', [], Conn, Optconn);
        annotation('textbox', [0.2 0.05 0.6 0.05], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'normalized', ...
                'String', sprintf('%s, %s, %s: %d all thresholded connections (%.2f)', strtype{it1}, bandlab{itb}, strmeth{itype}, length(IndselBands_cell{it1,it2,itb}), ThMat(it1,it2)), 'FontUnits', 'normalized', 'FitBoxToText', 'on');
    end
end
% % plot histograms
figure;
cnt = 0;
for it = 1:Nth
    cnt = cnt + 1;
    subplot(1,Nth,cnt);
    h = histogram(maxClusterExt_stat(2:end,it));
    h.Normalization = 'probability';
    h.BinWidth = 25;
    hold on; plot(ThrClusterSize(it), 0.7, '+'); hold off
%     axis([0 1.1*ThrClusterSize(it) 0 0.9]);
    axis([0 1.2*max(ThrClusterSize) 0 1]);
    set(gca, 'FontSize', 20);
    set(gca, 'YScale', 'log'); set(gca, 'YTick', [0.01 0.05 0.2 0.5 1])
end

% Opt=[]; Opt.nodesize=8; Opt.nodestyle='ok'; Opt.nodefacecolor='k'; Opt.edgesize=2; Opt.edgestyle='-'; Opt.edgecolor='b'; Opt.alpha=0.5;
% anatcol = [238 203 193]/255;
% % [ii,jj] = ind2sub([Ndip Ndip], indtriu(IndselBands_cell{1,2,itb})); Conn = sparse(ii,jj,1,Ndip,Ndip);
% [ii,jj] = ind2sub([Ndip Ndip], indtriu(Indsel_FClinks{1,2})); Conn = sparse(ii,jj,1,Ndip,Ndip);
% figure; plot_bsurfconn(surf, anatcol, Conn, Opt); axis off

% [TableAdj, nROIxROI, labelROIs] = connectivity_ROI(ii, jj, Ndip);
% flag = any(nROIxROI); conn = full(nROIxROI(flag,flag)), sparse(nROIxROI)