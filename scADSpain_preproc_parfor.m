% scspm_ADSpain_preproc_parfor
warning off

addpath('C:\WORK\MATLAB\Utiles');
if ~exist('spm.m', 'file')
    addpath('.\spm12')
    spm('defaults','EEG');
end
dirname = '.';
[ndata, fieldnames] = xlsread(fullfile(dirname, 'Resting_eyes_closed_30vs30.xlsx'));
id_subj = ndata(:, strcmp(fieldnames, 'ID_meg'));
timewin = [ndata(:,strncmp(fieldnames,'Start',5)) ndata(:,strncmp(fieldnames,'End',3))];
nsubj = length(id_subj);
for k = 41:nsubj
    close all;
    disp([k nsubj]);
    %% get data
    tag          = num2str(id_subj(k));
    dirnamesubj  = fullfile(dirname, 'RESULTS', tag);
    tmp          = dir(fullfile(dirnamesubj, '*.fif'));
    filename_dat = fullfile(dirnamesubj, tmp.name);

    % filename_out = fullfile(dirnamesubj, ['out_' tag '.mat']);
    filename_out = fullfile(dirnamesubj, ['obj_' tag '.mat']);
    
    filename_mri = fullfile(dirnamesubj, [tag '.nii']);
    filename_fid = fullfile(dirnamesubj, [tag '_fids.txt']);
    filename_mod_fid = fullfile(dirnamesubj, [tag '_mod_fids.xlsx']);
    
    %% Read MEG data and downsample
    S             = [];
    S.timewin     = timewin(k,:);
    S.dataset     = filename_dat;
    S.outfile     = filename_out;
    D             = spm_eeg_convert(S);
    S             = [];
    S.D           = D;
    S.fsample_new = 200;
    D             = spm_eeg_downsample(S);
    delete(S.D);
    %% Read anatomical image and compute meshes
    val                  = 1;
    D.inv{val}           = [];  % If want to be safe!
    D.val                = val;
    D.inv{val}.date      = char(date,datestr(now,15));
    D.inv{val}.comment   = 'Individual MRI';
    D.inv{val}.mesh.sMRI = filename_mri;
    D                    = spm_eeg_inv_mesh_ui(D, val, D.inv{val}.mesh.sMRI, 2);
    %% Data Reg
    fid = fopen(filename_fid);
    C = textscan(fid, '%d%d%d%s');
    fclose(fid);
    MRIfids           = [];
    MRIfids.fid.pnt   = double([C{1} C{2} C{3}]);
    MRIfids.fid.label = strrep(C{4}, 'NAS', 'Nasion');
    MRIfids.pnt       = D.inv{val}.mesh.fid.pnt;        % Scalp mesh points from MRI above
    MRIfids.unit      ='mm';
    MEGfids           = D.fiducials;
    
%     MEGfids           = D.fiducials;  % To remove nose points in next line
%     flag = (MEGfids.pnt(:,2) > 0) & (MEGfids.pnt(:,3) < 0);
%     MEGfids.pnt(flag,:) = [];
    
%     D = spm_eeg_inv_datareg_ui(D, val, MEGfids, MRIfids, false);
    
    % Manual coregistration
    D = spm_coreg(D, val, MEGfids, MRIfids);
    
    % save modified MRI fiducial points
    datareg = D.inv{val}.datareg;
    flag = strcmp({datareg.modality}, 'MEG');
    datareg = datareg(flag);
    mrifids = ft_transform_headshape(D.inv{val}.mesh.Affine\datareg.toMNI, datareg.fid_mri);
    tmp = cell(4,4);
    tmp(1,:) = {'Labels' 'X' 'Y' 'Z'};
    tmp(2:end,2:end) = num2cell(mrifids.fid.pnt);
    xlswrite(filename_mod_fid, tmp);

    %% Forward modelling   
    D.inv{val}.forward = struct([]);
    for ind = 1:length(D.inv{val}.datareg)
        if strcmp(D.inv{val}.datareg(ind).modality,'EEG')
            D.inv{val}.datareg(ind) = [];
            break;
        end
    end
    D.inv{val}.forward(1).voltype = 'Single Shell';
    D = spm_eeg_inv_forward(D);
    D.save
end
rmpath(genpath('.\spm12'));