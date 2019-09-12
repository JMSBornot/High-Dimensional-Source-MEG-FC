warning off
% scspm_ADSpain_display

cols = [...
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

% workdir = 'H:\BACKUP_LAPTOP_DEC182017\E';
% addpath(fullfile(workdir, 'MATLAB', 'Pipelines'));
% addpath(genpath(fullfile(workdir, 'MATLAB', 'spm12')));
% dirname = fullfile(workdir, 'DATOS', 'Spain_MEEG');
% [ndata, fieldnames] = xlsread(fullfile(dirname, 'Resting_eyes_closed.xlsx'));

addpath('C:\WORK\MATLAB\Pipelines')
% dirname = '\\sceis_cl1fs\shared\PhD Students\JoseSanchez\';
dirname = '\\sceis_cl1fs\shared\Functional Brain Mapping\JoseSanchez';
if ~exist('spm.m', 'file')
    addpath(genpath(fullfile(dirname, 'spm12')));
end
[ndata, fieldnames] = xlsread(fullfile(dirname, 'Resting_eyes_closed_30vs30.xlsx'));

id_subj = ndata(:, strcmp(fieldnames, 'ID_meg'));
timewin = [ndata(:,strncmp(fieldnames,'Start',5)) ndata(:,strncmp(fieldnames,'End',3))];
nsubj = length(id_subj);
anatcol = [238 203 193]/255;
twinlength = 60; % sec
Ns = 306; % # mag + plannar channels
Ndip = 1000;
foi = [4 13]; % Hz
ivert = cell(6,1);
for k = 1:nsubj
    close all;
    disp([k nsubj]);
    %% get data
    tag          = num2str(id_subj(k));
    dirnamesubj  = fullfile(dirname, 'RESULTS', tag);
    tmp          = dir(fullfile(dirnamesubj, 'dobj*.dat'));
    filename_dat = fullfile(dirnamesubj, tmp.name);
    %     filename_mri = fullfile(dirname, tag, [tag '.nii']);
    D = spm_eeg_load(filename_dat);
    Fs = D.fsample;
    %% update field names
    listfd = {'sMRI' 'def' 'tess_ctx' 'tess_scalp' 'tess_oskull' 'tess_iskull'};
    for it = 1:length(listfd)
        % D.inv{1}.mesh.(listfd{it}) = strrep(D.inv{1}.mesh.(listfd{it}), 'E:\DATOS', fullfile(workdir, 'DATOS'));
        D.inv{1}.mesh.(listfd{it}) = strrep(D.inv{1}.mesh.(listfd{it}), '.\RESULTS', fullfile(dirname, 'RESULTS'));
    end
    %% check coregistration
    tmp = D.inv{1}.datareg;
    flag = strcmp({tmp.modality}, 'MEG');
    datareg = D.inv{1}.datareg(flag);
    mesh = spm_eeg_inv_transform_mesh(datareg.fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
    grad = D.sensors('MEG');
    figure; ft_plot_sens(grad);
    hold on
    patch('vertices', mesh.tess_ctx.vert, 'faces', mesh.tess_ctx.face, ...
        'EdgeColor', cols(1,:), 'FaceColor', cols(1,:));
    patch('vertices', mesh.tess_iskull.vert, 'faces', mesh.tess_iskull.face, ...
        'EdgeColor', cols(2,:), 'FaceColor', 'none');
    patch('vertices', mesh.tess_scalp.vert, 'faces', mesh.tess_scalp.face, ...
        'EdgeColor', cols(3,:), 'FaceColor', 'none');
    meegfid  = datareg.fid_eeg;
    mrifid = datareg.fid_mri;
    
    flag = (meegfid.pnt(:,2) > 0) & (meegfid.pnt(:,3) < 0);
    
    plot3(meegfid.pnt(~flag,1), meegfid.pnt(~flag,2), meegfid.pnt(~flag,3), 'o', ...
        'MarkerFaceColor', cols(6,:), 'MarkerSize', 4, 'MarkerEdgeColor', 'k');
    plot3(meegfid.pnt(flag,1), meegfid.pnt(flag,2), meegfid.pnt(flag,3), 'o', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 4, 'MarkerEdgeColor', 'k');
    
    
    plot3(meegfid.fid.pnt(:,1), meegfid.fid.pnt(:,2), meegfid.fid.pnt(:,3), 'o', ...
        'MarkerFaceColor', cols(6,:), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
    plot3(mrifid.fid.pnt(:,1), mrifid.fid.pnt(:,2), mrifid.fid.pnt(:,3), 'd', ...
        'MarkerFaceColor', cols(4,:), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
    hold off
    axis image; view([1 0 0]);
    savefig(fullfile(dirnamesubj, [tag '_correg']));
    close all;
end