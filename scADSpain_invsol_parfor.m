% scADSpain_invsol_parfor
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
twinlength = 2; % sec
%% compute invsol
parfor k = 1:nsubj
    format compact
    close all;
    disp([k nsubj]);
    % --- get data
    tag          = num2str(id_subj(k));
    dirnamesubj  = fullfile(dirname, 'RESULTS', tag);
    tmp          = dir(fullfile(dirnamesubj, 'dobj_*.dat'));
    filename_dat = fullfile(dirnamesubj, tmp.name);
    %     filename_mri = fullfile(dirname, tag, [tag '.nii']);
    D = spm_eeg_load(filename_dat);
    Fs = D.fsample;
    % --- passband filtering
    freqband = [0.5 48];
    S = [];
    S.D = D;
    S.freq = freqband;
    S.type  = 'butterworth';
    S.order = 5;
    S.band  = 'bandpass'; %'low';
    S.dir   = 'twopass';
    D = spm_eeg_filter(S);
    % --- Do SPM inverse solutions
    inv_conditions = {'Undefined'};
    inv_modalities = {{'MEG' 'MEGPLANAR'}};
    inv_type =       {       'COH'       };
    inv_fboi = freqband;
    tmax = max(D.time);
    VxFxT = cell(1,400);
    val = 1;
    tonset = D.timeonset;
    D.inv{val}.inverse = [];   % Clear to be safe!
    D.val = val;
    D.inv{val}.comment = {sprintf('Ind: Val%d: Mod:%s Inv:%s',val,cat(2,inv_modalities{val}{:}),inv_type{val})};
    D.inv{val}.inverse.trials   = inv_conditions;
    D.inv{val}.inverse.Han      = 1;
    %         D.inv{val}.inverse.Han      = 0;
    D.inv{val}.inverse.lpf      = inv_fboi(1);
    D.inv{val}.inverse.hpf      = inv_fboi(2);
    D.inv{val}.inverse.type     = inv_type{val};
    D.inv{val}.inverse.modality = inv_modalities{val};    
    Ntr = 0;
    while(tonset + twinlength < tmax)
        inv_twin = 1000*(tonset + [0 twinlength]);
        D.inv{val}.inverse.woi = inv_twin;
        D = spm_eeg_invert(D);
        
        % tonset = tonset + 0.5; % sliding window of 0.5 sec overlapping
        tonset = tonset + twinlength; % non-overlapped windows
        
        Jsol = D.inv{val}.inverse.J{1}*D.inv{val}.inverse.T';
        Ntr = Ntr + 1;
        if (Ntr == 1)
            NFFT = size(Jsol,2);
            freq = Fs/2*linspace(0,1,NFFT/2+1);
            indsel = find(freq >= freqband(1) & freq <= freqband(2));
            freq = freq(indsel);
        end
        disp(Ntr);
        tmp = fft(Jsol,NFFT,2);
        VxFxT{Ntr} = tmp(:,indsel);
    end
    VxFxT = cell2mat(VxFxT(1:Ntr));
    VxFxT = reshape(VxFxT, size(VxFxT,1), [], Ntr);
    VxFxT = single(VxFxT);
    filename = fullfile(dirnamesubj, ['VxFxT_' inv_type{val} '.mat']);
    save_mat_invsol(filename, VxFxT, freq, NFFT, freqband);
    % save(['VxFxT_' inv_type{val} '.mat'], 'VxFxT', 'freq', 'NFFT', 'freqband');
    D.delete();
end