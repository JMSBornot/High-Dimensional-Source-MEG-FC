% scADSpain_sourceFC_parfor
warning off
% dirname = '\\sceis_cl1fs\shared\PhD Students\JoseSanchez';
dirname = '\\sceis_cl1fs\shared\Functional Brain Mapping\JoseSanchez';
cd(dirname);

addpath('C:\WORK\MATLAB\Utiles');
if ~exist('spm.m', 'file')
    addpath('.\spm12')
    spm('defaults','EEG');
end
[ndata, fieldnames] = xlsread(fullfile(dirname, 'Resting_eyes_closed_30vs30.xlsx'));
id_subj = ndata(:, strcmp(fieldnames, 'ID_meg'));
timewin = [ndata(:,strncmp(fieldnames,'Start',5)) ndata(:,strncmp(fieldnames,'End',3))];
nsubj = length(id_subj);
% parfor k = 2:40 %nsubj
% parfor k = 46:60 %nsubj
parfor k = 1:nsubj
    format compact;
    disp([k nsubj id_subj(k)]); drawnow;
    val = 1;
    inv_type = {'COH'};
    tag = num2str(id_subj(k));
    fprintf('subject: %s\n', tag);
    dirnamesubj = fullfile(dirname, 'RESULTS', tag);
    filename = fullfile(dirnamesubj, ['VxFxT_' inv_type{val} '.mat']);
    [VxFxT, freq, NFFT, freqband] = load_mat_invsol(filename);
%     bandlim = [0.5 4.0 6.0 8.0 10.5 13.0 20.0 30.0 48];
%     bandlab = {'delta' 'thetaI' 'thetaII' 'alphaI' 'alphaII' 'betaI' 'betaII' 'gamma'};
%     nband = length(bandlim)-1;
    isize = [500*ones(1,15) 696];
    nblock = length(isize);
    indprev = find(triu(ones(isize(1)),1));
    indlast = find(triu(ones(isize(end)),1));
    [Nd, Nf, Ntr] = size(VxFxT);
    FxVxT = permute(VxFxT,[2 1 3]);
    VxFxT = [];
    Sp = mean(abs(FxVxT).^2,3);
    FxVxT = mat2cell(FxVxT,Nf,isize,Ntr);
    Sp = mat2cell(Sp,Nf,isize);
    for itr = 1:nblock
        disp([itr nblock]); drawnow;
        for itc = itr:nblock
            fname = fullfile(dirname, 'RESULTS', tag, sprintf('EIC_rn%dcn%d_%04.1fHz.mat',itr,itc,freq(end)));
            if exist(fname, 'file'), continue; end
            % we will add x trial by trial to estimates the mean.
            iCOH = 0;
            den = 0;
            for it = 1:Ntr
                tmp = bsxfun(@times, FxVxT{itr}(:,:,it), reshape(conj(FxVxT{itc}(:,:,it)),[Nf 1 isize(itc)])); % imaginary xspectrum
                tmp = reshape(tmp, [Nf isize(itr)*isize(itc)]);
                % select only block upper-triangular indices
                if (itr == itc)
                    if (itr < nblock)
                        tmp = tmp(:,indprev);
                    else
                        tmp = tmp(:,indlast);
                    end
                end
                iCOH = iCOH + tmp; % imaginary coherence
                % den = den + abs(hilbert(imag(tmp))); % denominator for EIC normalization
                den = den + abs(hilbert(imag(tmp))).^2;
            end
            EIC = abs(hilbert(imag(iCOH/Ntr)./sqrt(den/Ntr))).^2;
            EIC(isnan(EIC)) = 0;
            
%             den = bsxfun(@times, Sp{itr}, reshape(Sp{itc},[Nf 1 isize(itc)])); % denominator for COH normalization
%             den = reshape(den, [Nf isize(itr)*isize(itc)]);
%             if (itr == itc)
%                 if (itr < nblock)
%                     den = den(:,indprev);
%                 else
%                     den = den(:,indlast);
%                 end
%             end
%             iCOH = imag(iCOH/Ntr)./sqrt(den);
            
            for it = 1:Nf
                block = EIC(it,:);
                fname = fullfile(dirname, 'RESULTS', tag, sprintf('EIC_rn%dcn%d_%04.1fHz.mat',itr,itc,freq(it)));
                save_data(fname, block);
            end
            
%             % summarize results by frequency band
%             for it = 1:nband
%                 block = mean(EIC(freq>=bandlim(it) & freq<bandlim(it+1),:));
%                 fname = fullfile(dirname, 'RESULTS', tag, sprintf('EIC_rn%dcn%d_%s.mat',itr,itc,bandlab{it}));
%                 save_data(fname, block);
%                 block = mean(iCOH(freq>=bandlim(it) & freq<bandlim(it+1),:));
%                 fname = fullfile(dirname, 'RESULTS', tag, sprintf('iCOH_rn%dcn%d_%s.mat',itr,itc,bandlab{it}));
%                 save_data(fname, block);
%             end
        end
    end
    disp('done!'); drawnow;
end