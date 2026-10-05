%% *ELL7420: Advanced Digital Signal Processing*
%% Design Assignment 2
%% Pitch Detection - Autocorrelation (time domain) vs Cepstrum (frequency domain)
% Pitch (fundamental frequency, F0) is estimated every 10 ms for four speech
% recordings with two detectors:
%%
% # *Method A - Autocorrelation (time domain).* The frame is low-pass filtered
% and centre clipped, and the lag of the largest autocorrelation peak is the
% pitch period (Rabiner, 1977).
% # *Method B - Cepstrum (frequency domain).* The inverse transform of the log
% magnitude spectrum has a peak at the quefrency equal to the pitch period
% (Noll, 1967).
%%
% s1.wav and s2.wav come with reference pitch tracks (p1.mat, p2.mat), so both
% detectors are scored against them. s3.wav and s4.wav have no reference; their
% pitch tracks are the output of this assignment.
%
% Needs the Signal Processing Toolbox (firpm, hamming, medfilt1, spectrogram).

clear; clc; close all;
rng(0);   % fixes the random numbers so the noise test gives the same answer every run
%% 1. Data Preparation
% speech_samples holds four recordings (8 kHz, 3 s, 16 bit) and two reference
% tracks:
%%
% * s1..s4 - 24000 x 1 speech signals.
% * p1, p2 - 300 x 1 reference pitch in Hz, one value per 10 ms. Zero means
% unvoiced or silence.

if isfolder('speech_samples')
    data_dir = 'speech_samples';
elseif isfolder(fullfile('..','speech_samples'))
    data_dir = fullfile('..','speech_samples');
else
    error('speech_samples folder not found next to, or one level above, this script.');
end
out_dir = 'assets';
if ~isfolder(out_dir), mkdir(out_dir); end
log_file = fullfile(out_dir,'results_log.txt');
if isfile(log_file), delete(log_file); end
diary(log_file);                                      % every number in the report comes from this log

% One colour per entity, used in every figure
C.ref   = [0.13 0.13 0.13];                           % reference track
C.acf   = [31 95 168]/255;                            % autocorrelation method
C.cep   = [176 31 82]/255;                            % cepstrum method
C.muted = [0.60 0.60 0.60];
set(groot,'defaultAxesFontSize',8,'defaultAxesBox','off','defaultAxesTickDir','out', ...
    'defaultAxesXGrid','on','defaultAxesYGrid','on','defaultAxesGridAlpha',0.12, ...
    'defaultLineLineWidth',1.1,'defaultTextFontSize',8,'defaultAxesTitleFontSizeMultiplier',1);

names      = {'s1','s2','s3','s4'};
num_files  = numel(names);
speech     = cell(1,num_files);
for i = 1:num_files
    [speech{i}, fs] = audioread(fullfile(data_dir,[names{i} '.wav']));
end
tmp = load(fullfile(data_dir,'p1.mat'));  reference{1} = tmp.p1(:);
tmp = load(fullfile(data_dir,'p2.mat'));  reference{2} = tmp.p2(:);

fprintf('Loaded %d recordings, %d samples each, fs = %d Hz.\n\n', ...
    num_files, numel(speech{1}), fs);
%% 2. Analysis parameters
% *2.A Framing (common to both methods)*
%
% Frame step = 10 ms (matches the reference tracks)
%
% Pitch search range = 70 Hz to 400 Hz (lags 20 to 115 samples)
%
% *2.B Autocorrelation method*
%
% Window = 30 ms rectangular, low-pass 0-900 Hz, centre clipping at 68 %
%
% *2.C Cepstrum method*
%
% Window = 40 ms Hamming, 1024-point FFT
%
% *2.D Voicing decision and smoothing (common)*
%
% Voiced when the peak clears the method threshold and the frame is not silent.
% A 5-point median smoother removes isolated errors.

P.fs              = fs;
P.hop             = round(0.010*fs);     % 10 ms frame step -> 80 samples
P.num_frames      = floor(numel(speech{1})/P.hop);   % 300 frames
P.offset          = round(P.hop/2);      % window centre sits mid-way through each 10 ms segment
P.f0_min          = 70;                  % Hz - lowest pitch searched
P.f0_max          = 400;                 % Hz - highest pitch searched
P.lag_min         = floor(fs/P.f0_max);  % 20 samples
P.lag_max         = ceil(fs/P.f0_min);   % 115 samples

P.acf_win         = round(0.030*fs);     % 240 samples
P.lpf             = firpm(40,[0 900 1700 fs/2]/(fs/2),[1 1 0 0]);  % Rabiner's 0-900 Hz low-pass
P.clip_frac       = 0.68;                % centre clipping level, fraction of the frame peak
P.acf_thr         = 0.50;                % taper-corrected autocorrelation peak needed for "voiced"

P.cep_win         = round(0.040*fs);     % 320 samples
P.nfft            = 1024;                % FFT length used for the cepstrum
P.cep_thr         = 0.12;                % cepstral peak needed for "voiced"

P.energy_win      = round(0.030*fs);     % 240 samples - window of the silence gate
P.energy_floor_dB = -30;                 % frames this far below the loudest frame are silence
P.med_len         = 5;                   % median smoother length (frames)
P.gross_tol       = 0.20;                % error above 20 % of the reference is a gross error

frame_time = ((0:P.num_frames-1)' + 0.5) * P.hop / fs;   % centre of each 10 ms segment, s
%% 3. Sanity check on a synthetic signal
% A harmonic signal of known pitch must come back within 0.6 % from both
% detectors. This fails loudly if a lag or quefrency index is off by one.

for f0_true = [90 150 220 330]
    t_syn = (0:numel(speech{1})-1)'/fs;
    x_syn = zeros(size(t_syn));
    for h = 1:floor(0.45*fs/f0_true)
        x_syn = x_syn + cos(2*pi*h*f0_true*t_syn)/h;
    end
    f_a = pitch_acf(x_syn,P);  f_c = pitch_cep(x_syn,P);
    mid = 50:250;
    assert(all(abs(f_a(mid)-f0_true) < 0.006*f0_true), 'ACF sanity check failed at %d Hz', f0_true);
    assert(all(abs(f_c(mid)-f0_true) < 0.006*f0_true), 'Cepstrum sanity check failed at %d Hz', f0_true);
end
fprintf('Sanity check passed: both detectors recover 90, 150, 220 and 330 Hz within 0.6 %%.\n\n');
%% 4. Run both detectors on all four recordings
% Each detector returns, for every frame, a raw pitch candidate and the height
% of the peak it came from. The voicing decision and the median smoother then
% turn the candidates into the final track.

for i = 1:num_files
    [raw_acf{i}, peak_acf{i}, map_acf{i}, detail_acf{i}] = pitch_acf(speech{i},P);
    [raw_cep{i}, peak_cep{i}, map_cep{i}, detail_cep{i}] = pitch_cep(speech{i},P);
    energy{i}    = frame_energy_dB(speech{i},P);
    track_acf{i} = make_track(raw_acf{i}, peak_acf{i}, energy{i}, P.acf_thr, P);
    track_cep{i} = make_track(raw_cep{i}, peak_cep{i}, energy{i}, P.cep_thr, P);
end
%% 5. Score against the reference (s1, s2)
% Error measures follow Rabiner et al. (1976):
%%
% * V->UV, UV->V - voicing errors, % of reference voiced / unvoiced frames.
% * VDE - voicing decision error, % of all frames with the wrong voicing.
% * GPE - gross pitch error, % of frames voiced in both tracks whose pitch is
% off by more than 20 %.
% * Fine error - mean, standard deviation and mean absolute error over the
% remaining frames.
% * FFE - F0 frame error, % of all frames with either a voicing or a gross error.

P_raw = P;  P_raw.med_len = 1;                         % same pipeline without the smoother
ref_all = [reference{1}; reference{2}];
fprintf('%-6s %-9s %-9s %7s %7s %7s %7s %9s %8s %8s %7s\n', ...
    'file','method','smoothing','V->UV%','UV->V%','VDE%','GPE%','mean Hz','std Hz','MAE Hz','FFE%');
score_rows = {};
for i = 1:3                                            % 3 = s1 and s2 pooled
    if i < 3, sel = i; label = names{i}; else, sel = 1:2; label = 's1+s2'; end
    ref_i = vertcat(reference{sel});
    for smooth = [false true]
        if smooth, Q = P; tag = 'median-5'; else, Q = P_raw; tag = 'none'; end
        est_a = [];  est_c = [];
        for j = sel
            est_a = [est_a; make_track(raw_acf{j},peak_acf{j},energy{j},P.acf_thr,Q)]; %#ok<AGROW>
            est_c = [est_c; make_track(raw_cep{j},peak_cep{j},energy{j},P.cep_thr,Q)]; %#ok<AGROW>
        end
        Sa = score_track(est_a, ref_i, P.gross_tol);  Sc = score_track(est_c, ref_i, P.gross_tol);
        print_score(label,'ACF',tag,Sa);  print_score(label,'Cepstrum',tag,Sc);
        score_rows(end+1,:) = [{label,'ACF',tag}, struct2cell(Sa)']; %#ok<SAGROW>
        score_rows(end+1,:) = [{label,'Cepstrum',tag}, struct2cell(Sc)']; %#ok<SAGROW>
        if smooth, score_acf(i) = Sa; score_cep(i) = Sc; end
    end
end
writetable(cell2table(score_rows,'VariableNames', ...
    [{'file','method','smoothing'}, fieldnames(score_acf(1))']), fullfile(out_dir,'metrics_s1_s2.csv'));

% Raw candidates on reference-voiced frames: no voicing decision, no smoothing.
% This isolates the pitch estimator from the voicing detector.
fprintf('\nRaw pitch candidates on reference-voiced frames (oracle voicing):\n');
fprintf('%-6s %-9s %7s %8s %9s %9s %7s\n','file','method','frames','gross %','halving','doubling','other');
for i = 1:2
    v = reference{i} > 0;
    for m = 1:2
        if m == 1, est = raw_acf{i}(v); mname = 'ACF'; else, est = raw_cep{i}(v); mname = 'Cepstrum'; end
        ratio   = est ./ reference{i}(v);
        gross   = abs(ratio-1) > P.gross_tol;
        halved  = gross & ratio > 0.4 & ratio < 0.6;   % detector picked twice the true period
        doubled = gross & ratio > 1.8 & ratio < 2.2;   % detector picked half the true period
        oracle_gross(i,m) = 100*mean(gross);
        fprintf('%-6s %-9s %7d %8.1f %9d %9d %7d\n', names{i}, mname, sum(v), oracle_gross(i,m), ...
            sum(halved), sum(doubled), sum(gross & ~halved & ~doubled));
        voiced_frames = find(v);
        gross_frames{i,m} = voiced_frames(gross)'; %#ok<SAGROW>
    end
end
for i = 1:2
    fprintf('  %s frames with a raw gross error: ACF [%s], cepstrum [%s]\n', names{i}, ...
        num2str(gross_frames{i,1}), num2str(gross_frames{i,2}));
end
%% 6. Voicing threshold sweep
% The threshold on the peak height trades missed voiced frames (V->UV) against
% false voiced frames (UV->V). The sweep shows where the chosen thresholds sit
% and how well each peak height separates voiced from unvoiced speech.

thr_acf_grid = 0.05:0.01:1.00;
thr_cep_grid = 0.02:0.005:0.30;
sweep_acf = sweep_threshold(raw_acf, peak_acf, energy, reference, thr_acf_grid, P);
sweep_cep = sweep_threshold(raw_cep, peak_cep, energy, reference, thr_cep_grid, P);
[~,ia] = min(sweep_acf.ffe(:,3));  [~,ic] = min(sweep_cep.ffe(:,3));
fprintf('\nThreshold sweep, s1+s2 pooled: ACF lowest FFE %.2f %% at %.2f, cepstrum lowest FFE %.2f %% at %.3f\n', ...
    sweep_acf.ffe(ia,3), thr_acf_grid(ia), sweep_cep.ffe(ic,3), thr_cep_grid(ic));
fprintf('Chosen thresholds: ACF %.2f (FFE %.2f %%), cepstrum %.2f (FFE %.2f %%)\n', ...
    P.acf_thr, score_acf(3).ffe, P.cep_thr, score_cep(3).ffe);
for i = 1:2                                            % tune on one file, test on the other
    [~,ja] = min(sweep_acf.ffe(:,i));  [~,jc] = min(sweep_cep.ffe(:,i));
    fprintf('  tuned on %s only: ACF thr %.2f -> FFE on %s %.2f %% | cepstrum thr %.3f -> FFE on %s %.2f %%\n', ...
        names{i}, thr_acf_grid(ja), names{3-i}, sweep_acf.ffe(ja,3-i), ...
        thr_cep_grid(jc), names{3-i}, sweep_cep.ffe(jc,3-i));
end
%% 7. Frame alignment check
% The reference gives one value per 10 ms but not where its analysis window
% sits. The window centre is moved either side of the mid-segment position
% used here to see how much the scores depend on that choice.

fprintf('\nFrame alignment (window centre relative to the middle of the segment), pooled FFE %%:\n');
for shift = [-80 -60 -40 -20 0 20 40 80 120]
    Q = P;  Q.offset = P.offset + shift;  out = run_pair(speech, reference, Q);
    fprintf('  %+5.1f ms: ACF %.2f, cepstrum %.2f\n', 1000*shift/fs, out.ffe(3,:));
end
%% 8. Pitch estimates for s3 and s4
% No reference exists for these two recordings. The two tracks are saved, and
% their agreement with each other is the only internal check available.

p3_acf = track_acf{3};  p3_cep = track_cep{3};
p4_acf = track_acf{4};  p4_cep = track_cep{4};
save('p3.mat','p3_acf','p3_cep');
save('p4.mat','p4_acf','p4_cep');

track_table = table((1:P.num_frames)', frame_time, reference{1}, track_acf{1}, track_cep{1}, ...
    reference{2}, track_acf{2}, track_cep{2}, p3_acf, p3_cep, p4_acf, p4_cep, 'VariableNames', ...
    {'frame','time_s','s1_ref','s1_acf','s1_cep','s2_ref','s2_acf','s2_cep','s3_acf','s3_cep','s4_acf','s4_cep'});
track_table{:,3:end} = round(track_table{:,3:end},1);
writetable(track_table, fullfile(out_dir,'pitch_tracks.csv'));

fprintf('\nPitch tracks without reference:\n');
fprintf('%-4s %-9s %7s %8s %8s %8s %8s\n','file','method','voiced','mean Hz','median','min','max');
for i = 3:4
    for m = 1:2
        if m == 1, tr = track_acf{i}; mname = 'ACF'; else, tr = track_cep{i}; mname = 'Cepstrum'; end
        fprintf('%-4s %-9s %7d %8.1f %8.1f %8.1f %8.1f\n', names{i}, mname, sum(tr>0), ...
            mean(tr(tr>0)), median(tr(tr>0)), min(tr(tr>0)), max(tr(tr>0)));
    end
end
fprintf('\nAgreement between the two methods (all four recordings):\n');
fprintf('%-4s %12s %12s %10s %10s %12s\n','file','same voicing','both voiced','within 5%','within 20%','mean |diff|');
for i = 1:num_files
    a = track_acf{i};  c = track_cep{i};  both = a > 0 & c > 0;
    rel = abs(a(both)-c(both)) ./ (0.5*(a(both)+c(both)));
    agree(i,:) = [100*mean((a>0) == (c>0)), sum(both), 100*mean(rel<0.05), 100*mean(rel<0.20), ...
        mean(abs(a(both)-c(both)))]; %#ok<SAGROW>
    fprintf('%-4s %11.1f%% %12d %9.1f%% %9.1f%% %9.2f Hz\n', names{i}, agree(i,:));
end
fprintf('\nFrames of s3 and s4 where the methods differ in voicing or by more than 20 %% in pitch:\n');
fprintf('%-4s %5s %7s %8s %8s %9s %9s %9s\n','file','frame','time s','ACF Hz','cep Hz','energy dB','ACF peak','cep peak');
for i = 3:4
    a = track_acf{i};  c = track_cep{i};
    differ = find((a>0) ~= (c>0) | (a>0 & c>0 & abs(a-c) > P.gross_tol*0.5*(a+c)));
    for k = differ'
        fprintf('%-4s %5d %7.3f %8.1f %8.1f %9.1f %9.2f %9.3f\n', names{i}, k, frame_time(k), a(k), c(k), ...
            energy{i}(k), peak_acf{i}(k), peak_cep{i}(k));
    end
end
%% 9. Noise test
% White Gaussian noise is added to s1 and s2 at a set signal-to-noise ratio
% (measured over the whole recording) and both detectors are scored again with
% the thresholds left untouched.

snr_levels_dB = [0 5 10 15 20 30];
num_trials    = 20;
noise_raw = zeros(numel(snr_levels_dB),2,num_trials);  noise_vde = noise_raw;  noise_ffe = noise_raw;
noise_raw_nolpf = zeros(numel(snr_levels_dB),num_trials);
P_nolpf = P;  P_nolpf.lpf = 1;                         % autocorrelation detector without its low-pass filter
for s = 1:numel(snr_levels_dB)
    for trial = 1:num_trials
        for i = 1:2
            noise_std = sqrt(mean(speech{i}.^2) / 10^(snr_levels_dB(s)/10));
            noisy{i}  = speech{i} + noise_std*randn(size(speech{i}));
        end
        out = run_pair(noisy, reference, P);
        noise_raw(s,:,trial) = out.raw_gross(3,:);
        noise_vde(s,:,trial) = out.vde(3,:);
        noise_ffe(s,:,trial) = out.ffe(3,:);
        out = run_pair(noisy, reference, P_nolpf);
        noise_raw_nolpf(s,trial) = out.raw_gross(3,1);
    end
end
clean = run_pair(speech, reference, P);
clean_nolpf = run_pair(speech, reference, P_nolpf);
fprintf('\nNoise test, s1+s2 pooled, mean of %d trials (ACF | cepstrum):\n', num_trials);
fprintf('%8s %19s %19s %19s %14s\n','SNR dB','raw gross error %','VDE %','FFE %','ACF no LPF raw');
for s = 1:numel(snr_levels_dB)
    fprintf('%8d %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %14.1f\n', snr_levels_dB(s), ...
        mean(noise_raw(s,:,:),3), mean(noise_vde(s,:,:),3), mean(noise_ffe(s,:,:),3), mean(noise_raw_nolpf(s,:)));
end
fprintf('%8s %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %14.1f\n','clean', clean.raw_gross(3,:), clean.vde(3,:), ...
    clean.ffe(3,:), clean_nolpf.raw_gross(3,1));
for i = 1:2
    active = energy{i} >= P.energy_floor_dB;  F = make_frames(speech{i}, P.energy_win, P);
    fprintf('  %s: power of frames above the silence gate is %.1f dB above the whole-file power\n', names{i}, ...
        10*log10(mean(mean(F(active,:).^2,2)) / mean(speech{i}.^2)));
end
fprintf('  the low-pass filter passes %.1f %% of white-noise power (%.1f dB)\n', 100*sum(P.lpf.^2), 10*log10(sum(P.lpf.^2)));
%% 10. Sensitivity to window length and clipping level
% Both detectors are re-run with one design parameter changed at a time (the
% window sweep gives both the same length). The raw gross error needs no
% threshold, so it shows the effect on the pitch estimator alone.

win_ms_grid = [20 25 30 35 40 50 60];
fprintf('\nWindow length sweep, raw gross error %% (s1 ACF, s1 cep, s2 ACF, s2 cep) and pooled FFE %% (ACF, cep):\n');
for j = 1:numel(win_ms_grid)
    Q = P;  Q.acf_win = round(win_ms_grid(j)*fs/1000);  Q.cep_win = Q.acf_win;
    out = run_pair(speech, reference, Q);
    win_raw(j,:) = [out.raw_gross(1,:) out.raw_gross(2,:)]; %#ok<SAGROW>
    win_ffe(j,:) = out.ffe(3,:); %#ok<SAGROW>
    fprintf('  %2d ms: %5.1f %5.1f %5.1f %5.1f | %5.2f %5.2f\n', win_ms_grid(j), win_raw(j,:), win_ffe(j,:));
end
clip_grid = [0 0.2 0.3 0.4 0.5 0.6 0.68 0.8];
fprintf('\nCentre clipping level sweep, ACF raw gross error %% (s1, s2):\n');
for j = 1:numel(clip_grid)
    Q = P;  Q.clip_frac = clip_grid(j);  out = run_pair(speech, reference, Q);
    clip_raw(j,:) = out.raw_gross(1:2,1)'; %#ok<SAGROW>
    fprintf('  clip %.2f: %5.1f %5.1f\n', clip_grid(j), clip_raw(j,:));
end
Q = P;  Q.lpf = 1;  Q.clip_frac = 0;  out = run_pair(speech, reference, Q);
plain_raw = out.raw_gross(1:2,1)';
fprintf('  plain autocorrelation (no low-pass, no clipping): %5.1f %5.1f\n', plain_raw);
%% 11. Run time
% Time to process one 3 s recording (300 frames), measured with timeit.

time_acf = timeit(@() pitch_acf(speech{1},P));
time_cep = timeit(@() pitch_cep(speech{1},P));
fprintf('\nRun time for one 3 s recording: ACF %.2f ms (%.1f us/frame), cepstrum %.2f ms (%.1f us/frame)\n', ...
    1e3*time_acf, 1e6*time_acf/P.num_frames, 1e3*time_cep, 1e6*time_cep/P.num_frames);
%% 12. Figures
% All figures are written to assets/ as PNG files.

lag_ms    = (0:P.lag_max+1)/fs*1000;                   % lag / quefrency axis in ms
search_ms = [P.lag_min P.lag_max]/fs*1000;
method_name = {'autocorrelation','cepstrum'};  method_colour = {C.acf, C.cep};

% ---- Fig 2: block diagrams of the two detectors ----
fig = new_figure(17,5.2);  ax = axes(fig,'Position',[0.01 0.02 0.98 0.96]);  hold(ax,'on');  axis(ax,'off');
draw_chain(ax, 2.1, {sprintf('Low-pass\n0-900 Hz'), sprintf('30 ms frame\nrectangular'), ...
    sprintf('Centre clip\n68 %% of peak'), sprintf('Autocorrelation\nr(k)'), sprintf('Peak search\n2.5-14.4 ms'), ...
    sprintf('Voiced if peak\n\\geq %.2f and\nnot silent',P.acf_thr), sprintf('Median\nsmoothing')}, C.acf, ...
    'Method A - autocorrelation (time domain)', 'F_0 = f_s / k*');
draw_chain(ax, 0.6, {sprintf('40 ms frame\nHamming'), sprintf('FFT\n1024 points'), sprintf('Logarithm\nlog |X(k)|'), ...
    sprintf('Inverse FFT\nc(q)'), sprintf('Peak search\n2.5-14.4 ms'), ...
    sprintf('Voiced if peak\n\\geq %.2f and\nnot silent',P.cep_thr), sprintf('Median\nsmoothing')}, C.cep, ...
    'Method B - cepstrum (frequency domain)', 'F_0 = f_s / q*');
xlim(ax,[-0.8 12.2]);  ylim(ax,[0 3]);
save_figure(fig,out_dir,'fig02_block_diagrams');

% ---- Fig 1: the four recordings ----
fig = new_figure(17,10);  tl = tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
t = (0:numel(speech{1})-1)'/fs;
for i = 1:num_files
    ax = nexttile(tl);  hold(ax,'on');
    if i <= 2                                           % shade the reference voiced segments
        edges = diff([0; reference{i} > 0; 0]);
        on = find(edges == 1);  off = find(edges == -1) - 1;
        for k = 1:numel(on)
            xregion(ax,(on(k)-1)*P.hop/fs, off(k)*P.hop/fs,'FaceColor',C.muted,'FaceAlpha',0.25,'EdgeColor','none');
        end
        label = sprintf('%s.wav - shaded: voiced in the reference (%d of %d frames)', names{i}, ...
            sum(reference{i}>0), P.num_frames);
    else
        label = sprintf('%s.wav - no reference', names{i});
    end
    plot(ax,t,speech{i},'Color',C.ref,'LineWidth',0.4);
    ylim(ax,[-1 1]);  xlim(ax,[0 3]);  ylabel(ax,'amplitude');
    title(ax,label,'FontWeight','normal');
    if i < num_files, ax.XTickLabel = []; else, xlabel(ax,'time (s)'); end
end
save_figure(fig,out_dir,'fig01_recordings');

% ---- Fig 3 and 4: one voiced and one unvoiced frame of s2 through each detector ----
k_voiced   = 70;                                        % steady vowel, reference 133.3 Hz
isolated   = conv(double(reference{2} > 0), ones(7,1), 'same') == 0;   % no voiced frame within 30 ms
candidates = find(isolated & track_acf{2}==0 & track_cep{2}==0);
[~,j]      = max(energy{2}(candidates));
k_unvoiced = candidates(j);                             % loudest frame well inside an unvoiced stretch
frame_ids  = [k_voiced k_unvoiced];  frame_kind = {'voiced','unvoiced'};
fprintf('\nExample frames of s2: voiced %d (%.3f s, reference %.1f Hz), unvoiced %d (%.3f s, %.1f dB)\n', ...
    k_voiced, frame_time(k_voiced), reference{2}(k_voiced), k_unvoiced, frame_time(k_unvoiced), energy{2}(k_unvoiced));
fprintf('  voiced frame: ACF %.1f Hz (r = %.2f, taper-corrected %.2f), cepstrum %.1f Hz (peak %.3f)\n', ...
    raw_acf{2}(k_voiced), detail_acf{2}.r_peak(k_voiced), peak_acf{2}(k_voiced), raw_cep{2}(k_voiced), peak_cep{2}(k_voiced));
fprintf('  unvoiced frame: ACF r = %.2f, taper-corrected %.2f, cepstrum peak %.3f\n', ...
    detail_acf{2}.r_peak(k_unvoiced), peak_acf{2}(k_unvoiced), peak_cep{2}(k_unvoiced));

fig = new_figure(17,9.5);  tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
t_ms = (0:P.acf_win-1)/fs*1000;
for r = 1:2
    k = frame_ids(r);  D = detail_acf{2};
    ax = nexttile(tl);  hold(ax,'on');
    plot(ax,t_ms,D.F(k,:),'Color',C.ref);
    yline(ax,D.CL(k)*[1 -1],'--','Color',C.muted,'LineWidth',0.9);
    title(ax,sprintf('%s frame, low-passed',frame_kind{r}),'FontWeight','normal');
    xlabel(ax,'time in frame (ms)');  ylabel(ax,'amplitude');  xlim(ax,[0 30]);
    ax = nexttile(tl);
    plot(ax,t_ms,D.Fc(k,:),'Color',C.acf);
    title(ax,'after centre clipping','FontWeight','normal');
    xlabel(ax,'time in frame (ms)');  xlim(ax,[0 30]);
    ax = nexttile(tl);  hold(ax,'on');
    r_plain = real(ifft(abs(fft(D.F(k,:),512)).^2));  r_plain = r_plain(1:P.lag_max+2)/r_plain(1);
    xregion(ax,search_ms(1),search_ms(2),'FaceColor',C.muted,'FaceAlpha',0.15,'EdgeColor','none');
    h1 = plot(ax,lag_ms,r_plain,'Color',C.muted);
    h2 = plot(ax,lag_ms,map_acf{2}(k,:),'Color',C.acf);
    k_search = P.lag_min:P.lag_max;                     % r(k*) N/(N-k*) >= threshold, drawn on r(k)
    plot(ax,k_search/fs*1000,P.acf_thr*(P.acf_win-k_search)/P.acf_win,':','Color',C.ref,'LineWidth',0.9);
    plot(ax,1000/raw_acf{2}(k),D.r_peak(k),'o','MarkerSize',6,'MarkerFaceColor',C.acf,'MarkerEdgeColor','w');
    text(ax,1000/raw_acf{2}(k)+0.4,min(D.r_peak(k)+0.14,0.95),sprintf('%.2f ms = %.1f Hz', ...
        1000/raw_acf{2}(k),raw_acf{2}(k)),'FontSize',7);
    title(ax,'normalised autocorrelation','FontWeight','normal');
    xlabel(ax,'lag (ms)');  ylabel(ax,'r(k)');  xlim(ax,[0 lag_ms(end)]);  ylim(ax,[-1 1.05]);
    if r == 1, legend(ax,[h1 h2],{'without clipping','with clipping'},'Location','southeast','Box','off'); end
end
save_figure(fig,out_dir,'fig03_acf_frames');

fig = new_figure(17,9.5);  tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
t_ms = (0:P.cep_win-1)/fs*1000;  f_axis = (0:P.nfft/2)*fs/P.nfft;
for r = 1:2
    k = frame_ids(r);  D = detail_cep{2};
    ax = nexttile(tl);
    plot(ax,t_ms,D.F(k,:),'Color',C.ref);
    title(ax,sprintf('%s frame, Hamming window',frame_kind{r}),'FontWeight','normal');
    xlabel(ax,'time in frame (ms)');  ylabel(ax,'amplitude');  xlim(ax,[0 40]);
    ax = nexttile(tl);  hold(ax,'on');
    log_spec = log(abs(fft(D.F(k,:),P.nfft)) + eps);
    cep_full = real(ifft(log_spec));
    lifter   = zeros(1,P.nfft);  lifter([1:P.lag_min, P.nfft-P.lag_min+2:P.nfft]) = 1;
    envelope = real(fft(cep_full.*lifter));             % low-quefrency part = vocal tract envelope
    h1 = plot(ax,f_axis,20/log(10)*log_spec(1:P.nfft/2+1),'Color',C.muted,'LineWidth',0.7);
    h2 = plot(ax,f_axis,20/log(10)*envelope(1:P.nfft/2+1),'Color',C.cep,'LineWidth',1.3);
    title(ax,'log magnitude spectrum','FontWeight','normal');
    xlabel(ax,'frequency (Hz)');  ylabel(ax,'dB');  xlim(ax,[0 fs/2]);
    yl = ylim(ax);  ylim(ax,[yl(1) yl(2)+22]);
    if r == 1, legend(ax,[h1 h2],{'log |X(k)|','envelope'},'Location','northeast','Box','off'); end
    ax = nexttile(tl);  hold(ax,'on');
    xregion(ax,search_ms(1),search_ms(2),'FaceColor',C.muted,'FaceAlpha',0.15,'EdgeColor','none');
    plot(ax,lag_ms,map_cep{2}(k,:),'Color',C.cep);
    yline(ax,P.cep_thr,':','Color',C.ref,'LineWidth',0.9);
    plot(ax,1000/raw_cep{2}(k),peak_cep{2}(k),'o','MarkerSize',6,'MarkerFaceColor',C.cep,'MarkerEdgeColor','w');
    text(ax,1000/raw_cep{2}(k)+0.4,peak_cep{2}(k)+0.06,sprintf('%.2f ms = %.1f Hz', ...
        1000/raw_cep{2}(k),raw_cep{2}(k)),'FontSize',7,'BackgroundColor','w','Margin',0.5);
    title(ax,'real cepstrum','FontWeight','normal');
    xlabel(ax,'quefrency (ms)');  ylabel(ax,'c(q)');  xlim(ax,[0 lag_ms(end)]);  ylim(ax,[-0.15 0.55]);
end
save_figure(fig,out_dir,'fig04_cepstrum_frames');

% ---- Fig 5: what each detector searches - autocorrelation and cepstrum over time ----
fig = new_figure(17,11);  tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ramp = linspace(0,1,128)';
for i = 1:2
    for m = 1:2
        ax = nexttile(tl);  hold(ax,'on');
        if m == 1
            imagesc(ax,frame_time,lag_ms,map_acf{i}');  clim(ax,[0 1]);
            colormap(ax,1 - ramp*(1-C.acf));  what = 'autocorrelation r(k)';
        else
            imagesc(ax,frame_time,lag_ms,map_cep{i}');  clim(ax,[0 0.5]);
            colormap(ax,1 - ramp*(1-C.cep));  what = 'cepstrum c(q)';
        end
        period_ms = 1000 ./ reference{i};  period_ms(reference{i} == 0) = NaN;
        plot(ax,frame_time,period_ms,'.','Color',C.ref,'MarkerSize',5);
        axis(ax,'xy');  xlim(ax,[0 3]);  ylim(ax,search_ms);  grid(ax,'off');
        title(ax,sprintf('%s: %s, dots = reference period',names{i},what),'FontWeight','normal');
        if m == 1, ylabel(ax,'lag / quefrency (ms)'); end
        if i == 2, xlabel(ax,'time (s)'); end
        colorbar(ax);
    end
end
save_figure(fig,out_dir,'fig05_pitch_maps');

% ---- Fig 6, 7, 11, 12: pitch tracking on each recording ----
fig_no = [6 7 11 12];
for i = 1:num_files
    if i <= 2, ref_i = reference{i}; else, ref_i = []; end
    fig = plot_tracking(speech{i}, names{i}, ref_i, track_acf{i}, track_cep{i}, frame_time, P, C);
    save_figure(fig,out_dir,sprintf('fig%02d_tracking_%s',fig_no(i),names{i}));
end

% ---- Fig 8: raw candidates against the reference ----
fig = new_figure(17,8.2);  tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
marker = {'o','^'};
for m = 1:2
    ax = nexttile(tl);  hold(ax,'on');
    f_line = [60 420];
    patch(ax,[f_line fliplr(f_line)],[f_line*(1-P.gross_tol) fliplr(f_line)*(1+P.gross_tol)], ...
        C.muted,'FaceAlpha',0.18,'EdgeColor','none');
    plot(ax,f_line,f_line,'-','Color',C.muted,'LineWidth',0.8);
    plot(ax,f_line,2*f_line,':','Color',C.muted,'LineWidth',0.8);
    plot(ax,f_line,0.5*f_line,':','Color',C.muted,'LineWidth',0.8);
    for i = 1:2
        v = reference{i} > 0;
        if m == 1, est = raw_acf{i}(v); else, est = raw_cep{i}(v); end
        hs(i) = scatter(ax,reference{i}(v),est,14,marker{i},'MarkerFaceColor',method_colour{m}, ...
            'MarkerEdgeColor','w','MarkerFaceAlpha',0.75,'LineWidth',0.4); %#ok<SAGROW>
    end
    text(ax,150,330,'2 x reference','FontSize',7,'Color',C.ref,'Rotation',38);
    text(ax,250,112,'reference / 2','FontSize',7,'Color',C.ref,'Rotation',38);
    set(ax,'XScale','log','YScale','log');  xlim(ax,f_line);  ylim(ax,f_line);
    xticks(ax,[70 100 150 200 300 400]);  yticks(ax,[70 100 150 200 300 400]);
    xlabel(ax,'reference pitch (Hz)');  ylabel(ax,'raw candidate (Hz)');
    title(ax,sprintf('%s: gross errors s1 %.1f %%, s2 %.1f %%',method_name{m},oracle_gross(1,m),oracle_gross(2,m)), ...
        'FontWeight','normal');
    legend(ax,hs,{'s1','s2'},'Location','northwest','Box','off');
end
save_figure(fig,out_dir,'fig08_scatter_reference');

% ---- Fig 9: fine error distribution ----
fig = new_figure(17,6.5);  tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
est_all = {[track_acf{1}; track_acf{2}], [track_cep{1}; track_cep{2}]};
fprintf('\nFine pitch error, s1+s2 pooled, %% of reference:\n');
for m = 1:2
    ax = nexttile(tl);  hold(ax,'on');
    both = ref_all > 0 & est_all{m} > 0;
    rel  = 100*(est_all{m}(both) - ref_all(both)) ./ ref_all(both);
    rel  = rel(abs(rel) <= 100*P.gross_tol);
    histogram(ax,rel,-20.5:1:20.5,'FaceColor',method_colour{m},'EdgeColor','w','FaceAlpha',0.9);
    xline(ax,0,'-','Color',C.ref,'LineWidth',0.8);
    xlabel(ax,'pitch error (% of reference)');  ylabel(ax,'frames');  xlim(ax,[-20 20]);
    title(ax,sprintf('%s: mean %+.2f %%, std %.2f %%, %d frames',method_name{m},mean(rel),std(rel),numel(rel)), ...
        'FontWeight','normal');
    fprintf('  %-16s mean %+.2f, std %.2f, mean |error| %.2f, within 5 %%: %.1f %% of %d frames\n', ...
        method_name{m}, mean(rel), std(rel), mean(abs(rel)), 100*mean(abs(rel)<5), numel(rel));
end
save_figure(fig,out_dir,'fig09_fine_error');

% ---- Fig 10: voicing decision ----
fig = new_figure(17,6.8);  tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
peak_all = {[peak_acf{1}; peak_acf{2}], [peak_cep{1}; peak_cep{2}]};
loud     = [energy{1}; energy{2}] >= P.energy_floor_dB;
fprintf('\nVoicing histograms: %d of %d frames pass the silence gate (%d reference voiced, %d unvoiced); largest peaks %.2f and %.2f\n', ...
    sum(loud), numel(loud), sum(loud & ref_all > 0), sum(loud & ref_all == 0), max(peak_all{1}(loud)), max(peak_all{2}(loud)));
bin_edges = {0:0.05:1.25, 0:0.04:1.12};  thr_now = [P.acf_thr P.cep_thr];
axis_name = {'taper-corrected autocorrelation peak','cepstral peak c(q*)'};
for m = 1:2
    ax = nexttile(tl);  hold(ax,'on');
    h1 = histogram(ax,peak_all{m}(loud & ref_all == 0),bin_edges{m},'FaceColor',C.muted,'EdgeColor','w','FaceAlpha',0.9);
    h2 = histogram(ax,peak_all{m}(loud & ref_all > 0),bin_edges{m},'FaceColor',method_colour{m},'EdgeColor','w','FaceAlpha',0.75);
    xline(ax,thr_now(m),'--','Color',C.ref,'LineWidth',0.9);
    xlabel(ax,axis_name{m});  ylabel(ax,'frames');
    yl = ylim(ax);  ylim(ax,[0 1.35*yl(2)]);
    title(ax,sprintf('threshold %.2f (dashed)',thr_now(m)),'FontWeight','normal');
    legend(ax,[h1 h2],{'reference unvoiced','reference voiced'},'Location','northeast','Box','off');
end
ax = nexttile(tl);  hold(ax,'on');
h1 = plot(ax,sweep_acf.uv_v(:,3),sweep_acf.v_uv(:,3),'-','Color',C.acf);
h2 = plot(ax,sweep_cep.uv_v(:,3),sweep_cep.v_uv(:,3),'-','Color',C.cep);
plot(ax,score_acf(3).uv_v,score_acf(3).v_uv,'o','MarkerSize',6,'MarkerFaceColor',C.acf,'MarkerEdgeColor','w');
plot(ax,score_cep(3).uv_v,score_cep(3).v_uv,'o','MarkerSize',6,'MarkerFaceColor',C.cep,'MarkerEdgeColor','w');
xlabel(ax,'UV->V error (%)');  ylabel(ax,'V->UV error (%)');  xlim(ax,[0 25]);  ylim(ax,[0 25]);
title(ax,'threshold sweep, dot = chosen','FontWeight','normal');
legend(ax,[h1 h2],method_name,'Location','northeast','Box','off');
save_figure(fig,out_dir,'fig10_voicing');

% ---- Fig 13: agreement between the two methods on s3 and s4 ----
fig = new_figure(17,8.2);  tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for i = 3:4
    ax = nexttile(tl);  hold(ax,'on');
    a = track_acf{i};  c = track_cep{i};  both = a > 0 & c > 0;
    f_line = [50 450];
    patch(ax,[f_line fliplr(f_line)],[f_line*0.95 fliplr(f_line)*1.05],C.muted,'FaceAlpha',0.25,'EdgeColor','none');
    plot(ax,f_line,f_line,'-','Color',C.muted,'LineWidth',0.8);
    scatter(ax,a(both),c(both),14,'o','MarkerFaceColor',C.ref,'MarkerEdgeColor','w','MarkerFaceAlpha',0.7,'LineWidth',0.4);
    lo = 0.9*min([a(both); c(both)]);  hi = 1.1*max([a(both); c(both)]);
    xlim(ax,[lo hi]);  ylim(ax,[lo hi]);
    xlabel(ax,'autocorrelation estimate (Hz)');  ylabel(ax,'cepstrum estimate (Hz)');
    title(ax,sprintf('%s: %d frames voiced in both, %.0f %% within 5 %% (band)',names{i},sum(both),agree(i,3)), ...
        'FontWeight','normal');
end
save_figure(fig,out_dir,'fig13_agreement_s3_s4');

% ---- Fig 14: noise test ----
fig = new_figure(17,6.8);  tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
noise_data = {noise_raw, noise_vde, noise_ffe};
clean_data = {clean.raw_gross(3,:), clean.vde(3,:), clean.ffe(3,:)};
panel_name = {'raw gross pitch error (%)','voicing decision error (%)','F0 frame error (%)'};
x_pos = 1:numel(snr_levels_dB)+1;
x_lab = [arrayfun(@num2str,snr_levels_dB,'UniformOutput',false), {'clean'}];
for q = 1:3
    ax = nexttile(tl);  hold(ax,'on');
    for m = 1:2
        mu = [mean(noise_data{q}(:,m,:),3); clean_data{q}(m)];
        lo = [min(noise_data{q}(:,m,:),[],3); clean_data{q}(m)];
        hi = [max(noise_data{q}(:,m,:),[],3); clean_data{q}(m)];
        patch(ax,[x_pos fliplr(x_pos)],[lo' fliplr(hi')],method_colour{m},'FaceAlpha',0.15,'EdgeColor','none');
        hn(m) = plot(ax,x_pos,mu,'-o','Color',method_colour{m},'MarkerSize',4, ...
            'MarkerFaceColor',method_colour{m},'MarkerEdgeColor','w'); %#ok<SAGROW>
    end
    if q == 1
        hn(3) = plot(ax,x_pos,[mean(noise_raw_nolpf,2); clean_nolpf.raw_gross(3,1)],'--^','Color',C.acf, ...
            'MarkerSize',4,'MarkerFaceColor',C.acf,'MarkerEdgeColor','w');
    end
    xticks(ax,x_pos);  xticklabels(ax,x_lab);  xlim(ax,[0.7 x_pos(end)+0.3]);
    yl = ylim(ax);  ylim(ax,[0 yl(2)]);
    xlabel(ax,'SNR (dB)');  title(ax,panel_name{q},'FontWeight','normal');
    if q == 1, legend(ax,hn,[method_name {'autocorr., no low-pass'}],'Location','northeast','Box','off'); end
end
save_figure(fig,out_dir,'fig14_noise');

% ---- Fig 15: sensitivity to window length and clipping level ----
fig = new_figure(17,6.8);  tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);  hold(ax,'on');
h1 = plot(ax,win_ms_grid,win_raw(:,1),'-o','Color',C.acf,'MarkerSize',4,'MarkerFaceColor',C.acf,'MarkerEdgeColor','w');
h2 = plot(ax,win_ms_grid,win_raw(:,3),'--^','Color',C.acf,'MarkerSize',4,'MarkerFaceColor',C.acf,'MarkerEdgeColor','w');
h3 = plot(ax,win_ms_grid,win_raw(:,2),'-o','Color',C.cep,'MarkerSize',4,'MarkerFaceColor',C.cep,'MarkerEdgeColor','w');
h4 = plot(ax,win_ms_grid,win_raw(:,4),'--^','Color',C.cep,'MarkerSize',4,'MarkerFaceColor',C.cep,'MarkerEdgeColor','w');
xlabel(ax,'window length (ms)');  ylabel(ax,'raw gross error (%)');  ylim(ax,[0 max(win_raw(:))+5]);
title(ax,'window length: pitch estimator','FontWeight','normal');
legend(ax,[h1 h2 h3 h4],{'ACF, s1','ACF, s2','cepstrum, s1','cepstrum, s2'},'Location','northeast','Box','off');
ax = nexttile(tl);  hold(ax,'on');
h1 = plot(ax,win_ms_grid,win_ffe(:,1),'-o','Color',C.acf,'MarkerSize',4,'MarkerFaceColor',C.acf,'MarkerEdgeColor','w');
h2 = plot(ax,win_ms_grid,win_ffe(:,2),'-o','Color',C.cep,'MarkerSize',4,'MarkerFaceColor',C.cep,'MarkerEdgeColor','w');
xlabel(ax,'window length (ms)');  ylabel(ax,'F0 frame error, s1+s2 (%)');  ylim(ax,[0 max(win_ffe(:))+3]);
title(ax,'window length: whole detector','FontWeight','normal');
legend(ax,[h1 h2],method_name,'Location','north','Box','off');
ax = nexttile(tl);  hold(ax,'on');
h1 = plot(ax,100*clip_grid,clip_raw(:,1),'-o','Color',C.acf,'MarkerSize',4,'MarkerFaceColor',C.acf,'MarkerEdgeColor','w');
h2 = plot(ax,100*clip_grid,clip_raw(:,2),'--^','Color',C.acf,'MarkerSize',4,'MarkerFaceColor',C.acf,'MarkerEdgeColor','w');
xline(ax,100*P.clip_frac,':','Color',C.ref,'LineWidth',0.9);
xlabel(ax,'centre clipping level (% of peak)');  ylabel(ax,'raw gross error (%)');  ylim(ax,[0 max(clip_raw(:))+3]);
title(ax,'ACF: clipping level (dotted = used)','FontWeight','normal');
legend(ax,[h1 h2],{'s1','s2'},'Location','north','Box','off');
save_figure(fig,out_dir,'fig15_sensitivity');

fprintf('\nFigures written to %s\n', out_dir);
diary off;
reset(groot);                                          % give back the default plot settings
%% Local functions

function F = make_frames(x, win_len, P)
% Row k holds the win_len samples centred on the k-th 10 ms segment.
half   = floor(win_len/2);
xp     = [zeros(half+P.hop,1); x(:); zeros(win_len+P.hop,1)];   % one extra hop of zeros each side
centre = (0:P.num_frames-1)'*P.hop + P.offset + P.hop;
F      = xp(centre + (1:win_len));
end

function E = frame_energy_dB(x, P)
% Frame energy in dB relative to the loudest frame of the recording.
F = make_frames(x, P.energy_win, P);
E = 10*log10(mean(F.^2,2) + eps);
E = E - max(E);
end

function d = parabolic_offset(Y, col)
% Sub-sample position of a peak from the three points around column col.
row = (1:size(Y,1))';
a = Y(sub2ind(size(Y),row,col-1));
b = Y(sub2ind(size(Y),row,col));
c = Y(sub2ind(size(Y),row,col+1));
d = 0.5*(a-c)./(a-2*b+c);
d(~isfinite(d) | abs(d) > 0.5) = 0;    % not a true local maximum (edge of the search range)
end

function [f0, peak, R, D] = pitch_acf(x, P)
% Autocorrelation pitch detector: low-pass, centre clip, autocorrelate, pick peak.
x_lp  = conv(x, P.lpf, 'same');                          % odd symmetric FIR -> zero delay
F     = make_frames(x_lp, P.acf_win, P);
third = floor(P.acf_win/3);
CL    = P.clip_frac * min(max(abs(F(:,1:third)),[],2), max(abs(F(:,end-third+1:end)),[],2));
Fc    = sign(F) .* max(abs(F)-CL, 0);                    % centre clipper
nfft  = 2^nextpow2(2*P.acf_win);
R     = real(ifft(abs(fft(Fc,nfft,2)).^2, [], 2));       % autocorrelation via the power spectrum
R     = R(:,1:P.lag_max+2) ./ max(R(:,1),eps);           % lags 0..lag_max+1, r(0) = 1
[r_peak,k] = max(R(:,P.lag_min+1:P.lag_max+1), [], 2);
lag   = k + P.lag_min - 1;
peak  = r_peak .* P.acf_win ./ (P.acf_win - lag);        % undo the (N-k)/N taper: a periodic frame scores 1 at any pitch
f0    = P.fs ./ (lag + parabolic_offset(R, lag+1));
D     = struct('F',F,'Fc',Fc,'CL',CL,'r_peak',r_peak);   % intermediate signals, for the figures
end

function [f0, peak, C, D] = pitch_cep(x, P)
% Cepstrum pitch detector: window, log magnitude spectrum, inverse FFT, pick peak.
F  = make_frames(x, P.cep_win, P) .* hamming(P.cep_win)';
C  = real(ifft(log(abs(fft(F,P.nfft,2)) + eps), [], 2));    % real cepstrum
C  = C(:,1:P.lag_max+2);                                    % quefrency 0..lag_max+1 samples
[peak,k] = max(C(:,P.lag_min+1:P.lag_max+1), [], 2);
q  = k + P.lag_min - 1;
q  = q + parabolic_offset(C, q+1);
f0 = P.fs ./ q;
D  = struct('F',F);
end

function f0 = make_track(f0_raw, peak, E, thr, P)
% Voicing decision followed by median smoothing.
voiced = (peak >= thr) & (E >= P.energy_floor_dB);
voiced_s = medfilt1(double(voiced), P.med_len) > 0.5;       % majority vote on voicing
f0 = f0_raw;  f0(~voiced) = NaN;
half = (P.med_len-1)/2;  n = numel(f0);
W  = [nan(half,1); f0; nan(half,1)];
W  = sort(W((1:n)' + (0:P.med_len-1)), 2);                  % each row: the neighbours in order, NaN last
m  = max(ceil(sum(~isnan(W),2)/2), 1);                      % middle voiced value; the lower one if the count is even,
f0 = W(sub2ind(size(W),(1:n)',m));                          % so the output is always a measured candidate
f0(~voiced_s | isnan(f0)) = 0;
end

function S = score_track(est, ref, tol)
% Voicing and pitch error measures of an estimated track against a reference.
v_ref = ref > 0;  v_est = est > 0;  both = v_ref & v_est;
err   = est(both) - ref(both);
gross = abs(err) > tol*ref(both);
S.v_uv      = 100*sum(v_ref & ~v_est)/sum(v_ref);
S.uv_v      = 100*sum(~v_ref & v_est)/sum(~v_ref);
S.vde       = 100*mean(v_ref ~= v_est);
S.gpe       = 100*sum(gross)/max(sum(both),1);
S.fine_mean = mean(err(~gross));
S.fine_std  = std(err(~gross));
S.fine_mae  = mean(abs(err(~gross)));
S.ffe       = 100*(sum(v_ref ~= v_est) + sum(gross))/numel(ref);
S.n_both    = sum(both);
S.n_gross   = sum(gross);
end

function print_score(file, method, tag, S)
fprintf('%-6s %-9s %-9s %7.1f %7.1f %7.1f %7.1f %9.2f %8.2f %8.2f %7.1f\n', file, method, tag, ...
    S.v_uv, S.uv_v, S.vde, S.gpe, S.fine_mean, S.fine_std, S.fine_mae, S.ffe);
end

function W = sweep_threshold(raw, peak, E, ref, grid, P)
% Voicing and frame errors as the voicing threshold is varied.
% Columns: s1, s2, s1+s2 pooled.
for j = 1:numel(grid)
    est = {make_track(raw{1},peak{1},E{1},grid(j),P), make_track(raw{2},peak{2},E{2},grid(j),P)};
    for i = 1:3
        if i < 3, S = score_track(est{i}, ref{i}, P.gross_tol);
        else,     S = score_track([est{1}; est{2}], [ref{1}; ref{2}], P.gross_tol); end
        W.v_uv(j,i) = S.v_uv;  W.uv_v(j,i) = S.uv_v;  W.ffe(j,i) = S.ffe;
    end
end
end

function out = run_pair(speech, reference, P)
% Runs both detectors on s1 and s2 with the settings in P.
% Every output is 3 x 2: rows s1, s2, pooled; columns ACF, cepstrum.
est_a = cell(1,3);  est_c = cell(1,3);  n_gross = zeros(2,2);  n_voiced = zeros(2,1);
for i = 1:2
    E = frame_energy_dB(speech{i},P);
    [fa,pa] = pitch_acf(speech{i},P);  [fc,pc] = pitch_cep(speech{i},P);
    v = reference{i} > 0;  tol = P.gross_tol*reference{i}(v);
    n_gross(i,:) = [sum(abs(fa(v)-reference{i}(v)) > tol), sum(abs(fc(v)-reference{i}(v)) > tol)];
    n_voiced(i)  = sum(v);
    est_a{i} = make_track(fa,pa,E,P.acf_thr,P);  est_c{i} = make_track(fc,pc,E,P.cep_thr,P);
end
est_a{3} = [est_a{1}; est_a{2}];  est_c{3} = [est_c{1}; est_c{2}];
ref = {reference{1}, reference{2}, [reference{1}; reference{2}]};
out.raw_gross = 100*[n_gross; sum(n_gross,1)] ./ [n_voiced; sum(n_voiced)];
for i = 1:3
    Sa = score_track(est_a{i}, ref{i}, P.gross_tol);  Sc = score_track(est_c{i}, ref{i}, P.gross_tol);
    out.vde(i,:) = [Sa.vde Sc.vde];  out.ffe(i,:) = [Sa.ffe Sc.ffe];
end
end

function fig = new_figure(width_cm, height_cm)
fig = figure('Visible','off','Color','w','Units','centimeters','Position',[1 1 width_cm height_cm]);
if isprop(fig,'Theme'), fig.Theme = 'light'; end
end

function save_figure(fig, out_dir, name)
exportgraphics(fig, fullfile(out_dir,[name '.png']), 'Resolution', 220);
close(fig);
end

function v = gaps(v)
% Zero (unvoiced) becomes NaN so that plotted lines break there.
v(v == 0) = NaN;
end

function draw_chain(ax, y, labels, colour, heading, result)
% One row of processing blocks joined by arrows.
w = 1.30;  gap = 0.25;  h = 0.80;  ink = [0.25 0.25 0.25];
text(ax,-0.32,y,'x(n)','HorizontalAlignment','right','FontSize',8);
for k = 0:numel(labels)
    x_end = k*(w+gap);                                  % arrow that ends at block k+1
    plot(ax,[x_end-gap+0.03 x_end-0.07],[y y],'Color',ink,'LineWidth',0.9);
    patch(ax,x_end-[0.11 0.11 0.01],y+[0.07 -0.07 0],ink,'EdgeColor','none');
    if k < numel(labels)
        rectangle(ax,'Position',[x_end y-h/2 w h],'Curvature',0.18,'EdgeColor',colour,'FaceColor','w','LineWidth',1.2);
        text(ax,x_end+w/2,y,labels{k+1},'HorizontalAlignment','center','FontSize',6.5);
    end
end
text(ax,numel(labels)*(w+gap)+0.04,y,result,'FontSize',8);
text(ax,0,y+h/2+0.12,heading,'FontWeight','bold','FontSize',8.5,'VerticalAlignment','bottom');
end

function fig = plot_tracking(x, name, ref, est_acf, est_cep, frame_time, P, C)
% Waveform, narrow-band spectrogram, pitch tracks and voicing for one recording.
fs = P.fs;  t = (0:numel(x)-1)'/fs;
fig = new_figure(17,15);  tl = tiledlayout(fig,8,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl,[1 1]);
plot(ax1,t,x,'Color',C.ref,'LineWidth',0.4);
ylim(ax1,[-1 1]);  ylabel(ax1,'amplitude');  ax1.XTickLabel = [];
title(ax1,sprintf('%s.wav',name),'FontWeight','normal');

ax2 = nexttile(tl,[3 1]);  hold(ax2,'on');
[S,f,ts] = spectrogram(x, hamming(P.cep_win), P.cep_win-P.hop/2, P.nfft, fs);
S_dB = 20*log10(abs(S)+eps);
imagesc(ax2,ts,f,S_dB);  axis(ax2,'xy');  grid(ax2,'off');
clim(ax2,max(S_dB(:))+[-55 0]);  colormap(ax2,1-0.8*linspace(0,1,128)'*[1 1 1]);
h1 = plot(ax2,frame_time,gaps(est_acf),'o','Color',C.acf,'MarkerSize',3.2,'LineWidth',0.7);
h2 = plot(ax2,frame_time,gaps(est_cep),'.','Color',C.cep,'MarkerSize',6);
ylim(ax2,[0 1000]);  ylabel(ax2,'frequency (Hz)');  ax2.XTickLabel = [];
legend(ax2,[h1 h2],{'autocorrelation','cepstrum'},'Location','northeast','Orientation','horizontal','Box','on');
title(ax2,'narrow-band spectrogram (40 ms window) with both estimates on the first harmonic','FontWeight','normal');

ax3 = nexttile(tl,[3 1]);  hold(ax3,'on');
hh = gobjects(0);  leg = {};
if ~isempty(ref)
    hh(end+1) = plot(ax3,frame_time,gaps(ref),'-','Color',[0.78 0.78 0.78],'LineWidth',4.5);  leg{end+1} = 'reference';
end
hh(end+1) = plot(ax3,frame_time,gaps(est_acf),'-o','Color',C.acf,'MarkerSize',2.6,'LineWidth',0.8);  leg{end+1} = 'autocorrelation';
hh(end+1) = plot(ax3,frame_time,gaps(est_cep),'-','Color',C.cep,'Marker','.','MarkerSize',6,'LineWidth',0.8);  leg{end+1} = 'cepstrum';
all_f0 = [est_acf(est_acf>0); est_cep(est_cep>0); ref(ref>0)];
ylim(ax3,[max(50,0.85*min(all_f0)) 1.15*max(all_f0)]);  ylabel(ax3,'pitch (Hz)');  ax3.XTickLabel = [];
legend(ax3,hh,leg,'Location','northeast','Orientation','horizontal','Box','on');
title(ax3,'pitch track, one value every 10 ms','FontWeight','normal');

ax4 = nexttile(tl,[1 1]);
rows = {est_cep > 0, est_acf > 0};  row_colour = {C.cep, C.acf};  row_name = {'cepstrum','ACF'};
if ~isempty(ref), rows{end+1} = ref > 0;  row_colour{end+1} = [0.45 0.45 0.45];  row_name{end+1} = 'reference'; end
strip = ones(numel(rows),numel(frame_time),3);
for r = 1:numel(rows)
    for ch = 1:3, strip(r,rows{r},ch) = row_colour{r}(ch); end
end
image(ax4,frame_time,1:numel(rows),strip);  axis(ax4,'xy');  grid(ax4,'off');
yticks(ax4,1:numel(rows));  yticklabels(ax4,row_name);  ylim(ax4,[0.5 numel(rows)+0.5]);
xlabel(ax4,'time (s)');  title(ax4,'frames declared voiced','FontWeight','normal');

linkaxes([ax1 ax2 ax3 ax4],'x');  xlim(ax1,[0 3]);
end
