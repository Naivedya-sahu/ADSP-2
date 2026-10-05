%% ELL7420 Design Assignment 2: Pitch detection by autocorrelation
% This is the time-domain part of our pitch detection assignment.
% Run this file with the 'speech samples' folder beside it. It reads s1-s4,
% checks s1 and s2 against the given pitch values, and saves the plots/CSV files.

%% 1. Find the recordings and prepare the results folder
base = fileparts(mfilename('fullpath'));
% Check the common folder names next to the Live Script and one level up.
parents = {base, fullfile(base,'..')};
names = {'speech samples','speech_samples'};
data = '';
for p = 1:numel(parents)
    for n = 1:numel(names)
        folder = fullfile(parents{p},names{n});
        options = {folder,fullfile(folder,'speech_samples'), ...
                   fullfile(folder,'speech samples')};
        for q = 1:numel(options)
            if isfile(fullfile(options{q},'s1.wav'))
                data = options{q};
                break;
            end
        end
        if ~isempty(data), break; end
    end
    if ~isempty(data), break; end
end
if isempty(data)
    % If the folder is elsewhere, select the folder that has s1.wav inside.
    data = uigetdir(base,'Select the folder containing s1.wav');
    assert(~isequal(data,0),'No speech samples folder was selected');
end
assert(isfile(fullfile(data,'s1.wav')), ...
       'Please select the folder that directly contains s1.wav');
fprintf('Reading speech files from: %s\n',data);
out = fullfile(base,'results');
plots = fullfile(out,'report_plots');
if ~isfolder(plots), mkdir(plots); end

%% 2. Estimate pitch for each speech sample
for s = 1:4
    [x,fs] = audioread(fullfile(data,sprintf('s%d.wav',s)));
    x = mean(x,2); % Average channels if the recording is stereo.
    if s == 1, show_method(x,fs,plots); end
    [time,pitch] = track_pitch(x,fs);
    reference = nan(size(pitch));

    % Only s1 and s2 have a supplied pitch track. Zero means unvoiced.
    if s <= 2
        key = sprintf('p%d',s);
        saved = load(fullfile(data,[key '.mat']));
        reference = saved.(key); reference = reference(:);
        assert(numel(reference)==numel(pitch),'Reference length does not match');
        actual = reference > 0; detected = pitch > 0;
        tp = sum(actual & detected); fp = sum(~actual & detected);
        fn = sum(actual & ~detected); tn = sum(~actual & ~detected);
        matched = actual & detected; % Compare pitch error where both say voiced.
        mae = mean(abs(pitch(matched)-reference(matched)));
        within20 = sum(abs(pitch(matched)-reference(matched)) <= ...
                       .2*reference(matched));
        fprintf(['s%d: TP=%d FP=%d FN=%d TN=%d | precision=%.1f%% ' ...
                 'recall=%.1f%% | voiced MAE=%.2f Hz | within 20%%=%d/%d\n'], ...
                 s,tp,fp,fn,tn,100*tp/max(tp+fp,1), ...
                 100*tp/max(tp+fn,1),mae,within20,sum(actual));
    else
        % There is no ground truth for s3/s4, so we can only report our estimates.
        fprintf('s%d: %d voiced of %d frames; no reference supplied\n', ...
                s,sum(pitch>0),numel(pitch));
    end

    % 3. Save frame-wise values and the pitch track.
    % The reference column is NaN for s3/s4 because no reference was given.
    tableOut = table(time,pitch,reference,'VariableNames', ...
                     {'time_s','estimated_hz','reference_hz'});
    writetable(tableOut,fullfile(out,sprintf('s%d_autocorrelation.csv',s)));
    fig = figure('Color','w','Name',sprintf('%d. s%d pitch track',s+2,s));
    plot(time,pitch,'b','LineWidth',1.2); hold on;
    if s <= 2
        plot(time,reference,'r--');
        legend('Autocorrelation','Provided reference','Location','best');
    end
    xlabel('Time (s)'); ylabel('Pitch (Hz)');
    title(sprintf('%d. Pitch track for s%d',s+2,s));
    xlim([0 numel(x)/fs]); ylim([0 450]); grid on;
    saveas(fig,fullfile(plots,sprintf('%02d_s%d_pitch_track.png',s+2,s)));
end

function [time,pitch] = track_pitch(x,fs)
    % Work every 10 ms so our estimates line up with p1/p2.
    % Use 0 Hz when a frame does not pass the voiced tests.
    hop = round(.01*fs);
    count = ceil(numel(x)/hop);
    time = ((0:count-1)'*hop+floor(hop/2))/fs;
    pitch = zeros(count,1);
    for i = 1:count
        [~,lags,r,energy] = frame_correlation(x,fs,i);
        if energy < .008, continue; end % Too quiet to give a useful pitch.
        [period,strength] = best_period(lags,r);
        % A regular speech cycle should give a strong correlation peak.
        if strength >= .50, pitch(i) = fs/period; end
    end
end

function [frame,lags,r,energy] = frame_correlation(x,fs,i)
    % Take 40 ms around this frame's time point, padding at the ends.
    L = round(.04*fs); hop = round(.01*fs);
    index = (i-1)*hop+floor(hop/2)-floor(L/2)+(0:L-1)';
    valid = index>=0 & index<numel(x);
    frame = zeros(L,1); frame(valid) = x(index(valid)+1);
    % Remove the mean and taper the edges before comparing delayed copies.
    window = .54-.46*cos(2*pi*(0:L-1)'/(L-1));
    frame = (frame-mean(frame)).*window;
    energy = sqrt(mean(frame.^2));
    % Pitch = fs/lag, so these lags cover roughly 70 to 400 Hz.
    lags = (ceil(fs/400):min(L-2,floor(fs/70)))';
    r = zeros(size(lags));
    for j = 1:numel(lags)
        k = lags(j); a = frame(1:end-k); b = frame(k+1:end);
        % Normalise so loudness has less effect on the peak height.
        r(j) = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2)+1e-12);
    end
end

function [period,strength] = best_period(lags,r)
    % Find local peaks, then take the first one close to the strongest peak.
    % This favours the shortest convincing period over its multiples.
    peaks = find(r(2:end-1)>r(1:end-2) & r(2:end-1)>=r(3:end))+1;
    period = NaN; strength = 0;
    if isempty(peaks), return; end
    peaks = peaks(r(peaks)>=.82*max(r(peaks)));
    j = peaks(1); strength = r(j);
    % Fit around the peak to estimate a lag between integer samples.
    curvature = r(j-1)-2*r(j)+r(j+1);
    shift = 0;
    if abs(curvature)>1e-12
        shift = .5*(r(j-1)-r(j+1))/curvature;
        shift = max(-.5,min(.5,shift));
    end
    period = lags(j)+shift;
end

function show_method(x,fs,folder)
    % Two example frames from s1 show how we made the voiced decision.
    [voice,lags,r,voiceRms] = frame_correlation(x,fs,176); % Near 1.755 s
    [quiet,~,rq,quietRms] = frame_correlation(x,fs,21);   % Near 0.205 s
    [period,~] = best_period(lags,r);
    k = round(period); L = numel(voice);

    % Figure 1 follows the signal -> delayed copy -> correlation peak.
    fig = figure('Color','w','Name','1. Autocorrelation extraction');
    subplot(3,1,1);
    plot((0:numel(x)-1)/fs,x); hold on;
    xline(1.735,'r--'); xline(1.775,'r--');
    title('1. Select a 40 ms voiced frame in s1');
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;
    subplot(3,1,2);
    t = (0:L-k-1)/fs*1000;
    plot(t,voice(1:end-k)); hold on; plot(t,voice(k+1:end));
    title(sprintf('2. Compare with a copy delayed by %d samples',k));
    xlabel('Time in frame (ms)'); ylabel('Windowed amplitude');
    legend('Original','Delayed copy'); grid on;
    subplot(3,1,3);
    plot(lags,r); hold on; yline(.5,'m--');
    plot(k,r(lags==k),'ro','MarkerFaceColor','r');
    title(sprintf('3. Correlation peak: lag %.2f, pitch %.1f Hz',period,fs/period));
    xlabel('Lag (samples)'); ylabel('Normalized correlation'); grid on;
    saveas(fig,fullfile(folder,'01_autocorrelation_extraction_s1.png'));

    % Figure 2 compares a clear voiced frame with a quiet frame.
    fig = figure('Color','w','Name','2. Voiced and silent frames');
    frames = {voice,quiet}; curves = {r,rq};
    levels = [voiceRms quietRms]; names = {'Voiced','Silent'};
    for c = 1:2
        subplot(2,2,c);
        plot((0:L-1)/fs*1000,frames{c}); grid on;
        title(sprintf('%s frame: RMS %.4f',names{c},levels(c)));
        xlabel('Time (ms)'); ylabel('Amplitude');
        subplot(2,2,c+2);
        plot(lags,curves{c}); hold on; yline(.5,'m--'); grid on;
        xlabel('Lag (samples)'); ylabel('Normalized correlation');
    end
    sgtitle('Voiced frame receives pitch; silent frame receives 0 Hz');
    saveas(fig,fullfile(folder,'02_voiced_unvoiced_decision_s1.png'));
end
