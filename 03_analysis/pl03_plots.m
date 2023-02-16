cfg=[]; 
cfg.dataset = [inpath{p},indat(v).name];
cfg.trialdef.length = Inf; % Continuous data

cfg=ft_definetrial(cfg); 

% 2.2. Preprocessing
data_p = ft_preprocessing(cfg); % data_p is the preprocessed data

%%
data_app = data_pscr;
for t = 2:length(data_pscr.trial)
    data_app.trial{1} = [data_app.trial{1}, data_pscr.trial{t}];
    data_app.time{1} = [data_app.time{1}, data_pscr.time{t}];
end

%%
figure;
subplot(4,1,1),plot(data_p.time{1},data_p.trial{1}(1,:),'linewidth',2); 
    title('Raw Data')
    xlabel('Time (s)'); ylabel('mV')
subplot(4,1,2),plot(data_ps.time{1},data_ps.trial{1},'linewidth',2); 
    title('0.1 - 23 Hz Filtered Data')
    xlabel('Time (s)'); ylabel('mV')
subplot(4,1,3),plot(data_app.time{1},data_app.trial{1},'linewidth',2); 
    title('Cleaned Data')
    xlabel('Time (s)'); ylabel('mV')
subplot(4,1,4),plot(data_fft.freq,data_fft.powspctrm,'linewidth',2); 
    title('Frequency Transformed Data')
    xlim([data_fft.freq(1) data_fft.freq(end)]);
    ; xlabel('Frequency (Hz)'); ylabel('Amplitude')