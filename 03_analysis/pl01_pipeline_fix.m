%% Pipeline to process the Planarian Data
% Julian Keil 08.12.2021
%
% Steps:
% 0. Set Basics
% 1. Define Data in Folder
% 2. Loop Data
% 2.1. Define current dataset
% 2.2. Read in continuous data, filter and demean
% 2.3. Compute variance for each channel across time to remove empty
% channels
% 
% TO DO:
% Check if limits make sense
%% 0. Set Basics
% This script requires the FieldTrip Toolbox: https://www.fieldtriptoolbox.org/
clear all
close all
clc

% Set limits
minvar = 0.1; % Minimum required standard deviation to count as actual data

%% 1. Define Data in Folder
inpath = '/Plattwurm/roh/';
indat = dir([inpath,'*.bdf']);
outpath = '/Plattwurm/preproc/';

%% 2. Loop Data
doc = {1,1}
for v = 1:length(indat) % Loop around all participants
    %% 2.1. Define current dataset
    cfg=[]; 
    cfg.dataset = [inpath,indat(v).name];
    cfg.trialdef.length = Inf; % Continuous data

    cfg=ft_definetrial(cfg); 

    %% 2.2. Preprocessing
    %cfg.demean = 'yes';
    %cfg.detrend = 'yes';
    cfg.hpfilter = 'yes';
    cfg.hpfreq = .5;
    cfg.hpfilttype = 'firws';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 20; % Data is heavily contaminated by Line noise
    cfg.lpfilttype = 'firws';

    data_p = ft_preprocessing(cfg); % data_p is the preprocessed data

    %% 2.3. Compute variance across time to identify empty channels
    % 2.3.1. Compute variance
    var_l = zeros(1,length(data_p.label));
    for l = 1:length(data_p.label)
        var_l(l) = std(data_p.trial{1}(l,:)); % Check Standard Deviation
    end
    goodchan = var_l > minvar; % Good channels are those with a variance unequal to zero
    
    % 2.3.2. Keep only good channels and remove first and last 500ms to avoid edge
    % artifacts
    cfg = [];
    cfg.channel = data_p.label(goodchan);
    cfg.latency = [data_p.time{1}(251) data_p.time{1}(end-251)];
    
    data_ps = ft_selectdata(cfg,data_p); % data_ps is the, preprocessed, selected data
    
    %% 2.4. Cut into smaller segments
    % 2.4.1. Cut my data into pieces   
    cfg = [];
    cfg.length = 10; % 10 seconds
    cfg.overlap = .25; % 2.5s overlap to get a continuous change in the spectrum
    data_psc = ft_redefinetrial(cfg,data_ps);
    %% 0. Set the cfgs
    if exist('cfg')
        % Trials
        if isfield(cfg,'keeptrials')
            if strcmpi(cfg.keeptrials,'no')
                trial_flag = 1;
            elseif strcmpi(cfg.keeptrials,'yes')
                meth_flag = 2;
            end
        else
            trial_flag = 1; % Default: Don't keep trials
        end
        %Channels
        if isfield(cfg,'keepchannels')
            if strcmpi(cfg.keepchannels,'no')
                chan_flag = 1;
            elseif strcmpi(cfg.keepchannels,'yes')
                chan_flag = 2;
            end
        else
            chan_flag = 1; % Default: Don't keep channels
        end
        % Thresholds
        if isfield(cfg,'threshold')
            thresh = cfg.threshold;
        else
            thresh = 2.5; % Default is 2.5 STD above the mean
        end
        % Peaks
        if isfield(cfg,'peak')
            peak = cfg.peak;
        else
            peak = 1000; % Default is 1000mV
        end
    end % exist

    % To start, define all trials and channels as good.
    goodchannels = 1:length(data_psc.label);
    goodtrials = 1:length(data_psc.trial);
    %% 1. First. remove the Channels 
    % We'll restrict this to outliers based on 2.5 STD above the mean
    % First, compute the std across trials for each channel
    % Then, compute the mean and std of the stds across channels
    
    if chan_flag == 1
        % 1.1. STD across trials
        tmp_c = zeros(1,length(data_psc.trial)*length(data_psc.trial{1})); % Build an empty vector to save time
        for c = 1:length(data_psc.label) % Loop channels
            ti = 1; % Index of trials
            for t = 1:length(data_psc.trial) % Loop trials
                tmp_c(ti:ti+length(data_psc.trial{t})-1) = data_psc.trial{t}(c,:); % Stitch all trials together
                ti = ti+length(data_psc.trial{t}); % Move the trial index
            end
            std_c(c) = std(tmp_c); % And compute the std across trials
        end

        %1.1 Compute std for each trial
        m_std = 0;
        s_std = 0;
        trialstd = (1);
        trialmean = (1);
        for h = 1:length(data_psc.label)
            for i = 1:length(data_psc.trial)
                trialstd(i,h) = std(data_psc.trial{1,i}(h:length(data_psc.trial{1,1})));
                trialmean(i,h) = mean(data_psc.trial{1,i}(h:length(data_psc.trial{1,1})));
            end
            m_std(1,h) = mean(trialstd(1:length(data_psc.trial),h)); % Mean of STD across trials
            s_std(1,h) = std(trialstd(1:length(data_psc.trial),h)); % STD of STD across trials
            
            goodchannels(h) = std_c(h) <= m_std(h) + (thresh*s_std(h)); % Mean + 2.5*STD
            goodtrials(h,1:length(trialstd)) = transpose(trialstd(1:length(trialstd),h) <= m_std(h) + (thresh*s_std(h)));
        end
        
        fprintf('\n Good channels: %d \n\n',sum(goodchannels)); 

        % 3.1.1.4. Keep only good channels
        cfg = [];
        cfg.channel = data_psc.label(logical(goodchannels));

        data_psc = ft_selectdata(cfg,data_psc);

        % Keep only good trials

        
        data_psc.trial = data_psc.trial(logical(goodtrials(1,1:length(goodtrials))));
        data_psc.time = data_psc.time(logical(goodtrials(1,1:length(goodtrials))));
%         data_psc.sampleinfo = data_psc.sampleinfo(logical(goodtrials));
%         if it doesn't work, transpose
        infostruct = 1:length(data_psc.sampleinfo);
        for i = 1:width(data_psc.sampleinfo)
        infostruct(i,1:length(data_psc.sampleinfo)) = goodtrials(1, 1:length(data_psc.sampleinfo));
        end
        data_psc.sampleinfo = data_psc.sampleinfo(transpose(logical(infostruct)));
    end


    if goodchannels > 0 
        %% 2.5. Frequency Analysis
        % Set the parameters
        cfg = [];
        cfg.foi = [10.^(log10(.5):((log10(20)-log10(.5))/35):log10(20))];
        cfg.pad = 'nextpow2';
        cfg.tapsmofrq = 1;
        cfg.method = 'mtmfft';
        % Compute the exponent
        cfg.output = 'fooof_aperiodic';
        fractal = ft_freqanalysis(cfg, data_psc);
        % Get the exponent
        for c = 1:length(data_psc.label)
            data_expo(c) = fractal.fooofparams(c).aperiodic_params(2); % second putput is exponent
        end
        % Get the Power Spectrum
        cfg.output = 'pow';
        data_fft = ft_freqanalysis(cfg, data_psc);

        % Save the data
        save([outpath,indat(v).name,'.mat'],'data_expo','data_fft','-V7.3');

        % Compute statistics
        trialstats_1(1, v) = convertCharsToStrings(indat(v).name);
        trialstats_2(v) = length(goodtrials);
        trialstats_3(v) = sum(goodtrials(1,1:length(goodtrials)));
        trialpercent(v) = (trialstats_3/trialstats_2)*100;
    else
        fprintf('No good channels');
    end
end % Move down

trialstats = {trialstats_1, trialstats_2, trialstats_3, trialpercent};
% Write trial statistics

