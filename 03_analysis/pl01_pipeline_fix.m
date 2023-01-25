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
ft_defaults;

% Set limits
minvar = 0.1; % Minimum required standard deviation to count as actual data

%% 1. Define Data in Folder
inpath{1} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/01_Raw/01_Dunkel/'; 
inpath{2} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/01_Raw/02_Hell/'; 
inpath{3} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/01_Raw/03_Mllebend/'; 
inpath{4} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/01_Raw/04_Mltot/'; 

outpath{1} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/01_Dunkel/'; 
outpath{2} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/02_Hell/'; 
outpath{3} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/03_Mllebend/'; 
outpath{4} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/04_Mltot/'; 

% inpath{1} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/01_Raw/05_HellDunkelHell/'; 
% inpath{2} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/01_Raw/06_KLDunkel/'; 
% inpath{3} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/01_Raw/07_KLHell/'; 
% 
% outpath{1} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/05_HellDunkelHell/'; 
% outpath{2} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/06_KLDunkel/'; 
% outpath{3} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/07_KLHell/'; 
%% 2. Loop Data
for p = 1:length(inpath)
    %% 2.1 Get the list of worms
    indat = dir([inpath{p},'*.bdf']);
    
    % 2.2 Loop worms
    for v = 1:length(indat)
        % 2.2.1. Define current dataset
        cfg=[]; 
        cfg.dataset = [inpath{p},indat(v).name];
        cfg.trialdef.length = Inf; % Continuous data

        cfg=ft_definetrial(cfg); 

        % 2.2. Preprocessing
        cfg.hpfilter = 'yes';
        cfg.hpfreq = .1;
        cfg.hpfilttype = 'firws';
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 23; % Data is heavily contaminated by Line noise
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
    cfg = [];
    cfg.length = 6; % 10 seconds
    %cfg.overlap = .25; % 3.75 s overlap to get a continuous change in the spectrum
    data_psc = ft_redefinetrial(cfg,data_ps);
  
    %% 2.5. Automatic artefact removal
        % 2.5.1 Get mean within trials for peak threshold
        clear mean_tmp mean_t
        for c = 1:length(data_psc.label)
            for t = 1:length(data_psc.trial)
               mean_tmp(c,t) = mean(data_psc.trial{t}(c,:));
            end
        end
        mean_t = mean(mean_tmp,1);
    
        % 2.5.2 Automatically reject trials and channels
        cfg = [];
        cfg.keepchannels = 'no';
        cfg.keeptrials = 'no';
        cfg.threshold = 2.5; % Threshold for variance etc.
        cfg.peak = mean(mean_t) + (5*std(mean_t)); % Threshold for absolute peaks

        data_pscr = vt_autoreject(cfg,data_psc);

    %% 2.6. If we have enough data, continue
    if  length(data_pscr.trial) > 5 & size(data_pscr.trial{1},1) > 0 
        % 2.6.1 Frequency Analysis
        % Set the parameters
        cfg = [];
        cfg.foi = [10.^(log10(.5):((log10(20)-log10(.5))/22):log10(20))];
        cfg.pad = 'nextpow2';
        cfg.tapsmofrq = 1;
        cfg.method = 'mtmfft';
        % Compute the exponent
        cfg.output = 'fooof_aperiodic';
        
        fractal = ft_freqanalysis(cfg, data_pscr);
        % Extract the exponent
        for c = 1:length(data_pscr.label)
            data_expo(c) = fractal.fooofparams(c).aperiodic_params(2); % second putput is exponent
        end
        
        % Get the Power Spectrum
        cfg.output = 'pow';
        data_fft = ft_freqanalysis(cfg, data_pscr);

        % Save the data
        sprintf(indat(v).name)
        save([outpath{p},indat(v).name,'.mat'],'data_expo','data_fft','-V7.3');

        % Compute statistics
        trialstats{p}(v,1) = convertCharsToStrings(indat(v).name);
        trialstats{p}(v,2) = size(data_pscr.trial{1},1);
    else
        fprintf('No enough data :-( \n');
    end
    end % Move down
end

% Write trial statistics

