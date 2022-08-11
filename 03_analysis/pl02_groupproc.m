%% Script to read in individual files, and process groups
% Julian Keil 11.02.2022
%
% Steps:
% 0. Set Basics
% 1. Define Data in Folder
% 2. Loop Data
% 2.1. Define current dataset
% 2.2. Read in preprocessed data
% 2.3. Group average
% 2.3.1. Save FFT as matix
% channels

%% 0. Set Basics
% This script requires the FieldTrip Toolbox: https://www.fieldtriptoolbox.org/
clear all
close all
clc

ft_defaults;
%% 1. Define Data in Folder
inpath{1} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/01_Dunkel/'; 
inpath{2} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/02_Hell/'; 
inpath{3} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/03_Mllebend/'; 
inpath{4} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/04_Mltot/'; 

outpath = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/05_Output/';

for p = 1:length(inpath)
    folders = find((inpath{p} == '/'));
    names{p} = [inpath{p}(folders(end-1)+1:folders(end)-1),'.csv'];
end

%% 2. Loop Data to collect single recordings
cond = cell(1,4);
for p = 1:length(inpath)
    indat = dir([inpath{p},'*.mat']);
    cond{p} = cell(1,1);
    for v = 1:length(indat)
        load([inpath{p},indat(v).name]);
        for c = 1:length(data_fft.label)
            cond{p}{v,c}.name = indat(v).name(1:end-8);
            cond{p}{v,c}.label = {'Fp1'};
            cond{p}{v,c}.dimord = 'chan_freq';
            cond{p}{v,c}.freq = data_fft.freq;
            cond{p}{v,c}.powspctrm = squeeze(data_fft.powspctrm(c,:));
            cond{p}{v,c}.exponent = data_expo(c);
        end
    end
    cond{p} = reshape(cond{p},[1, numel(cond{p})]);
    cond{p} = cond{p}(~cellfun('isempty',cond{p}));
end

%% 3. Group Average
for p = 1:length(cond)
    cfg = [];
    cfg.keepindividual = 'no';

    GA{p} = ft_freqgrandaverage(cfg,cond{p}{:});
end

% %% 4. Write CSV
% for p = 1:length(GA)
%     % 4.1. Make Table
%     clear tbl
%     tbl = array2table(squeeze(GA{p}.powspctrm));
% 
%     % 4.2. Prepare Column Names = Frequencies
%     colNames = {};
%     for c = 1:length(GA{p}.freq)
%         tmp = num2str(GA{p}.freq(c));
%         tmp = replace(tmp,'.','_');
%         colNames{c} = ['f_',tmp];
%     end
%     tbl.Properties.VariableNames = colNames;
% 
%     % 4.3. Add Row Names = Dataset Names
%     for v = 1:size(tbl,1)
%         tbl.Properties.RowNames{v} = cond{p}{v}.name;
%     end
% 
%     % 4.4 Now save it 
%     cd(outpath);
%     writetable(tbl,names{p},'WriteVariableNames',1,'WriteRowNames',1);
% end
%  
%% 5. Stats
    % 5.1 Define the Contrasts
    contrasts = [1 2;... % 1 dunkel vs hell
                 1 3;... % 2 dunkel vs mllebend
                 1 4;... % 3 dunkel vs mltoto
                 2 3;... % 4 hell vs mllebend
                 2 4;... % 5 hell vs mltot
                 3 4]; % 6 mllebend vs. mltot

    % 5.2. Loop the conditions fro pairwise comparisons
    for p = 1:length(contrasts)
        % Compare the spectra
        cfg=[];
        cfg.method = 'montecarlo';
        cfg.statistic ='indepsamplesT';
        cfg.numrandomization = 10000;
        cfg.correctm = 'fdr';
        cfg.neighbours = [];
        cfg.alpha = 0.01;
        cfg.tail = 0;
        cfg.correcttail = 'alpha';
        cfg.ivar = 1;

        cfg.design=[ones(1,length(cond{contrasts(p,1)})),ones(1,length(cond{contrasts(p,2)}))*2];

        pow_stats{p} = ft_freqstatistics(cfg,cond{contrasts(p,1)}{:},cond{contrasts(p,2)}{:});
        
        % 5.2.1 Test the exponents
        clear exp1 exp2
        for v = 1:length(cond{contrasts(p,1)})
            exp1(v) = cond{contrasts(p,1)}{v}.exponent;
        end
        for v = 1:length(cond{contrasts(p,2)})
            exp2(v) = cond{contrasts(p,2)}{v}.exponent;
        end
        
        [H,P,CI,exp_stats{p}] = ttest2(exp1,exp2,'vartype','unequal');
        exp_stats{p}.p = P;
        exp_stats{p}.CI = CI;
    end

%% 6. Plots

figure; 
    % Dunkel   
    loglog(GA{1}.freq,squeeze(mean(GA{1}.powspctrm,1)),'linewidth',3,'color',[1,0,0]);
    hold on
    
    % Hell
    loglog(GA{1}.freq,squeeze(mean(GA{2}.powspctrm,1)),'linewidth',3,'color',[0,0,1]);
    
    % Mltot
    loglog(GA{1}.freq,squeeze(mean(GA{4}.powspctrm,1)),'linewidth',3,'color',[0 0 0]);
    
    plot(GA{1}.freq,pow_stats{1}.mask,'r*');

xlim([GA{1}.freq(1) GA{1}.freq(end)]);


figure; 
    % Dunkel
    loglog(GA{1}.freq,squeeze(GA{1}.powspctrm),'linewidth',1,'color',[1,0.5,0.5]);
    hold on
    loglog(GA{1}.freq,squeeze(mean(GA{1}.powspctrm,1)),'linewidth',3,'color',[1,0,0]);
    
    % Hell
    loglog(GA{1}.freq,squeeze(GA{2}.powspctrm),'linewidth',1,'color',[0.5,0.5,1]);
    hold on
    loglog(GA{1}.freq,squeeze(mean(GA{2}.powspctrm,1)),'linewidth',3,'color',[0,0,1]);
    
    % Mllebend
    loglog(GA{1}.freq,squeeze(GA{3}.powspctrm),'linewidth',1,'color',[0.5,1,0.5]);
    hold on
    loglog(GA{1}.freq,squeeze(mean(GA{3}.powspctrm,1)),'linewidth',3,'color',[0,1,0]);
    
    % Mltot
    loglog(GA{1}.freq,squeeze(GA{4}.powspctrm),'linewidth',1,'color',[0.25,0.25,0.25]);
    hold on
    loglog(GA{1}.freq,squeeze(mean(GA{4}.powspctrm,1)),'linewidth',3,'color',[0 0 0]);

xlim([GA{1}.freq(1) GA{1}.freq(end)]);


plot(GA{1}.freq,squeeze(GA{1}.powspctrm),'linewidth',1,'color',[1,0.5,0.5]);
hold on
plot(GA{1}.freq,squeeze(mean(GA{1}.powspctrm,1)),'linewidth',3,'color',[1,0,0]);

figure;
loglog(light_GA.freq,mean(squeeze(light_GA.powspctrm)),'linewidth',3,'color',[1,0,0]);
hold on
for v = 1:size(light_GA.powspctrm,1)
    loglog(light_GA.freq,squeeze(light_GA.powspctrm(v,:,:)),'linewidth',1,'color',[1,0.5,0.5]);
end
loglog(light_GA.freq,mean(squeeze(light_GA.powspctrm)),'linewidth',3,'color',[1,0,0]);

for v = 1:size(dark_GA.powspctrm,1)
    loglog(dark_GA.freq,squeeze(dark_GA.powspctrm(v,:,:)),'linewidth',1,'color',[0.5,0.5,1]);
end

loglog(dark_GA.freq,mean(squeeze(dark_GA.powspctrm)),'linewidth',3,'color',[0,0,1]);
loglog(stats.freq,stats.mask,'k*');
title(strcat('red =  ', firstn, ', blue =  ', secondn, ', stars = sig. diff'));

% Test the exponents
[H,P,CI,STATS] = ttest2(light_exp',dark_exp','vartype','unequal')
