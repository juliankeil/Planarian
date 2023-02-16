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
inpath{5} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/06_KLDunkel/'; 
inpath{6} = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/07_KLHell/';

outpath = '/Users/juliankeil/Documents/Arbeit/Kiel/Abschlussarbeiten/Lang/GitHub/Planarian/02_data/02_Preproc/101_Output/';

for p = 1:length(inpath)
    folders = find((inpath{p} == '/'));
    names{p} = [inpath{p}(folders(end-1)+1:folders(end)-1),'.csv'];
end

%% 2. Loop Data to collect single recordings
cond = cell(1,6);
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
    cfg.keepindividual = 'yes';

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
                 3 4;...% 6 mllebend vs. mltot
                 5 6;... % 7 Kopflos dunkel vs. kopflos hell
                 1 5;... % 8 dunkel vs. kopflow dunkel
                 2 6;... % 9 hell vs. kopflos hell
                 ]; %

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
        exp_stats{p}.means = [mean(exp1),mean(exp2)];
        exp_stats{p}.stds = [std(exp1),std(exp2)];
    end

%% 6. Plots
% Figure 2
figure;

loglog(GA{1}.freq,mean(squeeze(GA{1}.powspctrm)),'linewidth',3,'color',[1,0,0]);
hold on
loglog(GA{2}.freq,mean(squeeze(GA{2}.powspctrm)),'linewidth',3,'color',[0,0,1]);
loglog(GA{4}.freq,mean(squeeze(GA{4}.powspctrm)),'linewidth',3,'color',[0,0,0]);


% Dark
for v = 1:size(GA{1}.powspctrm,1)
    loglog(GA{1}.freq,squeeze(GA{1}.powspctrm(v,:,:)),'linewidth',1,'color',[1,0.5,0.5]);
end

% Light
for v = 1:size(GA{2}.powspctrm,1)
    loglog(GA{2}.freq,squeeze(GA{2}.powspctrm(v,:,:)),'linewidth',1,'color',[0.5,0.5,1]);
end

% MLtot
for v = [2:5,7:18]%size(GA{4}.powspctrm,1)
    loglog(GA{4}.freq,squeeze(GA{4}.powspctrm(v,:,:)),'linewidth',1,'color',[0.5,0.5,0.5]);
end

loglog(GA{1}.freq,mean(squeeze(GA{1}.powspctrm)),'linewidth',3,'color',[1,0,0]);
loglog(GA{2}.freq,mean(squeeze(GA{2}.powspctrm)),'linewidth',3,'color',[0,0,1]);
loglog(GA{4}.freq,mean(squeeze(GA{4}.powspctrm)),'linewidth',3,'color',[0,0,0]);

plot(pow_stats{1}.freq,pow_stats{1}.mask.*10^-2.5,'g--','linewidth',3);

xlim([GA{1}.freq(1) GA{1}.freq(end)]);
legend('Dark','Light','Larvae');
xlabel('log Frequency (Hz)');
ylabel('log Amplitude');

% Figure 3
% Dark
figure;
subplot(1,2,1)
loglog(GA{1}.freq,mean(squeeze(GA{1}.powspctrm)),'linewidth',3,'color',[1,0,0]);
hold on
loglog(GA{5}.freq,mean(squeeze(GA{5}.powspctrm)),'linewidth',3,'color',[0,0,1]);


% Dark
for v = 1:size(GA{1}.powspctrm,1)
    loglog(GA{1}.freq,squeeze(GA{1}.powspctrm(v,:,:)),'linewidth',1,'color',[1,0.5,0.5]);
end

% Dark Decap
for v = 1:size(GA{5}.powspctrm,1)
    loglog(GA{5}.freq,squeeze(GA{5}.powspctrm(v,:,:)),'linewidth',1,'color',[0.5,0.5,1]);
end

plot(pow_stats{8}.freq,pow_stats{8}.mask.*min(GA{1}.powspctrm(:)),'g--','linewidth',3);

loglog(GA{1}.freq,mean(squeeze(GA{1}.powspctrm)),'linewidth',3,'color',[1,0,0]);
loglog(GA{5}.freq,mean(squeeze(GA{5}.powspctrm)),'linewidth',3,'color',[0,0,1]);
xlim([GA{1}.freq(1) GA{1}.freq(end)]);
legend('Dark intact','Dark decapitated');
xlabel('log Frequency (Hz)');
ylabel('log Amplitude');

% Light
subplot(1,2,2)
loglog(GA{2}.freq,mean(squeeze(GA{2}.powspctrm)),'linewidth',3,'color',[1,0,0]);
hold on
loglog(GA{6}.freq,mean(squeeze(GA{6}.powspctrm)),'linewidth',3,'color',[0,0,1]);

% Dark
for v = 1:size(GA{2}.powspctrm,1)
    loglog(GA{2}.freq,squeeze(GA{2}.powspctrm(v,:,:)),'linewidth',1,'color',[1,0.5,0.5]);
end

% Dark Decap
for v = 1:size(GA{6}.powspctrm,1)
    loglog(GA{6}.freq,squeeze(GA{6}.powspctrm(v,:,:)),'linewidth',1,'color',[0.5,0.5,1]);
end

plot(pow_stats{9}.freq,pow_stats{9}.mask.*min(GA{1}.powspctrm(:)),'g--','linewidth',3);

loglog(GA{2}.freq,mean(squeeze(GA{2}.powspctrm)),'linewidth',3,'color',[1,0,0]);
loglog(GA{6}.freq,mean(squeeze(GA{6}.powspctrm)),'linewidth',3,'color',[0,0,1]);
xlim([GA{2}.freq(1) GA{2}.freq(end)]);
legend('Light intact','Light decapitated');
xlabel('log Frequency (Hz)');
ylabel('log Amplitude');


