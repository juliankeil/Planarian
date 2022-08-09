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

%% 1. Define Data in Folder
inpath{1} = '/Plattwurm/preproc/mltot/'; % Erster Datensatz
inpath{2} = '/Plattwurm/preproc/leer/'; % Zweiter Datensatz
outpath = '/Plattwurm/Output/';

nposfirst = find((inpath{1} == '/')==1)
firstn = inpath{1}((nposfirst(length(nposfirst)-1)+1):(nposfirst(length(nposfirst))-1))
npossecond = find((inpath{2} == '/')==1)
secondn = inpath{2}((npossecond(length(npossecond)-1)+1):(npossecond(length(npossecond))-1))

tblnamefirst = strcat(firstn, '.csv');
tblnamesecond = strcat(secondn, '.csv');

%% 2. Loop Data

% 2.1. Light
p = 1; % Light
clear indat
indat = dir([inpath{p},'*.mat']);
% 2.1.1. Read in
light_doc = {1;1};
for v = 1:length(indat)
    tmp = load([inpath{p},indat(v).name],'data_fft');
    light_fft{1} = 1;
    light_exp(1) = 1;
    if length(tmp.data_fft.label) == 1  % Wenn nur ein Kanal: Alles in Ordnung
        light_fft{length(light_fft)+1} = tmp.data_fft;
        tmp2 = load([inpath{p},indat(v).name],'data_expo');
        light_exp(length(light_exp)+1) = tmp2.data_expo;
        light_doc(1,length(light_fft)) = {indat(v).name};
        light_doc(2,length(light_fft)) = tmp.data_fft.label;
    else % Wenn mehr als ein Kanal: Jeder Kanal wird als eigener Wurm behandelt
        for j = 1:length(tmp.data_fft.label)
        tmp2 = tmp;
        tmp2.data_fft.label = tmp.data_fft.label(j);
        tmp2.data_fft.powspctrm = tmp.data_fft.powspctrm(j,1:length(tmp.data_fft.powspctrm));
        light_fft{length(light_fft)+1} = tmp2.data_fft;
        tmp3 = load([inpath{p},indat(v).name],'data_expo');
        light_exp(length(light_exp)+1) = tmp3.data_expo(j);
        light_doc(1,length(light_fft)) = {indat(v).name};
        light_doc{2,length(light_fft)} = tmp.data_fft.label(j);
        end
    end
end

% Clean up
light_fft = light_fft(2:length(light_fft));
light_exp = light_exp(2:length(light_exp));
light_doc = light_doc(1:2,2:length(light_doc));


% Renaming for Group average
channelname = {'Fp1'};
for i = 1:length(light_fft)
    light_fft{i}.label = channelname;
end

% 2.1.2. Group Average
cfg = [];
cfg.keepindividual = 'yes';

light_GA = ft_freqgrandaverage(cfg,light_fft{:});

% 2.1.3. Write CSV
    % 2.1.3.1. Make Table
    clear tbl
    tbl = array2table(squeeze(light_GA.powspctrm));

    % 2.1.3.2. Prepare Column Names = Frequencies
    colNames = {};
    for c = 1:length(light_GA.freq)
        tmp = num2str(light_GA.freq(c));
        tmp = replace(tmp,'.','_');
        colNames{c} = tmp;
    end

    % 2.1.3.3. Add Column Names
    for v = 1:size(tbl,2)
        tbl.Properties.VariableNames{v} = ['f_',colNames{v}];
    end

    % 2.1.3.4. Add Row Names = Dataset Names
    for v = 1:size(tbl,1)
        tmpname = char(light_doc(1,v));
        tmpname = tmpname(find(tmpname~='_'));
        tmpname2 = char(light_doc{2,v}(1));
        tmpname_final = strcat(tmpname(1:(length(tmpname)-8)),tmpname2);
        tbl.Properties.RowNames{v} = tmpname_final;
    end

    % 2.1.3.5 Now save it 
    cd(outpath);
    writetable(tbl,tblnamefirst,'WriteVariableNames',1,'WriteRowNames',1);
    
% 2.2. Dark
p = 2; % Dark
clear indat
indat = dir([inpath{p},'*.mat']);
% 2.2.1.
dark_doc = {1;1};
for v = 1:length(indat)
    tmp = load([inpath{p},indat(v).name],'data_fft');
    dark_fft{1} = 1;
    dark_exp(1) = 1;
    if length(tmp.data_fft.label) == 1  % Wenn nur ein Kanal: Alles in Ordnung
        dark_fft{length(dark_fft)+1} = tmp.data_fft;
        tmp2 = load([inpath{p},indat(v).name],'data_expo');
        dark_exp(length(dark_exp)+1) = tmp2.data_expo;
        dark_doc(1,length(dark_fft)) = {indat(v).name};
        dark_doc(2,length(dark_fft)) = tmp.data_fft.label;
    else % Wenn mehr als ein Kanal: Jeder Kanal wird als eigener Wurm behandelt
        for j = 1:length(tmp.data_fft.label)
        tmp2 = tmp;
        tmp2.data_fft.label = tmp.data_fft.label(j);
        tmp2.data_fft.powspctrm = tmp.data_fft.powspctrm(j,1:length(tmp.data_fft.powspctrm));
        dark_fft{length(dark_fft)+1} = tmp2.data_fft;
        tmp3 = load([inpath{p},indat(v).name],'data_expo');
        dark_exp(length(dark_exp)+1) = tmp3.data_expo(j);
        dark_doc(1,length(dark_fft)) = {indat(v).name};
        dark_doc{2,length(dark_fft)} = tmp.data_fft.label(j);
        end
    end
end

% Clean up
dark_fft = dark_fft(2:length(dark_fft));
dark_exp = dark_exp(2:length(dark_exp));
dark_doc = dark_doc(1:2,2:length(dark_doc));


% 2.2.2.1 Renaming for Group Average
for i = 1:length(dark_fft)
    dark_fft{i}.label = channelname;
end


% 2.2.2.2 Actual Group Average
cfg = [];
cfg.keepindividual = 'yes';
dark_GA = ft_freqgrandaverage(cfg,dark_fft{:});

% 2.1.3. Write CSV
    % 2.1.3.1. Make Table
    clear tbl
    tbl = array2table(squeeze(dark_GA.powspctrm));

    % 2.1.3.2. Prepare Column Names = Frequencies
    colNames = {};
    for c = 1:length(dark_GA.freq)
        tmp = num2str(dark_GA.freq(c));
        tmp = replace(tmp,'.','_');
        colNames{c} = tmp;
    end

    % 2.1.3.3. Add Column Names
    for v = 1:size(tbl,2)
        tbl.Properties.VariableNames{v} = ['f_',colNames{v}];
    end

    % 2.1.3.4. Add Row Names = Dataset Names
    for v = 1:size(tbl,1)
        tmpname = char(dark_doc(1,v));
        tmpname = tmpname(find(tmpname~='_'));
        tmpname2 = char(dark_doc{2,v}(1));
        tmpname_final = strcat(tmpname(1:(length(tmpname)-8)),tmpname2);
        tbl.Properties.RowNames{v} = tmpname_final;
    end

    % 2.1.3.5 Now save it 
    cd(outpath);
    writetable(tbl,tblnamesecond,'WriteVariableNames',1,'WriteRowNames',1);

%% Stats

% Compare the spectra
cfg=[];
cfg.method = 'montecarlo';
cfg.statistic ='indepsamplesT';
cfg.numrandomization = 1000;
cfg.correctm = 'cluster';
cfg.neighbours = [];
cfg.alpha = 0.05;
cfg.tail = 0;
cfg.correcttail = 'alpha';
cfg.ivar = 1;

% % if cfg.statistic == 'depsamplesT'
% cfg.design=[1:length(light_fft),1:length(dark_fft);ones(1,length(light_fft)),ones(1,length(dark_fft))*2];
% cfg.uvar   = 1;
% cfg.ivar   = 2;
% % else
cfg.design=[ones(1,length(light_fft)),ones(1,length(dark_fft))*2];
% end

stats=ft_freqstatistics(cfg,light_fft{:},dark_fft{:})

% Make a group plot
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
