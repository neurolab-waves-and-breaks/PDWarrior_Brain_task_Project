function [] = load_EEG_data(Const, Paths)


DIODE = 38;
C3    = 15;
C4    = 17;
Oz    = 31;
Fz    = 6; % good result
F8    = 8;
F7    = 4;
Pz    = 26;
FPz    = 2;
F5    = 4;
Ri_APB  = 33;
Ri_FDI  = 34;


%% INFORMATION FOR P1
% Force task was only performed with the right hand (plus one rec. with the
% left at the very end)
% REST COMPARISON NOT GOOD (loDopa dopa PRE vs. hiDopa dopa POST)
% fname.P1.preExerc_loDopa.Rest = ['D>\PD_warrior_tests\data\P1\EEG\PD_Warrior_PreExerc_OFF_Rest_2023-03-23_11-49-15.cnt'];  % contains REST 
fname.P1.preExerc_loDopa.Force.file = ['D:\PD_warrior_tests\data\P1\EEG\PD_Warrior_PreExerc_OFF_FORCE_right_2023-03-23_11-52-48.cnt'];  % contains FORCE 
fname.P1.preExerc_loDopa.Force.time = [9, 138];  % contains FORCE 
fname.P1.preExerc_hiDopa.Force.file = ['D:\PD_warrior_tests\data\P1\EEG\PD_Warrior_PreExerc_ON_FORCE_2023-03-23_12-59-59.cnt'];  % contains FORCE 
fname.P1.preExerc_hiDopa.Force.time = [26, 155];  % contains FORCE 

% fname.P1.preExerc_loDopa.Rest = ['D>\PD_warrior_tests\data\P1\EEG\PD_Warrior_PostExerc_ON_REVSTOP_REST_AT_07min40_2023-03-23_14-55-54.cnt'];  % contains REST AND REVSTOP

fname.P1.postExerc_loDopa.Force.file = ['D:\PD_warrior_tests\data\P1\EEG\PD_Warrior_PostExerc_OFF_FORCE_2023-03-23_14-46-25.cnt'];  % contains FORCE
fname.P1.postExerc_loDopa.Force.time = [312, 438];  % contains FORCE 
fname.P1.postExerc_hiDopa.Force.file = ['D:\PD_warrior_tests\data\P1\EEG\PD_Warrior_PostExerc_ON_FORCE_2023-03-23_16-16-08.cnt'];  % contains FORCE
fname.P1.postExerc_hiDopa.Force.time = [62, 190];  % contains FORCE 



%% INFORMATION FOR P3
% Force task was performed first with the right hand, and then with the
% left hand, FOCUS on right hand task only, as this is where we had EMG
% UNFORTUNATELY RIGHT HAND performance in preExerc_loDopa seems very noisy
% in EEG

fname.P3.preExerc_loDopa.Rest.file = ['D:\PD_warrior_tests\data\P3\EEG\PDW3_PreExerc_LowDopa_2023-04-24_10-27-55.cnt'];  % contains REST and REVSTOP
fname.P3.preExerc_loDopa.Rest.time = [263, 382]; 
fname.P3.preExerc_loDopa.Force.file = ['D:\PD_warrior_tests\data\P3\EEG\PDW3_PreExerc_LowDopa_2023-04-24_10-44-26.cnt'];  % contains FORCE
fname.P3.preExerc_loDopa.Force.time = [8, 137]; 

fname.P3.preExerc_hiDopa.Rest.file = ['D:\PD_warrior_tests\data\P3\EEG\PDW3_PreExerc_HighDopa_2023-04-24_11-58-18.cnt'];  % contains REST ONLY
fname.P3.preExerc_hiDopa.Rest.time = [0.1, 120];  % contains REST AND REVSTOP
fname.P3.preExerc_hiDopa.Force.file = ['D:\PD_warrior_tests\data\P3\EEG\PDW3_PreExerc_HighDopa_2023-04-24_11-50-19.cnt'];  % contains FORCE
fname.P3.preExerc_hiDopa.Force.time = [8, 138]; 

fname.P3.postExerc_hiDopa.Rest.file = ['D:\PD_warrior_tests\data\P3\EEG\PDW3_PostExerc_HighDopa_2023-04-24_13-28-13.cnt'];  % contains REST ONLY
fname.P3.postExerc_hiDopa.Rest.time = [0.1, 122];  
fname.P3.postExerc_hiDopa.Force.file = ['D:\PD_warrior_tests\data\P3\EEG\PDW3_PostExerc_HighDopa_2023-04-24_13-19-30.cnt'];  % contains FORCE?
fname.P3.postExerc_hiDopa.Force.time = [20, 138]; 

fname.P3.postExerc_loDopa.Rest.file = ['D:\PD_warrior_tests\data\P3\EEG\PDW3_PostExerc_LowDopa_2023-04-24_14-38-13.cnt'];  % contains REST AND REVSTOP
fname.P3.postExerc_loDopa.Rest.time = [403, 523];  
fname.P3.postExerc_loDopa.Force.file = ['D:\PD_warrior_tests\data\P3\EEG\PDW3_PostExerc_LowDopa_2023-04-24_14-28-30.cnt'];  % contains FORCE?
fname.P3.postExerc_loDopa.Force.time = [26, 153];  


curr_fname = fname.P3.postExerc_hiDopa.Force.file;
info = eepv4_read_info(curr_fname);
% info.triggers
data = eepv4_read(curr_fname, 1, info.sample_count);
SR   = info.sample_rate;


% % PLOT LIGHT AND EEG channel so that you coan identify what .time to
% % specify
% figure; 
% time = (1:size(data.samples,2)) / SR;
% subplot(2,1,1)
% plot(time, data.samples(DIODE, :))
% subplot(2,1,2)
% plot(time, data.samples(C3, :))





%% COMPUTE COHERENCE AND POWER

participants = {'P1', 'P3'};
participants = {'P3'};
conds = {'preExerc_loDopa', 'postExerc_loDopa', 'preExerc_hiDopa', 'postExerc_hiDopa'};
conds = {'preExerc_loDopa', 'preExerc_hiDopa', 'postExerc_loDopa'};

Const.allCols = [Const.cols.lightorange; Const.cols.lightblue; Const.cols.darkorange ; ];

task = 'Force';

POW_SPECTRUM_FREQS = 5:2:80;

chans_to_plot = {'C3'}; %, 'C4', 'Oz', 'Ri_APB', 'Ri_FDI'};
% chans_to_plot_EMG = {'C4'}; %, 'Ri_FDI'};
chans_to_plot_EMG = {'Ri_APB'}; %, 'Ri_FDI'};
% chans_to_plot_EMG = {'C4'}; %, 'Ri_FDI'};


for p = 1:numel(participants)
    subj = participants{p};
    for c = 1:numel(conds)
        
        for chan = 1:numel(chans_to_plot)
            
            currCond = conds{c};
            currChan = eval(chans_to_plot{chan});
            curr_fname =  fname.(subj).(currCond).(task).file;
            info = eepv4_read_info(curr_fname);
            SR   = info.sample_rate;
            data = eepv4_read(curr_fname, 1, info.sample_count);

            cutTime   = fname.(subj).(currCond).(task).time;
            time_idcs = (cutTime(1)*SR) : (cutTime(2)*SR);

            EEG_pre  = data.samples(currChan, time_idcs);
            EEG_mean = movmean(EEG_pre, SR*2);       
            EEG      = EEG_pre - EEG_mean;         

            
            EMG_pre = data.samples(eval(chans_to_plot_EMG{chan}), time_idcs);
            EMG_mean = movmean(EMG_pre, SR*2);       
            EMG      = EMG_pre - EMG_mean;         
            FLT_ORDER = 2;
            EMG = abs(butterworth_filter2(EMG',  [20, 200], FLT_ORDER, 'twopass', SR));  

            FREQS = 3:2:45;

            for fi = 1:numel(FREQS)               
                sig1 = butterworth_filter2(EEG',  [FREQS(fi)-1, FREQS(fi)+1], FLT_ORDER, 'twopass', SR);
                sig1  = hilbert(sig1);

%                 sig2 = ft_preproc_highpassfilter(EMG, SR, FREQS(fi)-1, 4, 'but','twopass'); % twopass
%                 sig2 = ft_preproc_lowpassfilter(sig2,  SR, FREQS(fi)+1, 4, 'but','twopass'); % twopass
                sig2 = butterworth_filter2(EMG,  [FREQS(fi)-1, FREQS(fi)+1], FLT_ORDER, 'twopass', SR);
                sig2 = hilbert(sig2);

                % Match the time
                if numel(sig1) > 120*SR
                    sig1 = sig1(1:120*SR)';
                    sig2 = sig2(1:120*SR)';
                else
                    sig1 = sig1';
                    sig2 = sig2';
                end

                spec1 = mean(abs(sig1).^2,2);
                spec2 = mean(abs(sig2).^2,2);
                specX = abs(mean( abs(sig1).*abs(sig2) .* exp(1i*(angle(sig1)-angle(sig2))) ,2)).^2;
                % compute spectral coherence
                spectcoher(fi) = specX./(spec1.*spec2);

                phaseDiffs = angle(sig1) - angle(sig2);
                PLV(fi) = abs(nanmean(exp(1i*phaseDiffs), 2));
            end
    
            coh.(subj).(currCond).(task) = spectcoher;
            allPLV.(subj).(currCond).(task) = PLV;
        end
    end   
end



%% PLOTTING 

for p = 1:numel(participants)
    subj = participants{p};
    fig_h = figure;
    for chan = 1:numel(chans_to_plot)
        for coh_PLV = 1:2
            subplot(2,3,coh_PLV)
            
            for c = 1:numel(conds)
                
                currCond = conds{c};
                if coh_PLV == 1
                    powSpect = coh.(subj).(currCond).(task);
                    ylabeling = 'Coherence';
                elseif coh_PLV == 2
                     powSpect = allPLV.(subj).(currCond).(task);
                    ylabeling = 'Phase coupling';                   
                end
                % Normalization step
%                 powSpect = (powSpect - mean(powSpect(3:end))) / mean(powSpect(4:end));
                if c == 1 || c == 2
                    h(c) = plot(smooth(powSpect, 5), '-.', 'LineWidth', 3, 'Color', Const.allCols(c,:)); hold on
                else
                    h(c) = plot(smooth(powSpect, 5), '-', 'LineWidth', 3, 'Color', Const.allCols(c,:)); hold on
                end
                xTicks = 1:5:numel(POW_SPECTRUM_FREQS);
                set(gca, 'XTick', xTicks);
                freqCenter = POW_SPECTRUM_FREQS(xTicks);
                xlim([1, numel(powSpect)])
                set(gca, 'XTickLabel', freqCenter);            
                yLims = get(gca, 'YLim');
        %         plot(fltFrq_idx*[1,1], yLims, '--', 'Color', 'k', 'LineWidth', 2)
                %set(gca, 'YLim', yLims)
                title(strrep([subj, ' ', chans_to_plot{chan}, ' ', chans_to_plot_EMG{chan}, task], '_', ' '))
                xlabel('Frequencies [Hz]')
                ylabel(ylabeling)
                ytickformat('%.1f')
                xlim([2,45/2])
                set(gca, 'box', 'off')
                set(gca, 'FontSize', 16)
            end
            legend(h, strrep(conds, '_', ' '))
        %     saveas(fig_h, [plots_path, 'perSubject_EMG\', participant_ID, ' ', chans_to_plot{chan}, ' ', chans_to_plot_EMG{chan}, task'_EMG_powSpectra.png'])     
        end
    end
end