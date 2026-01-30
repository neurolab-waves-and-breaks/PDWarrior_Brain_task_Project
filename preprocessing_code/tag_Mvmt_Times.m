function [Evt, fig_h] = tag_Mvmt_Times(condition, Data, Const, Paths, plotOn)

FullFileName = [Paths.matfiles, '/Evts_', condition, '.mat'];
% SR = Data_RS.SR
SR = Data.SR;
evts_present = false;
overshoot_start_idcs = [];
rts = [];
try  % Load it in case I want to modify single events (in L 45)
    load(FullFileName, 'Evt')
    evts_present = true;
    incr_Cue = Evt.incr_Cue;
    rts      = Evt.rts;
    middleFrcChange_idcs = round(Evt.middleFrcChange_time * SR);
    startChange_idcs = round(Evt.startChange_time * SR);
    end_Change_idcs  = round(Evt.endChange_time * SR);
    overshoot_start_idcs = round(Evt.startOvershoot_time * SR);
    overshoot_end_idcs   = round(Evt.endOvershoot_time * SR);
    overshoot_start = Evt.startOvershoot;
    overshoot_end   = Evt.endOvershoot;
end

smoothingSpan = Const.smoothingSpan;
%targetScale  = Data.TargetScale;
targetScale  = Data.targetScale;
force_smooth  = smooth(Data.cursorScale, smoothingSpan);

%difftargetScale  = [0; diff(targetScale)]
difftargetScale  = [0; diff(targetScale)'];   
cutEnd_And_Start = Const.cutEnd_And_StartSec*SR;
difftargetScale(1:cutEnd_And_Start) = nan;
difftargetScale(end-cutEnd_And_Start:end) = nan;

diffFrc = [0; diff(force_smooth)];                  
diffFrc(1:cutEnd_And_Start) = nan;
diffFrc(end-cutEnd_And_Start:end) = nan;

normCueChange = abs(difftargetScale)>0.2;
tagCueChange  = diff(normCueChange) > 0;

% Get the direction of the change 
ONLY_PLOT_ALL_EVTS = false;
if ~ONLY_PLOT_ALL_EVTS
    idcsCueChange      = find(tagCueChange);
    % some events were double-detected, focus only on the onset
    doubleEvts = [false; diff(idcsCueChange)/SR < 0.2];
    idcsCueChange = idcsCueChange(~doubleEvts);
    
    incr_Cue           = (targetScale(idcsCueChange+5))-(targetScale(idcsCueChange-5))>0;

    % Find maximum force change within a certain window 
    % clear maxFrcChange_idcs startChange_idcs end_Change_idcs rts overshoot_start_idcs overshoot_end_idcs overshoot_start overshoot_end;
    minCueDist = Const.minCueDist; %1.25; 
    toDo = 1:numel(idcsCueChange);
    
    for t = toDo
        currWin = (idcsCueChange(t)-SR*1.5): idcsCueChange(t)+3*SR;
        if currWin(end) > numel(diffFrc)
            currWin = (idcsCueChange(t)-SR*1.5): numel(diffFrc);
        end
        currFrcChange = diffFrc(currWin); % until 1.5 sec time
        if incr_Cue(t)
            [maxi, idcsMax] = max(currFrcChange);
        else  % if it was a decreasing cue
            [maxi, idcsMax] = min(currFrcChange);
        end
        maxFrcChange_idcs(t) = idcsCueChange(t) + idcsMax;
        
        if plotOn  % show the limited window for tagging events
            close all
            fig_h = figure; plot(targetScale(currWin)); hold on   
            plot(force_smooth(currWin)); 
            
            if numel(rts) >= t
                if evts_present && ~isnan(rts(t))  
                    plot(overshoot_start_idcs(t) - currWin(1), force_smooth(overshoot_start_idcs(t)), 'x', 'Color', 'k', 'LineWidth', 2);
                    plot(overshoot_end_idcs(t)  - currWin(1), force_smooth(overshoot_end_idcs(t)), 'x', 'Color', 'k', 'LineWidth', 2);
                    plot(startChange_idcs(t)  - currWin(1), force_smooth(startChange_idcs(t)), 'o', 'Color', 'r');
                    plot(end_Change_idcs(t)  - currWin(1), force_smooth(end_Change_idcs(t)), 'o', 'Color', 'g');
                    plot(middleFrcChange_idcs(t)  - currWin(1), force_smooth(middleFrcChange_idcs(t)), 'o', 'Color', 'b');                
                end
            end
            xlim([1,numel(currWin)])
            ylims = get(gca, 'YLim');
            xlims = get(gca, 'Xlim');
            %for inputs = 1:4
            [x,y] = ginput(1); % collect 4 points in a row
            if y(1) > ylims(2) % skip if you clicked above the top border of the plot
                continue
            end
            if x(1) < xlims(1) % tag as invalid if you clicked to the left of the plot
                fprintf('Trial discarded.\n');
                rts(t) = nan;
                middleFrcChange_idcs(t) = nan;
                overshoot_start_idcs(t) = nan;
                startChange_idcs(t) = nan;
                end_Change_idcs(t) = nan;
                overshoot_end_idcs(t)  = nan;
                overshoot_start(t) = nan;
                overshoot_end(t) = nan;
                continue
            end
            [x2,~] = ginput(3); % collect 4 points in a row
            x(2:4) = x2;
            x = x + currWin(1);
            % 4 inputs in case there is over and undershoot:
            % Always click on the beginning of the first movement, if there is
            % no clear undershoot click on the same spot twice
            % Then click on the end of the movement (peak of the curve=overshoot), and 
            % where the movement actually stops
            if (x(2)-x(1) < - 50) | (x(4)-x(3) < -50)  % click to the left 
                fprintf('Trial discarded.\n');
                rts(t) = nan;
                middleFrcChange_idcs(t) = nan;
                overshoot_start_idcs(t) = nan;
                startChange_idcs(t) = nan;
                end_Change_idcs(t) = nan;
                overshoot_end_idcs(t)  = nan;
                overshoot_start(t) = nan;
                overshoot_end(t) = nan;
                continue
            else
                fprintf('Event %i\n', t)
                overshoot_start_idcs(t) = round(x(1)) ;
                startChange_idcs(t) = round(x(2));
                end_Change_idcs(t)  = round(x(3));
                overshoot_end_idcs(t) = round(x(4));
                middleFrcChange_idcs(t) = round(startChange_idcs(t)+(end_Change_idcs(t)-startChange_idcs(t))/2);
                overshoot_start(t) = force_smooth(startChange_idcs(t)) - ...
                                     force_smooth(overshoot_start_idcs(t));
                overshoot_end(t) = force_smooth(end_Change_idcs(t)) - ...
                                   force_smooth(overshoot_end_idcs(t));
                fprintf('Overshoot start =%.1f, end=%.1f \n', overshoot_start(t), overshoot_end(t))
%                 close all
%                 fig_h = figure; plot(targetScale); hold on   
%                 plot(force_smooth); 
                plot(overshoot_start_idcs(t) - currWin(1), force_smooth(overshoot_start_idcs(t)), 'x', 'Color', 'k');
                plot(overshoot_end_idcs(t)  - currWin(1), force_smooth(overshoot_end_idcs(t)), 'x', 'Color', 'k');
                plot(startChange_idcs(t)  - currWin(1), force_smooth(startChange_idcs(t)), 'o', 'Color', 'r');
                plot(end_Change_idcs(t)  - currWin(1), force_smooth(end_Change_idcs(t)), 'o', 'Color', 'g');
                plot(middleFrcChange_idcs(t)  - currWin(1), force_smooth(middleFrcChange_idcs(t)), 'o', 'Color', 'b');
                rts(t) = (overshoot_start_idcs(t) - idcsCueChange(t)) /SR;                
            end
        end        
        nans = isnan(rts) | (rts <= 0);
        rts(nans) = nan;
        startChange_idcs(nans) = nan;
        end_Change_idcs(nans) = nan;
        middleFrcChange_idcs(nans) = nan;
        overshoot_end_idcs(nans)   = nan;
        overshoot_start_idcs(nans) = nan;
        overshoot_end(nans)   = nan;
        overshoot_start(nans) = nan;

        Evt.targetChange_time    = idcsCueChange' / SR;  
        Evt.incr_Cue             = incr_Cue;
        Evt.rts                  = rts ;
        Evt.middleFrcChange_time = middleFrcChange_idcs /SR;
        Evt.startChange_time     = startChange_idcs /SR;
        Evt.endChange_time       = end_Change_idcs /SR;
        Evt.startOvershoot_time  = overshoot_start_idcs /SR;
        Evt.endOvershoot_time    = overshoot_end_idcs /SR;
        Evt.startOvershoot       = overshoot_start;
        Evt.endOvershoot         = overshoot_end;
        
        save(FullFileName, 'Evt', '-v7.3')                    
    end
end

if plotOn
    fig_h = figure; 
    %load([Paths.SaveMatFiles, 'conversionFactor/conversionFactor.mat'])
    conversionFactor = 1;
    
    subplot(2,1,1)
    plot(targetScale / conversionFactor, 'LineWidth', 1.3, 'Color', 0.7*[1,1,1]); hold on   
    increments = middleFrcChange_idcs(incr_Cue);
    decrements = middleFrcChange_idcs(~incr_Cue);
    decrements(isnan(decrements)) = [];
    increments(isnan(increments)) = [];
    plot(force_smooth / conversionFactor, 'LineWidth', 1.2);   
    startChange = startChange_idcs;
    startChange(isnan(startChange)) = [];
    endChange = end_Change_idcs;
    endChange(isnan(endChange)) = [];
    plot(startChange, force_smooth(startChange) / conversionFactor, 'x', 'Color', 'k', 'LineWidth', 1.3, 'MarkerSize', 4);
    plot(endChange, force_smooth(endChange) / conversionFactor, 'x', 'Color', 'k', 'LineWidth', 1.3, 'MarkerSize', 4);
    title(strrep(condition, '_', ' '))
    xTicks = get(gca, 'XTick');
    set(gca, 'XTickLabel', xTicks/1000)
%     ylim([0, 3.2])
%     set(gca, 'FontSize', 14) % Laptop
    set(gca, 'FontSize', 16) % Work PC
    ylabel('Force [N]')
    xlabel('Time [s]')
    
    if not(isfolder([Paths.plots, 'force_traces/']))
        mkdir([Paths.plots, 'force_traces/'])
    end
    saveas(fig_h, [Paths.plots, 'force_traces/Frc_trace_', condition, '.jpg'])       
    close all
end


