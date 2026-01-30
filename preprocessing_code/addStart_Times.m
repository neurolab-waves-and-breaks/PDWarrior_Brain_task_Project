function [] = addStart_Times(condition, Data, Const, Paths)

FullFileName = [Paths.matfiles, '/Evts_', condition, '.mat'];
SR = Data.SR;

% try
    load(FullFileName, 'Evt')


    smoothingSpan = Const.smoothingSpan;
    targetScale   = Data.targetScale;
    force_smooth  = smooth(Data.cursorScale, smoothingSpan);

    difftargetScale  = [0; diff(targetScale)'];   
    cutEnd_And_Start = Const.cutEnd_And_StartSec*SR;
    difftargetScale(1:cutEnd_And_Start) = nan;
    difftargetScale(end-cutEnd_And_Start:end) = nan;

    diffFrc = [0; diff(force_smooth)];                  
    diffFrc(1:cutEnd_And_Start) = nan;
    diffFrc(end-cutEnd_And_Start:end) = nan;

    normCueChange = abs(difftargetScale)>0.2;
    tagCueChange  = diff(normCueChange) > 0;

    idcsCueChange      = find(tagCueChange);
    % some events were double-detected, focus only on the onset
    doubleEvts = [false; diff(idcsCueChange)/SR < 0.2];
    idcsCueChange = idcsCueChange(~doubleEvts);   

    Evt.targetChange_time = idcsCueChange' / SR;  
    if numel(Evt.targetChange_time) > numel(Evt.startChange_time)
       Evt.targetChange_time(end) = []; 
    end
    Evt.RT_mainMvmt       = Evt.startChange_time - Evt.targetChange_time;

    Evt.oldTargetSize = targetScale(idcsCueChange - (0.1*SR));
    Evt.newTargetSize = targetScale(idcsCueChange + (0.1*SR));

    startMvmt = round(Evt.startOvershoot_time *SR);
    endMvmt   = round(Evt.endOvershoot_time * SR);
    nans      = isnan(startMvmt);

    Evt.main_mvmtDur   = Evt.endChange_time - Evt.startChange_time;
    Evt.total_mvmtDur  = Evt.endOvershoot_time - Evt.startOvershoot_time;

    start_mainMvmt = round(Evt.startChange_time *SR);
    end_mainMvmt   = round(Evt.endChange_time * SR);
    Evt.speed_mainMvmt(~nans) = abs(force_smooth(end_mainMvmt(~nans))' - force_smooth(start_mainMvmt(~nans))') ./ Evt.main_mvmtDur(~nans);
    forceChange = abs(force_smooth(end_mainMvmt(~nans))' - force_smooth(start_mainMvmt(~nans))');
    Evt.speed_mainMvmt(nans)  = nan;

    for i = 1:numel(Evt.rts)
        if isnan(Evt.rts(i))
            Evt.peakVel(i)     = nan;
            Evt.prctPlateau(i) =  nan;
            Evt.startError(i)  = nan;
            Evt.endError(i)    = nan;
        else
            currForce = force_smooth(start_mainMvmt(i):end_mainMvmt(i));
            if numel(currForce)== 1
                warning('Wrongly tagged in %s %i', FullFileName, i)
                Evt.peakVel(i)     = nan;
                Evt.prctPlateau(i) =  nan;
                Evt.startError(i)  = nan;
                Evt.endError(i)    = nan;
                Evt.rts(i)         = nan;
                continue
            end
            Evt.peakVel(i) = max(abs(diff(currForce)));
            noMvmtThresh = 0.06; 
            tag_plateau = abs(diff(currForce)) < noMvmtThresh;
            
%             figure; plot(currForce); hold on; plot(find(tag_plateau), currForce(tag_plateau), 'x', 'Color', 'r')
            Evt.prctPlateau(i) = sum(abs(diff(currForce)) < noMvmtThresh) / numel(currForce);
            
            % Take the mean of the error in a window ranging 0.1s before
            % and 0.1s after movement start and end respectivelz
            Evt.startError(i) = mean(abs(force_smooth((startMvmt(i)-0.1*SR) : startMvmt(i))' - Evt.oldTargetSize(i)));
            if (endMvmt(i)+0.1*SR) <= numel(force_smooth)
                Evt.endError(i)   = mean(abs(force_smooth(endMvmt(i) : (endMvmt(i)+0.1*SR))'- Evt.newTargetSize(i)));
            else
                Evt.endError(i)   = mean(abs(force_smooth(endMvmt(i) : end)'- Evt.newTargetSize(i)));                
            end           
        end
    end
    
    
%     fig_h = figure;
%     supertitle(strrep(condition, '_', ' '))
%     subplot(2,2,1); plot(forceChange, Evt.main_mvmtDur(~nans), 'o'); xlabel('Force change'); ylabel('Mvmt duration')
%     subplot(2,2,2); plot(forceChange,Evt.speed_mainMvmt(~nans), 'o'); xlabel('Force change'); ylabel('Speed')
%     subplot(2,2,3); plot(Evt.main_mvmtDur(~nans),Evt.speed_mainMvmt(~nans), 'o'); xlabel('Mvmt duration'); ylabel('Speed')
%     saveas(fig_h, [Paths.plots, 'FORCE_task_speed/Speed__', condition '.png'])


    save(FullFileName, 'Evt', '-v7.3')                    
% end