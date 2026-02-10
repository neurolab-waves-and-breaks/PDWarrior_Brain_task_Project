function [] = behav_analysis_BRAIN_Test_reordered(Const, Paths)
%% BEHAV_ANALYSIS_BRAIN_Test 

[Const] = storeConst();
participant_dirs = dir(Paths.Data);
participant_dirs = participant_dirs([participant_dirs.isdir]);
participant_dirs = participant_dirs(~ismember({participant_dirs.name}, {'.', '..'}));

participants = {participant_dirs.name};
fprintf('Found %d participants: %s\n', length(participants), strjoin(participants, ', '));

exerc = {'preE', 'postE'};
hand  = {'right', 'left'};

allTravel = [];
allDwell = [];
allV = [];
handrep = [];
exerep = [];
IDrep = [];
recrep = [];
IDrep1 = [];
exerep1 = [];
handrep1 = [];
recrep1 = [];
HANDP = [];
EXERCP = [];
IDP = [];
RECORP = [];
IDrep2 = [];
exerep2 = [];
handrep2 = [];
recrep2 = [];
metrics = [];

%% Extract participant characteristics
for p = 1:numel(participants)
    subj = participants{p};
    part_folder_info = dir([Paths.Data, subj]);
    all_fnames = {part_folder_info.name};
    
    for h = 1:numel(hand)
        for e = 1:numel(exerc)
            for rec = 1:2
                fileIdx = ~cellfun(@isempty, strfind(all_fnames, [exerc{e}, 'xerc_', hand{h}, '_v', num2str(rec)]));            
                if sum(fileIdx) == 0
                   continue
                end
            
                data = readtable(fullfile(Paths.Data, subj, all_fnames{fileIdx}));

                % Dysmetria score
                correctLetters = strcmp(data.letters, 's') + strcmp(data.letters, 'semicolon');
                adjacentKeys = (strcmp(data.letters, 'a') + strcmp(data.letters, 'w') + strcmp(data.letters, 'e') + strcmp(data.letters, 'd') + strcmp(data.letters, 'x') + strcmp(data.letters, 'z') + ...
                                strcmp(data.letters, 'l') + strcmp(data.letters, 'p') + strcmp(data.letters, 'leftbracket') + strcmp(data.letters, 'quote') + ...
                                strcmp(data.letters, 'slash') + strcmp(data.letters, 'period')) * 2;
                otherKeys = ~(correctLetters + adjacentKeys) * 3;      

                dysmetScores = correctLetters + adjacentKeys + otherKeys;
                    
                strokeItvl = diff(data.strokeOnset);
                strokeItvl = [0; strokeItvl];  

                dwellTime = data.dwellTime(2:end);
                dwellTime(dwellTime < 0.01) = [];
                [dwellTime_outliers] = medAbsDeviation(dwellTime);
                dwellTime = dwellTime(~dwellTime_outliers);
                                        
                %% Recreate keypairs to calculate distance    
                travelTime = strokeItvl(2:end) - data.dwellTime(1:end-1);
                
                keypairs = [];
                for key = 1:length(data.letters)-1 
                    keypairs{key} = strcat(data.letters{key}, data.letters{key+1});
                end
                keypairs = reshape(keypairs, [], 1);

                %% Define distance between pressed keys
                velocitydata = [];
                valid_keypair_mask = true(length(keypairs), 1);

                for pair = 1:length(keypairs)
                    if travelTime(pair) > 0.1 && travelTime(pair) < 1.5
                        if strcmp(keypairs{pair}, 'ssemicolon') || strcmp(keypairs{pair}, 'semicolons') || ...
                                strcmp(keypairs{pair}, 'xslash') || strcmp(keypairs{pair}, 'slashx') || ...
                                strcmp(keypairs{pair}, 'eleftbracket') || strcmp(keypairs{pair}, 'leftbrackete') || ...
                                strcmp(keypairs{pair}, 'zperiod') || strcmp(keypairs{pair}, 'periodz') || ...
                                strcmp(keypairs{pair}, 'al') || strcmp(keypairs{pair}, 'la') || strcmp(keypairs{pair}, 'wp') || ...
                                strcmp(keypairs{pair}, 'dquote') || strcmp(keypairs{pair}, 'quoted') || strcmp(keypairs{pair}, 'pw')
                            velocitydata(end+1, :) = [15, travelTime(pair)];

                        elseif strcmp(keypairs{pair}, 'ld') || strcmp(keypairs{pair}, 'dl')
                            velocitydata(end+1, :) = [11.5, travelTime(pair)];

                        elseif strcmp(keypairs{pair}, 'lx') || strcmp(keypairs{pair}, 'xl')
                            velocitydata(end+1, :) = [12.5, travelTime(pair)];

                        elseif strcmp(keypairs{pair}, 'sl') || strcmp(keypairs{pair}, 'ls') || ...
                                strcmp(keypairs{pair}, 'dsemicolon') || strcmp(keypairs{pair}, 'semicolond') || ...
                                strcmp(keypairs{pair}, 'ak') || strcmp(keypairs{pair}, 'ka') || ...
                                strcmp(keypairs{pair}, 'eo') || strcmp(keypairs{pair}, 'oe') || ...
                                strcmp(keypairs{pair}, 'fquote') || strcmp(keypairs{pair}, 'quotef') || ...
                                strcmp(keypairs{pair}, 'xperiod') || strcmp(keypairs{pair}, 'periodx')
                            velocitydata(end+1, :) = [13.5, travelTime(pair)];

                        elseif strcmp(keypairs{pair}, 'sp') || strcmp(keypairs{pair}, 'ps') || ...
                                strcmp(keypairs{pair}, 'zl') || strcmp(keypairs{pair}, 'lz') || ...
                                strcmp(keypairs{pair}, 'esemicolon') || strcmp(keypairs{pair}, 'semicolone') || ...
                                strcmp(keypairs{pair}, 'xsemicolon') || strcmp(keypairs{pair}, 'semicolonx') || ...
                                strcmp(keypairs{pair}, 'speriod') || strcmp(keypairs{pair}, 'periods') || ...
                                strcmp(keypairs{pair}, 'slashd') || strcmp(keypairs{pair}, 'dslash')
                            velocitydata(end+1, :) = [14.5, travelTime(pair)];

                        elseif strcmp(keypairs{pair}, 'xquote') || strcmp(keypairs{pair}, 'quotex') || ...
                                strcmp(keypairs{pair}, 'sslash') || strcmp(keypairs{pair}, 'slashs') || ...
                                strcmp(keypairs{pair}, 'wsemicolon') || strcmp(keypairs{pair}, 'semicolonw') || ...
                                strcmp(keypairs{pair}, 'zsemicolon') || strcmp(keypairs{pair}, 'semicolonz')
                            velocitydata(end+1, :) = [16, travelTime(pair)];

                        elseif strcmp(keypairs{pair}, 'sleftbracket') || strcmp(keypairs{pair}, 'leftbrackets') || ...
                                strcmp(keypairs{pair}, 'ap') || strcmp(keypairs{pair}, 'pa') || ...
                                strcmp(keypairs{pair}, 'speriod') || strcmp(keypairs{pair}, 'pa') || ...
                                strcmp(keypairs{pair}, 'qleftbracket') || strcmp(keypairs{pair}, 'leftbracketq')
                            velocitydata(end+1, :) = [16.5, travelTime(pair)];

                        elseif strcmp(keypairs{pair}, 'squote') || strcmp(keypairs{pair}, 'quotes') || ...
                                strcmp(keypairs{pair}, 'asemicolon') || strcmp(keypairs{pair}, 'semicolona') || ...
                                strcmp(keypairs{pair}, 'zslash') || strcmp(keypairs{pair}, 'slashz') || ...
                                strcmp(keypairs{pair}, 'qp') || strcmp(keypairs{pair}, 'pq') || ...
                                strcmp(keypairs{pair}, 'wleftbracket') || strcmp(keypairs{pair}, 'leftbracketw')
                            velocitydata(end+1, :) = [17, travelTime(pair)];

                        elseif strcmp(keypairs{pair}, 'slasha') || strcmp(keypairs{pair}, 'aslash') || ...
                                strcmp(keypairs{pair}, 'quotez') || strcmp(keypairs{pair}, 'zquote') || ...
                                strcmp(keypairs{pair}, 'zslash') || strcmp(keypairs{pair}, 'slashz') || ...
                                strcmp(keypairs{pair}, 'qp') || strcmp(keypairs{pair}, 'pq') || ...
                                strcmp(keypairs{pair}, 'wleftbracket') || strcmp(keypairs{pair}, 'leftbracketw')
                            velocitydata(end+1, :) = [18, travelTime(pair)];

                        elseif strcmp(keypairs{pair}, 'aquote') || strcmp(keypairs{pair}, 'quotea') || ...      
                                strcmp(keypairs{pair}, 'qleftbracket') || strcmp(keypairs{pair}, 'leftbracketq')
                            velocitydata(end+1, :) = [19, travelTime(pair)];                    

                        else
                            fprintf('Unrecognized keypair: %s\n', keypairs{pair});
                            valid_keypair_mask(pair) = false;
                        end
                    else
                        valid_keypair_mask(pair) = false;
                    end
                end

                dwellTime_filtered = dwellTime(valid_keypair_mask(1:min(length(valid_keypair_mask), length(dwellTime))));
                [dwellTime_outliers] = medAbsDeviation(dwellTime_filtered);
                dwellTime = dwellTime_filtered(~dwellTime_outliers);
                dwellTime(dwellTime <= 0) = [];
                dwellTime(dwellTime > 1) = []; 

                travelTime_filtered = travelTime(valid_keypair_mask(1:min(length(valid_keypair_mask), length(travelTime))));
                [travelTime_outliers] = medAbsDeviation(travelTime_filtered);
                travelTime = travelTime_filtered(~travelTime_outliers);
                travelTime(travelTime < 0.1) = [];
                travelTime(travelTime > 1.5) = []; 

                velocity = velocitydata(:, 1) ./ velocitydata(:, 2);
                [Velocity_outliers] = medAbsDeviation(velocity);
                velocity = velocity(~Velocity_outliers);

                IDrep = [IDrep; repmat(participants(p), size(travelTime))];
                handrep = [handrep; repmat(hand(h), size(travelTime))];
                exerep = [exerep; repmat(exerc(e), size(travelTime))];
                recrep = [recrep; repmat(rec, size(travelTime))];

                IDrep1 = [IDrep1; repmat(participants(p), size(dwellTime))];
                handrep1 = [handrep1; repmat(hand(h), size(dwellTime))];
                exerep1 = [exerep1; repmat(exerc(e), size(dwellTime))];
                recrep1 = [recrep1; repmat(rec, size(dwellTime))];

                IDrep2 = [IDrep2; repmat(participants(p), size(velocity))];
                handrep2 = [handrep2; repmat(hand(h), size(velocity))];
                exerep2 = [exerep2; repmat(exerc(e), size(velocity))];
                recrep2 = [recrep2; repmat(rec, size(velocity))];

                allTravel = [allTravel; travelTime];
                allDwell = [allDwell; dwellTime];
                allV = [allV; velocity];

                PCvelocity = (velocity - velocity(1)) / velocity(1) * 100;
                [PCvelocity_outliers] = medAbsDeviation(PCvelocity);
                PCvelocity = PCvelocity(~PCvelocity_outliers);
                time_steps = 1:numel(PCvelocity);                  
                coefficients = polyfit(time_steps, PCvelocity, 1);

                %% Define output matrix
                clabel = ['P', num2str(p), 'H', num2str(h), exerc{e}, 'Exercise', 'Rec', num2str(rec)];
                if ~isfield(metrics, clabel)
                    metrics.(clabel) = struct();
                end
                        
                metrics.(clabel).nrTaps = numel(dwellTime);
                metrics.(clabel).dysmetriaScore = nanmean(dysmetScores);
                metrics.(clabel).incoordScore = nanvar(travelTime);
                metrics.(clabel).nrErrors = nansum(~correctLetters);             
                metrics.(clabel).sequenceEffect = coefficients(1);
                metrics.(clabel).velocity = (nanmean(velocitydata(:, 1)) ./ (nanmean(velocitydata(:, 2))));
                metrics.(clabel).travelTimes = nanmean(travelTime);
                metrics.(clabel).dwellTimes = nanmean(dwellTime);
                
                IDP = [IDP; participants(p)];
                HANDP = [HANDP, hand(h)];
                RECORP = [RECORP, rec]; 
                EXERCP = [EXERCP, exerc(e)];
            end
        end
    end    
end

save('metrics.mat', 'metrics')
travelTable = table(IDrep, handrep, exerep, recrep, allTravel);
dwellTable = table(IDrep1, handrep1, exerep1, recrep1, allDwell);
velocitytable = table(IDrep2, handrep2, exerep2, recrep2, allV);

writetable(travelTable, 'all_travel_times.csv');
writetable(dwellTable, 'all_dwell_times.csv');
writetable(velocitytable, 'all_velocity.csv');

%% Save to metrics
fieldNames = fieldnames(metrics);
Subfieldnames = fieldnames(metrics.(fieldNames{1}));
cellArray = cell(length(fieldNames)+1, length(Subfieldnames)+1);
cellArray(1, 2:end) = Subfieldnames;
cellArray{1, 1} = 'participants';

for row = 1:numel(fieldNames)
    field = fieldNames{row};
    values = cell2mat(struct2cell(metrics.(field)))';
    [cellArray{row+1, 1}] = deal(field);
    for ce = 1:numel(values)
        [cellArray{row+1, ce+1}] = num2cell(values(ce));
    end
end

metricstable = cell2table(cellArray(1:end, :));
% csvFilePath = '\\yourpath_to_csv\metrics.csv';
% writetable(metricstable, csvFilePath, 'WriteRowNames', true);

%% Parse and organize metrics
field_names = fieldnames(metrics);

participant_ids = {};
hands = {};
exercises = {};
recordings = [];
nrTaps = [];
dysmetriaScore = [];
incoordScore = [];
nrErrors = [];
sequenceEffect = [];
velocity = [];
travelTimes = [];
dwellTimes = [];

for i = 1:length(field_names)
    field_name = field_names{i};
    
    % Parse field name: e.g., 'P1H1preEExerciseRec1'
    p_idx = strfind(field_name, 'P');
    h_idx = strfind(field_name, 'H');
    participant_id = field_name(p_idx+1:h_idx-1);
    
    hand_num = field_name(h_idx+1:h_idx+1);
    if strcmp(hand_num, '1')
        hand = 'right';
    else
        hand = 'left';
    end
    
    ex_idx = strfind(field_name, 'Exercise');
    exercise_start = h_idx + 2;
    exercise_end = ex_idx - 1;
    exercise = field_name(exercise_start:exercise_end);
    
    rec_idx = strfind(field_name, 'Rec');
    recording = str2double(field_name(rec_idx+3:end));
    
    participant_ids{end+1} = participant_id;
    hands{end+1} = hand;
    exercises{end+1} = exercise;
    recordings(end+1) = recording;
    
    nrTaps(end+1) = metrics.(field_name).nrTaps;
    dysmetriaScore(end+1) = metrics.(field_name).dysmetriaScore;
    incoordScore(end+1) = metrics.(field_name).incoordScore;
    nrErrors(end+1) = metrics.(field_name).nrErrors;
    sequenceEffect(end+1) = metrics.(field_name).sequenceEffect;
    velocity(end+1) = metrics.(field_name).velocity;
    travelTimes(end+1) = metrics.(field_name).travelTimes;
    dwellTimes(end+1) = metrics.(field_name).dwellTimes;
end

organized_table = table(participant_ids', hands', exercises', recordings', ...
                       nrTaps', dysmetriaScore', incoordScore', nrErrors', ...
                       sequenceEffect', velocity', travelTimes', dwellTimes', ...
                       'VariableNames', {'ParticipantID', 'Hand', 'Exercise', 'Recording', ...
                                       'nrTaps', 'dysmetriaScore', 'incoordScore', 'nrErrors', ...
                                       'sequenceEffect', 'velocity', 'travelTimes', 'dwellTimes'});

fprintf('Organized metrics table created with %d rows\n', height(organized_table));
disp(organized_table(1:10, :));

save('organized_metrics.mat', 'organized_table');
writetable(organized_table, 'organized_metrics.csv');


end