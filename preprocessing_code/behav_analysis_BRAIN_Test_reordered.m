function [] = behav_analysis_BRAIN_Test_reordered(Const, Paths)
%% BEHAV_ANALYSIS_BRAIN_Test 

[Const] = storeConst();
participant_dirs = dir(Paths.Data);
participant_dirs = participant_dirs([participant_dirs.isdir]); % Keep only directories
participant_dirs = participant_dirs(~ismember({participant_dirs.name}, {'.', '..'})); % Remove . and ..

% Extract participant IDs dynamically
participants = {participant_dirs.name};
fprintf('Found %d participants: %s\n', length(participants), strjoin(participants, ', '));



exerc = {'preE', 'postE'};
% meds  = { 'OFF', 'ON', 'ON-OFF'};
hand  = {'right', 'left'};
% PLOTTING
FLEXIBLE_YLIMS = true;

allTravel = [];
allDwell =[];
allV = [];
handrep=[];
exerep=[];
IDrep=[];
recrep=[];
IDrep1=[];
exerep1=[];
handrep1=[];
recrep1=[];
HANDP = [];
EXERCP=[];
IDP=[];
RECORP=[];
IDrep2=[];
exerep2=[];
handrep2=[];
recrep2=[];
metrics = [];
%% extract participant characteristics
for p = 1:numel(participants)
    subj = participants{p};
    part_folder_info = dir([Paths.Data, subj]);
    all_fnames       = {part_folder_info.name};
    
    for h = 1:numel(hand)
        conds = {};
         for e = 1:numel(exerc)
            
            
            for rec = 1:2
               
                TravelDistribution = [];
              
                DwellDistribution = [];
                fileIdx =~cellfun(@isempty, strfind(all_fnames, [exerc{e}, 'xerc_', hand{h}, '_v', num2str(rec)]));            
                if sum(fileIdx)==0
                   continue
                end
            
                conds     = {conds{:}, [exerc{e},num2str(rec)]};
                data = readtable(fullfile(Paths.Data,subj,all_fnames{fileIdx}));

               % Dysmetria score(DS; a measure of the average accuracy of key strikes where the central key
               % scores 1, adjacent keys are 2, and all other keys are 3).
                correctLetters = strcmp(data.letters, 's')  + strcmp(data.letters, 'semicolon');
                adjacentKeys   =  (strcmp(data.letters, 'a') + strcmp(data.letters, 'w') + strcmp(data.letters, 'e') + strcmp(data.letters, 'd') + strcmp(data.letters, 'x') + strcmp(data.letters, 'z') + ...
                                      strcmp(data.letters, 'l') + strcmp(data.letters, 'p') + strcmp(data.letters, 'leftbracket') + strcmp(data.letters, 'quote') + ...
                                      strcmp(data.letters, 'slash') + strcmp(data.letters, 'period')) * 2; % multiply by 2 to fill i with twos
                otherKeys      =  ~(correctLetters + adjacentKeys)  * 3;      

                dysmetScores = correctLetters + adjacentKeys + otherKeys;
                    
                strokeItvl = [];
                strokeItvl = diff(data.strokeOnset);
                strokeItvl = [0;strokeItvl];  

                dwellTime = [];
                dwellTime = data.dwellTime(2:end);
                dwellTime(dwellTime<0.01) = [];  % computer error of registering events that are too short, occurred in one condition in P2
                [dwellTime_outliers] = medAbsDeviation(dwellTime);
                dwellTime = dwellTime(~dwellTime_outliers);
                % in case dwellTime classified as outliers while not
                % matching with stroke Itvv

                                        
            %% recreate the keypairs to calculate the distance    
                travelTime = strokeItvl(2:end) - data.dwellTime(1:end-1);% minus the corresponding stroke onset. 
                % the last stroke interval only includes the dwell time for the second last stroke                               
                keypairs= [];
                for key = 1: length(data.letters)-1 
                    keypairs{key} = strcat(data.letters{key},data.letters{key+1});
                end
                
                keypairs =reshape(keypairs,[],1);

                %% define the distance between pressed keys

                velocitydata = [];

                valid_keypair_mask = true(length(keypairs), 1);  % Track which keypairs are valid

                for pair = 1: length(keypairs)
                    if travelTime(pair) >0.1 && travelTime(pair) < 1.5
                        if strcmp(keypairs{pair},'ssemicolon')|| strcmp(keypairs{pair},'semicolons')||...
                                strcmp(keypairs{pair},'xslash')|| strcmp(keypairs{pair},'slashx')||...
                                strcmp(keypairs{pair},'eleftbracket')|| strcmp(keypairs{pair},'leftbrackete')||...
                                strcmp(keypairs{pair},'zperiod')|| strcmp(keypairs{pair},'periodz')||...
                                strcmp(keypairs{pair},'al')||strcmp(keypairs{pair},'la')|| strcmp(keypairs{pair},'wp')||...  %
                                strcmp(keypairs{pair},'dquote')||strcmp(keypairs{pair},'quoted')|| strcmp(keypairs{pair},'pw')

                            velocitydata(end+1,:) = [15,travelTime(pair)];% distance is 15 cm between these key pairs

                        elseif strcmp(keypairs{pair},'ld')|| strcmp(keypairs{pair},'dl')

                            velocitydata(end+1,:) = [11.5,travelTime(pair)];
                        elseif strcmp(keypairs{pair},'lx')|| strcmp(keypairs{pair},'xl')

                            velocitydata(end+1,:) = [12.5,travelTime(pair)];

                        elseif strcmp(keypairs{pair},'sl')|| strcmp(keypairs{pair},'ls')|| ...
                                strcmp(keypairs{pair},'dsemicolon')|| strcmp(keypairs{pair},'semicolond')|| ...
                                strcmp(keypairs{pair},'ak')||strcmp(keypairs{pair},'ka')||...
                                 strcmp(keypairs{pair},'eo')||strcmp(keypairs{pair},'oe')||...
                                strcmp(keypairs{pair},'fquote')||strcmp(keypairs{pair},'quotef')||...
                                strcmp(keypairs{pair},'xperiod')||strcmp(keypairs{pair},'periodx')

                            velocitydata(end+1,:) = [13.5,travelTime(pair)];


                       

                        elseif strcmp(keypairs{pair},'sp')|| strcmp(keypairs{pair},'ps')||...
                                strcmp(keypairs{pair},'zl')|| strcmp(keypairs{pair},'lz')||...
                                 strcmp(keypairs{pair},'esemicolon')|| strcmp(keypairs{pair},'semicolone')||...
                                  strcmp(keypairs{pair},'xsemicolon')|| strcmp(keypairs{pair},'semicolonx')||...
                              strcmp(keypairs{pair},'speriod')|| strcmp(keypairs{pair},'periods')||...
                                strcmp(keypairs{pair},'slashd')|| strcmp(keypairs{pair},'dslash')

                            velocitydata(end+1,:) = [14.5,travelTime(pair)];

                        elseif strcmp(keypairs{pair},'xquote')|| strcmp(keypairs{pair},'quotex')|| ...
                                strcmp(keypairs{pair},'sslash')|| strcmp(keypairs{pair},'slashs')|| ...
                                 strcmp(keypairs{pair},'wsemicolon')|| strcmp(keypairs{pair},'semicolonw')|| ...
                                strcmp(keypairs{pair},'zsemicolon')|| strcmp(keypairs{pair},'semicolonz')
                            velocitydata(end+1,:) = [16,travelTime(pair)];

                        elseif strcmp(keypairs{pair},'sleftbracket')|| strcmp(keypairs{pair},'leftbrackets')|| ...
                                strcmp(keypairs{pair},'ap')|| strcmp(keypairs{pair},'pa')|| ...
                                strcmp(keypairs{pair},'speriod')|| strcmp(keypairs{pair},'pa')|| ...
                                strcmp(keypairs{pair},'qleftbracket')|| strcmp(keypairs{pair},'leftbracketq')
                            velocitydata(end+1,:) = [16.5,travelTime(pair)];

                        elseif strcmp(keypairs{pair},'squote')|| strcmp(keypairs{pair},'quotes')|| ...
                                strcmp(keypairs{pair},'asemicolon')|| strcmp(keypairs{pair},'semicolona')|| ...
                                strcmp(keypairs{pair},'zslash')|| strcmp(keypairs{pair},'slashz')|| ...
                                strcmp(keypairs{pair},'qp')|| strcmp(keypairs{pair},'pq') || ...
                                strcmp(keypairs{pair},'wleftbracket')|| strcmp(keypairs{pair},'leftbracketw')
                            velocitydata(end+1,:) = [17,travelTime(pair)];
                         elseif strcmp(keypairs{pair},'slasha')|| strcmp(keypairs{pair},'aslash')|| ...
                                strcmp(keypairs{pair},'quotez')|| strcmp(keypairs{pair},'zquote')|| ...
                                strcmp(keypairs{pair},'zslash')|| strcmp(keypairs{pair},'slashz')|| ...
                                strcmp(keypairs{pair},'qp')|| strcmp(keypairs{pair},'pq') || ...
                                strcmp(keypairs{pair},'wleftbracket')|| strcmp(keypairs{pair},'leftbracketw')
                            velocitydata(end+1,:) = [18,travelTime(pair)];
                         elseif strcmp(keypairs{pair},'aquote')|| strcmp(keypairs{pair},'quotea')|| ...      
                                strcmp(keypairs{pair},'qleftbracket')|| strcmp(keypairs{pair},'leftbracketq')
                            velocitydata(end+1,:) = [19,travelTime(pair)];                    

                        else
                             fprintf('Unrecognized keypair: %s\n', keypairs{pair});
                             valid_keypair_mask(pair) = false;
                        end
                    else
                        velocitydata = velocitydata;
                        continue             
                     
                    end

                end
                dwellTime_filtered = dwellTime(valid_keypair_mask(1:min(length(valid_keypair_mask), length(dwellTime))));
                [dwellTime_outliers] = medAbsDeviation(dwellTime_filtered);
                dwellTime =dwellTime_filtered(~dwellTime_outliers);
                dwellTime(dwellTime <=0) = [];
                dwellTime(dwellTime> 1) = []; 
                travelTime_filtered =travelTime(valid_keypair_mask(1:min(length(valid_keypair_mask), length(travelTime))));
                [travelTime_outliers] = medAbsDeviation( travelTime_filtered);
                travelTime = travelTime_filtered(~ travelTime_outliers);
                travelTime(travelTime<0.1) = [];
                travelTime(travelTime> 1.5) = []; 
                velocity = velocitydata(:,1)./velocitydata(:,2);
                [Velocity_outliers] = medAbsDeviation(velocity);
                velocity = velocity(~Velocity_outliers);



                IDrep = [IDrep;repmat(participants(p),size(travelTime))];
                handrep= [handrep;repmat(hand(h),size(travelTime))];
                exerep= [exerep;repmat(exerc(e),size(travelTime))];
                recrep = [recrep;repmat(rec,size(travelTime))];

                IDrep1 = [IDrep1;repmat(participants(p),size(dwellTime))];
                handrep1= [handrep1;repmat(hand(h),size(dwellTime))];
                exerep1= [exerep1;repmat(exerc(e),size(dwellTime))];
                recrep1 = [recrep1;repmat(rec,size(dwellTime))];
                IDrep2 = [IDrep2;repmat(participants(p),size(velocity))];
                handrep2= [handrep2;repmat(hand(h),size(velocity))];
                exerep2= [exerep2;repmat(exerc(e),size(velocity))];
                recrep2 = [recrep2;repmat(rec,size(velocity))];
                allTravel = [allTravel; travelTime];
                allDwell= [allDwell; dwellTime];
                allV = [allV;velocity];

               
                PCvelocity = (velocity-velocity(1))/velocity(1)*100;
                [PCvelocity_outliers] = medAbsDeviation(PCvelocity);
                PCvelocity = PCvelocity(~PCvelocity_outliers);
                time_steps = 1:numel(PCvelocity);                  
                coefficients = polyfit(time_steps,PCvelocity,1);

               %% define output matrix
              
                clabel = ['P',num2str(p),'H',num2str(h), exerc{e},'Exercise', 'Rec', num2str(rec)];
                if ~isfield(metrics,clabel)
                    metrics.(clabel) = struct();
                end
                        
                metrics.(clabel).nrTaps           = numel(dwellTime);%taps that are not neighbouring
                metrics.(clabel).dysmetriaScore =  nanmean(dysmetScores);
                % Incoordination score (IS; a measure of rhythm given by the variance
                    % in the traveling times between key presses)
                
                metrics.(clabel).incoordScore   = nanvar(travelTime);
                

              
                metrics.(clabel).nrErrors       = nansum(~correctLetters); % but this if accidently tapped twice also not considered as error             
                metrics.(clabel).sequenceEffect =  coefficients(1);
                metrics.(clabel).velocity     = (nanmean(velocitydata(:,1))./(nanmean(velocitydata(:,2))));

                metrics.(clabel).travelTimes    =nanmean (travelTime);
                metrics.(clabel).dwellTimes     = nanmean(dwellTime);

                
                %% 
                IDP = [IDP;participants(p)];
                HANDP = [HANDP,hand(h)];
                RECORP = [RECORP, rec]; 
                EXERCP = [EXERCP, exerc(e)];
                % dwellTime = [];

              
            end
        end
        
       
    end    
end
save('metrics.mat','metrics')
travelTable=table(IDrep,handrep,exerep,recrep,allTravel);
dwellTable = table(IDrep1,handrep1,exerep1,recrep1,allDwell);
velocitytable = table(IDrep2,handrep2,exerep2,recrep2,allV);


writetable(travelTable,'all_travel_times.csv');
writetable(dwellTable,'all_dwell_times.csv');
writetable (velocitytable,'all_velocity.csv');


     

                 
               


%% save to metrics

fieldNames = fieldnames(metrics);
Subfieldnames = fieldnames(metrics.(fieldNames{1}));
cellArray = cell(length(fieldNames)+1,length(Subfieldnames)+1);
cellArray (1,2:end) =Subfieldnames;
cellArray{1,1} = 'participants';

for row = 1:numel(fieldNames)
    field=fieldNames{row};
 
    values=cell2mat(struct2cell(metrics.(field)))';
    [cellArray{row+1,1}]=deal(field);
    for ce = 1: numel(values)
        [cellArray{row+1,ce+1}]=num2cell(values(ce));
    end
end

metricstable = cell2table(cellArray(1:end,:));
csvFilePath = '\\rdsfcifs.acrc.bris.ac.uk\Fischer_StudentProjects\PhD_Rui_Ni\metrics.csv';
writetable(metricstable,csvFilePath,'WriteRowNames',true);
%%
field_names = fieldnames(metrics);

% Initialize arrays to store table data
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

% Parse each field name and extract data
for i = 1:length(field_names)
    field_name = field_names{i};
    
    % Parse the field name: e.g., 'P1H1preEExerciseRec1'
    % Extract participant ID
    p_idx = strfind(field_name, 'P');
    h_idx = strfind(field_name, 'H');
    participant_id = field_name(p_idx+1:h_idx-1);
    
    % Extract hand
    ex_idx = strfind(field_name, 'Exercise');
    hand_num = field_name(h_idx+1:h_idx+1);
    if strcmp(hand_num, '1')
        hand = 'right';
    else
        hand = 'left';
    end
    
    % Extract exercise type
    exercise_start = h_idx + 2;
    exercise_end = ex_idx - 1;
    exercise = field_name(exercise_start:exercise_end);
    
    % Extract recording number
    rec_idx = strfind(field_name, 'Rec');
    recording = str2double(field_name(rec_idx+3:end));
    
    % Store parsed information
    participant_ids{end+1} = participant_id;
    hands{end+1} = hand;
    exercises{end+1} = exercise;
    recordings(end+1) = recording;
    
    % Extract metrics
    nrTaps(end+1) = metrics.(field_name).nrTaps;
    dysmetriaScore(end+1) = metrics.(field_name).dysmetriaScore;
    incoordScore(end+1) = metrics.(field_name).incoordScore;
    nrErrors(end+1) = metrics.(field_name).nrErrors;
    sequenceEffect(end+1) = metrics.(field_name).sequenceEffect;
    velocity(end+1) = metrics.(field_name).velocity;
    travelTimes(end+1) = metrics.(field_name).travelTimes;
    dwellTimes(end+1) = metrics.(field_name).dwellTimes;
end

% Create organized table
organized_table = table(participant_ids', hands', exercises', recordings', ...
                       nrTaps', dysmetriaScore', incoordScore', nrErrors', ...
                       sequenceEffect', velocity', travelTimes', dwellTimes', ...
                       'VariableNames', {'ParticipantID', 'Hand', 'Exercise', 'Recording', ...
                                       'nrTaps', 'dysmetriaScore', 'incoordScore', 'nrErrors', ...
                                       'sequenceEffect', 'velocity', 'travelTimes', 'dwellTimes'});

% Display the organized table
fprintf('Organized metrics table created with %d rows\n', height(organized_table));
disp(organized_table(1:10,:)); % Show first 10 rows

% Save the organized table
save('organized_metrics.mat', 'organized_table');
writetable(organized_table, 'organized_metrics.csv');

%% Calculate pre/post average behavior for right/left hand

% Get unique participants
unique_participants = unique(organized_table.ParticipantID);
unique_hands = {'right', 'left'};
behavior_measures = {'nrTaps', 'dysmetriaScore', 'incoordScore', 'nrErrors', ...
                    'sequenceEffect', 'velocity', 'travelTimes', 'dwellTimes'};

% Initialize results structure
pre_post_results = struct();

for p = 1:length(unique_participants)
    participant = unique_participants{p};
    
    for h = 1:length(unique_hands)
        hand = unique_hands{h};
        
        % Filter data for current participant and hand
        participant_hand_data = organized_table(strcmp(organized_table.ParticipantID, participant) & ...
                                              strcmp(organized_table.Hand, hand), :);
        
        if isempty(participant_hand_data)
            continue;
        end
        
        % Separate pre and post exercise data
        pre_data = participant_hand_data(strcmp(participant_hand_data.Exercise, 'preE'), :);
        post_data = participant_hand_data(strcmp(participant_hand_data.Exercise, 'postE'), :);
        
        % Create field name for results
        result_field = ['P', participant, '_', hand];
        pre_post_results.(result_field) = struct();
        
        % Calculate averages for each behavior measure
        for b = 1:length(behavior_measures)
            behavior = behavior_measures{b};
            
            % Pre-exercise average (across recordings)
            if ~isempty(pre_data)
                pre_values = pre_data.(behavior);
                pre_values = pre_values(~isnan(pre_values)); % Remove NaN values
                pre_post_results.(result_field).(['pre_', behavior]) = mean(pre_values);
                pre_post_results.(result_field).(['pre_', behavior, '_std']) = std(pre_values);
                pre_post_results.(result_field).(['pre_', behavior, '_n']) = length(pre_values);
            else
                pre_post_results.(result_field).(['pre_', behavior]) = NaN;
                pre_post_results.(result_field).(['pre_', behavior, '_std']) = NaN;
                pre_post_results.(result_field).(['pre_', behavior, '_n']) = 0;
            end
            
            % Post-exercise average (across recordings)
            if ~isempty(post_data)
                post_values = post_data.(behavior);
                post_values = post_values(~isnan(post_values)); % Remove NaN values
                pre_post_results.(result_field).(['post_', behavior]) = mean(post_values);
                pre_post_results.(result_field).(['post_', behavior, '_std']) = std(post_values);
                pre_post_results.(result_field).(['post_', behavior, '_n']) = length(post_values);
            else
                pre_post_results.(result_field).(['post_', behavior]) = NaN;
                pre_post_results.(result_field).(['post_', behavior, '_std']) = NaN;
                pre_post_results.(result_field).(['post_', behavior, '_n']) = 0;
            end
            
            % Calculate change (post - pre) and percent change
            if ~isnan(pre_post_results.(result_field).(['pre_', behavior])) && ...
               ~isnan(pre_post_results.(result_field).(['post_', behavior]))
                
                pre_val = pre_post_results.(result_field).(['pre_', behavior]);
                post_val = pre_post_results.(result_field).(['post_', behavior]);
                
                pre_post_results.(result_field).(['change_', behavior]) = post_val - pre_val;
                
                if pre_val ~= 0
                    pre_post_results.(result_field).(['percent_change_', behavior]) = ...
                        ((post_val - pre_val) / abs(pre_val)) * 100;
                else
                    pre_post_results.(result_field).(['percent_change_', behavior]) = NaN;
                end
            else
                pre_post_results.(result_field).(['change_', behavior]) = NaN;
                pre_post_results.(result_field).(['percent_change_', behavior]) = NaN;
            end
        end
    end
end

% Save pre/post results
save('pre_post_results.mat', 'pre_post_results');

%% Create summary table for pre/post analysis
summary_participant_ids = {};
summary_hands = {};
pre_nrTaps = []; post_nrTaps = []; change_nrTaps = [];
pre_dysmetriaScore = []; post_dysmetriaScore = []; change_dysmetriaScore = [];
pre_incoordScore = []; post_incoordScore = []; change_incoordScore = [];
pre_nrErrors = []; post_nrErrors = []; change_nrErrors = [];
pre_sequenceEffect = []; post_sequenceEffect = []; change_sequenceEffect = [];
pre_velocity = []; post_velocity = []; change_velocity = [];
pre_travelTimes = []; post_travelTimes = []; change_travelTimes = [];
pre_dwellTimes = []; post_dwellTimes = []; change_dwellTimes = [];

result_fields = fieldnames(pre_post_results);
for i = 1:length(result_fields)
    field_name = result_fields{i};
    
    % Parse participant and hand from field name
    parts = strsplit(field_name, '_');
    participant = parts{1}(2:end); % Remove 'P' prefix
    hand = parts{2};
    
    summary_participant_ids{end+1} = participant;
    summary_hands{end+1} = hand;
    
    % Extract pre/post values for each measure
    pre_nrTaps(end+1) = pre_post_results.(field_name).pre_nrTaps;
    post_nrTaps(end+1) = pre_post_results.(field_name).post_nrTaps;
    change_nrTaps(end+1) = pre_post_results.(field_name).change_nrTaps;
    
    pre_dysmetriaScore(end+1) = pre_post_results.(field_name).pre_dysmetriaScore;
    post_dysmetriaScore(end+1) = pre_post_results.(field_name).post_dysmetriaScore;
    change_dysmetriaScore(end+1) = pre_post_results.(field_name).change_dysmetriaScore;
    
    pre_incoordScore(end+1) = pre_post_results.(field_name).pre_incoordScore;
    post_incoordScore(end+1) = pre_post_results.(field_name).post_incoordScore;
    change_incoordScore(end+1) = pre_post_results.(field_name).change_incoordScore;
    
    pre_nrErrors(end+1) = pre_post_results.(field_name).pre_nrErrors;
    post_nrErrors(end+1) = pre_post_results.(field_name).post_nrErrors;
    change_nrErrors(end+1) = pre_post_results.(field_name).change_nrErrors;
    
    pre_sequenceEffect(end+1) = pre_post_results.(field_name).pre_sequenceEffect;
    post_sequenceEffect(end+1) = pre_post_results.(field_name).post_sequenceEffect;
    change_sequenceEffect(end+1) = pre_post_results.(field_name).change_sequenceEffect;
    
    pre_velocity(end+1) = pre_post_results.(field_name).pre_velocity;
    post_velocity(end+1) = pre_post_results.(field_name).post_velocity;
    change_velocity(end+1) = pre_post_results.(field_name).change_velocity;
    
    pre_travelTimes(end+1) = pre_post_results.(field_name).pre_travelTimes;
    post_travelTimes(end+1) = pre_post_results.(field_name).post_travelTimes;
    change_travelTimes(end+1) = pre_post_results.(field_name).change_travelTimes;
    
    pre_dwellTimes(end+1) = pre_post_results.(field_name).pre_dwellTimes;
    post_dwellTimes(end+1) = pre_post_results.(field_name).post_dwellTimes;
    change_dwellTimes(end+1) = pre_post_results.(field_name).change_dwellTimes;
end

% Create summary table
summary_table = table(summary_participant_ids', summary_hands', ...
                     pre_nrTaps', post_nrTaps', change_nrTaps', ...
                     pre_dysmetriaScore', post_dysmetriaScore', change_dysmetriaScore', ...
                     pre_incoordScore', post_incoordScore', change_incoordScore', ...
                     pre_nrErrors', post_nrErrors', change_nrErrors', ...
                     pre_sequenceEffect', post_sequenceEffect', change_sequenceEffect', ...
                     pre_velocity', post_velocity', change_velocity', ...
                     pre_travelTimes', post_travelTimes', change_travelTimes', ...
                     pre_dwellTimes', post_dwellTimes', change_dwellTimes', ...
                     'VariableNames', {'ParticipantID', 'Hand', ...
                                     'Pre_nrTaps', 'Post_nrTaps', 'Change_nrTaps', ...
                                     'Pre_dysmetriaScore', 'Post_dysmetriaScore', 'Change_dysmetriaScore', ...
                                     'Pre_incoordScore', 'Post_incoordScore', 'Change_incoordScore', ...
                                     'Pre_nrErrors', 'Post_nrErrors', 'Change_nrErrors', ...
                                     'Pre_sequenceEffect', 'Post_sequenceEffect', 'Change_sequenceEffect', ...
                                     'Pre_velocity', 'Post_velocity', 'Change_velocity', ...
                                     'Pre_travelTimes', 'Post_travelTimes', 'Change_travelTimes', ...
                                     'Pre_dwellTimes', 'Post_dwellTimes', 'Change_dwellTimes'});

% Display summary
fprintf('\nPre/Post Exercise Summary Table:\n');
disp(summary_table);

% Save summary table
save('pre_post_summary.mat', 'summary_table');
writetable(summary_table, 'pre_post_summary.csv');

%% Display overall statistics by hand
fprintf('\n=== OVERALL STATISTICS BY HAND ===\n');

for h = 1:length(unique_hands)
    hand = unique_hands{h};
    hand_data = summary_table(strcmp(summary_table.Hand, hand), :);
    
    fprintf('\n--- %s Hand ---\n', upper(hand));
    
    for b = 1:length(behavior_measures)
        behavior = behavior_measures{b};
        
        pre_col = ['Pre_', behavior];
        post_col = ['Post_', behavior];
        change_col = ['Change_', behavior];
        
        pre_values = hand_data.(pre_col);
        post_values = hand_data.(post_col);
        change_values = hand_data.(change_col);
        
        % Remove NaN values
        pre_values = pre_values(~isnan(pre_values));
        post_values = post_values(~isnan(post_values));
        change_values = change_values(~isnan(change_values));
        
        fprintf('%s:\n', behavior);
        fprintf('  Pre:  Mean=%.3f, SD=%.3f, N=%d\n', mean(pre_values), std(pre_values), length(pre_values));
        fprintf('  Post: Mean=%.3f, SD=%.3f, N=%d\n', mean(post_values), std(post_values), length(post_values));
        fprintf('  Change: Mean=%.3f, SD=%.3f, N=%d\n', mean(change_values), std(change_values), length(change_values));
        
        % Perform paired t-test if we have paired data
        if length(pre_values) == length(post_values) && length(pre_values) > 1
            [h_test, p_val] = ttest(pre_values, post_values);
            if p_val < 0.05
                fprintf('Paired t-test: p=%.4f (significant)\n', p_val);
            else
                fprintf('Paired t-test: p=%.4f (not significant)\n', p_val);
            end
        end
        fprintf('\n');
    end
end

fprintf('Analysis complete! Files saved:\n');
fprintf('- organized_metrics.mat & .csv\n');
fprintf('- pre_post_results.mat\n');
fprintf('- pre_post_summary.mat & .csv\n');
        %%
        fprintf('\nGenerating per-file summary table...\n');

% Step 1: Prepare table from all metrics
summary_rows = {};
fieldNames = fieldnames(metrics);

for i = 1:numel(fieldNames)
    field = fieldNames{i};  % e.g., 'C1P2H1'
    entry = metrics.(field);

    % Parse participant/hand/exercise/recording from field label
    pid = IDP{i};
    hand_ = HANDP{i};
    rec_ = RECORP(i);
    exer_ = EXERCP{i};

    % Build row
    row = {pid, hand_, exer_, rec_, ...
           entry.travelTimes, ...
           entry.dwellTimes, ...
           entry.velocity, ...
           entry.nrErrors};

    summary_rows(end+1, :) = row;
end

% Step 2: Create table
summary_table = cell2table(summary_rows, ...
    'VariableNames', {'Participant', 'Hand', 'Exercise', 'Recording', ...
                      'AvgTravelTime', 'AvgDwellTime', 'AvgVelocity', 'NumErrors'});

% Step 3: Write to CSV
writetable(summary_table, 'summary_per_file_metrics.csv');
fprintf('✅ Summary written to: summary_per_file_metrics.csv\n');

                
end