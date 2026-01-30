%% Load and process the demographic data
% Load the Excel file
demographic_data = readtable('fulltable1.xlsx');

% Display the structure
fprintf('Demographic data loaded with %d rows and %d columns\n', height(demographic_data), width(demographic_data));
fprintf('Column names: %s\n', strjoin(demographic_data.Properties.VariableNames, ', '));

% Clean column names (remove spaces and special characters)
demographic_data.Properties.VariableNames = matlab.lang.makeValidName(demographic_data.Properties.VariableNames);

% Display first few rows
fprintf('\nFirst 5 rows of demographic data:\n');
disp(demographic_data(1:5,:));

%% Calculate Disease Progression using PCA
variables_for_pca = { 'dailyMedicationDose', 'diseaseDuration'};

% Check if the specified columns exist
missing_vars = setdiff(variables_for_pca, demographic_data.Properties.VariableNames);
if ~isempty(missing_vars)
    error('Missing variables in the table: %s', strjoin(missing_vars, ', '));
end

% Extract data and remove rows with any missing values
id_col = 'patientsId';



% Compute mean values per participant
collapsed_data = varfun(@mean, demographic_data, ...
    'InputVariables', variables_for_pca, ...
    'GroupingVariables', id_col);

% Extract the matrix for PCA (skip first column: ParticipantID)
X = collapsed_data{:, 2:end};

% Remove rows with NaN
valid_rows = all(~isnan(X), 2);
X_clean = X(valid_rows, :);
std_vals = std(X_clean);
disp('STD per column:');
disp(std_vals);

% Keep only columns with non-zero variance
X_clean = X_clean(:, std_vals > 1e-5);
X_clean_z = zscore(X_clean);  % mean 0, std 1 per column
[coeff, score, latent, tsquared, explained, mu] = pca(X_clean_z);
% Perform PCA
%[coeff, score, latent, tsquared, explained, mu] = pca(X_clean);

% Display the results
fprintf('\nPCA completed on %d participants.\n', size(X_clean, 1));
disp('Explained variance by each component (%):');
disp(explained);

% Optionally, plot explained variance
figure;
pareto(explained);
title('PCA on Trait-Level Demographics');
xlabel('Principal Component');
ylabel('Variance Explained (%)');
%% Get unique participants from demographic data and prepare for merging
% Get unique participants (since demographic data might have multiple rows per participant)
unique_participants = unique(demographic_data.patientsId);
fprintf('\nUnique participants in demographic data: %s\n', mat2str(unique_participants));

% Create a summary demographic table with one row per participant
summary = table();
for i = 1:length(unique_participants)
    participant_id = unique_participants(i);
    participant_rows = demographic_data(demographic_data.patientsId == participant_id, :);
    
    % Take the first row for demographic info (should be same across rows for same participant)
    demo_row = participant_rows(1, :);
    
    % For disease progression, take the mean if multiple values exist
    % if sum(~isnan(participant_rows.diseaseProgression)) > 0
    %     demo_row.diseaseProgression = nanmean(participant_rows.diseaseProgression);
    % end
    
    summary = [summary; demo_row];
end

fprintf('Demo summary table created with %d participants\n', height(summary));

%% Load organized table (assuming it exists from previous code)
% If organized_table doesn't exist, load it
if ~exist('organized_table', 'var')
    if exist('organized_metrics.mat', 'file')
        load('organized_metrics.mat', 'organized_table');
        fprintf('Loaded organized_table from file\n');
    else
        error('organized_table not found. Please run the previous analysis first.');
    end
end

fprintf('Organized table has %d rows\n', height(organized_table));
fprintf('Unique participants in organized table: %s\n', strjoin(unique(organized_table.ParticipantID), ', '));

%% Filter organized_table to include only participants in demographic data
% Convert participant IDs to same format for comparison
organized_table.ParticipantID_num = str2double(organized_table.ParticipantID);
participants_in_both = intersect(unique(organized_table.ParticipantID_num), unique_participants);
fprintf('\nParticipants present in both tables: %s\n', mat2str(participants_in_both));

% Filter organized table
filtered_organized_table = organized_table(ismember(organized_table.ParticipantID_num, participants_in_both), :);
fprintf('Filtered organized table: %d rows (from %d original rows)\n', height(filtered_organized_table), height(organized_table));

%% Create the combined table with "Recording with more affected hand" column
% Initialize the new column
recording_with_more_affected = NaN(height(filtered_organized_table), 1);

for i = 1:height(filtered_organized_table)
    participant_id = filtered_organized_table.ParticipantID_num(i);
    hand_in_recording = filtered_organized_table.Hand{i};
    
    % Find the corresponding demographic data
    demo_row = summary(summary.patientsId == participant_id, :);
    
    if ~isempty(demo_row)
        most_affected_side = strtrim(demo_row.mostAffectedSide{1}); % Remove trailing spaces
        
        % Compare hand in recording with most affected side
        % Convert to consistent format for comparison
        if strcmpi(hand_in_recording, 'right') && strcmpi(most_affected_side, 'right')
            recording_with_more_affected(i) = 1;
        elseif strcmpi(hand_in_recording, 'left') && strcmpi(most_affected_side, 'left')
            recording_with_more_affected(i) = 1;
        elseif strcmpi(most_affected_side, 'na')
            recording_with_more_affected(i) = NaN; % Set to NaN when most affected side is 'na'
        else
            recording_with_more_affected(i) = 0;
        end
        
        % Debug output for first few cases
        if i <= 5
            fprintf('P%d: Recording hand=%s, Most affected=%s, More affected recording=%g\n', ...
                    participant_id, hand_in_recording, most_affected_side, recording_with_more_affected(i));
        end
    else
        fprintf('Warning: No demographic data found for participant %d\n', participant_id);
        recording_with_more_affected(i) = NaN;
    end
end

% Add the new column to the filtered table
filtered_organized_table.RecordingWithMoreAffectedHand = recording_with_more_affected;

%% Merge demographic information into the organized table
% Add demographic columns to the organized table
age_merged = NaN(height(filtered_organized_table), 1);
sex_merged = cell(height(filtered_organized_table), 1);
disease_duration_merged = NaN(height(filtered_organized_table), 1);
daily_medication_dose_merged = NaN(height(filtered_organized_table), 1);
most_affected_side_merged = cell(height(filtered_organized_table), 1);
dominant_hand_merged = cell(height(filtered_organized_table), 1);
relative_plasma_levodopa_merged = NaN(height(filtered_organized_table), 1);
plasma_levodopa_level_merged = NaN(height(filtered_organized_table), 1);
execise_frequency_merged = NaN(height(filtered_organized_table),1);
mood_merged = NaN(height(filtered_organized_table),1);
selfReport_merged = NaN(height(filtered_organized_table),1);
for i = 1:height(filtered_organized_table)
    participant_id = filtered_organized_table.ParticipantID_num(i);
    exercise_condition = filtered_organized_table.Exercise{i}; % 'preE' or 'postE'
    
    % Map exercise condition to the Exercise column in demographic data
    % preE = 0 (pre-exercise), postE = 1 (post-exercise)
    if strcmp(exercise_condition, 'preE')
        exercise_code = 0;
    elseif strcmp(exercise_condition, 'postE')
        exercise_code = 1;
    else
        exercise_code = NaN;
    end
    
    % Find the specific row that matches both participant ID and exercise condition
    matching_rows = demographic_data(demographic_data.patientsId == participant_id & ...
                                   demographic_data.Exercise == exercise_code, :);
    
    if ~isempty(matching_rows)
        % Use the first matching row 
        demo_row = matching_rows(1, :);
        
        age_merged(i) = demo_row.age_yearsOld_;
        sex_merged{i} = demo_row.sex_F_M_{1};
        disease_duration_merged(i) = demo_row.diseaseDuration;
        daily_medication_dose_merged(i) = demo_row.dailyMedicationDose;      
        most_affected_side_merged{i} = strtrim(demo_row.mostAffectedSide{1});
        dominant_hand_merged{i} = strtrim(demo_row.dominantHand{1});
        relative_plasma_levodopa_merged(i) = demo_row.relativePlasmaLevodopaLevel;
        plasma_levodopa_level_merged(i) = demo_row.plasmaLevodopaLevel_mg_;
        execise_frequency_merged(i) = demo_row.exerciseFrequency;
        mood_merged(i) = demo_row.mood;
        selfReport_merged(i) = demo_row.selfReportPdWarriorImprovements;


        % Debug output for first few cases
        if i <= 5
            fprintf('P%d %s: Relative levodopa=%.3f, Absolute levodopa=%.1f mg\n', ...
                    participant_id, exercise_condition, ...
                    demo_row.relativePlasmaLevodopaLevel, demo_row.plasmaLevodopaLevel_mg_);
        end
    else
        % No matching row found - fill with demographic summary for basic info
        demo_summary_row = demo_summary(demo_summary.patientsId == participant_id, :);
        if ~isempty(demo_summary_row)
            age_merged(i) = demo_summary_row.age_years_old_;
            sex_merged{i} = demo_summary_row.sex_F_M_{1};
            disease_duration_merged(i) = demo_summary_row.diseaseDuration;
            daily_medication_dose_merged(i) = demo_summary_row.dailyMedicationDose;
            most_affected_side_merged{i} = strtrim(demo_summary_row.mostAffectedSide{1});
            dominant_hand_merged{i} = strtrim(demo_summary_row.dominantHand{1});
        else
            sex_merged{i} = '';
            most_affected_side_merged{i} = '';
            dominant_hand_merged{i} = '';
        end
        
        fprintf('Warning: No levodopa data found for P%d %s condition\n', participant_id, exercise_condition);
        relative_plasma_levodopa_merged(i) = NaN;
        plasma_levodopa_level_merged(i) = NaN;
    end
end

% Create the final combined table
combined_table = [filtered_organized_table, ...
                  table(age_merged, sex_merged, disease_duration_merged, ...
                        daily_medication_dose_merged,  relative_plasma_levodopa_merged, ...
                        most_affected_side_merged, dominant_hand_merged, execise_frequency_merged, mood_merged, selfReport_merged,...
                        'VariableNames', {'Age', 'Gender', 'DiseaseDuration', ...
                                        'DailyMedicationDose', 'relativePlasmaLevodopa', ...
                                        'MostAffectedSide', 'DominantHand','exercise_frequency','mood','self_report_improvements'})];

% Display the combined table
fprintf('\n=== COMBINED DEMOGRAPHIC AND BEHAVIORAL TABLE ===\n');
fprintf('Combined table created with %d rows and %d columns\n', height(combined_table), width(combined_table));
disp(combined_table(1:5,:));
save('combined_table');
writetable(combined_table, 'combined_demographic_behavioral.csv');
%%
% Load demographic data if not already loaded
% Load demographic data if not already loaded
if ~exist('demographic_data', 'var')
    demographic_data = readtable('fulltable1.xlsx');
    fprintf('Demographic data loaded with %d rows and %d columns\n', height(demographic_data), width(demographic_data));
    
    % Clean column names
    demographic_data.Properties.VariableNames = matlab.lang.makeValidName(demographic_data.Properties.VariableNames);
end

% Get unique participants from demographic data
unique_demo_participants = unique(demographic_data.patientsId);
fprintf('Unique participants in demographic data: %s\n', mat2str(unique_demo_participants));

%% Create demographic summary table (one row per participant for basic info)
demo_summary = table();
for i = 1:length(unique_demo_participants)
    participant_id = unique_demo_participants(i);
    participant_rows = demographic_data(demographic_data.patientsId == participant_id, :);
    
    % Take the first row for demographic info (should be same across rows for same participant)
    demo_row = participant_rows(1, :);
    demo_summary = [demo_summary; demo_row];
end

fprintf('Demo summary table created with %d participants\n', height(demo_summary));



%% Load the individual tables if they don't exist
if ~exist('travelTable', 'var')
    if exist('all_travel_times.csv', 'file')
        travelTable = readtable('all_travel_times.csv');
        fprintf('Loaded travelTable from CSV file\n');
    else
        error('travelTable not found. Please run the first analysis script first.');
    end
end

if ~exist('dwellTable', 'var')
    if exist('all_dwell_times.csv', 'file')
        dwellTable = readtable('all_dwell_times.csv');
        fprintf('Loaded dwellTable from CSV file\n');
    else
        error('dwellTable not found. Please run the first analysis script first.');
    end
end

if ~exist('velocitytable', 'var')
    if exist('all_velocity.csv', 'file')
        velocitytable = readtable('all_velocity.csv');
        fprintf('Loaded velocitytable from CSV file\n');
    else
        error('velocitytable not found. Please run the first analysis script first.');
    end
end

%% Merge demographic data with each table
fprintf('\n=== MERGING TRAVEL TABLE ===\n');
combined_travel_table = merge_demographic_with_table(travelTable, 'IDrep', 'exerep', demographic_data, demo_summary);
fprintf('Combined travel table created with %d rows and %d columns\n', height(combined_travel_table), width(combined_travel_table));

fprintf('\n=== MERGING DWELL TABLE ===\n');
combined_dwell_table = merge_demographic_with_table(dwellTable, 'IDrep1', 'exerep1', demographic_data, demo_summary);
fprintf('Combined dwell table created with %d rows and %d columns\n', height(combined_dwell_table), width(combined_dwell_table));

fprintf('\n=== MERGING VELOCITY TABLE ===\n');
combined_velocity_table = merge_demographic_with_table(velocitytable, 'IDrep2', 'exerep2', demographic_data, demo_summary);
fprintf('Combined velocity table created with %d rows and %d columns\n', height(combined_velocity_table), width(combined_velocity_table));

%% Display sample data from each merged table
fprintf('\n=== SAMPLE DATA FROM MERGED TABLES ===\n');

fprintf('\nCombined Travel Table (first 5 rows):\n');
disp(combined_travel_table(1:5,:));

fprintf('\nCombined Dwell Table (first 5 rows):\n');
disp(combined_dwell_table(1:5,:));

fprintf('\nCombined Velocity Table (first 5 rows):\n');
disp(combined_velocity_table(1:5,:));

%% Save the merged tables
save('combined_travel_table.mat', 'combined_travel_table');
save('combined_dwell_table.mat', 'combined_dwell_table');
save('combined_velocity_table.mat', 'combined_velocity_table');

writetable(combined_travel_table, 'combined_travel_demographic.csv');
writetable(combined_dwell_table, 'combined_dwell_demographic.csv');
writetable(combined_velocity_table, 'combined_velocity_demographic.csv');

%% Summary statistics
unique_participants = unique(combined_table.ParticipantID);
behavior_measures = {'velocity', 'travelTimes', 'dwellTimes', 'nrErrors','nrTaps', 'dysmetriaScore','relativePlasmaLevodopa'};  % adjust based on your real column names
behavior_var_map = struct( ...
    'velocity',      'Velocity', ...
    'travelTimes',   'TravelTime', ...
    'dwellTimes',    'DwellTime', ...
    'nrErrors',      'NumErrors', ...
    'nrTaps',        'NumTaps', ...
    'dysmetriaScore', 'DS', ...
    'relativePlasmaLevodopa','rpl' ...
);



% Initialize results structure
pre_post_results_by_affection = struct();

% Group codes: 1 = More Affected Hand, 0 = Less Affected Hand
group_codes = [1, 0];
group_labels = {'more_affected', 'less_affected'};

for p = 1:length(unique_participants)
    participant = unique_participants{p};
    
    for g = 1:length(group_codes)
        code = group_codes(g);
        hand_type = group_labels{g};
        
        % Filter data for current participant and hand affection level
        participant_hand_data = combined_table(strcmp(combined_table.ParticipantID, participant) & ...
                                             combined_table.RecordingWithMoreAffectedHand == code, :);
        
        if isempty(participant_hand_data)
            continue;
        end
        
        % Separate pre and post exercise data
        pre_data = participant_hand_data(strcmp(participant_hand_data.Exercise, 'preE'), :);
        post_data = participant_hand_data(strcmp(participant_hand_data.Exercise, 'postE'), :);
        
        % Create field name for results
        result_field = ['P', participant, '_', hand_type];
        pre_post_results_by_affection.(result_field) = struct();
        
        % Calculate averages for each behavior measure
        for b = 1:length(behavior_measures)
            behavior = behavior_measures{b};
            
            % Pre-exercise average (across recordings)
            if ~isempty(pre_data) && any(~isnan(pre_data.(behavior)))
                pre_values = pre_data.(behavior);
                pre_values = pre_values(~isnan(pre_values)); % Remove NaN values
                pre_post_results_by_affection.(result_field).(['pre_', behavior]) = mean(pre_values);
                pre_post_results_by_affection.(result_field).(['pre_', behavior, '_std']) = std(pre_values);
                pre_post_results_by_affection.(result_field).(['pre_', behavior, '_n']) = length(pre_values);
            else
                pre_post_results_by_affection.(result_field).(['pre_', behavior]) = NaN;
                pre_post_results_by_affection.(result_field).(['pre_', behavior, '_std']) = NaN;
                pre_post_results_by_affection.(result_field).(['pre_', behavior, '_n']) = 0;
            end
            
            % Post-exercise average (across recordings)
            if ~isempty(post_data) && any(~isnan(post_data.(behavior)))
                post_values = post_data.(behavior);
                post_values = post_values(~isnan(post_values)); % Remove NaN values
                pre_post_results_by_affection.(result_field).(['post_', behavior]) = mean(post_values);
                pre_post_results_by_affection.(result_field).(['post_', behavior, '_std']) = std(post_values);
                pre_post_results_by_affection.(result_field).(['post_', behavior, '_n']) = length(post_values);
            else
                pre_post_results_by_affection.(result_field).(['post_', behavior]) = NaN;
                pre_post_results_by_affection.(result_field).(['post_', behavior, '_std']) = NaN;
                pre_post_results_by_affection.(result_field).(['post_', behavior, '_n']) = 0;
            end
            
            % Calculate change (post - pre) and percent change
            if ~isnan(pre_post_results_by_affection.(result_field).(['pre_', behavior])) && ...
               ~isnan(pre_post_results_by_affection.(result_field).(['post_', behavior]))
                
                pre_val = pre_post_results_by_affection.(result_field).(['pre_', behavior]);
                post_val = pre_post_results_by_affection.(result_field).(['post_', behavior]);
                
                pre_post_results_by_affection.(result_field).(['change_', behavior]) = post_val - pre_val;
                
                if pre_val ~= 0
                    pre_post_results_by_affection.(result_field).(['percent_change_', behavior]) = ...
                        ((post_val - pre_val) / abs(pre_val)) * 100;
                else
                    pre_post_results_by_affection.(result_field).(['percent_change_', behavior]) = NaN;
                end
            else
                pre_post_results_by_affection.(result_field).(['change_', behavior]) = NaN;
                pre_post_results_by_affection.(result_field).(['percent_change_', behavior]) = NaN;
            end
        end
    end
end

% Save pre/post results
save('pre_post_results_by_affection.mat', 'pre_post_results_by_affection');

%% Create summary table for pre/post analysis by hand affection
summary_participant_ids = {};
summary_hand_types = {};
pre_Velocity = []; post_Velocity = []; change_Velocity = []; percentChange_Velocity = [];
pre_TravelTime = []; post_TravelTime = []; change_TravelTime = []; percentChange_TravelTime =[];
pre_DwellTime = []; post_DwellTime = []; change_DwellTime = []; percentChange_DwellTime = [];
pre_NumErrors = []; post_NumErrors = []; change_NumErrors = []; percentChange_NumErrors =[];
pre_NumTaps = []; post_NumTaps = []; change_NumTaps = [];  percentChange_NumTaps =[];
pre_DS = []; post_DS = []; change_DS = [];    percentChange_DS =[];
pre_rpl  = []; post_rpl = []; change_rpl = []; percentChange_rpl =[];
 
result_fields = fieldnames(pre_post_results_by_affection);
for i = 1:length(result_fields)
    field_name = result_fields{i};
    
    % Parse participant and hand type from field name
    parts = strsplit(field_name, '_');
    participant = parts{1}(2:end); % Remove 'P' prefix
    hand_type = strjoin(parts(2:end), '_'); % Rejoin in case of multi-word labels
    
    summary_participant_ids{end+1} = participant;
    summary_hand_types{end+1} = hand_type;
    
    % Extract pre/post values for each measure
    pre_Velocity(end+1) = pre_post_results_by_affection.(field_name).pre_velocity;
    post_Velocity(end+1) = pre_post_results_by_affection.(field_name).post_velocity;
    change_Velocity(end+1) = pre_post_results_by_affection.(field_name).change_velocity;
    
    pre_TravelTime(end+1) = pre_post_results_by_affection.(field_name).pre_travelTimes;
    post_TravelTime(end+1) = pre_post_results_by_affection.(field_name).post_travelTimes;
    change_TravelTime(end+1) = pre_post_results_by_affection.(field_name).change_travelTimes;
    
    pre_DwellTime(end+1) = pre_post_results_by_affection.(field_name).pre_dwellTimes;
    post_DwellTime(end+1) = pre_post_results_by_affection.(field_name).post_dwellTimes;
    change_DwellTime(end+1) = pre_post_results_by_affection.(field_name).change_dwellTimes;
    
    pre_NumErrors(end+1) = pre_post_results_by_affection.(field_name).pre_nrErrors;
    post_NumErrors(end+1) = pre_post_results_by_affection.(field_name).post_nrErrors;
    change_NumErrors(end+1) = pre_post_results_by_affection.(field_name).change_nrErrors;

    pre_NumTaps(end+1) = pre_post_results_by_affection.(field_name).pre_nrTaps;
    post_NumTaps(end+1) = pre_post_results_by_affection.(field_name).post_nrTaps;
    change_NumTaps(end+1) = pre_post_results_by_affection.(field_name).change_nrTaps;

    pre_DS(end+1) = pre_post_results_by_affection.(field_name).pre_dysmetriaScore;
    post_DS(end+1) = pre_post_results_by_affection.(field_name).post_dysmetriaScore;
    change_DS(end+1) = pre_post_results_by_affection.(field_name).change_dysmetriaScore;
    
    pre_rpl(end+1) = pre_post_results_by_affection.(field_name).pre_relativePlasmaLevodopa;
    post_rpl(end+1) = pre_post_results_by_affection.(field_name).post_relativePlasmaLevodopa;
    change_rpl(end+1) = pre_post_results_by_affection.(field_name).change_relativePlasmaLevodopa;

    calc_percent_change = @(pre, post) ...
        (isnan(pre) || pre == 0 || isnan(post)) * NaN + ...
        (~isnan(pre) && pre ~= 0 && ~isnan(post)) * ((post - pre) / abs(pre) * 100);

    percentChange_Velocity(end+1)      = calc_percent_change(pre_Velocity(end), post_Velocity(end));
    percentChange_TravelTime(end+1)    = calc_percent_change(pre_TravelTime(end), post_TravelTime(end));
    percentChange_DwellTime(end+1)     = calc_percent_change(pre_DwellTime(end), post_DwellTime(end));
    percentChange_NumErrors(end+1)     = calc_percent_change(pre_NumErrors(end), post_NumErrors(end));
    percentChange_NumTaps(end+1)       = calc_percent_change(pre_NumTaps(end), post_NumTaps(end));
    percentChange_DS(end+1)            = calc_percent_change(pre_DS(end), post_DS(end));
    percentChange_rpl(end+1)           = calc_percent_change(pre_rpl(end), post_rpl(end));
end

% Create summary table
summary_table_by_affection = table(summary_participant_ids', summary_hand_types', ...
                     pre_Velocity', post_Velocity', change_Velocity',  percentChange_Velocity',...
                     pre_TravelTime', post_TravelTime', change_TravelTime',  percentChange_TravelTime', ...
                     pre_DwellTime', post_DwellTime', change_DwellTime',   percentChange_DwellTime',...
                     pre_NumErrors', post_NumErrors', change_NumErrors', percentChange_NumErrors',...
                     pre_NumTaps', post_NumTaps', change_NumTaps', percentChange_NumTaps',...
                     pre_DS', post_DS', change_DS', percentChange_DS',...
                     pre_rpl', post_rpl', change_rpl',percentChange_rpl', ...
                     'VariableNames', {'ParticipantID', 'HandType', ...
                                     'Pre_Velocity', 'Post_Velocity', 'Change_Velocity','percent_change_Velocity', ...
                                     'Pre_TravelTime', 'Post_TravelTime', 'Change_TravelTime', 'percent_change_TravelTime', ...
                                     'Pre_DwellTime', 'Post_DwellTime', 'Change_DwellTime', 'percent_change_DwellTime',...
                                     'Pre_NumErrors', 'Post_NumErrors', 'Change_NumErrors', 'percent_change_NumErrors',...
                                     'Pre_NumTaps', 'Post_NumTaps', 'Change_NumTaps' , 'percent_change_NumTaps',...
                                     'Pre_DS', 'Post_DS', 'Change_DS', 'percent_change_DS',...
                                     'Pre_rpl', 'Post_rpl', 'Change_rpl','percent_change_rpl'});

% Display summary
fprintf('\nPre/Post Exercise Summary Table by Hand Affection:\n');
disp(summary_table_by_affection);

% Save summary table
save('pre_post_summary_by_affection.mat', 'summary_table_by_affection');
writetable(summary_table_by_affection, 'pre_post_summary_by_affection.csv');

%% Display overall statistics by hand affection (WITH OUTLIER REMOVAL)
fprintf('\n=== STATISTICS SEPARATED BY AFFECTED HAND (WITH OUTLIER REMOVAL) ===\n');
group_labels = {'More Affected Hand', 'Less Affected Hand'};
group_codes  = [1, 0];

% Store outlier information for reporting
outlier_report = struct();

for g = 1:2
    code = group_codes(g);
    group_name = group_labels{g};
    hand_type_filter = group_labels{g};
    if g == 1
        hand_type_filter = 'more_affected';
    else
        hand_type_filter = 'less_affected';
    end

    % Filter summary table by hand type
    group_data = summary_table_by_affection(strcmp(summary_table_by_affection.HandType, hand_type_filter), :);

    fprintf('\n--- %s ---\n', group_name);
    outlier_report.(hand_type_filter) = struct();

    for b = 1:length(behavior_measures)
        behavior = behavior_measures{b};
        behavior_label = behavior_var_map.(behavior);  % Correct casing

        pre_col    = ['Pre_', behavior_label];
        post_col   = ['Post_', behavior_label];
        change_col = ['Change_', behavior_label];

        pre_values = group_data.(pre_col);
        post_values = group_data.(post_col);
        change_values = group_data.(change_col);

        % Remove outliers from each dataset
        [pre_values_clean, pre_outliers] = remove_outliers_mad(pre_values);
        [post_values_clean, post_outliers] = remove_outliers_mad(post_values);
        [change_values_clean, change_outliers] = remove_outliers_mad(change_values);

        % Store outlier information
        outlier_report.(hand_type_filter).(behavior) = struct();
        outlier_report.(hand_type_filter).(behavior).pre_outliers_removed = sum(pre_outliers);
        outlier_report.(hand_type_filter).(behavior).post_outliers_removed = sum(post_outliers);
        outlier_report.(hand_type_filter).(behavior).change_outliers_removed = sum(change_outliers);
        outlier_report.(hand_type_filter).(behavior).original_n = length(pre_values);

        fprintf('%s:\n', behavior);
        fprintf('  Pre:  Mean=%.3f, SD=%.3f, N=%d (removed %d outliers)\n', ...
                mean(pre_values_clean), std(pre_values_clean), length(pre_values_clean), sum(pre_outliers));
        fprintf('  Post: Mean=%.3f, SD=%.3f, N=%d (removed %d outliers)\n', ...
                mean(post_values_clean), std(post_values_clean), length(post_values_clean), sum(post_outliers));
        fprintf('  Change: Mean=%.3f, SD=%.3f, N=%d (removed %d outliers)\n', ...
                mean(change_values_clean), std(change_values_clean), length(change_values_clean), sum(change_outliers));

        % For paired tests, we need to ensure we're comparing the same participants
        % Get indices of participants that have both pre and post data (after outlier removal)
        pre_valid = ~isnan(pre_values) & ~pre_outliers;
        post_valid = ~isnan(post_values) & ~post_outliers;
        both_valid = pre_valid & post_valid;
        
        if sum(both_valid) > 1
            pre_paired = pre_values(both_valid);
            post_paired = post_values(both_valid);
            differences = post_paired - pre_paired;
            
            % Remove outliers from differences as well
            [differences_clean, diff_outliers] = remove_outliers_mad(differences);
            
            % Get the corresponding pre and post values for the clean differences
            clean_indices = both_valid;
            clean_indices(both_valid) = ~diff_outliers; % Update indices based on difference outliers
            
            pre_final = pre_values(clean_indices);
            post_final = post_values(clean_indices);
            
            fprintf('  Paired analysis: N=%d (removed %d additional outliers from differences)\n', ...
                    length(differences_clean), sum(diff_outliers));

            % Visual + statistical check
            is_ttest_ok = true;

            if length(differences_clean) < 10
                is_ttest_ok = false;  % too small for t-test reliability
            else
                % Normality test (less strict)
                if kstest(differences_clean)  % p < 0.05 means NOT normal
                    if abs(skewness(differences_clean)) > 1  % strong skew = non-normal
                        is_ttest_ok = false;
                    end
                end
            end

            if is_ttest_ok
                % Paired t-test
                [~, p_val, ci, stats] = ttest(pre_final, post_final);
                test_type = 'Paired t-test';
                t_or_z_val = stats.tstat;
                df = stats.df;
                ci_str = sprintf('95%% CI: [%.3f, %.3f]', ci(1), ci(2));
            else
                % Wilcoxon signed-rank test
                [p_val, ~, stats] = signrank(pre_final, post_final);
                test_type = 'Wilcoxon signed-rank test';
                t_or_z_val = stats.signedrank;
                df = NaN;
                ci_str = 'CI not available for non-parametric test';
            end

            sig_label = ternary(p_val < 0.05, '(significant)', '(not significant)');
            mean_diff = mean(differences_clean);

            if isnan(df)
                fprintf('  %s: signed-rank = %d, p = %.4f %s\n', test_type, t_or_z_val, p_val, sig_label);
            else
                fprintf('  %s: t(%d) = %.3f, p = %.4f %s\n', test_type, df, t_or_z_val, p_val, sig_label);
            end
            fprintf('    Mean difference (Post - Pre): %.3f\n', mean_diff);
            fprintf('    %s\n', ci_str);
        end
        fprintf('\n');
    end
end

%% Compare between more and less affected hands (WITH OUTLIER REMOVAL)
fprintf('\n=== COMPARISON BETWEEN MORE AND LESS AFFECTED HANDS (WITH OUTLIER REMOVAL) ===\n');

more_affected_data = summary_table_by_affection(strcmp(summary_table_by_affection.HandType, 'more_affected'), :);
less_affected_data = summary_table_by_affection(strcmp(summary_table_by_affection.HandType, 'less_affected'), :);

for b = 1:length(behavior_measures)
    behavior = behavior_measures{b};
    
    fprintf('\n--- %s ---\n', behavior);
    behavior_label = behavior_var_map.(behavior);
    pre_col  = ['Pre_', behavior_label];
    post_col = ['Post_', behavior_label];
    change_col = ['Change_', behavior_label];
    
    % Pre-exercise comparison
    more_affected_pre = more_affected_data.(pre_col);
    less_affected_pre = less_affected_data.(pre_col);
    
    [more_affected_pre_clean, more_pre_outliers] = remove_outliers_mad(more_affected_pre);
    [less_affected_pre_clean, less_pre_outliers] = remove_outliers_mad(less_affected_pre);
    
    if length(more_affected_pre_clean) > 1 && length(less_affected_pre_clean) > 1
        [h_pre, p_pre] = ttest2(more_affected_pre_clean, less_affected_pre_clean);
        fprintf('Pre-exercise: More affected (M=%.3f, SD=%.3f, N=%d, %d outliers removed) vs Less affected (M=%.3f, SD=%.3f, N=%d, %d outliers removed), p=%.4f\n', ...
                mean(more_affected_pre_clean), std(more_affected_pre_clean), length(more_affected_pre_clean), sum(more_pre_outliers), ...
                mean(less_affected_pre_clean), std(less_affected_pre_clean), length(less_affected_pre_clean), sum(less_pre_outliers), p_pre);
    end
    
    % Post-exercise comparison
    more_affected_post = more_affected_data.(post_col);
    less_affected_post = less_affected_data.(post_col);
    
    [more_affected_post_clean, more_post_outliers] = remove_outliers_mad(more_affected_post);
    [less_affected_post_clean, less_post_outliers] = remove_outliers_mad(less_affected_post);
    
    if length(more_affected_post_clean) > 1 && length(less_affected_post_clean) > 1
        [h_post, p_post] = ttest2(more_affected_post_clean, less_affected_post_clean);
        fprintf('Post-exercise: More affected (M=%.3f, SD=%.3f, N=%d, %d outliers removed) vs Less affected (M=%.3f, SD=%.3f, N=%d, %d outliers removed), p=%.4f\n', ...
                mean(more_affected_post_clean), std(more_affected_post_clean), length(more_affected_post_clean), sum(more_post_outliers), ...
                mean(less_affected_post_clean), std(less_affected_post_clean), length(less_affected_post_clean), sum(less_post_outliers), p_post);
    end
    
    % Change comparison
    more_affected_change = more_affected_data.(change_col);
    less_affected_change = less_affected_data.(change_col);
    
    [more_affected_change_clean, more_change_outliers] = remove_outliers_mad(more_affected_change);
    [less_affected_change_clean, less_change_outliers] = remove_outliers_mad(less_affected_change);
    
    if length(more_affected_change_clean) > 1 && length(less_affected_change_clean) > 1
        [h_change, p_change] = ttest2(more_affected_change_clean, less_affected_change_clean);
        fprintf('Change: More affected (M=%.3f, SD=%.3f, N=%d, %d outliers removed) vs Less affected (M=%.3f, SD=%.3f, N=%d, %d outliers removed), p=%.4f\n', ...
                mean(more_affected_change_clean), std(more_affected_change_clean), length(more_affected_change_clean), sum(more_change_outliers), ...
                mean(less_affected_change_clean), std(less_affected_change_clean), length(less_affected_change_clean), sum(less_change_outliers), p_change);
    end
    
    % Percent change comparison
    pct_col = ['percent_change_',  behavior];
    
    % Use raw structure again to retrieve percent changes
    more_affected_pct = [];
    less_affected_pct = [];
    for i = 1:length(result_fields)
        fname = result_fields{i};
        if contains(fname, 'more_affected') && isfield(pre_post_results_by_affection.(fname), pct_col)
            more_affected_pct(end+1) = pre_post_results_by_affection.(fname).(pct_col);
        elseif contains(fname, 'less_affected') && isfield(pre_post_results_by_affection.(fname), pct_col)
            less_affected_pct(end+1) = pre_post_results_by_affection.(fname).(pct_col);
        end
    end

    [more_affected_pct_clean, more_pct_outliers] = remove_outliers_mad(more_affected_pct);
    [less_affected_pct_clean, less_pct_outliers] = remove_outliers_mad(less_affected_pct);

    if length(more_affected_pct_clean) > 1 && length(less_affected_pct_clean) > 1
        [~, p_pct] = ttest2(more_affected_pct_clean, less_affected_pct_clean);
        fprintf('Percent Change: More affected (M=%.2f%%, SD=%.2f, N=%d, %d outliers removed) vs Less affected (M=%.2f%%, SD=%.2f, N=%d, %d outliers removed), p=%.4f\n', ...
                mean(more_affected_pct_clean), std(more_affected_pct_clean), length(more_affected_pct_clean), sum(more_pct_outliers), ...
                mean(less_affected_pct_clean), std(less_affected_pct_clean), length(less_affected_pct_clean), sum(less_pct_outliers), p_pct);
    end
    
end

% Save outlier report
save('outlier_report.mat', 'outlier_report');

% Display outlier summary
fprintf('\n=== OUTLIER REMOVAL SUMMARY ===\n');
for g = 1:2
    if g == 1
        hand_type_filter = 'more_affected';
        group_name = 'More Affected Hand';
    else
        hand_type_filter = 'less_affected';
        group_name = 'Less Affected Hand';
    end
    
    fprintf('\n--- %s ---\n', group_name);
    for b = 1:length(behavior_measures)
        behavior = behavior_measures{b};
        if isfield(outlier_report, hand_type_filter) && isfield(outlier_report.(hand_type_filter), behavior)
            orig_n = outlier_report.(hand_type_filter).(behavior).original_n;
            pre_removed = outlier_report.(hand_type_filter).(behavior).pre_outliers_removed;
            post_removed = outlier_report.(hand_type_filter).(behavior).post_outliers_removed;
            change_removed = outlier_report.(hand_type_filter).(behavior).change_outliers_removed;
            
            fprintf('%s: Original N=%d, Pre outliers removed=%d, Post outliers removed=%d, Change outliers removed=%d\n', ...
                    behavior, orig_n, pre_removed, post_removed, change_removed);
        end
    end
end

fprintf('\nAnalysis by hand affection complete! Files saved:\n');
fprintf('- pre_post_results_by_affection.mat\n');
fprintf('- pre_post_summary_by_affection.mat & .csv\n');
fprintf('- outlier_report.mat\n');
fprintf('\nNote: Statistics were computed after removing outliers using the MAD method (3*MAD threshold).\n');
%%
function merged_table = merge_demographic_with_table(data_table, participant_col, exercise_col, demographic_data, demo_summary)
    
    % Convert participant IDs to numeric format for comparison
    if iscell(data_table.(participant_col))
        % Handle string format like 'P01', 'P02', etc.
        participant_ids_cell = data_table.(participant_col);
        participant_ids_num = zeros(length(participant_ids_cell), 1);
        for j = 1:length(participant_ids_cell)
            pid_str = participant_ids_cell{j};
            if startsWith(pid_str, 'P')
                % Remove 'P' and convert to number
                participant_ids_num(j) = str2double(pid_str(2:end));
            else
                participant_ids_num(j) = str2double(pid_str);
            end
        end
    else
        participant_ids_num = data_table.(participant_col);
    end
    
    % Filter to include only participants in demographic data
    unique_demo_participants = unique(demographic_data.patientsId);
    participants_in_both = intersect(unique(participant_ids_num), unique_demo_participants);
    fprintf('Participants present in both tables: %s\n', mat2str(participants_in_both));
    
    % Filter the data table
    filtered_table = data_table(ismember(participant_ids_num, participants_in_both), :);
    fprintf('Filtered table: %d rows (from %d original rows)\n', height(filtered_table), height(data_table));
    
    % Initialize demographic columns
    n_rows = height(filtered_table);
    age_merged = NaN(n_rows, 1);
    sex_merged = cell(n_rows, 1);
    disease_duration_merged = NaN(n_rows, 1);
    daily_medication_dose_merged = NaN(n_rows, 1);
    most_affected_side_merged = cell(n_rows, 1);
    dominant_hand_merged = cell(n_rows, 1);
    relative_plasma_levodopa_merged = NaN(n_rows, 1);
    plasma_levodopa_level_merged = NaN(n_rows, 1);
    recording_with_more_affected = NaN(n_rows, 1);
    
    % Convert participant IDs in filtered table to numeric
    if iscell(filtered_table.(participant_col))
        % Handle string format like 'P01', 'P02', etc.
        participant_ids_cell = filtered_table.(participant_col);
        filtered_participant_ids = zeros(length(participant_ids_cell), 1);
        for j = 1:length(participant_ids_cell)
            pid_str = participant_ids_cell{j};
            if startsWith(pid_str, 'P')
                % Remove 'P' and convert to number
                filtered_participant_ids(j) = str2double(pid_str(2:end));
            else
                filtered_participant_ids(j) = str2double(pid_str);
            end
        end
    else
        filtered_participant_ids = filtered_table.(participant_col);
    end
    
    % Merge demographic information
    for i = 1:n_rows
        participant_id = filtered_participant_ids(i);
        exercise_condition = filtered_table.(exercise_col){i}; % 'preE' or 'postE'
        
        % Map exercise condition to numeric code
        if strcmp(exercise_condition, 'preE')
            exercise_code = 0;
        elseif strcmp(exercise_condition, 'postE')
            exercise_code = 1;
        else
            exercise_code = NaN;
        end
        
        % Find matching row in demographic data (specific exercise condition)
        matching_rows = demographic_data(demographic_data.patientsId == participant_id & ...
                                       demographic_data.Exercise == exercise_code, :);
        
        if ~isempty(matching_rows)
            % Use the first matching row
            demo_row = matching_rows(1, :);
            
            age_merged(i) = demo_row.age_yearsOld_;
            sex_merged{i} = demo_row.sex_F_M_{1};
            disease_duration_merged(i) = demo_row.diseaseDuration;
            daily_medication_dose_merged(i) = demo_row.dailyMedicationDose;
            most_affected_side_merged{i} = strtrim(demo_row.mostAffectedSide{1});
            dominant_hand_merged{i} = strtrim(demo_row.dominantHand{1});
            relative_plasma_levodopa_merged(i) = demo_row.relativePlasmaLevodopaLevel;
            plasma_levodopa_level_merged(i) = demo_row.plasmaLevodopaLevel_mg_;
        else
            % No specific exercise condition found - use summary for basic info
            demo_summary_row = demo_summary(demo_summary.patientsId == participant_id, :);
            if ~isempty(demo_summary_row)
                age_merged(i) = demo_summary_row.age_yearsOld_;
                sex_merged{i} = demo_summary_row.sex_F_M_{1};
                disease_duration_merged(i) = demo_summary_row.diseaseDuration;
                daily_medication_dose_merged(i) = demo_summary_row.dailyMedicationDose;
                most_affected_side_merged{i} = strtrim(demo_summary_row.mostAffectedSide{1});
                dominant_hand_merged{i} = strtrim(demo_summary_row.dominantHand{1});
            else
                sex_merged{i} = '';
                most_affected_side_merged{i} = '';
                dominant_hand_merged{i} = '';
            end
            
            fprintf('Warning: No levodopa data found for P%d %s condition\n', participant_id, exercise_condition);
            relative_plasma_levodopa_merged(i) = NaN;
            plasma_levodopa_level_merged(i) = NaN;
        end
        
        % Calculate "Recording with more affected hand" column
        hand_in_recording = '';
        if any(strcmp(filtered_table.Properties.VariableNames, 'handrep'))
            hand_in_recording = filtered_table.handrep{i};
        elseif any(strcmp(filtered_table.Properties.VariableNames, 'handrep1'))
            hand_in_recording = filtered_table.handrep1{i};
        elseif any(strcmp(filtered_table.Properties.VariableNames, 'handrep2'))
            hand_in_recording = filtered_table.handrep2{i};
        end
        
        if ~isempty(most_affected_side_merged{i}) && ~isempty(hand_in_recording)
            most_affected_side = most_affected_side_merged{i};
            
            if strcmpi(hand_in_recording, 'right') && strcmpi(most_affected_side, 'right')
                recording_with_more_affected(i) = 1;
            elseif strcmpi(hand_in_recording, 'left') && strcmpi(most_affected_side, 'left')
                recording_with_more_affected(i) = 1;
            elseif strcmpi(most_affected_side, 'na')
                recording_with_more_affected(i) = NaN;
            else
                recording_with_more_affected(i) = 0;
            end
        else
            recording_with_more_affected(i) = NaN;
        end
    end
    
    % Create the merged table
    demographic_columns = table(age_merged, sex_merged, disease_duration_merged, ...
                               daily_medication_dose_merged, relative_plasma_levodopa_merged, ...
                               plasma_levodopa_level_merged, most_affected_side_merged, ...
                               dominant_hand_merged, recording_with_more_affected, ...
                               'VariableNames', {'Age', 'Gender', 'DiseaseDuration', ...
                                               'DailyMedicationDose', 'RelativePlasmaLevodopa', ...
                                               'PlasmaLevodopaLevel', 'MostAffectedSide', ...
                                               'DominantHand', 'RecordingWithMoreAffectedHand'});
    
    merged_table = [filtered_table, demographic_columns];
end
function result = ternary(cond, valTrue, valFalse)
    if cond
        result = valTrue;
    else
        result = valFalse;
    end
end
% Function to remove outliers using MAD method (Median Absolute Deviation)
function [clean_data, outlier_indices] = remove_outliers_mad(data, multiplier)
    if nargin < 2
        multiplier = 3; % Standard MAD multiplier (3*MAD is common threshold)
    end
    
    % Remove NaN values first
    valid_indices = ~isnan(data);
    valid_data = data(valid_indices);
    
    if length(valid_data) < 3 % Need at least 3 points for meaningful MAD calculation
        clean_data = data;
        outlier_indices = false(size(data));
        return;
    end
    
    % Calculate median and MAD
    data_median = median(valid_data);
    absolute_deviations = abs(valid_data - data_median);
    mad_value = median(absolute_deviations);
    
    % Handle case where MAD is 0 (all values are the same)
    if mad_value == 0
        clean_data = data;
        outlier_indices = false(size(data));
        return;
    end
    
    % Scale MAD to be consistent with standard deviation for normal distributions
    % This factor (1.4826) makes MAD comparable to standard deviation
    scaled_mad = mad_value * 1.4826;
    
    % Define outlier bounds using scaled MAD
    threshold = multiplier * scaled_mad;
    
    % Identify outliers in the original data
    outlier_indices = (abs(data - data_median) > threshold) | isnan(data);
    clean_data = data(~outlier_indices);
end