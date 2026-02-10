%% Load and process the demographic data
% Load the Excel file
demographic_data = readtable('demographics.xlsx');

% Display the structure
fprintf('Demographic data loaded with %d rows and %d columns\n', height(demographic_data), width(demographic_data));
fprintf('Column names: %s\n', strjoin(demographic_data.Properties.VariableNames, ', '));

% Clean column names (remove spaces and special characters)
demographic_data.Properties.VariableNames = matlab.lang.makeValidName(demographic_data.Properties.VariableNames);



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
if ~exist('demographic_data', 'var')
    demographic_data = readtable('demographics.xlsx');
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


%% Save the merged tables
save('combined_travel_table.mat', 'combined_travel_table');
save('combined_dwell_table.mat', 'combined_dwell_table');
save('combined_velocity_table.mat', 'combined_velocity_table');

writetable(combined_travel_table, 'combined_travel_demographic.csv');
writetable(combined_dwell_table, 'combined_dwell_demographic.csv');
writetable(combined_velocity_table, 'combined_velocity_demographic.csv');


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