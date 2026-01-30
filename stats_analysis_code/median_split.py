import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# =====================
# Load dataset
# =====================
df = pd.read_csv('combined_demographic_behavioral.csv')
df = df.dropna(subset=['RecordingWithMoreAffectedHand'])

# =====================
# Settings
# =====================
factors_to_split = ['Age', 'DiseaseDuration', 'relativePlasmaLevodopa']
behavior_vars = ['velocity', 'dwellTimes']

# Detect recording column
possible_rec_cols = ['RecordingNumber','Recording','rec','Rec','recording','RecordingIndex','RecordingID']
REC_COL = ''
for c in possible_rec_cols:
    if c in df.columns:
        REC_COL = c
        break
if REC_COL == '':
    raise ValueError("Couldn't find a recording column. Please add one of these columns with values 1/2: " + ", ".join(possible_rec_cols))

# =====================
# Median split setup
# =====================
def create_median_split(df, column):
    median_val = df[column].median()
    df[f'{column}_group'] = df[column].apply(lambda x: 'Low' if x <= median_val else 'High')
    return df, median_val

median_values = {}
for factor in factors_to_split:
    if factor in df.columns:
        df, median_val = create_median_split(df, factor)
        median_values[factor] = median_val

# =====================
# Compute group/recording-level means and SDs
# =====================
def calculate_means(data, description, factor_name=None, group=None, rec=None):
    d = data.copy()
    result = {'Description': description}
    if factor_name and group:
        result['Factor'] = factor_name
        result['Group'] = group
    if rec is not None:
        d = d[d[REC_COL] == rec]
        result['Recording'] = rec

    for rep, rep_name in [(1, 'MoreAffected'), (0, 'LessAffected')]:
        rep_data = d[d['RecordingWithMoreAffectedHand'] == rep]
        for condition in ['preE', 'postE']:
            cond_data = rep_data[rep_data['Exercise'] == condition]
            for behavior in behavior_vars:
                col_base = f'{rep_name}_{condition}_{behavior}'
                if len(cond_data) > 0 and behavior in cond_data.columns:
                    vals = cond_data[behavior].dropna()
                    result[f'{col_base}_mean'] = vals.mean()
                    result[f'{col_base}_std'] = vals.std()
                else:
                    result[f'{col_base}_mean'] = np.nan
                    result[f'{col_base}_std'] = np.nan
    return result

# =====================
# Generate mean tables by recording and group
# =====================
results_by_rec = []
for rec in [1, 2]:
    results_by_rec.append(calculate_means(df, f'Overall (Recording {rec})', rec=rec))
    for factor in factors_to_split:
        if f'{factor}_group' in df.columns:
            for group_label in ['Low', 'High']:
                sub = df[df[f'{factor}_group'] == group_label]
                results_by_rec.append(calculate_means(sub, f'{factor} {group_label}', factor, group_label, rec))

results_by_rec_df = pd.DataFrame(results_by_rec)
results_by_rec_df.to_csv('median_split_results_by_recording.csv', index=False)

# =====================
# Compute Rec2 − Rec1 differences (means and stds)
# =====================
rec1 = results_by_rec_df[results_by_rec_df['Recording'] == 1].copy()
rec2 = results_by_rec_df[results_by_rec_df['Recording'] == 2].copy()
merged = pd.merge(rec1, rec2, on=['Description', 'Factor', 'Group'], how='outer', suffixes=('_Rec1', '_Rec2'))

diff_rows = []
for _, row in merged.iterrows():
    out = {
        'Description': row.get('Description', ''),
        'Factor': row.get('Factor', ''),
        'Group': row.get('Group', '')
    }
    for hand in ['MoreAffected', 'LessAffected']:
        for cond in ['preE', 'postE']:
            for beh in behavior_vars:
                base_mean = f'{hand}_{cond}_{beh}_mean'
                base_std = f'{hand}_{cond}_{beh}_std'
                val1_mean = row.get(base_mean + '_Rec1', np.nan)
                val2_mean = row.get(base_mean + '_Rec2', np.nan)
                val1_std = row.get(base_std + '_Rec1', np.nan)
                val2_std = row.get(base_std + '_Rec2', np.nan)
                out[base_mean + '_Rec2_minus_Rec1'] = val2_mean - val1_mean if pd.notna(val1_mean) and pd.notna(val2_mean) else np.nan
                out[base_std + '_Rec2_minus_Rec1'] = val2_std - val1_std if pd.notna(val1_std) and pd.notna(val2_std) else np.nan
    diff_rows.append(out)

diff_df = pd.DataFrame(diff_rows)
diff_df.to_csv('median_split_differences_rec2_minus_rec1.csv', index=False)

# =====================
# Compute Pre–Post differences (within each recording) including stds
# =====================
prepost_rows = []
for _, row in results_by_rec_df.iterrows():
    out = {
        'Description': row.get('Description', ''),
        'Factor': row.get('Factor', ''),
        'Group': row.get('Group', ''),
        'Recording': row.get('Recording', np.nan)
    }
    for hand in ['MoreAffected','LessAffected']:
        for beh in behavior_vars:
            pre_mean = row.get(f'{hand}_preE_{beh}_mean', np.nan)
            post_mean = row.get(f'{hand}_postE_{beh}_mean', np.nan)
            pre_std = row.get(f'{hand}_preE_{beh}_std', np.nan)
            post_std = row.get(f'{hand}_postE_{beh}_std', np.nan)
            if pd.notna(pre_mean) and pd.notna(post_mean):
                diff = post_mean - pre_mean
                pct = (diff / pre_mean) * 100 if pre_mean != 0 else np.nan
            else:
                diff, pct = np.nan, np.nan
            if pd.notna(pre_std) and pd.notna(post_std):
                diff_std = post_std - pre_std
            else:
                diff_std = np.nan
            out[f'{hand}_{beh}_PrePostDiff'] = diff
            out[f'{hand}_{beh}_PrePostPct'] = pct
            out[f'{hand}_{beh}_PrePostStdDiff'] = diff_std
    prepost_rows.append(out)

prepost_df = pd.DataFrame(prepost_rows)
prepost_df.to_csv('median_split_differences_prepost_by_recording.csv', index=False)

# =====================
# Combine Rec2–Rec1 and Pre–Post results
# =====================
final_merge = pd.merge(diff_df, prepost_df, on=['Description','Factor','Group'], how='outer')
final_merge.to_csv('median_split_combined_differences.csv', index=False)

print('All exports complete:')
print(' - median_split_results_by_recording.csv')
print(' - median_split_differences_rec2_minus_rec1.csv')
print(' - median_split_differences_prepost_by_recording.csv')
print(' - median_split_combined_differences.csv')
