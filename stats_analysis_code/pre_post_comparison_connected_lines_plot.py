import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import ttest_rel, wilcoxon

IQR_MULTIPLIER = 1.5

# ----------------------------
#  Prepare averaged dataframe
# ----------------------------
data = pd.read_csv('combined_demographic_behavioral.csv')

grouped = (
    data.loc[data['Recording'].isin([1, 2])]
        .groupby(['ParticipantID_num', 'Exercise', 'RecordingWithMoreAffectedHand'], as_index=False)
        .mean(numeric_only=True)
)

# Map exercise labels to 0/1
grouped['Exercise'] = grouped['Exercise'].map({'preE': 0, 'postE': 1}).astype(int)

# Merge redundent gender
if 'Gender' in data.columns:
    gender_map = data[['ParticipantID_num', 'Gender']].drop_duplicates()
    grouped = grouped.merge(gender_map, on='ParticipantID_num', how='left')
    grouped['gender(female or not)'] = grouped['Gender'].map({'male': 0, 'female': 1})

# Rename columns
column_mapping = {
    'ParticipantID_num': 'patients id',
    'RecordingWithMoreAffectedHand': 'Record with most affected hand(1) or not(0)',
    'velocity': 'Speed (cm/s)',
    'dwellTimes': 'Dwell Time (s)',
    'nrErrors': 'Number of Errors'
}
grouped = grouped.rename(columns={k: v for k, v in column_mapping.items() if k in grouped.columns})

df = grouped.copy()

# ----------------------------
# stats helper
# ----------------------------
def paired_stats(df_hand: pd.DataFrame, metric: str, iqr_multiplier: float = 1.5):
    """Return paired stats for pre(0) vs post(1). Outlier removal for non-Error metrics only."""
    paired = []
    for pid in df_hand['patients id'].unique():
        tmp = df_hand[df_hand['patients id'] == pid]
        if tmp['Exercise'].nunique() == 2:
            pre = tmp.loc[tmp['Exercise'] == 0, metric].iloc[0]
            post = tmp.loc[tmp['Exercise'] == 1, metric].iloc[0]
            paired.append((pre, post))

    if len(paired) < 2:
        return {
            "N": len(paired),
            "N_removed": 0,
            "test": None,
            "p": np.nan,
            "paired_df": pd.DataFrame(paired, columns=['pre', 'post'])
        }

    paired_df = pd.DataFrame(paired, columns=['pre', 'post'])

    # if number of errors, wilcoxon
    if metric == 'Number of Errors':
        stat, p = wilcoxon(paired_df['pre'], paired_df['post'])
        return {"N": len(paired_df), "N_removed": 0, "test": "Wilcoxon", "p": p, "paired_df": paired_df}

    # Outlier removal 
    pre_q1, pre_q3 = paired_df['pre'].quantile([0.25, 0.75])
    post_q1, post_q3 = paired_df['post'].quantile([0.25, 0.75])
    pre_iqr = pre_q3 - pre_q1
    post_iqr = post_q3 - post_q1

    pre_ok = (paired_df['pre'] >= pre_q1 - iqr_multiplier * pre_iqr) & (paired_df['pre'] <= pre_q3 + iqr_multiplier * pre_iqr)
    post_ok = (paired_df['post'] >= post_q1 - iqr_multiplier * post_iqr) & (paired_df['post'] <= post_q3 + iqr_multiplier * post_iqr)

    clean = paired_df[pre_ok & post_ok]
    n_removed = len(paired_df) - len(clean)

    if len(clean) < 2:
        return {"N": len(paired_df), "N_removed": n_removed, "test": "paired t-test", "p": np.nan, "paired_df": clean}

    t, p = ttest_rel(clean['pre'], clean['post'])
    return {"N": len(paired_df), "N_removed": n_removed, "test": "paired t-test", "p": p, "paired_df": clean}

def format_p(p, test_name=None):
    if np.isnan(p):
        return "Stats N/A"

    # base p-text
    if p < 0.001:
        p_text = "p < 0.001"
    elif p < 0.01:
        p_text = "p < 0.01"
    elif p < 0.05:
        p_text = "p < 0.05"
    else:
        p_text = f"p = {p:.3f}"

    # only append for Wilcoxon
    if test_name == "Wilcoxon":
        p_text += " (Wilcoxon)"
    return p_text


# ----------------------------
# Plot function
# ----------------------------
def plot_panel(df, ax, metric, affected_hand, column_index):
    df_hand = df[df['Record with most affected hand(1) or not(0)'] == affected_hand]

    colors = ['skyblue', 'lightcoral']  # pre, post

    # Half violins
    for ex in [0, 1]:
        vals = df_hand.loc[df_hand['Exercise'] == ex, metric].dropna().values
        if len(vals) == 0:
            continue
        parts = ax.violinplot(vals, points=100, positions=[ex], widths=0.25,
                              showmeans=False, showextrema=False, showmedians=False)
        for b in parts['bodies']:
            m = np.mean(b.get_paths()[0].vertices[:, 0])
            if ex == 0:
                b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
            else:
                b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
            b.set_color(colors[ex])
            b.set_alpha(0.4)

    # Boxplot overlay
    sns.boxplot(x='Exercise', y=metric, data=df_hand, width=0.1, ax=ax,
                hue='Exercise', palette=colors, legend=False, zorder=10)

    # Paired lines/points
    for pid in df_hand['patients id'].unique():
        tmp = df_hand[df_hand['patients id'] == pid]
        if tmp['Exercise'].nunique() == 2:
           pre = tmp.loc[tmp['Exercise'] == 0, metric].iloc[0]
           post = tmp.loc[tmp['Exercise'] == 1, metric].iloc[0]

           female = tmp['gender(female or not)'].iloc[0] if 'gender(female or not)' in tmp.columns else 0
           marker = '^' if female == 1 else 's'

           ax.plot([0.3, 0.7], [pre, post],
                color='gray', alpha=0.6, linewidth=0.8,
                linestyle='--')

           ax.scatter([0.3], [pre], color='gray', marker=marker, s=20, alpha=0.6)
           ax.scatter([0.7], [post], color='gray', marker=marker, s=20, alpha=0.6)


    # Stats text 
    stats = paired_stats(df_hand, metric, IQR_MULTIPLIER)
    ax.text(0.5, 0.95, format_p(stats["p"], stats["test"]),
            ha='center', va='top', transform=ax.transAxes, fontsize=12, fontweight='bold')

    # Axis styling
    ax.set_ylabel(metric, fontsize=14, fontweight='bold')
    ax.set_xlabel("")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Pre Exercise', 'Post Exercise'], fontsize=14, fontweight='bold')

    if column_index == 1:
        ax.yaxis.set_visible(False)
        ax.spines['left'].set_visible(False)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ----------------------------
# Create figure
# ----------------------------
fig, axs = plt.subplots(3, 2, figsize=(18, 24),
                        gridspec_kw={'width_ratios': [1, 1], 'hspace': 0.4, 'wspace': 0.3})

metrics = ['Speed (cm/s)', 'Dwell Time (s)', 'Number of Errors']

for r, metric in enumerate(metrics):
    for affected_hand in [1, 0]:
        c = 1 - affected_hand
        plot_panel(df, axs[r, c], metric, affected_hand, c)

fig.text(0.25, 0.95, 'More Affected Hand', ha='center', va='center', fontsize=16, fontweight='bold')
fig.text(0.75, 0.95, 'Less Affected Hand', ha='center', va='center', fontsize=16, fontweight='bold')

# Legend
lines = [
    plt.Line2D([0], [0], color='gray', marker='s', linestyle='-', markersize=12, alpha=0.6),
    plt.Line2D([0], [0], color='gray', marker='^', linestyle='--', markersize=12, alpha=0.6),
]
fig.legend(lines, ['Male', 'Female'], loc='upper center', bbox_to_anchor=(0.52, 0.9), ncol=1, fontsize=14)

plt.tight_layout(pad=3.0, rect=[0, 0, 1, 0.88])
plt.savefig('bw_full_comparison.svg', bbox_inches='tight')

plt.show()

# ----------------------------
#  Print stats summary 
# ----------------------------
print("\n" + "="*60)
print("STATISTICAL RESULTS")
print("="*60)

for metric in metrics:
    if metric not in df.columns:
        print(f"\n{metric}: column not found")
        continue

    print(f"\n{metric}\n" + "-"*40)
    for affected_hand in [1, 0]:
        label = "More Affected Hand" if affected_hand == 1 else "Less Affected Hand"
        df_hand = df[df['Record with most affected hand(1) or not(0)'] == affected_hand]
        stats = paired_stats(df_hand, metric, IQR_MULTIPLIER)

        print(f"{label}: N={stats['N']}, removed={stats['N_removed']}, test={stats['test']}, p={stats['p'] if not np.isnan(stats['p']) else 'NA'}")
