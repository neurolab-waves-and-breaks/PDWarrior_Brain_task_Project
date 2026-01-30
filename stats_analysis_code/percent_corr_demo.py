import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
import warnings
import matplotlib.lines as mlines

warnings.filterwarnings('ignore')


# Load data

df = pd.read_csv('combined_demographic_behavioral.csv')
df = df.dropna(subset=['RecordingWithMoreAffectedHand'])


# Settings

selected_demos = ['self_report_improvements', 'Age', 'DiseaseDuration', 'DailyMedicationDose', 'exercise_frequency']
categorical_vars = ['self_report_improvements', 'mood', 'exercise_frequency']

BEH_SPEED = 'velocity'    
BEH_DWELL = 'dwellTimes'    

axis_labels = {
    'self_report_improvements': 'Self-report improvements',
    'exercise_frequency': 'Exercise frequency',
    'DailyMedicationDose': 'Daily medication dose',
    'DiseaseDuration': 'Disease duration',
    'Age': 'Age (years)',
    'mood': 'Mood',
    'relativePlasmaLevodopa': 'Relative plasma levodopa',
}

# Locate participant ID column

possible_id_cols = ['ParticipantID', 'pid', 'PID', 'Participant', 'subj', 'Subject', 'ID', 'id', 'participant_id', 'PartID']
ID_COL = ''
for _c in possible_id_cols:
    if _c in df.columns:
        ID_COL = _c
        break
if ID_COL == '':
    raise ValueError(
        'No participant ID column found. Add one of these columns to your CSV: ' + ', '.join(possible_id_cols)
    )


# Helpers


def get_ticks(var, values):
    """Consistent x ticks (matching your earlier code)."""
    if len(values) == 0:
        return [], []
    x_min, x_max = (min(values), max(values))

    if var == 'Age':
        return [50, 60, 70, 80], ['50', '60', '70', '80']
    if var == 'relativePlasmaLevodopa':
        return [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0']
    if var == 'self_report_improvements':
        return [0, 1, 2, 3], ['None', '10-30%', '40-60%', '60-80%']
    if var == 'mood':
        return [-1, 0, 1], ['Bad', 'Average', 'Good']
    if var == 'exercise_frequency':
        return [2, 3, 4, 5], ['1', '2', '3', '4']

    if x_max == x_min:
        return [x_min], [str(x_min)]

    step = max(1, int((x_max - x_min) // 5) or 1)
    ticks = list(range(int(x_min), int(x_max) + 1, int(step)))
    return ticks, [str(t) for t in ticks]


# Compute percent change per participant, per hand
# % change = (post - pre) / pre * 100

work = df[df['Exercise'].isin(['preE', 'postE'])].copy()


def one_demo_per_id(demo):
    post = work[work['Exercise'] == 'postE'].groupby(ID_COL)[demo].first()
    pre = work[work['Exercise'] == 'preE'].groupby(ID_COL)[demo].first()
    return post.combine_first(pre)


# gender map
if 'Gender' in df.columns:
    gender_map = df[[ID_COL, 'Gender']].dropna().drop_duplicates(subset=[ID_COL])
else:
    gender_map = None


def percent_change_by_hand(hand_val, beh_cols):
    dd = work[work['RecordingWithMoreAffectedHand'] == hand_val].copy()
    out = {}
    for beh in beh_cols:
        if beh not in dd.columns:
            continue
        pv = dd.pivot_table(index=ID_COL, columns='Exercise', values=beh, aggfunc='first')
        pre = pv.get('preE')
        post = pv.get('postE')
        out[beh] = (post - pre) / pre * 100.0
    out_df = pd.DataFrame(out)
    out_df = out_df.replace([np.inf, -np.inf], np.nan)
    return out_df


ma = percent_change_by_hand(1, [BEH_SPEED, BEH_DWELL])  # more affected
la = percent_change_by_hand(0, [BEH_SPEED, BEH_DWELL])  # less affected

# Build a wide table keyed by ID
wide = pd.DataFrame(index=ma.index.union(la.index))

# speed
if BEH_SPEED in ma.columns:
    wide['pct_ma_vel'] = ma[BEH_SPEED]
if BEH_SPEED in la.columns:
    wide['pct_la_vel'] = la[BEH_SPEED]

# dwell
if BEH_DWELL in ma.columns:
    wide['pct_ma_dwell'] = ma[BEH_DWELL]
if BEH_DWELL in la.columns:
    wide['pct_la_dwell'] = la[BEH_DWELL]

# demo columns
for demo in selected_demos:
    if demo in work.columns:
        wide[demo] = one_demo_per_id(demo)

# gender
wide = wide.reset_index().rename(columns={'index': ID_COL})
if gender_map is not None:
    wide = wide.merge(gender_map, on=ID_COL, how='left')


# Plotting: 2x2 figures


def _panel(ax, data, xcol, ycol, title, ylabel, color='orange', is_categorical=False):
    """Single subplot: reg line + CI band + grey points with gender markers + r/p annotation."""
    sub = data[[xcol, ycol] + (['Gender'] if 'Gender' in data.columns else [])].dropna(subset=[xcol, ycol]).copy()

    if sub.shape[0] < 3:
        ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes, ha='center', va='center')
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel(axis_labels.get(xcol, xcol))
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        return

    # Regression line (no scatter) with CI
    sns.regplot(x=xcol, y=ycol, data=sub, ax=ax, scatter=False, ci=95, color=color)

    # Scatter points (grey; marker by gender)
    if 'Gender' in sub.columns:
        for gender in sub['Gender'].dropna().unique():
            marker = {'M': 's', 'F': '^'}.get(str(gender)[0].upper(), 'o')
            gd = sub[sub['Gender'] == gender]
            ax.scatter(gd[xcol], gd[ycol], marker=marker, s=28, alpha=0.85, color='grey')
    else:
        ax.scatter(sub[xcol], sub[ycol], s=28, alpha=0.85, color='grey')

    # r / p
    try:
        if is_categorical:
            r, p = spearmanr(sub[xcol], sub[ycol])
        else:
            r, p = pearsonr(sub[xcol], sub[ycol])
        ax.text(0.68, 0.92, f"r={r:.2f}\np={p:.3f}", transform=ax.transAxes, va='top')
    except Exception:
        pass

    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel(axis_labels.get(xcol, xcol))
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def plot_2x2_vs_demo(wide_df, demo, color='orange', save_prefix=None):
    """2x2 grid with the SAME spacing logic as your reference plot:
    - x labels/ticks only on bottom row
    - y labels/ticks only on left column
    - explicit margins so y-labels never clip when saved
    """

    needed = [demo, 'pct_ma_vel', 'pct_la_vel', 'pct_ma_dwell', 'pct_la_dwell']
    missing = [c for c in needed if c not in wide_df.columns]
    if missing:
        print(f"[WARN] Skip {demo}: missing columns {missing}")
        return None

    fig, axes = plt.subplots(
    2, 2,
    figsize=(7, 6),
    sharey='row' 
    )

    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.14, top=0.90, wspace=0.25, hspace=0.28)

    x_vals = wide_df[demo].dropna().tolist()
    x_ticks, x_labels = get_ticks(demo, x_vals)
    is_cat = demo in categorical_vars

    def _format_x(ax):
        if x_ticks:
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels, fontsize=11)

        if demo in ['self_report_improvements', 'mood', 'exercise_frequency'] and x_ticks:
            ax.set_xlim(min(x_ticks) - 0.5, max(x_ticks) + 0.5)
        elif demo == 'relativePlasmaLevodopa':
            ax.set_xlim(0.0, 1.05)
        elif demo in ['Age', 'DiseaseDuration', 'DailyMedicationDose'] and len(x_vals) > 1:
            pad = 0.05 * (max(x_vals) - min(x_vals))
            ax.set_xlim(min(x_vals) - pad, max(x_vals) + pad)
        elif len(x_vals) > 1:
            ax.set_xlim(min(x_vals), max(x_vals))

    # ---- 4 panels ----
    _panel(axes[0, 0], wide_df, demo, 'pct_ma_vel',
           title='More affected hand', ylabel='% change of Speed',
           color=color, is_categorical=is_cat)

    _panel(axes[0, 1], wide_df, demo, 'pct_la_vel',
           title='Less affected hand', ylabel='% change of Speed',
           color=color, is_categorical=is_cat)

    _panel(axes[1, 0], wide_df, demo, 'pct_ma_dwell',
           title='', ylabel='% change of dwell time',
           color=color, is_categorical=is_cat)

    _panel(axes[1, 1], wide_df, demo, 'pct_la_dwell',
           title='', ylabel='% change of dwell time',
           color=color, is_categorical=is_cat)

    # Top row: no x labels/ticklabels
    for ax in axes[0, :]:
        ax.set_xlabel('')
        ax.set_xticklabels([])

    # Bottom row: x label
    for ax in axes[1, :]:
        ax.set_xlabel(axis_labels.get(demo, demo), fontsize=13, fontweight='bold')
        ax.tick_params(axis='x', labelsize=11)

    # Left column: y labels + y ticks
    for ax in axes[:, 0]:
        ax.tick_params(axis='y', labelleft=True, labelsize=11)
        ax.yaxis.label.set_size(13)
        ax.yaxis.label.set_weight('bold')
        ax.yaxis.labelpad = 10

    # Right column: hide y tick labels ONLY 
    for ax in axes[:, 1]:
        ax.set_ylabel('')
        ax.tick_params(axis='y', left=True, labelleft=False)


    for ax in axes.flat:
        _format_x(ax)


    

    if save_prefix is None:
        save_prefix = f'{demo}_percent_change_vs_demo_2x2'

    fig.savefig(f'{save_prefix}.svg', bbox_inches='tight', pad_inches=0.2)

    return fig


# Generate one 2x2 plot per demo

for demo in selected_demos:
    if demo not in wide.columns:
        continue
    fig = plot_2x2_vs_demo(wide, demo, color='orange')
    if fig is not None:
        plt.show()


