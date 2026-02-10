import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import spearmanr, pearsonr, shapiro
import warnings
import matplotlib.lines as mlines
warnings.filterwarnings('ignore')

# Load the dataset
df = pd.read_csv('combined_demographic_behavioral.csv')
df = df.dropna(subset=['RecordingWithMoreAffectedHand'])

# Global variables
demographic_vars = ['Age', 'DiseaseDuration','DailyMedicationDose', 'self_report_improvements', 'mood', 'exercise_frequency', 'relativePlasmaLevodopa']
categorical_vars = ['self_report_improvements', 'mood', 'exercise_frequency']
behavior_vars = ['velocity', 'dwellTimes']
conditions = [('preE', 1, 'More Affected Pre'), ('postE', 1, 'More Affected Post'),
              ('preE', 0, 'Less Affected Pre'), ('postE', 0, 'Less Affected Post')]
y_lims = {'velocity': (0, 80), 'dwellTimes': (0, 0.28)}
y_ticks = {'velocity': [0, 20, 40, 60, 80], 'dwellTimes': [0, 0.05, 0.1, 0.15, 0.2, 0.25]}

# Label mapping with units
behavior_labels = {
    'velocity': 'Speed (cm/s)',
    'dwellTimes': 'Dwell time (s)'
}

demo_labels = {
    'Age': 'Age (years)',
    'DiseaseDuration': 'Disease duration (years)',
    'DailyMedicationDose': 'LEDD (mg)',
    'self_report_improvements': 'Self-reported improvement (%)',
    'mood': 'Mood',
    'exercise_frequency': 'Exercise frequency (days/week)',
    'relativePlasmaLevodopa': 'Relative medication level',
}

def get_ticks(var, values):
    x_min, x_max = min(values), max(values)
    if var == 'Age': return [50, 60, 70, 80], ['50', '60', '70', '80']
    if var == 'relativePlasmaLevodopa': return [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0']
    if var == 'self_report_improvements': return [0,1,2,3], ['None','10-30%','40-60%','60-80%']
    if var == 'mood': return [-1,0,1], ['Bad','Average','Good']
    if var == 'exercise_frequency':
        return [2, 3, 4, 5], ['1', '2', '3', '4']  
    step = max(1, int((x_max - x_min) // 5) or 1)
    ticks = list(range(int(x_min), int(x_max)+1, int(step)))
    return ticks, [str(t) for t in ticks]


def test_normality(data, alpha=0.05):
    """
    Test normality using Shapiro-Wilk test.    
    
    Returns:
    --------
    is_normal : bool
        True if data appears normally distributed
    p_value : float
        P-value from Shapiro-Wilk test
    """
    # Remove NaN values
    clean_data = data[~np.isnan(data)]
    
    # Need at least 3 observations for Shapiro-Wilk
    if len(clean_data) < 3:
        return False, np.nan
    
    # Shapiro-Wilk test
    try:
        stat, p_value = shapiro(clean_data)
        is_normal = p_value > alpha
        return is_normal, p_value
    except Exception as e:
        print(f"Normality test failed: {e}")
        return False, np.nan


def choose_correlation_method(x_data, y_data, var_name, is_predefined_categorical=False):
    """
    Choose appropriate correlation method based on:
    1. Whether variable is predefined as categorical/ordinal
    2. Normality of both variables (using Shapiro-Wilk test)
    
    Returns:
    --------
    method : str
        'spearman' or 'pearson'
    reason : str
        Explanation for the choice
    x_p : float
        P-value from normality test for x
    y_p : float
        P-value from normality test for y
    """
    # If predefined as categorical, always use Spearman
    if is_predefined_categorical:
        return 'spearman', 'predefined categorical/ordinal variable', np.nan, np.nan
    
    # Test normality of both variables
    x_normal, x_p = test_normality(x_data)
    y_normal, y_p = test_normality(y_data)
    
    # If either variable is not normal, use Spearman
    if not x_normal or not y_normal:
        reason = f'non-normal distribution'
        return 'spearman', reason, x_p, y_p
    
    # Both variables are normal, use Pearson
    reason = f'both normally distributed'
    return 'pearson', reason, x_p, y_p


def plot_demo_vs_behavior(df, demo):
    fig, axes = plt.subplots(2, 4, figsize=(18, 9.5), constrained_layout=False)
    fig.subplots_adjust(top=0.85, right=0.80, wspace=0.3, hspace=0.25)

    x_vals = []
    for cond in conditions:
        subset = df[(df['Exercise'] == cond[0]) & (df['RecordingWithMoreAffectedHand'] == cond[1])]
        if demo in subset.columns:
            x_vals.extend(subset[demo].dropna().tolist())
    if len(x_vals) == 0:
        x_vals = [0, 1]
    x_ticks, x_labels = get_ticks(demo, x_vals)

    color_map = {'preE': 'skyblue', 'postE': 'lightcoral'}
    
    # Check if this demographic is categorical
    is_categorical = demo in categorical_vars
    
    print(f"\n{'='*60}")
    print(f"Analyzing correlations for: {demo}")
    print(f"{'='*60}\n")

    for i, behavior in enumerate(behavior_vars):
        for j, (exercise, hand, title) in enumerate(conditions):
            ax = axes[i][j]
            data = df[(df['Exercise'] == exercise) & (df['RecordingWithMoreAffectedHand'] == hand)]
            if demo in data.columns and behavior in data.columns:
                use_cols = [demo, behavior]
                if 'Gender' in data.columns:
                    use_cols.append('Gender')
                data = data[use_cols].dropna()
                if len(data) >= 3:
                    x, y = data[demo], data[behavior]
                    sns.regplot(x=x, y=y, ax=ax, scatter=False, ci=95, color=color_map[exercise])
                    if 'Gender' in data.columns:
                        for gender in data['Gender'].dropna().unique():
                            marker = {'M': 's', 'F': '^'}.get(str(gender)[0].upper(), 'o')
                            gd = data[data['Gender'] == gender]
                            ax.scatter(gd[demo], gd[behavior], color='gray', marker=marker, s=25, alpha=0.7)
                    else:
                        ax.scatter(x, y, color='gray', s=25, alpha=0.7)
                    
                    # Choose correlation method based on normality
                    try:
                        x_data = x.values
                        y_data = y.values
                        
                        method, reason, x_p, y_p = choose_correlation_method(x_data, y_data, demo, is_categorical)
                        
                        if method == 'spearman':
                            corr, p = spearmanr(x_data, y_data)
                            method_label = "ρ"  # Spearman's rho
                        else:
                            corr, p = pearsonr(x_data, y_data)
                            method_label = "r"  # Pearson's r
                        
                        ax.text(0.68, 0.93, f'{method_label}={corr:.2f}\np={p:.3f}', 
                                transform=ax.transAxes, va='top', fontsize=12)
                        
                        # Print diagnostic information
                        hand_label = "More affected" if hand == 1 else "Less affected"
                        exercise_label = "Pre-exercise" if exercise == 'preE' else "Post-exercise"
                        print(f"{hand_label} - {exercise_label} - {behavior}:")
                        print(f"  Method: {method.capitalize()}")
                        print(f"  Reason: {reason}")
                        if not np.isnan(x_p) and not np.isnan(y_p):
                            print(f"  Normality tests: x_p={x_p:.3f}, y_p={y_p:.3f}")
                        print(f"  {method_label}={corr:.2f}, p={p:.3f}")
                        print()
                        
                    except Exception as e:
                        print(f"Correlation failed for {hand_label} - {exercise_label} - {behavior}: {e}")
                        pass
                else:
                    ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes, ha='center')

            if demo in ['self_report_improvements', 'mood', 'exercise_frequency']:
                ax.set_xlim(min(x_ticks) - 0.5, max(x_ticks) + 0.5)
            elif demo == 'relativePlasmaLevodopa':
                ax.set_xlim(0.0, 1.05)
            elif demo in ['Age', 'DiseaseDuration', 'DailyMedicationDose']:
                if x_vals:
                    pad = 0.05 * (max(x_vals) - min(x_vals))
                    ax.set_xlim(min(x_vals) - pad, max(x_vals) + pad)
            else:
                ax.set_xlim(min(x_vals), max(x_vals))

            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels, fontsize=11)

            if i == len(behavior_vars) - 1:
                ax.set_xlabel(demo_labels.get(demo, demo), fontsize=13, fontweight='bold')
            else:
                ax.set_xlabel('')
                ax.set_xticklabels([])

            ax.set_ylim(y_lims[behavior])
            ax.set_yticks(y_ticks[behavior])
            ax.tick_params(axis='y', labelsize=11)

            label = behavior_labels.get(behavior, behavior)
            if j == 0:
                ax.set_ylabel(label, fontsize=13, fontweight='bold')
            else:
                ax.set_ylabel('')
                ax.set_yticklabels([])

            if i == 0:
                ax.set_title(('Pre-exercise' if exercise == 'preE' else 'Post-exercise'), fontsize=13, fontweight='bold')
            else:
                ax.set_title('')

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

    triangle = mlines.Line2D([], [], color='gray', marker='^', linestyle='None', markersize=8, label='Female')
    square = mlines.Line2D([], [], color='gray', marker='s', linestyle='None', markersize=8, label='Male')
    fig.legend(handles=[square, triangle], loc='upper right', bbox_to_anchor=(0.96, 0.93), frameon=False, fontsize=12)

    top_axes = axes[0]
    x_left = (top_axes[0].get_position().x0 + top_axes[1].get_position().x1) / 2
    x_right = (top_axes[2].get_position().x0 + top_axes[3].get_position().x1) / 2
    y_header = top_axes[0].get_position().y1 + 0.04

    fig.text(x_left,  y_header, 'More affected hand', ha='center', va='bottom', fontsize=14, fontweight='bold', transform=fig.transFigure)
    fig.text(x_right, y_header, 'Less affected hand', ha='center', va='bottom', fontsize=14, fontweight='bold', transform=fig.transFigure)

    return fig

for demo_var in demographic_vars:
    if demo_var in df.columns:
        fig = plot_demo_vs_behavior(df, demo_var)
        plt.show()
        fig.savefig(f'{demo_var}_vs_behaviors.svg', bbox_inches="tight")
        fig.savefig(f'{demo_var}_vs_behaviors.png', bbox_inches="tight", dpi=300)

print("\n=== All plots saved ===")