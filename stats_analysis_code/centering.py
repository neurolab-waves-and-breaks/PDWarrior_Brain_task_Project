import pandas as pd

#  Load CSV file
df = pd.read_csv('combined_dwell_demographic.csv')  # <-- change filename if needed

#  Variables to center
vars_to_center = ['Age', 'DiseaseDuration', 'RelativePlasmaLevodopa']

#  Create centered versions (subtract the mean)
for var in vars_to_center:
    centered_name = var + '_c'
    df[centered_name] = df[var] - df[var].mean()

#  Save to a new CSV to use in GLMM
df.to_csv('combined_dwell_demographic_centered.csv', index=False)
