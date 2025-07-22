
# GET NUMBERS FROM count_subsets.py

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib
import os
matplotlib.rcParams["svg.fonttype"] = "none" # do not embed directly, instead the downstream program will view with system fonts
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams["text.usetex"] = False 


# Change items here
# SPECIES = "human" # or "mouse"
# identical_percentages =  [28, 75, 87, 63, 67] # human numbers
# subset_data = [85, 98, 82, 88] # human numbers
# results_dir = "/home/drc/Downloads/refgenomes_pics_test/21jul2025/fix_the_density_maps/"
SPECIES = "mouse"
identical_percentages =  [22, 81, 94, 81, 78] # human numbers
subset_data = [81, 97, 86, 81] # human numbers
results_dir = "/home/drc/Downloads/refgenomes_pics_test/21jul2025/fix_the_density_maps/"

# Data for "Identical to" / "Duplicates"
identical_data = {
    'Comparison Type': [
        'Duplicates',
        'Names',
        'Lengths',
        'Sequences',
        'Name-len-pairs'
    ],
    'Percentage': identical_percentages
}
df_identical = pd.DataFrame(identical_data)

# Data for "Subsets"
subset_data = {
    'Comparison Type': [
        'Name-len-pairs',
        'Lengths',
        'Sequences',
        'Names'
    ],
    'Percentage':  subset_data
    }
df_subset = pd.DataFrame(subset_data)

# --- Visualization ---

# Set a nice style for the plots
sns.set_theme(style="whitegrid")

# Create a figure with two subplots (one for each type of comparison)
fig, axes = plt.subplots(1, 2, figsize=(8, 2)) # 1 row, 2 columns

# Plot for "Identical to" / "Duplicates"
sns.barplot(
    x='Percentage',
    y='Comparison Type',
    data=df_identical.sort_values(by='Percentage', ascending=False), # Sort for better readability
    palette='viridis_r', # A nice color palette
    ax=axes[0]
)
axes[0].set_title('Duplicates & identical to another reference')
axes[0].set_xlabel('Percentage (%)')
axes[0].set_ylabel('') # Clear y-axis label, categories are self-explanatory
axes[0].set_xlim(0, 100) # Ensure x-axis goes from 0 to 100

# Add percentage labels on the bars
for p in axes[0].patches:
    width = p.get_width()
    axes[0].text(width + 2, # x-position (a bit to the right of the bar)
                 p.get_y() + p.get_height() / 2, # y-position (center of the bar)
                 f'{width:.0f}%', # Label text
                 va='center')

# Plot for "Subsets"
sns.barplot(
    x='Percentage',
    y='Comparison Type',
    data=df_subset.sort_values(by='Percentage', ascending=False), # Sort for better readability
    palette='viridis_r', # Another nice color palette
    ax=axes[1]
)
axes[1].set_title('Subsets of another reference')
axes[1].set_xlabel('Percentage (%)')
axes[1].set_ylabel('')
axes[1].set_xlim(0, 100) # Ensure x-axis goes from 0 to 100

# Add percentage labels on the bars
for p in axes[1].patches:
    width = p.get_width()
    axes[1].text(width + 2,
                 p.get_y() + p.get_height() / 2,
                 f'{width:.0f}%',
                 va='center')

plt.tight_layout() # Adjust layout to prevent labels from overlapping
#plt.show()
output_path = os.path.join(results_dir, f'percentages_bargraph_{SPECIES}.svg')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

# --- Alternative: Single Plot with a Grouping Variable (less ideal for your specific text) ---
# If you wanted to combine them, you would need a 'Category' column
# and then use a grouped bar chart, but the narrative flow suggests
# separate comparisons.
# For example:
# combined_data = {
#     'Comparison Type': [
#         'Duplicates (Human Genomes)', 'Identical (Names)', 'Identical (Lengths)', 'Identical (Sequences)', 'Identical (Name-Len Pairs)',
#         'Subset (Name-Length-Pairs)', 'Subset (Lengths)', 'Subset (Sequences)', 'Subset (Names)'
#     ],
#     'Percentage': [28, 75, 87, 63, 67, 85, 98, 82, 88],
#     'Group': [
#         'Identical/Duplicates', 'Identical/Duplicates', 'Identical/Duplicates', 'Identical/Duplicates', 'Identical/Duplicates',
#         'Subsets', 'Subsets', 'Subsets', 'Subsets'
#     ]
# }
# df_combined = pd.DataFrame(combined_data)
#
# plt.figure(figsize=(10, 8))
# sns.barplot(x='Percentage', y='Comparison Type', hue='Group', data=df_combined.sort_values(by=['Group', 'Percentage'], ascending=[True, False]), palette='Set2')
# plt.title('Comparison of Reference Relationships by Type')
# plt.xlabel('Percentage (%)')
# plt.ylabel('')
# plt.xlim(0, 100)
# plt.legend(title='Relationship Type')
# plt.tight_layout()
# plt.show()