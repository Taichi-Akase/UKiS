# This program calculates the length distribution of
# the first introns and all introns based on the exon
# information of the transcript obtained from Refseq database.
# After measurement, the results are output as a boxplot.

import pandas as pd
import matplotlib.pyplot as plt

data_path = "./Source_Data/RefSeq.csv"
temp_df = pd.read_csv(data_path)

# ID, exon coordinates and strand information are retrieved and stored in list.
temp_ex_li = []
for i in temp_df.index:
    temp_id = str(temp_df["RefSeqID"][i])
    ex_con = int(temp_df["exonCount"][i])
    ex_st = [int(i2) for i2 in temp_df["exonStarts"][i].rstrip(",").split(",")]
    ex_en = [int(i2) for i2 in temp_df["exonEnds"][i].rstrip(",").split(",")]
    ex_dir = str(temp_df["strand"][i])

    temp_ex_li.append([temp_id,ex_con,ex_st,ex_en,ex_dir])

# Calculate introns length from exon coordinates.
temp_int_li = []
for i in temp_ex_li:
    temp_subst = []
    intron_num = int(i[1] -1)
    intron_st = i[2]
    intron_en = i[3]

    for num in range(intron_num):
        intron_len = intron_st[num+1] - intron_en[num] + 1
        temp_subst.append(intron_len)
    temp_int_li.append(i+[temp_subst])

# Summarize the total intron length
total_int_num = []
for i in temp_int_li:
    for j in i[5]:
        total_int_num.append(j)

# Summarize the first intron length
first_int_num = []
for i in temp_int_li:
    if i[4] == "+": first_int = i[5][0]
    elif i[4] == "-": first_int = i[5][-1]
    first_int_num.append(first_int)


# Drawing box plot
plt.rcParams["font.size"] = 25
fig, ax = plt.subplots(figsize=(8, 16))

bp = ax.boxplot([first_int_num,total_int_num],
    whis=[0,100], widths=0.5, # Setting of plot area and box width
    patch_artist = True,
    medianprops=dict(color='black', linewidth=2),  # Plot the median line
    boxprops=dict(facecolor='#d3d3d3', color='black', linewidth=2), #
)

# Plot area settings
spines = 2
ax.spines["top"].set_linewidth(spines)
ax.spines["left"].set_linewidth(spines)
ax.spines["bottom"].set_linewidth(spines)
ax.spines["right"].set_linewidth(spines)

# Visualizing the line of TP53 1st intron length
TP53 = 10754 # Length of TP53 1st intron
bp = ax.hlines(TP53, 0, 3, "black", linestyles='dashed')
ax.set_xticklabels(['1st introns','all introns'])

# Setting of X axis
plt.xlim([0,3])

# Setting of Y axis
plt.ylim([0.5,3000000])
plt.yscale("log")
plt.ylabel('Intron length (bp)')

# Plot and save
plt.savefig("IntronSize.png",dpi=300)
plt.show()
