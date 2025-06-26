import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns

macs2_peak = ["Result/MACS/KO_1_peaks.bed","Result/MACS/KO_2_peaks.bed","Result/MACS/WT_1_peaks.bed","Result/MACS/WT_2_peaks.bed"]
hues = ["KO1","KO2","WT1","WT2"]


peaks=[]
names=[]
for x,y in zip(hues,macs2_peak):
    tb = pd.read_csv(y, sep='\t', header=None)
    peak_widths = list(tb[2] - tb[1])
    peaks.extend(peak_widths)
    names.extend([x] * len(peak_widths))

df=pd.DataFrame(dict(name=names,peak=peaks))

sns.kdeplot(data = df, x = "peak", hue='name')
plt.savefig('Report/peak_widths.pdf',bbox_inches='tight')

kde = gaussian_kde(peak_widths)
x = np.linspace(min(peak_widths), max(peak_widths), 1000)
y = kde(x)
plt.plot(x, y, color='purple')
plt.xlabel("Peak Width")
plt.ylabel("Density")
plt.title("Distribution of Peak Widths")
plt.savefig('Report/peak_widths.pdf',bbox_inches='tight')
plt.close('all')

