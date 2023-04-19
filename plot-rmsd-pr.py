import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('rmsd-pr.dat', delim_whitespace=True)
df = df.apply(lambda x: x*0.00005 if x.name in df.columns[:1] else x)
plt.style.use('classic')
# bold = '#7F3C8D,#11A579,#3969AC,#F2B701,#E73F74,#80BA5A,#E68310,#008695,#CF1C90,#f97b72,#4b4b8f,#A5AA99'.split(',')
bold = ['#fcff5d', '#7dfc00', '#0ec434', '#228c68', '#8ad8e8', '#235b54', '#29bdab', '#3998f5', '#37294f', '#277da7', '#3750db', '#f22020', '#991919', '#ffcba5', '#e68f66', '#c56133', '#96341c', '#632819', '#ffc413', '#f47a22', '#2f2aa0', '#b732cc', '#772b9d', '#f07cab', '#d30b94', '#edeff3', '#c3a5b4', '#946aa2', '#5d4c86', '#201923']
count = 0
for traj in df.columns[1:]:
    plt.plot(df[df.columns[0]], df[traj], label = traj, alpha=0.8, marker='.', linestyle='None', color = bold[count], markersize=4)
    count += 1
# plt.legend(markerscale=5)
plt.xlabel('Time (μs)', fontsize=16)
plt.ylabel('RMSD (Å)', fontsize=16)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.ylim(0, 20)
plt.show()
plt.savefig('rmsd30-pr.png', dpi=300, bbox_inches='tight')        

