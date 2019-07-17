import matplotlib.pyplot as plt 
import pandas as pd 

#set directory 
df =  pd.read_csv('mixing.index', sep='\s+')
#set plot 

plt.plot(df['Time'],df['Liquid1'])

#labels 
plt.xlabel('time (sec)')
plt.ylabel('Mixing Index')
plt.title('mixing index of liquid with respect to time ')

plt.legend()
plt.show()
