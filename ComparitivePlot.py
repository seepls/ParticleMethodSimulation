import matplotlib.pyplot as plt 
import pandas as pd 

#plot for DoubleHelix
#set directory 
df =  pd.read_csv('InterfaceD.area', sep='\s+',)
#set plot 

plt.plot(df['Time'],df['LiquidWallTh'] , color= "red")

#Plot for single helix 
#set directory 
df =  pd.read_csv('InterfaceS.area', sep='\s+',)
#set plot 

plt.plot(df['Time'],df['LiquidWallTh'],color = "blue")


# plot for Tripple Paddle 
#set directory 
df =  pd.read_csv('InterfaceT.area', sep='\s+',)

#set plot 

plt.plot(df['Time'],df['LiquidWallTh'],color = "black")

#labels 
plt.xlabel('time [s]')
plt.ylabel('InterfaceArea')
plt.title('InterfaceArea between the 2 liquid ')
#plt.title('Red= DoubleHelical Black= TripplePaddle Blue= SingleHelical')

plt.legend()
plt.show()
