import numpy as np
import pandas as pd
from datetime import datetime
import  numpy as np
def clean_df(path):
    data= pd.read_csv(path,encoding = "ISO-8859-8")
    a=data.describe()
    data.columns
    ##clean
    data.rename(columns={'נפח בתחילת יום (L)':'volume_i','נפח לאחר השקיה (L)':'volume_f','דשן (mL)':'fertilizer_adding(ml)',
                         'קציר (L)':'harvest(l)','לכלוך':'dirt',"מקור\סיבה":'source','פעולה':'action',"dens (cell/ml)":'Cell density (cell/ml)','מקור.סיבה':'source'},inplace=True)
    data.columns =map(str.lower,data.columns)
    #volume
    data["volume_i"]=data["volume_i"].fillna(0)
    data["volume_f"]=data["volume_f"].fillna(data["volume_i"])
    #date
    data.dropna(subset=["date"],inplace=True)
    data.reset_index(inplace=True,drop=True)
    #change data to standard datatype
    for i in range(len(data)):
        if type(data["date"][i]) != type(datetime.now()):
            data.loc[i,"date"]=datetime.strptime(data.loc[i,'date'], "%d/%m/%Y")
    # col area
    data["area"]=data["area"].str.replace(' ', '')
    data['area+reactor']=data["area"]+'_'+(data["reactor"].astype(int).astype(str))
    a=data["area+reactor"].value_counts()
    #no3
    data.groupby(by="no3 (ppm)").count()
    data['no3 (ppm)'].replace('>500',500,inplace=True)
    data['no3 (ppm)'].replace('500+++',500,inplace=True)
    data['no3 (ppm)'].replace('100-250',175,inplace=True)
    data['no3 (ppm)'].replace('100-50',75,inplace=True)
    data['no3 (ppm)'].replace('25-50',37,inplace=True)
    data['no3 (ppm)'].replace('250-100',175,inplace=True)
    data['no3 (ppm)'].replace('250-500',375,inplace=True)
    data['no3 (ppm)'].replace('50-100',75,inplace=True)
    data['no3 (ppm)']=data['no3 (ppm)'].apply(lambda x: float(x))

    #fertilizer adding-
    data['fertilizer_adding(ml)']=data['fertilizer_adding(ml)'].apply(lambda x: float(x))

    #so4
    if "so4 (g)" in data.columns:
        data['so4 (g)']=data['so4 (g)'].apply(lambda x: float(x))

    #harvest
    data['harvest(l)']=data['harvest(l)'].apply(lambda x: float(x))
    data["harvest(l)"].fillna(0,inplace=True)


    #dirt and organizsms
    if "foreign algae" in data.columns:
        for i in ['foreign algae','bact','ciliate/zoo','dirt','fungi']:
            data[i].replace('0-1', 0.5, inplace=True)
            data[i].replace('1-0', 0.5, inplace=True)
            data[i].replace('1-2', 1.5, inplace=True)
            data[i].replace('2-3', 2.5, inplace=True)
            data[i].replace('01-פבר', 1.5, inplace=True)# months are resulted becouse the encoding
            data[i].replace('02-מרץ', 2.5, inplace=True)# months are resulted becouse the encoding
            data[i]=data[i].apply(lambda x: float(x))

    ##manipulation
    data["cell density (cell/ml)"].replace('1.53 E7',1.53e7,inplace=True)
    data["cell density (cell/ml)"]=data["cell density (cell/ml)"].astype(float)
    #create total cell col
    data["total_cell"]=data["cell density (cell/ml)"]*data["volume_i"]
    #create col of father,nodeand harvest , eliminated,
    data["parent"]=np.where(data["action"]=='נזרע מ',data["source"],np.nan)
    data["nodes"]=np.where((data["action"]== 'זרע את ' )| (data["action"]=='זרע את'),data["source"],np.nan)
    data["harvest_alert"]=np.where((data["action"]=='קציר')|(data["action"]=='קציר וחיסול')|(data["action"]=='קציר חיסול'),1,np.nan)
    data["cancel"]=np.where((data["action"]=='חוסל')|(data["action"]=='קציר וחיסול')|(data["action"]=='קציר חיסול')|(data["action"]=='חיסול'),1,np.nan)
    data['divide']=np.where((data["action"]=='פוצל ל'),1,np.nan)
    data.info()
    a=data.groupby(by='action').count()
    data["parent"].value_counts()

    data.info()
    #create doubling time
    #reference http://textbookofbacteriology.net/growth_3.html
    data.reset_index(drop=True,inplace=True)
    data["doubling_time"]=np.nan

    flag=0 # flag for the first result
    list_reactors=list(data["area+reactor"].value_counts().index)
    for reactor in list_reactors: # loop for each reactor
        for i in range(len(data)): #loop for each data in the reactor
            if (data.loc[i,"area+reactor"]==reactor)&(data.loc[i,"cell density (cell/ml)"]>0)&(flag==0):# find the first result of the database
                j=i
                flag=1
                continue
            if (data.loc[i,"area+reactor"]==reactor) & (data.loc[i, "cell density (cell/ml)"] > 0) & (flag == 1): # find the second,third and so onof the data base
                if data.loc[j,"harvest(l)"]>0: #in case of harvest
                    total_cell_j=data.loc[j, "cell density (cell/ml)"]*(data.loc[j, "volume_i"]-data.loc[j,"harvest(l)"])
                if data.loc[j,"harvest(l)"]==0: #in case of harvest
                    total_cell_j = data.loc[j, "total_cell"]
                delta_hr = (data.loc[i, "date"] - data.loc[j, "date"]).total_seconds() / 86400 #compute delta time
                delta_count = float(data.loc[i, "total_cell"] /(data.loc[j, "total_cell"]+0.01))
                j=i # define new variable
                if delta_count != 1: # log1=0 so avoid it
                    data.loc[i, "Doubling time"] = delta_hr * np.log(2) / (np.log(delta_count))
    return data


temp=clean_df(path="Database2022.csv")
temp1=clean_df(path="Database2022_construct.csv")

temp=temp[(temp["Doubling time"]>0)&(temp["Doubling time"]<6)]
temp["Doubling time"].hist(bins=40)

temp1=temp1[(temp1["Doubling time"]>0)&(temp1["Doubling time"]<6)]
temp1["Doubling time"].hist(bins=40)
class MR():
    def __init__(self,name):
        self.name=name
        self.reactor=[1,2,3,4,5,6,7,8,9,10,11,12]

class TPU():
    # constructor method
    def __init__(self,name):
        self.name=name
        self.lr=[1,2,3,4,5,6,7,8,9,10,11,12]
        self.sr=[21,22,23]
tpu_list=[None]
mr_list=[None]
for i in range(1,4):
    tpu_list.append(TPU(i))
    mr_list.append(MR(i))




