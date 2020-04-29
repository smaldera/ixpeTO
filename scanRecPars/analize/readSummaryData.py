import numpy as np









class readSummaryData_1d:
 

   def __init__(self,fileName): 
       self.fileName=fileName

       self.readFile()


   def readFile(self):
       print("readSummaryData:  opneing file ",self.fileName) 
       infile=open(self.fileName,'r')
       npArrayList=[]

       for line in infile:
           #print(line[:-1])
           splittedLine=line[:-1].split(" ")
           #floatList=list(map(float, splittedLine))
           #print ("flaotList = " ,floatList )
           #floatList2=[float(i) for i in splittedLine ]
           #print ("float2=",floatList2 )

           npArrayList.append(np.array(splittedLine,dtype=np.float32))
           #print ("npArray=",npArray )


       self.energy=npArrayList[0]
       self.mod1=npArrayList[1]
       self.mod1_err=npArrayList[2]
       self.mod2=npArrayList[3]
       self.mod2_err=npArrayList[4]
       self.mod1std=npArrayList[5]
       self.mod1std_err=npArrayList[6]
       self.mod2std=npArrayList[7]
       self.mod2std_err=npArrayList[8]
       self.best_val=npArrayList[9]
       self.best_val_up=npArrayList[10]
       self.best_val_low=npArrayList[11]




   def find_bestParameter(self,currentEnergy):
        """
        returns the best value of  parameter for Polarized data for a GIVEN ENERGY
        """
        if currentEnergy==5.89:
           currentEnergy=6.4
           print("WARNING !!!!!!!! E=5.89 KeV ==>> uso best value trovato a 6.4 KeV !!!!!")
        
        
        index_summary=1e6
        try:
           index_summary=np.where(  np.logical_and ( self.energy<(float(currentEnergy)+0.05), self.energy >(float(currentEnergy)-0.05) ) )[0][0]
           print ("readSummaryData:  energia  trovata! index = ",index_summary)
        except:
            print ("readSummaryData:  energia *NON* trovata nello scan ploarizzato")

            
        bestPar=1e6 
        if ( index_summary<1000):
            bestPar=self.best_val[index_summary]
        
        return bestPar   
   





       
###############################################



class readSummaryData_unpolDFF_1d:
 

   def __init__(self,fileName): 
       self.fileName=fileName

       self.readFile()


   def readFile(self):
       print("readSummaryData:  opneing file ",self.fileName) 
       infile=open(self.fileName,'r')
       npArrayList=[]

       for line in infile:
           #print(line[:-1])
           splittedLine=line[:-1].split(" ")
           #floatList=list(map(float, splittedLine))
           #print ("flaotList = " ,floatList )
           #floatList2=[float(i) for i in splittedLine ]
           #print ("float2=",floatList2 )

           npArrayList.append(np.array(splittedLine,dtype=np.float32))
           #print ("npArray=",npArray )


       self.energy=npArrayList[0]
       self.mod1=npArrayList[1]
       self.mod1_err=npArrayList[2]
       self.mod2=npArrayList[3]
       self.mod2_err=npArrayList[4]
       self.mod1std=npArrayList[5]
       self.mod1std_err=npArrayList[6]
       self.mod2std=npArrayList[7]
       self.mod2std_err=npArrayList[8]
       self.best_val=npArrayList[9]
       self.best_val_up=npArrayList[10]
       self.best_val_low=npArrayList[11]

       self.mod_phi1_bestPol=npArrayList[12]
       self.mod_phi2_bestPol=npArrayList[13]
       self.mod_phi1_bestPolErr=npArrayList[14]
       self.mod_phi2_bestPolErr=npArrayList[14]
       self.bestParPol=npArrayList[15]




       
class plotData2d:


   def __init__(self,fileName): 
       self.fileName=fileName

       self.readFile()


   def readFile(self):   
       infile=open(self.fileName,'r')
       npArrayList=[]

       for line in infile:
           #print(line[:-1])
           splittedLine=line[:-1].split(" ")
           #floatList=list(map(float, splittedLine))
           #print ("flaotList = " ,floatList )
           #floatList2=[float(i) for i in splittedLine ]
           #print ("float2=",floatList2 )

           npArrayList.append(np.array(splittedLine,dtype=np.float32))
           #print ("npArray=",npArray )


       self.energy=npArrayList[0]
       self.mod1=npArrayList[1]
       self.mod1_err=npArrayList[2]
       self.mod2=npArrayList[3]
       self.mod2_err=npArrayList[4]
       self.mod1std=npArrayList[5]
       self.mod1std_err=npArrayList[6]
       self.mod2std=npArrayList[7]
       self.mod2std_err=npArrayList[8]

       self.best_var1=npArrayList[9]
       self.best_var1_up=npArrayList[10]
       self.best_var1_low=npArrayList[11]

       self.best_var2=npArrayList[12]
       self.best_var2_up=npArrayList[13]
       self.best_var2_low=npArrayList[14]



#######################
