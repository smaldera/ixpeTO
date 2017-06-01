
import math
import ROOT
import numpy


class filtra_segnale:

  
  def  init_filtro(self):

          
        print "inizio init filtro,... "

        #self.M = 240 #lunghezza del filtro

        #self.M = 8 #lunghezza del filtro DEVE ESSERE PARI!!!
        #self.FreqCamp = 1e+9     #frequenza di campionamento in Sample/s
        #self.FreqCut =2e8  #frequenza di taglio in Hz

        #OK!!!!
        #self.M = 80 #lunghezza del filtro DEVE ESSERE PARI!!!
        #self.FreqCamp = 1e+9     #frequenza di campionamento in Sample/s
        #self.FreqCut =1e7  #frequenza di taglio in Hz

        # ok new
        self.M = 18 #lunghezza del filtro DEVE ESSERE PARI!!!
        self.FreqCamp = 1e+9     #frequenza di campionamento in Sample/s
        self.FreqCut =7e7  #frequenza di taglio in Hz
        

        #self.M = 8 #lunghezza del filtro DEVE ESSERE PARI!!!
        #self.FreqCamp = 1e+9     #frequenza di campionamento in Sample/s
        #self.FreqCut =2e8  #frequenza di taglio in Hz
        
        
        
        self.Fc = self.FreqCut/self.FreqCamp #taglio
        self.Pi = 3.14159265
        self.H=numpy.array([0]*self.M,dtype=float)   
        self.Win=numpy.array([0]*self.M,dtype=float) 


        self.H[self.M/2]=2.*self.Pi*self.Fc

        for  i in range (0, self.M ):        
	        
               self.Win[i]=(2.*self.Pi*i)/self.M	
               self.Win[i]=(0.42-0.5*math.cos(self.Win[i])+0.08*math.cos(2*self.Win[i])) #Blackman window
    
               if (i<self.M/2):
                  self.H[i]=(math.sin(2.*self.Pi*self.Fc*(i-self.M/2.)))/(i-self.M/2.)
             
               if (i>self.M/2):
                  self.H[i]=(math.sin(2.*self.Pi*self.Fc*(i-self.M/2.)))/(i-self.M/2.) 

               self.H[i]=self.H[i]*self.Win[i]

    
  
 
        sum=0
        for  i in range (0, self.M):
        
               sum=sum+self.H[i]

  
        for  i in range (0, self.M):
               self.H[i]=self.H[i]/sum

        return 1




  def  filtra(self, X):

        #aggiungo  M/2 bin a zero in x
        x0=numpy.zeros(self.M/2)
        Xfake=numpy.concatenate( (x0,X), axis=0)

        NBIN=len(Xfake);
        Yfake=numpy.array([0]*(NBIN),dtype=float) 
        print "lenx = ",len(Xfake)
        if NBIN<self.M:
          print "histogram too short!!!!!"
          return -1
        
    
         
        for  j in range (self.M, NBIN): 
           Yfake[j]=0  # ridondante!!!!
           for  i in range (0,self.M): 
		  #Y[j]=Y[j]+X[j-i]*self.H[i]         # ogni X viene moltiplicato per gli M punti della H; l'aumentare della lunghezza del filtro permette di calcolare una larghezza di banda + stretta M=4/BW
		  Yfake[j]=Yfake[j]+Xfake[j-i]*self.H[i]         # ogni X viene moltiplicato per gli M punti della H; l'aumentare della lunghezza del filtro permette di calcolare una larghezza di banda + stretta M=4/BW	 

        Y0=Yfake[self.M:]
        y1=numpy.zeros(self.M/2)
        Y=numpy.concatenate((Y0,y1))

                           
        return Y


                  



#===========================================================================

if __name__ == '__main__':

  
   h_traccia = ROOT.TH1F("h_traccia","",525,0,525)                 	       
   h_tracciaF = ROOT.TH1F("h_tracciaF","",525,0,525)                 	       


   print "creo oggetto..."
   simoFilter=filtra_segnale()

   print "init filtro.. "
   aa=simoFilter.init_filtro();
   #print aa
       	  
 
   x=numpy.array([0]*525,dtype=float)
   A=10
	   
   #	 //  0.3 - 0.003
   k=10
   omega = 0.003*k # //==> 3Mhz - 300 Mhz
   F=omega*1000.  # freq in MH
   for  jj in range (0,525):

      t=jj
      x[jj]=0
      if jj==300:
        x[jj]=10

     # if jj==299:
     #   x[jj]=11

      if jj==524:
        x[jj]=10
        
      h_traccia.Fill(jj,x[jj])                 	       

  

   
   y=simoFilter.filtra(x)
   
   print y  	       

   for  jj in range (0, 525):
      h_tracciaF.Fill(jj,y[jj])
     

   c1=ROOT.TCanvas ("c1","",0)
   c1.Divide(1,2)
   c1.cd(1)
   h_traccia.Draw()
   c1.cd(2)
   h_tracciaF.Draw()
   



