'''
Versione commentata del programma scritto da Konrad Mielke per leggere i file ascii con le tracce di ixpe (x pix, y pix, carica pix) e i corrispondenti valori MC (energy; x,y di absorption; theta e phi del pe). Es: simIxpe_5Kev.txt, prodotto da Simo, a partire dal LV0b.fits.
Questo programma produce due file in formato binario hdf5 contenenti le immagini degli eventi (trasformate con pixel quadrati) e le "labels" per ciascun evento (i valori MC)


'''

import ROOT

#usa pandas, un pacchetto di tools che fornisce data-structures (strutture dei dati, tipo numpy) e la data-analysis
#in questo progamma e` usato per leggere il file ascii. Sostanzialemente gestisce mega matricione con i dati  
import pandas as pd

#usa numpy, tools per la gestione degli array e di funzioni matematiche
import numpy as np

#usa h5py, che e` un interfaccia al formato di dati binario HDF5, che e` spesso usato per le immagini. Qui e` usato perche` e`  il formato che si deve dare a tensorflow per la creazione e l'uso della Neural Network. Anche ctapipe usa questo formato.
import h5py as h5

#usa tqdm, un tool per seguire l'evoluzione del programma python da terminale con comode sbarre di evoluzione (progress bar)
from tqdm import tqdm

#usa matplotlib per la creazione dei plot con le immagini degli eventi
import matplotlib.pyplot as plt

#i file ascii di input hanno per ogni evento:
#una prima linea con le info MC sull'evento, compreso il numero di pixel della traccia (npix)
#e npix linee con x, y e carica di ciascun pixel della traccia

#crea le liste in cui andranno a finire i numeri di linea per la linee iniziali di ciascun evento (info MC, s1) e delle linee con le info pixle-wise (s2)
s1 = list()
s2 = list()

#apre il file di input, con il metodo fico "with open(...) as ..." che serve per gestire gli errori, le exceptions e la chiusura finale del file
with open("simIxpe_5Kev_short.txt", "r") as f:
#readlines() legge tutte le linee fino a EOF e restituisce, lines, una list di string (credo di string...)
    lines = f.readlines()
	 #scrive quanti elementi ha la lista lines, quindi il numero totale di linee del file ascii	 
    print(len(lines))
	 
	 #scorre tutte le linee e salva il numero di linea corrsipondente alle linee "iniziali" dell'evento (valori MC, linee lunghe) in s1, lista di interi.
	 #usa anche tqdm, che fa apparire sul terminale una sbarra di scorrimento del loop
	 #perche` nel loop dentro with ? 
    for i in tqdm(range(0,len(lines))):
        if (len(lines[i]) > 100):
            s1.append(i)

#scorre una seconda volta tutte le linee e prima mette tutti i numeri di linea...
for i in tqdm(range(0,len(lines))):
    s2.append(i)
#...poi toglie i numeri di linea corrispondenti alle linee iniziali degli eventi (quelli in s1)...non troppo efficiente...
for i in tqdm((s1)):
    s2.remove(i)

'''
pixels e labels sono delle matrici, pixels (npixels_total x 3), labels (nevents x 14), che contengono tutti i dati del file 
vengono riempite con il read_csv (Read CSV (comma-separated) file into DataFrame) di pandas.
gli oggetti pixels e labels hanno certe funzioni, tipo .size, .values, e altre, tra cui .x, .y, .int per accedere alle variabili, chiamate cosi quando si e` costruito l'oggetto
'''  
pixels = pd.read_csv("simIxpe_5Kev_short.txt", skiprows=s1, header=None, delimiter = "  ", names = ["x","y","int"])
labels = pd.read_csv("simIxpe_5Kev_short.txt", skiprows=s2, header=None, delimiter = " ", names = ["garbage1","ID","garbage2","nPixels","garbage3","Energy","garbage4","x","garbage5","y","garbage6","phi","garbage7","theta"])

pixels_orig = pd.read_csv("simIxpe_5Kev_short.txt", skiprows=s1, header=None, delimiter = "  ", names = ["x","y","int"])

#non fare: 
#pixels_orig = pixels
#altrimenti pixels_orig e pixels puntano sempre agli stessi valori, cambi uno, cambia l'altro

number_of_all_pixels = int(pixels.size/3)
number_of_all_events = int(labels.size/14)

#questo non mi e` chiaro...moltiplica tutte le x per 40 e tutte le y per 40/3.464 ??? perche`???
#passa da x e y in mm a una nuova scala in cui Delta x = 1 e Delta y = 0.5
pixels.x *= 40
pixels.y *= 40/3.464 #3.464 = sqrt(3)*2

'''
Vuole trovare il valore massimo e minimo di x e y, di tutti i pixel di tutti gli eventi
'''
x_min = 0
y_min = 0
x_max = 0
y_max = 0

x_min_orig = 0
y_min_orig = 0
x_max_orig = 0
y_max_orig = 0

#che differenza c'e` tra pixels.x.values[i] e pixels.x[i], sara` il solito "valore singolo" o "vettore con tutti i valori" ? o un baco? R. un baco, innoquo.
for i in tqdm(range(0,number_of_all_pixels)):
	print(pixels_orig.x.values[i],pixels.x.values[i],pixels.x[i],pixels_orig.y.values[i],pixels.y.values[i],pixels.y[i])
	if pixels.x.values[i] < x_min:
		x_min = pixels.x[i]
	if pixels.y.values[i] < y_min:
		y_min = pixels.y[i]
	if pixels_orig.x.values[i] < x_min_orig:
		x_min_orig = pixels_orig.x[i]
	if pixels_orig.y.values[i] < y_min_orig:
		y_min_orig = pixels_orig.y[i]

		
#sottrae ad ogni pixel.x e .y i valori minimi...ma a che serve sottrarre il valore minimi assoluto, che non c'entra con il singolo evento ? inoltre poi calcola anche y_max e x_max ma poi non li usa...
pixels.x -= x_min
pixels.y -= y_min

pixels_orig.x -= x_min_orig
pixels_orig.y -= y_min_orig

for i in tqdm(range(0,number_of_all_pixels)):
	if pixels.x.values[i] > x_max:
		x_max = pixels.x[i]
	if pixels.y.values[i] + 0.5*pixels.x.values[i] > y_max:
		y_max = pixels.y[i] + 0.5*pixels.x.values[i]
	if pixels_orig.x.values[i] > x_max_orig:
		x_max_orig = pixels_orig.x[i]
	if pixels_orig.y.values[i] > y_max:
		y_max_orig = pixels_orig.y[i]
		
sum_pixels = 0
num_pixels = 0

max_dif_x = 0
max_dif_y = 0
max_dif_x_orig = 0
max_dif_y_orig = 0

#prepara e inizializza a 0 una matrice con 4 valori per ogni evento (size_x,size_y, xmin, ymin)
picture_size_list = np.zeros((number_of_all_events, 4), dtype = np.float16)
picture_size_list_orig = np.zeros((number_of_all_events, 4), dtype = np.float16)

#Fa un primo loop su tutti gli eventi per trovare le dimensioni di ciascuna immagine
for i in tqdm(range(0,int(number_of_all_events))):
    
# fa un loop su tutti i pixel dell evento, per trovare per ciascun evento la dimensione dell'immagine 
# facendo temp_pixels.size (numero totale di elementi) bisogna ricordarsi che e` una matrice npixels_totali x 3   	   
	 num_pixels = labels.nPixels[i]
	 sum_pixels += num_pixels
#prende solo i pixel relativi all'evento: da (ultimo pixel - i pixel di questo evento) a (ultimo pixel)
	 temp_pixels = pixels[sum_pixels-num_pixels:sum_pixels]
	 temp_pixels_orig = pixels_orig[sum_pixels-num_pixels:sum_pixels]
	 #anche temp_pixels e` una matrice num_pixels x 3 (x,y,carica)
	 temp_x_min = 500  
	 temp_x_max = 0
	 temp_y_min = 500
	 temp_y_max = 0
	 temp_x_min_orig = 500  
	 temp_x_max_orig = 0
	 temp_y_min_orig = 500
	 temp_y_max_orig = 0    
	 
	 for k in range(0,int(temp_pixels.size/3)):
	 	x = int(temp_pixels.x[sum_pixels-num_pixels+k])
		#i valori y vengono shifati (passaggio da pixel esagonali a quadrati) e sui nuovi valori y cerca max e min da cui la size_y
		y = int(temp_pixels.y[sum_pixels-num_pixels+k]+0.5*temp_pixels.x[sum_pixels-num_pixels+k])
		if (x < temp_x_min):
			temp_x_min = x
		if (x > temp_x_max):
			temp_x_max = x
		if (y < temp_y_min):
			temp_y_min = y
		if (y > temp_y_max):
			temp_y_max = y
		size_x = temp_x_max - temp_x_min + 1
		size_y = temp_y_max - temp_y_min + 1
		picture_size_list[i,0] = size_x
		picture_size_list[i,1] = size_y
		picture_size_list[i,2] = temp_x_min
		picture_size_list[i,3] = temp_y_min
		if size_x > max_dif_x:
			max_dif_x = size_x
		if size_y > max_dif_y:
			max_dif_y = size_y
			
		x_orig = int(temp_pixels_orig.x[sum_pixels-num_pixels+k])
		y_orig = int(temp_pixels_orig.y[sum_pixels-num_pixels+k])
		if (x_orig < temp_x_min_orig):
			temp_x_min_orig = x_orig
		if (x_orig > temp_x_max_orig):
			temp_x_max_orig = x_orig
		if (y_orig < temp_y_min_orig):
			temp_y_min_orig = y
		if (y_orig > temp_y_max_orig):
			temp_y_max_orig = y_orig
		size_x_orig = temp_x_max_orig - temp_x_min_orig + 1
		size_y_orig = temp_y_max_orig - temp_y_min_orig + 1
		picture_size_list_orig[i,0] = size_x_orig
		picture_size_list_orig[i,1] = size_y_orig
		picture_size_list_orig[i,2] = temp_x_min_orig
		picture_size_list_orig[i,3] = temp_y_min_orig
		if size_x_orig > max_dif_x_orig:
			max_dif_x_orig = size_x_orig
		if size_y_orig > max_dif_y_orig:
			max_dif_y_orig = size_y_orig  

#crea le matricione finali con gli eventi (immagini) e labels (i valori MC)
#usa solo la max_dif_x e non la max_dif_y...o deve essere per forza quadrata (e perche` non prende la piu` grande delle 2 nel caso, ma la max_dif_x)? o ha semplicemente sbagliato... mi sa di bachetto...anche se forse e` voluto visto che usa max_dif_x anche quando shifta l'immagine lungo y...
#le immagini vengono shiftate in modo da risultare sempre piu` o meno centrate in un quadrato che parte da 0,0, ma secondo me sbaglia anche qui...
#...quell'aggiunta di (max_dif_x-picture_size_list[i,0])/2...
#Inoltre shifta le y di 0.5 * x (per passare da griglia con pixel esagonali a griglia con pixel quadrati) che btw implica anche da 6 vicini a 9 vicini...

sum_pixels = 0
num_pixels = 0

event_list = np.zeros((number_of_all_events, max_dif_x, max_dif_x, 1), dtype = np.float16)
event_list_orig = np.zeros((number_of_all_events, max_dif_x_orig, max_dif_x_orig, 1), dtype = np.float16)
label_list = np.zeros((number_of_all_events, 7))

for i in tqdm(range(0,int(number_of_all_events))):
	num_pixels = labels.nPixels[i]
	sum_pixels += num_pixels
	temp_pixels = pixels[sum_pixels-num_pixels:sum_pixels]
	temp_pixels_orig = pixels_orig[sum_pixels-num_pixels:sum_pixels]
	for k in range(0,int(temp_pixels.size/3)):
		event_list[i, int(temp_pixels.x[sum_pixels-num_pixels+k]-picture_size_list[i,2] + (max_dif_x-picture_size_list[i,0])/2),int(temp_pixels.y[sum_pixels-num_pixels+k]+0.5*temp_pixels.x[sum_pixels-num_pixels+k]-picture_size_list[i,3] + (max_dif_x-picture_size_list[i,1])/2), 0] = temp_pixels.int[sum_pixels-num_pixels+k]
		event_list_orig[i, int(temp_pixels_orig.x[sum_pixels-num_pixels+k]-picture_size_list_orig[i,2] + (max_dif_x_orig-picture_size_list_orig[i,0])/2),int(temp_pixels_orig.y[sum_pixels-num_pixels+k]-picture_size_list_orig[i,3] + (max_dif_x_orig-picture_size_list_orig[i,1])/2), 0] = temp_pixels_orig.int[sum_pixels-num_pixels+k]
				
for i in tqdm(range(0,int(number_of_all_events))):
    label_list[i,0] = labels.ID[i]
    label_list[i,1] = labels.nPixels[i]
    label_list[i,2] = labels.Energy[i]
    label_list[i,3] = labels.x[i]
    label_list[i,4] = labels.y[i]
    label_list[i,5] = labels.phi[i]
    label_list[i,6] = labels.theta[i]

'''
Scrive gli output in formato hdf5
che e` abbastanza comodo visto che prende direttamente il matricione numpy array
'''
outfile = h5.File("pixels.hdf5", "w")
outfile.create_dataset("data", data=event_list)
outfile.close()

outfile = h5.File("labels.hdf5", "w")
outfile.create_dataset("data", data=label_list)
outfile.close()

'''
Plotta 11 immagini 0, 100, 200,...,1000
forse se metto il plt.show fuori dal for li plotta tutti insieme
'''

for i in range(0,10):
	plt.imshow(event_list[i,:,:,0], interpolation="nearest")
	plt.colorbar()
	#plt.show()
	
fout = ROOT.TFile("prova.root","recreate")
h2hex=ROOT.TH2Poly("h2poly","",100,0,100,100,0,100)
h2hex.Honeycomb(0.,0.,0.5, 100 ,100);
num_pixels_root = 0
sum_pixels_root = 0
for i in range(0,10):	
	num_pixels_root = labels.nPixels[i]
	sum_pixels_root += num_pixels_root
	temp_pixels_orig2 = pixels_orig[sum_pixels_root-num_pixels_root:sum_pixels_root]
	temp_pixels_orig2 = pixels_orig[sum_pixels_root-num_pixels_root:sum_pixels_root]
	#temp_pixels_orig2.x *= 20
	#temp_pixels_orig2.y *= 20
	name = "h2poly"+str(i)
	h2hex.SetName(name)
	for k in range(0,int(temp_pixels_orig2.size/3)):                    
		h2hex.Fill(temp_pixels_orig2.x[k],temp_pixels_orig2.y[k],temp_pixels_orig2.int[k])
		print(k,temp_pixels_orig2.x[k],temp_pixels_orig2.y[k])
	h2hex.Write()	
fout.Close()	

	 
'''
Dubbi finali

line 65: perche` moltiplica le x e le y per quei fattori?
R. passa da mm e giglia esagonale a scala di interi e semiinteri schiacciando i pixel lungo y (non sono piu` esagoni regolari ma fatti di triangoli rettangoli 90 deg al centro)

line 77: pixels.x[i] o pixels.x.values[i]...che differenza c'e`? 
R. nessuna dal mio print

line 84: perche` sottrae il valori x,y minimi assoluti (tutti gli eventi)
R. credo sia inutile tanto poi immagine per immagine shifta tutto a (0,0)

line 150: forse metterei un max_dif_y al posto del secondo max_dif_x
R. forse vuole immagini quadrate, e cmq usa max_dif_x anche per shiftare lungo y le immagini

line 158: non sono sicuro che l'aggiunta di (max_dif_x-picture_size_list[i,0])/2 funzioni.
R. dai plot funziona, ma controllare la matematica
'''
