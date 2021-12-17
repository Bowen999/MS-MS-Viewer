import xml.etree.ElementTree as ET
import sys
from base64 import b64decode
from array import array
import matplotlib.pyplot as plt
import gzip

#1. xml File Analysis
if len(sys.argv)<4:
    print("Please input the file name, scan number, and peptide sequence")
    sys.exit(0)

#1.1 unzip and read the xml file
try:
    zipfile = gzip.open(sys.argv[1])
except FileNotFoundError:
    print("No such file.")
    sys.exit(0)
file = ET.parse(zipfile)
root = file.getroot()
ns = "{http://sashimi.sourceforge.net/schema/}"

#1.2 Find the scan num and peaks
scan_num = sys.argv[2]
try:
    for scans in root.findall(ns+"scan"):
        if scans.attrib["num"] == scan_num:
            scan = scans
    peakselt = scan.find(ns+"peaks")

#1.3 Read the peaks
    peaks = array('f', b64decode(peakselt.text))
    if sys.byteorder != 'big':
        peaks.byteswap()
        mzs = peaks[::2]
        ints = peaks[1::2]
    dict_peak=dict(zip(mzs,ints))
    threshold_ints = (max(ints))*0.05    #Set up a threshold to remove insignificant values
except NameError:
    print("Sorry, the scan number you input was not found.")
    sys.exit(0)

#1.4 Plot
plt.stem(mzs,ints,"grey",markerfmt=' ', label="peaks", basefmt="k")            



#2. Calculation of b_ions and y_ions
#2.1 Input the peptide sequence and molecular weight
pep_seq = sys.argv[3].upper() 
mw = {'A': 71.04,  'C': 103.01, 'D': 115.03, 'E': 129.04,\
      'F': 147.07, 'G': 57.02,  'H': 137.06, 'I': 113.08,\
      'K': 128.09, 'L': 113.08, 'M': 131.04, 'N': 114.04,\
      'P': 97.05,  'Q': 128.06, 'R': 156.10, 'S': 87.03,\
      'T': 101.05, 'V': 99.07,  'W': 186.08, 'Y': 163.06}
#Check input sequence
for word in pep_seq:
    if mw.get(word) == None:
        print("Please input the correct peptide sequence.")
        sys.exit(0)

# 2.2 Split and Restructure Sequence   
b_ions=[]
y_ions_rev=[]    
pep = list(pep_seq)    #split sequence    

#bi = from 1 to i, yi = Count i from back to front
for i in range(len(pep)):
    b_ions.append(pep[0:i+1])
    y_ions_rev.append(pep[i:])
    y_ions = y_ions_rev[::-1]  
#print(b_ions,"\n",y_ions)


# 2.3 Calculate the sum of the molecular weights of each b and y 
#     and output to the dict
b_list = []
y_list = []
dict_b = {}
dict_y = {}
n = 1
for b in b_ions:
    b_values = 1
    for ion in b:
        b_values += mw[ion]
    dict_b["b"+str(n)] = b_values
    n += 1
    b_values = 0
    
n = 1
for y in y_ions:
    y_values = 19
    for ion in y:
        y_values += mw[ion]
    dict_y["y"+str(n)] = y_values
    n += 1
    y_values = 0


#3. Filter the Matching Values         
matched_b = {}
matched_y = {}
#use ints of each ions to minus ints of peaks
#if Difference/ints < 0.0001, and ints > threshold, put it in the dict 
for mz in mzs:
    ints_value = dict_peak.get(mz)
    for b in dict_b.keys():
        diff1 = abs((dict_b[b] - mz)/mz)
        if diff1 < 0.001:
            if dict_peak[mz]>threshold_ints:
                matched_b[b]=mz,dict_peak[mz]
                print(b,":","mz is",str(mz)+"; ints is",ints_value)

    for y in dict_y.keys():
        diff2 = abs((dict_y[y] - mz)/mz)
        if diff2 < 0.0001:
            if dict_peak[mz]>threshold_ints:
                matched_y[y]=mz,dict_peak[mz]
                print(y,":","mz is",str(mz)+"; ints is",ints_value)
print("*******SUMMARY*******")
print("The scan number is", "\033[1;35m %s \033[0m"%str(scan_num)+\
      "\nThe peptide sequence is","\033[1;35m %s \033[0m"%pep_seq)
print("There are\033[1;35m %d \033[0m"%len(matched_b),"matched b ions, and\033[1;35m %d \033[0m"%len(matched_y),"matched y ions.")



#4. Plot Figures 
#4.1 Set plot rules
def plotstem(ions,c,m,n):
    xaxis=list(item[1][0] for item in ions.items()) 
    yaxis=list(item[1][1] for item in ions.items())
    plt.stem(xaxis, yaxis, linefmt=c, markerfmt=m, label=n, basefmt="k")
    for ion in ions.keys():
        plt.annotate(text=ion, xy =(ions.get(ion)[0],ions.get(ion)[1]),\
                     xytext=(ions.get(ion)[0],ions.get(ion)[1]+(max(ints)*0.03)),color=c)# Set the coordinates of the annotation according to the maximum ints

#.4.2 Plot
try:
    plotstem(matched_b,"tomato", " ", "b ions")
except ValueError:
    print("No matched b ion.")
try:
    plotstem(matched_y,"forestgreen", " ", "y ions")
except ValueError:
    print("No matched y ion.")

plt.title("MS/MS Viewer", fontweight="bold", color="red")
plt.xlabel(r"$m/z$")
plt.ylabel(r"$Intensity$")
plt.legend()
plt.show()
