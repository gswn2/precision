import os
import matplotlib.pyplot as plt
from statistics import mean
import scipy.interpolate as i
import numpy as np
import glob

#Recuperation des données
tailles=[]
val=800
while (val <= 12800):
    tailles.append(val)
    val*=2

results=[]

vectorisations=["AVX512","AVX2","AVX","SSE4_2"]
precisions=["s","d"]
operations=["potrf","trsm","syrk","gemm"]

coeurs=["1"]
fichiers = os.listdir("./")
for c in coeurs:
    resultats={"s":{"potrf":{"AVX512":[],"AVX2":[],"AVX":[],"SSE4_2":[]},
                "trsm":{"AVX512":[],"AVX2":[],"AVX":[],"SSE4_2":[]},
                "gemm":{"AVX512":[],"AVX2":[],"AVX":[],"SSE4_2":[]},
                "syrk":{"AVX512":[],"AVX2":[],"AVX":[],"SSE4_2":[]}},
           "d":{"potrf":{"AVX512":[],"AVX2":[],"AVX":[],"SSE4_2":[]},
                "trsm":{"AVX512":[],"AVX2":[],"AVX":[],"SSE4_2":[]},
                "gemm":{"AVX512":[],"AVX2":[],"AVX":[],"SSE4_2":[]},
                "syrk":{"AVX512":[],"AVX2":[],"AVX":[],"SSE4_2":[]}}}
    for p in precisions:
        for op in operations:
            for v in vectorisations:
                pattern='test_chameleon_'+p+'testing_'+v+'_'+op+'*'+c+"_20"
                for f in glob.glob(pattern):
                    file=open(f,'r')
                    size=f.split("_")[-3]
                    lignes=file.readlines()
                    t=[part.strip() for part in lignes[56].strip("| ").strip().split("|") if part.strip()][1]
                    time = float(t)
                    ligne_GFLOPS = mean([float(l.split()[-1]) * float(l.split()[-2]) for l in lignes[7:26]])         
                    if (int(c) <=18):
                        lignePKG = float([part.strip() for part in lignes[77].strip("| ").strip().split("|") if part.strip()][3])
                        ligneDRAM = float([part.strip() for part in lignes[81].strip("| ").strip().split("|") if part.strip()][3])
                    else:
                        lignePKG = float([part.strip() for part in lignes[77].strip("| ").strip().split("|") if part.strip()][1])
                        ligneDRAM = float([part.strip() for part in lignes[81].strip("| ").strip().split("|") if part.strip()][1])
                    energie = lignePKG+ligneDRAM
                    W = energie/time
                    resultats[p][op][v].append([lignePKG,ligneDRAM,ligne_GFLOPS,
                                                ligne_GFLOPS*20/energie,energie,int(size),W,time])

                resultats[p][op][v]=sorted(resultats[p][op][v],key=lambda x: x[5])
    results.append(resultats)

i=0
for resultats in results:
    for op in operations:
        p='simple'
        plt.figure(figsize=(10, 5))
        plt.plot(tailles, [ligne[3] for ligne in resultats['s'][op]['AVX512']], label='AVX512')
        plt.plot(tailles, [ligne[3] for ligne in resultats['s'][op]['AVX2']], label='AVX2')
        plt.plot(tailles, [ligne[3] for ligne in resultats['s'][op]['AVX']], label='AVX')
        plt.plot(tailles, [ligne[3] for ligne in resultats['s'][op]['SSE4_2']], label='SSE4_2')
        plt.xlabel('Tailles')
        plt.ylabel('Gflops/J')
        plt.title(op+':GFLOPS par J en fonction de la taille de la matrice:'+p+' precision')
        plt.legend()
        plt.grid(True)
        plt.savefig('graphs/'+op+":graphe_"+p+"_Gflop_J_"+coeurs[i])
        plt.close()

        p='double'
        plt.figure(figsize=(10, 5))
        plt.plot(tailles, [ligne[3] for ligne in resultats['d'][op]['AVX512']], label='AVX512')
        plt.plot(tailles, [ligne[3] for ligne in resultats['d'][op]['AVX2']], label='AVX2')
        plt.plot(tailles, [ligne[3] for ligne in resultats['d'][op]['AVX']], label='AVX')
        plt.plot(tailles, [ligne[3] for ligne in resultats['d'][op]['SSE4_2']], label='SSE4_2')
        plt.xlabel('Tailles')
        plt.ylabel('GFLOPS/J')
        plt.title(op+':GFLOPS par J en fonction de la taille de la matrice:'+p+' precision')
        plt.legend()
        plt.grid(True)
        plt.savefig('graphs/'+op+":graphe_"+p+"_Gflop_J"+coeurs[i])
        plt.close()

        ##graphe de J par seconde
        p='simple'
        plt.figure(figsize=(10, 5))
        plt.plot(tailles, [ligne[-2] for ligne in resultats['s'][op]['AVX2']], label='AVX2')
        plt.plot(tailles, [ligne[-2] for ligne in resultats['s'][op]['AVX']], label='AVX')
        plt.plot(tailles, [ligne[-2] for ligne in resultats['s'][op]['AVX512']], label='AVX512')
        plt.plot(tailles, [ligne[-2] for ligne in resultats['s'][op]['SSE4_2']], label='SSE4_2')
        plt.xlabel('Tailles')
        plt.ylabel('Puissance: W')
        plt.title(op+':W en fonction de la taille de la matrice:'+p+' precision')
        plt.legend()
        plt.grid(True)
        plt.savefig('graphs/'+op+":graphe_"+p+"_W"+coeurs[i])
        plt.close()
    #
        p='double'
        plt.figure(figsize=(10, 5))
        plt.plot(tailles, [ligne[-2] for ligne in resultats['d'][op]['AVX512']], label='AVX512')
        plt.plot(tailles, [ligne[-2] for ligne in resultats['d'][op]['AVX2']], label='AVX2')
        plt.plot(tailles, [ligne[-2] for ligne in resultats['d'][op]['AVX']], label='AVX')
        plt.plot(tailles, [ligne[-2] for ligne in resultats['d'][op]['SSE4_2']], label='SSE4_2')
        plt.xlabel('Tailles')
        plt.ylabel('Puissance: W')
        plt.title(op+': W  en fonction de la taille de la matrice:'+p+' precision')
        plt.legend()
        plt.grid(True)
        plt.savefig('graphs/'+op+":graphe_"+p+"_W"+coeurs[i])
        plt.close()
      
         ##graphe conso en J
        p='simple'
        plt.figure(figsize=(10, 5))
        plt.plot(tailles, [ligne[4] for ligne in resultats['s'][op]['AVX2']], label='AVX2')
        plt.plot(tailles, [ligne[4] for ligne in resultats['s'][op]['AVX']], label='AVX')
        plt.plot(tailles, [ligne[4] for ligne in resultats['s'][op]['AVX512']], label='AVX512')
        plt.plot(tailles, [ligne[4] for ligne in resultats['s'][op]['SSE4_2']], label='SSE4_2')
        plt.xlabel('Tailles')
        plt.ylabel('J')
        plt.title(op+':J consommé en fonction de la taille de la matrice:'+p+' precision')
        plt.legend()
        plt.grid(True)
        plt.savefig('graphs/'+op+":graphe_"+p+"_J"+coeurs[i])
        plt.close()
    #
        p='double'
        plt.figure(figsize=(10, 5))
        plt.plot(tailles, [ligne[4] for ligne in resultats['d'][op]['AVX512']], label='AVX512')
        plt.plot(tailles, [ligne[4] for ligne in resultats['d'][op]['AVX2']], label='AVX2')
        plt.plot(tailles, [ligne[4] for ligne in resultats['d'][op]['AVX']], label='AVX')
        plt.plot(tailles, [ligne[4] for ligne in resultats['d'][op]['SSE4_2']], label='SSE4_2')
        plt.xlabel('Tailles')
        plt.ylabel('J')
        plt.title(op+': J consommé en fonction de la taille de la matrice:'+p+' precision')
        plt.legend()
        plt.grid(True)
        plt.savefig('graphs/'+op+":graphe_"+p+"_J"+coeurs[i])
        plt.close()

        plt.figure(figsize=(10, 5))
        plt.plot(tailles, [ligne[-2] for ligne in resultats['s']['gemm']['AVX512']], label='simple precision')
        plt.plot(tailles, [ligne[-2] for ligne in resultats['d']['gemm']['AVX512']], label='double precision')
        plt.xlabel('Tailles')
        plt.ylabel('Puisance (W)')
        plt.title('gemm: : Puisance (W) en fonction de la taille de la matrice: simple precision')
        plt.legend()
        plt.grid(True)
        plt.savefig('graphs/Double vs Simple')
        plt.close()

        plt.figure(figsize=(10, 5))
        plt.plot(tailles, [ligne[3] for ligne in resultats['s']['gemm']['AVX512']], label='simple precision')
        plt.plot(tailles, [ligne[3] for ligne in resultats['d']['gemm']['AVX512']], label='double precision')
        plt.xlabel('Tailles')
        plt.ylabel('Puisance (W)')
        plt.title('gemm: : Puisance (W) en fonction de la taille de la matrice: simple precision')
        plt.legend()
        plt.grid(True)
        plt.savefig('graphs/double vs simple')
        plt.close()
        
    i+=1

#Puissance :













