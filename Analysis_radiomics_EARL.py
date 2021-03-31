#coding: utf-8


from scipy import stats
import pandas as pad
import matplotlib.pyplot as plt
import seaborn as sns
import codecs
import time
import math
import statsmodels.formula.api as smf
import statsmodels.api as sm
import numpy as np
from scipy.stats import pearsonr

## POUR CHANGER COULEUR sns.set_palette("bright")
#custom_palette = [sns.xkcd_rgb["medium green"], sns.xkcd_rgb["pale red"], "blue", "orange", "blue","yellow", "purple"]
custom_paletteGeneral = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["pale red"], sns.xkcd_rgb["medium green"], "orange", "blue","yellow", "purple"]
sns.set_palette(custom_paletteGeneral)

##########################
### Analyse effet Algo ###
##########################


df = pad.read_excel('//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/AnalyseRadiomic-EARL_29032021.xlsx', header=0)
df = pad.read_excel('C:/Users/corrau/Documents/Etude IA Radiomics TEP/AnalyseRadiomic-IA_17032021.xlsx', header=0)
df = pad.read_excel('/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/AnalyseRadiomic-EARL_17032021.xlsx', header=0)

#df = pad.read_excel('//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/AnalyseTexture.xls', header=0)
#df = pad.read_excel('//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/AnalyseTexture2.xls', header=0)



##################################################
###  	Test distribution loi normale	   	   ###
##################################################

### avec test D’Agostino
y=3 
for i in range(118): 
	NormOrigin = stats.normaltest(df[df.columns[y]][:int(len(df)/2)])
	NormIA = stats.normaltest(df[df.columns[y]][int(len(df)/2):])
	if NormOrigin[1]> 0.05:
		print("La variable " + str(df.columns[y] + " pour les données Origin, ne suit pas une loi normale (p = " + str(NormOrigin[1])) + ")")#p value
	if NormIA[1]> 0.05:
		print("La variable " + str(df.columns[y] + " pour les données IA, ne suit pas une loi normale (p = " + str(NormIA[1])) + ")")#p value 
	y+=1

### avec test D’Agostino

### avec test shapiro
y=3 
for i in range(118): 
	NormOrigin = stats.shapiro(df[df.columns[y]][:int(len(df)/2)])
	NormIA = stats.shapiro(df[df.columns[y]][int(len(df)/2):])
	if NormOrigin[1]> 0.05:
		print("La variable " + str(df.columns[y] + " pour les données Origin, ne suit pas une loi normale (p = " + str(NormOrigin[1])) + ")")#p value
	if NormIA[1]> 0.05:
		print("La variable " + str(df.columns[y] + " pour les données IA, ne suit pas une loi normale (p = " + str(NormIA[1])) + ")")#p value 
	y+=1
### avec test shapiro




##################################################
### 	     	Toutes loc confondues	   	   ###
##################################################

####
y=3 
for i in range(118): 
	results = stats.ttest_rel(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])  ### pour test non app = stats.ttest_ind
	ax = sns.boxplot(x='Groupe', y=df[df.columns[y]], data=df, showfliers = False) 
	if results[1] < 0.05:
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Toutes loc/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/toutes_loc/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Toutes loc/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df.columns[y] + " la valeur de p concernant la différence entre les deux images est de " + str(results[1]))) #p value 
	figure.clear()
	y+=1



#### pour avoir les paramètres qui ne bougent pas ###
y=3 
for i in range(118): 
	results = stats.ttest_rel(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])  ### pour test non app = stats.ttest_ind
	ax = sns.boxplot(x='Groupe', y=df[df.columns[y]], data=df, showfliers = False) 
	if results[1] > 0.05:
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Toutes loc/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/toutes_loc/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Toutes loc/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df.columns[y] + " la valeur de p est supérieur à 0.05, et est égale à : " + str(results[1]))) #p value 
	figure.clear()
	y+=1

##################################################
### 	     		Loc par Loc	   	   		   ###
##################################################

df_liver = df[df[" labels "] == "Liver"]
df_Lung = df[df[" labels "] == "Lung"]
df_BloodPool = df[df[" labels "] == "BloodPool"]
df_Muscle = df[df[" labels "] == "Muscle"]
df_Lesion = df[df[" labels "] == "Lesion"]



##########################################
### 	    	LIVER	   	   	       ###
##########################################
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_liver[df_liver.columns[y]][:int(len(df_liver)/2)], df_liver[df_liver.columns[y]][int(len(df_liver)/2):])  ### pour test non app = stats.ttest_ind
	ax = sns.boxplot(x='Groupe', y=df_liver[df_liver.columns[y]], data=df_liver, showfliers = False) 
	if results[1] < 0.05:
		figure = ax.get_figure()
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/Figures_Analyse_Stat/liver/t-test__" + str(df.columns[y]), dpi=400)
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Liver/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Liver/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_liver.columns[y] + " la valeur de p concernant la différence entre les deux images est de " + str(results[1]))) #p value 
	figure.clear()
	y+=1

#### pour avoir les paramètres qui ne bougent pas ###
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_liver[df_liver.columns[y]][:int(len(df_liver)/2)], df_liver[df_liver.columns[y]][int(len(df_liver)/2):])  ### pour test non app = stats.ttest_ind
	if results[1] > 0.05:
		ax = sns.boxplot(x='Groupe', y=df_liver[df_liver.columns[y]], data=df_liver, showfliers = False) 
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Liver/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/Figures_Analyse_Stat/liver/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Liver/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_liver.columns[y] + " la valeur de p est supérieur à 0.05, et est égale à : " + str(results[1])))
	figure.clear()
	y+=1


##########################################
### 	    	LUNG	   	   	       ###
##########################################
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_Lung[df_Lung.columns[y]][:int(len(df_Lung)/2)], df_Lung[df_Lung.columns[y]][int(len(df_Lung)/2):])  ### pour test non app = stats.ttest_ind
	ax = sns.boxplot(x='Groupe', y=df_Lung[df_Lung.columns[y]], data=df_Lung, showfliers = False) 
	if results[1] < 0.05:
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Lung/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/lung/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Lung/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_Lung.columns[y] + " la valeur de p concernant la différence entre les deux images est de " + str(results[1]))) #p value 
	figure.clear()
	y+=1

#### pour avoir les paramètres qui ne bougent pas ###
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_Lung[df_Lung.columns[y]][:int(len(df_Lung)/2)], df_Lung[df_Lung.columns[y]][int(len(df_Lung)/2):])  ### pour test non app = stats.ttest_ind
	if results[1] > 0.05:
		ax = sns.boxplot(x='Groupe', y=df_Lung[df_Lung.columns[y]], data=df_Lung, showfliers = False) 
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Lung/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/lung/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Lung/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_Lung.columns[y] + " la valeur de p est supérieur à 0.05, et est égale à : " + str(results[1])))
	figure.clear()
	y+=1


##########################################
### 	    	Muscle	   	   	       ###
##########################################
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_Muscle[df_Muscle.columns[y]][:int(len(df_Muscle)/2)], df_Muscle[df_Muscle.columns[y]][int(len(df_Muscle)/2):])  ### pour test non app = stats.ttest_ind
	ax = sns.boxplot(x='Groupe', y=df_Muscle[df_Muscle.columns[y]], data=df_Muscle, showfliers = False) 
	if results[1] < 0.05:
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Muscle/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/Figures_Analyse_Stat/Muscle/t-test__" + str(df_Muscle.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/muscle/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Muscle/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_Muscle.columns[y] + " la valeur de p concernant la différence entre les deux images est de " + str(results[1]))) #p value 
	figure.clear()
	y+=1

#### pour avoir les paramètres qui ne bougent pas ###
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_Muscle[df_Muscle.columns[y]][:int(len(df_Muscle)/2)], df_Muscle[df_Muscle.columns[y]][int(len(df_Muscle)/2):])  ### pour test non app = stats.ttest_ind
	if results[1] > 0.05:
		ax = sns.boxplot(x='Groupe', y=df_Muscle[df_Muscle.columns[y]], data=df_Muscle, showfliers = False) 
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Muscle/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/muscle/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Muscle/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_Muscle.columns[y] + " la valeur de p est supérieur à 0.05, et est égale à : " + str(results[1])))
	figure.clear()
	y+=1



##########################################
### 	    	BloodPool  	   	       ###
##########################################
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_BloodPool[df_BloodPool.columns[y]][:int(len(df_BloodPool)/2)], df_BloodPool[df_BloodPool.columns[y]][int(len(df_BloodPool)/2):])  ### pour test non app = stats.ttest_ind
	ax = sns.boxplot(x='Groupe', y=df_BloodPool[df_BloodPool.columns[y]], data=df_BloodPool, showfliers = False) 
	if results[1] < 0.05:
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/BloodPool/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/bloodpool/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/BloodPool/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_BloodPool.columns[y] + " la valeur de p concernant la différence entre les deux images est de " + str(results[1]))) #p value 
	figure.clear()
	y+=1

#### pour avoir les paramètres qui ne bougent pas ###
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_BloodPool[df_BloodPool.columns[y]][:int(len(df_BloodPool)/2)], df_BloodPool[df_BloodPool.columns[y]][int(len(df_BloodPool)/2):])  ### pour test non app = stats.ttest_ind
	if results[1] > 0.05:
		ax = sns.boxplot(x='Groupe', y=df_BloodPool[df_BloodPool.columns[y]], data=df_BloodPool, showfliers = False) 
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/BloodPool/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/bloodpool/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/BloodPool/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_BloodPool.columns[y] + " la valeur de p est supérieur à 0.05, et est égale à : " + str(results[1])))
	figure.clear()
	y+=1



##########################################
### 	    	Lesion  	   	       ###
##########################################
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_Lesion[df_Lesion.columns[y]][:int(len(df_Lesion)/2)], df_Lesion[df_Lesion.columns[y]][int(len(df_Lesion)/2):])  ### pour test non app = stats.ttest_ind
	ax = sns.boxplot(x='Groupe', y=df_Lesion[df_Lesion.columns[y]], data=df_Lesion, showfliers = False) 
	if results[1] < 0.05:
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Lesion/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/lesion/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Lesion/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_Lesion.columns[y] + " la valeur de p concernant la différence entre les deux images est de " + str(results[1]))) #p value 
	figure.clear()
	y+=1

#### pour avoir les paramètres qui ne bougent pas ###
y=3 
for i in range(118): 
	results = stats.ttest_rel(df_Lesion[df_Lesion.columns[y]][:int(len(df_Lesion)/2)], df_Lesion[df_Lesion.columns[y]][int(len(df_Lesion)/2):])  ### pour test non app = stats.ttest_ind
	if results[1] > 0.05:
		ax = sns.boxplot(x='Groupe', y=df_Lesion[df_Lesion.columns[y]], data=df_Lesion, showfliers = False)
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Figures_Analyse_Stat/Lesion/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("C:/Users/corrau/Documents/Etude IA Radiomics TEP/Figures_ttest_Origin_IA/lesion/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_Analyse_Stat/Lesion/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)	
		print("Pour la variable " + str(df_Lesion.columns[y] + " la valeur de p est supérieur à 0.05, et est égale à : " + str(results[1])))
	figure.clear()
	y+=1


#####################################################################################
### Etude corrélation entre facteurs avant et après application de l'aglo d'IA    ###
### 	est-ce que le paramètre de radiomic et stable après IA        		   	  ###
#####################################################################################

df_liver = df[df[" labels "] == "Liver"]
df_Lung = df[df[" labels "] == "Lung"]
df_BloodPool = df[df[" labels "] == "BloodPool"]
df_Muscle = df[df[" labels "] == "Muscle"]
df_Lesion = df[df[" labels "] == "Lesion"]
### Toutes loc ###

y=3 
for i in range(118): 
	results = pearsonr(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])
	if results[1] < 0.05:
		ax = sns.regplot(x=df[df.columns[y]][:int(len(df)/2)], y=df[df.columns[y]][int(len(df)/2):], data=df)
		ax.set(xlabel="Origin " + str(df.columns[y]), ylabel="EARL " + str(df.columns[y]))
		ax.text(0.05, 0.9, "R²=" + str(round(results[0]*results[0],3)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Corrélation features/AllLoc/Effet_algo__" + str(df.columns[y]), dpi=400)
		#figure.savefig("//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/Figures_stabilité_radiomics_toutes_loc/Effet_algo__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_stabilité_radiomics_toutes_loc/Effet_algo__" + str(df.columns[y]), dpi=400)
		print("Pour la variable " + str(df.columns[y] + " la valeur de p est inférieur à 0.05, et R² est égal à : " + str(round(results[0]*results[0],3))))
	figure.clear()
	y+=1
### Toutes loc ###

### Liver ###
y=3 
for i in range(118): 
	results = pearsonr(df_liver[df_liver.columns[y]][:int(len(df_liver)/2)], df_liver[df_liver.columns[y]][int(len(df_liver)/2):])
	if results[1] < 0.05:
		ax = sns.regplot(x=df_liver[df_liver.columns[y]][:int(len(df_liver)/2)], y=df_liver[df_liver.columns[y]][int(len(df_liver)/2):], data=df)
		ax.set(xlabel="Origin " + str(df.columns[y]), ylabel="EARL " + str(df.columns[y]))
		ax.text(0.05, 0.9, "R²=" + str(round(results[0]*results[0],3)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Corrélation features/Liver/Effet_algo__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_stabilité_radiomics_liver/Effet_algo__" + str(df.columns[y]), dpi=400)
		print("Pour la variable " + str(df.columns[y] + " la valeur de p est inférieur à 0.05, et R² est égal à : " + str(round(results[0]*results[0],3))))
	figure.clear()
	y+=1
### Liver ###

### Lung ###
y=3 
for i in range(118): 
	results = pearsonr(df_Lung[df_Lung.columns[y]][:int(len(df_Lung)/2)], df_Lung[df_Lung.columns[y]][int(len(df_Lung)/2):])  ### pour test non app = stats.ttest_ind
	if results[1] < 0.05:
		ax = sns.regplot(x=df_Lung[df_Lung.columns[y]][:int(len(df_Lung)/2)], y=df_Lung[df_Lung.columns[y]][int(len(df_Lung)/2):], data=df)
		ax.set(xlabel="Origin " + str(df.columns[y]), ylabel="EARL " + str(df.columns[y]))
		ax.text(0.05, 0.9, "R²=" + str(round(results[0]*results[0],3)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Corrélation features/Lung/Effet_algo__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_stabilité_radiomics_lung/Effet_algo__" + str(df.columns[y]), dpi=400)
		print("Pour la variable " + str(df.columns[y] + " la valeur de p est inférieur à 0.05, et R² est égal à : " + str(round(results[0]*results[0],3))))
	figure.clear()
	y+=1
### Lung ###

### Muscle ###
y=3 
for i in range(118): 
	results = pearsonr(df_Muscle[df_Muscle.columns[y]][:int(len(df_Muscle)/2)], df_Muscle[df_Muscle.columns[y]][int(len(df_Muscle)/2):])  ### pour test non app = stats.ttest_ind
	if results[1] < 0.05:
		ax = sns.regplot(x=df_Muscle[df_Muscle.columns[y]][:int(len(df_Muscle)/2)], y=df_Muscle[df_Muscle.columns[y]][int(len(df_Muscle)/2):], data=df)
		ax.set(xlabel="Origin " + str(df.columns[y]), ylabel="EARL " + str(df.columns[y]))
		ax.text(0.05, 0.9, "R²=" + str(round(results[0]*results[0],3)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Corrélation features/Muscle/Effet_algo__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_stabilité_radiomics_muscle/Effet_algo__" + str(df.columns[y]), dpi=400)
		print("Pour la variable " + str(df.columns[y] + " la valeur de p est inférieur à 0.05, et R² est égal à : " + str(round(results[0]*results[0],3))))
	figure.clear()
	y+=1
### Muscle ###

### Lesion ###
y=3 
for i in range(118): 
	results = pearsonr(df_Lesion[df_Lesion.columns[y]][:int(len(df_Lesion)/2)], df_Lesion[df_Lesion.columns[y]][int(len(df_Lesion)/2):])  ### pour test non app = stats.ttest_ind
	if results[1] < 0.05:
		ax = sns.regplot(x=df_Lesion[df_Lesion.columns[y]][:int(len(df_Lesion)/2)], y=df_Lesion[df_Lesion.columns[y]][int(len(df_Lesion)/2):], data=df)
		ax.set(xlabel="Origin " + str(df.columns[y]), ylabel="EARL " + str(df.columns[y]))
		ax.text(0.05, 0.9, "R²=" + str(round(results[0]*results[0],3)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Corrélation features/Lesion/Effet_algo__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_stabilité_radiomics_lesion/Effet_algo__" + str(df.columns[y]), dpi=400)
		print("Pour la variable " + str(df.columns[y] + " la valeur de p est inférieur à 0.05, et R² est égal à : " + str(round(results[0]*results[0],3))))
	figure.clear()
	y+=1
### Lesion ###

### BloodPool ###
y=3 
for i in range(118): 
	results = pearsonr(df_BloodPool[df_BloodPool.columns[y]][:int(len(df_BloodPool)/2)], df_BloodPool[df_BloodPool.columns[y]][int(len(df_BloodPool)/2):])  ### pour test non app = stats.ttest_ind
	if results[1] < 0.05:
		ax = sns.regplot(x=df_BloodPool[df_BloodPool.columns[y]][:int(len(df_BloodPool)/2)], y=df_BloodPool[df_BloodPool.columns[y]][int(len(df_BloodPool)/2):], data=df)
		ax.set(xlabel="Origin " + str(df.columns[y]), ylabel="EARL " + str(df.columns[y]))
		ax.text(0.05, 0.9, "R²=" + str(round(results[0]*results[0],3)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
		figure = ax.get_figure()
		figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/Corrélation features/BloodPool/Effet_algo__" + str(df.columns[y]), dpi=400)
		#figure.savefig("/home/aureliencd/Documents/Baclesse_ACD/IA Radiomics/Figures_stabilité_radiomics_bloodpool/Effet_algo__" + str(df.columns[y]), dpi=400)
		print("Pour la variable " + str(df.columns[y] + " la valeur de p est inférieur à 0.05, et R² est égal à : " + str(round(results[0]*results[0],3))))
	figure.clear()
	y+=1
### BloodPool ###




















###########################################################################################
### 																   		   	        ###
### Etude Concordance Correlation Coefficient (CCC) cf article Philippe Lambin et al    ###
### 																   		   	        ###
###########################################################################################


def ccc(x,y):
	sxy = np.sum((x - x.mean())*(y - y.mean()))/x.shape[0]      
	rhoc = 2*sxy / (np.var(x) + np.var(y) + (x.mean() - y.mean())**2)      
	return rhoc  



### Toutes loc ###
### Stockage des valeurs de CCC dans des listes ###
FeatursAllLocList = []
CCCAllLocValue = []

y=3 
for i in range(118): 
	result = ccc(df[df.columns[y]][:int(len(df)/2)].to_numpy(), df[df.columns[y]][int(len(df)/2):].to_numpy())  ### pour test non app = stats.ttest_ind
	#print("Pour la variable " + str(df.columns[y]) + " le CCC est de : " + str(round(result,3)))
	FeatursAllLocList.append(str(df.columns[y]))
	CCCAllLocValue.append(round(result,3))
	if result > 0.85:
		print("Pour la variable " + str(df.columns[y]) + " le CCC est de : " + str(round(result,3)))
	y+=1

plt.figure(figsize=(5,25))
plt.axvline(0.85, 0,1, linewidth=3, color='b')
ax = sns.barplot(x=CCCAllLocValue, y=FeatursAllLocList)
ax.set(xlabel="CCC ")
plt.show()
### Toutes loc ###


### Lesion ###
FeatursLesion = []
CCCLesion = []

y=3 
for i in range(118): 
	result = ccc(df_Lesion[df_Lesion.columns[y]][:int(len(df_Lesion)/2)].to_numpy(), df_Lesion[df_Lesion.columns[y]][int(len(df_Lesion)/2):].to_numpy())  ### pour test non app = stats.ttest_ind
	#print("Pour la variable " + str(df.columns[y]) + " le CCC est de : " + str(round(result,3)))
	FeatursLesion.append(str(df_Lesion.columns[y]))
	CCCLesion.append(round(result,3))
	y+=1

plt.figure(figsize=(5,25))
plt.axvline(0.85, 0,1, linewidth=3, color='b')
ax = sns.barplot(x=CCCLesion, y=FeatursLesion)
ax.set(xlabel="CCC ")
plt.show()
### Lesion ###


### Liver ###
FeatursLiver = []
CCCLiver = []

y=3 
for i in range(118): 
	result = ccc(df_liver[df_liver.columns[y]][:int(len(df_liver)/2)].to_numpy(), df_liver[df_liver.columns[y]][int(len(df_liver)/2):].to_numpy())  ### pour test non app = stats.ttest_ind
	#print("Pour la variable " + str(df.columns[y]) + " le CCC est de : " + str(round(result,3)))
	FeatursLiver.append(str(df_liver.columns[y]))
	CCCLiver.append(round(result,3))
	y+=1

plt.figure(figsize=(5,25))
plt.axvline(0.85, 0,1, linewidth=3, color='b')
ax = sns.barplot(x=CCCLiver, y=FeatursLiver)
ax.set(xlabel="CCC ")
plt.show()
### Liver ###


### Lung ###
FeatursLung = []
CCCLung = []

y=3 
for i in range(118): 
	result = ccc(df_Lung[df_Lung.columns[y]][:int(len(df_Lung)/2)].to_numpy(), df_Lung[df_Lung.columns[y]][int(len(df_Lung)/2):].to_numpy())  ### pour test non app = stats.ttest_ind
	#print("Pour la variable " + str(df.columns[y]) + " le CCC est de : " + str(round(result,3)))
	FeatursLung.append(str(df_Lung.columns[y]))
	CCCLung.append(round(result,3))
	y+=1

plt.figure(figsize=(5,25))
plt.axvline(0.85, 0,1, linewidth=3, color='b')
ax = sns.barplot(x=CCCLung, y=FeatursLung)
ax.set(xlabel="CCC ")
plt.show()
### Lung ###


### Muscle ###
FeatursMuscle = []
CCCMuscle = []

y=3 
for i in range(118): 
	result = ccc(df_Muscle[df_Muscle.columns[y]][:int(len(df_Muscle)/2)].to_numpy(), df_Muscle[df_Muscle.columns[y]][int(len(df_Muscle)/2):].to_numpy())  ### pour test non app = stats.ttest_ind
	#print("Pour la variable " + str(df.columns[y]) + " le CCC est de : " + str(round(result,3)))
	FeatursMuscle.append(str(df_Muscle.columns[y]))
	CCCMuscle.append(round(result,3))
	y+=1

plt.figure(figsize=(5,25))
plt.axvline(0.85, 0,1, linewidth=3, color='b')
ax = sns.barplot(x=CCCMuscle, y=FeatursMuscle)
ax.set(xlabel="CCC ")
plt.show()
### Muscle ###


### BloodPool ###
FeatursBloodPool = []
CCCBloodPool = []

y=3 
for i in range(118): 
	result = ccc(df_BloodPool[df_BloodPool.columns[y]][:int(len(df_BloodPool)/2)].to_numpy(), df_BloodPool[df_BloodPool.columns[y]][int(len(df_BloodPool)/2):].to_numpy())  ### pour test non app = stats.ttest_ind
	#print("Pour la variable " + str(df.columns[y]) + " le CCC est de : " + str(round(result,3)))
	FeatursBloodPool.append(str(df_BloodPool.columns[y]))
	CCCBloodPool.append(round(result,3))
	y+=1

plt.figure(figsize=(5,25))
plt.axvline(0.85, 0,1, linewidth=3, color='b')
ax = sns.barplot(x=CCCBloodPool, y=FeatursBloodPool)
ax.set(xlabel="CCC ")
plt.show()
### BloodPool ###



### Discrimination en fonction des classes de radiomics ###

### Lésion ###
CCCLesionIntensity = CCCLesion[:9]
for elm in CCCLesion[23:41]:
	CCCLesionIntensity.append(elm)
FeatursLesionIntensity = FeatursLesion[:9]
for elm in FeatursLesion[23:41]:
	FeatursLesionIntensity.append(elm)
FeatursLesionGLRLM = FeatursLesion[79:95]
CCCLesionShape = CCCLesion[9:23]
FeatursLesionShape = FeatursLesion[9:23]
CCCLesionGLCM = CCCLesion[41:65]
FeatursLesionGLCM = FeatursLesion[41:65]
CCCLesionGLDM = CCCLesion[65:79]
FeatursLesionGLDM = FeatursLesion[65:79]
CCCLesionGLRLM = CCCLesion[79:95]
FeatursLesionGLSZM = FeatursLesion[95:111]
CCCLesionGLSZM = CCCLesion[95:111]
CCCLesionNGTDM = CCCLesion[111:-2]
FeatursLesionNGTDM = FeatursLesion[111:-2]
CCCLesionIQwavelet = CCCLesion[-2:]
FeatursLesionIQwavelet = FeatursLesion[-2:]
### Lésion ###

### Liver ###
CCCLiverIntensity = CCCLiver[:9]	
for elm in CCCLiver[23:41]:	
	CCCLiverIntensity.append(elm)
FeatursLiverIntensity = FeatursLiver[:9]	
for elm in FeatursLiver[23:41]:	
	FeatursLiverIntensity.append(elm)
FeatursLiverGLRLM = FeatursLiver[79:95]	
CCCLiverShape = CCCLiver[9:23]	
FeatursLiverShape = FeatursLiver[9:23]	
CCCLiverGLCM = CCCLiver[41:65]	
FeatursLiverGLCM = FeatursLiver[41:65]	
CCCLiverGLDM = CCCLiver[65:79]	
FeatursLiverGLDM = FeatursLiver[65:79]	
CCCLiverGLRLM = CCCLiver[79:95]	
FeatursLiverGLSZM = FeatursLiver[95:111]	
CCCLiverGLSZM = CCCLiver[95:111]	
CCCLiverNGTDM = CCCLiver[111:-2]	
FeatursLiverNGTDM = FeatursLiver[111:-2]	
CCCLiverIQwavelet = CCCLiver[-2:]	
FeatursLiverIQwavelet = FeatursLiver[-2:]	
### Liver ###

### Lung ###
CCCLungIntensity = CCCLung[:9]	
for elm in CCCLung[23:41]:	
	CCCLungIntensity.append(elm)
FeatursLungIntensity = FeatursLung[:9]	
for elm in FeatursLung[23:41]:	
	FeatursLungIntensity.append(elm)
FeatursLungGLRLM = FeatursLung[79:95]	
CCCLungShape = CCCLung[9:23]	
FeatursLungShape = FeatursLung[9:23]	
CCCLungGLCM = CCCLung[41:65]	
FeatursLungGLCM = FeatursLung[41:65]	
CCCLungGLDM = CCCLung[65:79]	
FeatursLungGLDM = FeatursLung[65:79]	
CCCLungGLRLM = CCCLung[79:95]	
FeatursLungGLSZM = FeatursLung[95:111]	
CCCLungGLSZM = CCCLung[95:111]	
CCCLungNGTDM = CCCLung[111:-2]	
FeatursLungNGTDM = FeatursLung[111:-2]	
CCCLungIQwavelet = CCCLung[-2:]	
FeatursLungIQwavelet = FeatursLung[-2:]	
### Lung ###

def saveCCCRadiomicClass(CCC, Featurs, Class, Loc):
	df_CustomPalette = pad.DataFrame(CCC, columns = ['CCC'])
	df_CustomPalette['Featurs'] = Featurs
	custom_palette = {}
	for q in set(df_CustomPalette.Featurs):
		avr = df_CustomPalette[df_CustomPalette.Featurs == q].CCC
		if avr.values < 0.85:
			custom_palette[q] = sns.xkcd_rgb["pale red"]
		else:
			custom_palette[q] = sns.xkcd_rgb["windows blue"]
	plt.figure(figsize=(5,0.36*len(Featurs)))
	plt.axvline(0.85, 0,1, linewidth=3, color='b')
	ax = sns.barplot(x=CCC, y=Featurs, palette=custom_palette)
	ax.set(xlabel="CCC ")
	ax.set_title(Class)
	figure = ax.get_figure()
	figure.savefig("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/EARL___64bins/CCC_EARL/"+ str(Loc) + "/CCC_" + str(Loc) + str(Class), dpi=400, bbox_inches='tight')
	plt.show()


saveCCCRadiomicClass(CCCLesionIntensity, FeatursLesionIntensity, "Intensity", "lesion")
saveCCCRadiomicClass(CCCLesionShape, FeatursLesionShape, "Shape", "lesion")
saveCCCRadiomicClass(CCCLesionGLCM, FeatursLesionGLCM, "GLCM", "lesion")
saveCCCRadiomicClass(CCCLesionGLDM, FeatursLesionGLDM, "GLDM", "lesion")
saveCCCRadiomicClass(CCCLesionGLSZM, FeatursLesionGLSZM, "GLSZM", "lesion")
saveCCCRadiomicClass(CCCLesionNGTDM, FeatursLesionNGTDM, "NGTDM", "lesion")
saveCCCRadiomicClass(CCCLesionIQwavelet, FeatursLesionIQwavelet, "IQwavelet", "lesion")

saveCCCRadiomicClass(CCCLiverIntensity, FeatursLiverIntensity, "Intensity", "Liver")
saveCCCRadiomicClass(CCCLiverShape, FeatursLiverShape, "Shape", "Liver")
saveCCCRadiomicClass(CCCLiverGLCM, FeatursLiverGLCM, "GLCM", "Liver")
saveCCCRadiomicClass(CCCLiverGLDM, FeatursLiverGLDM, "GLDM", "Liver")
saveCCCRadiomicClass(CCCLiverGLSZM, FeatursLiverGLSZM, "GLSZM", "Liver")
saveCCCRadiomicClass(CCCLiverNGTDM, FeatursLiverNGTDM, "NGTDM", "Liver")
saveCCCRadiomicClass(CCCLiverIQwavelet, FeatursLiverIQwavelet, "IQwavelet", "Liver")

saveCCCRadiomicClass(CCCLungIntensity, FeatursLungIntensity, "Intensity", "Lung")
saveCCCRadiomicClass(CCCLungShape, FeatursLungShape, "Shape", "Lung")
saveCCCRadiomicClass(CCCLungGLCM, FeatursLungGLCM, "GLCM", "Lung")
saveCCCRadiomicClass(CCCLungGLDM, FeatursLungGLDM, "GLDM", "Lung")
saveCCCRadiomicClass(CCCLungGLSZM, FeatursLungGLSZM, "GLSZM", "Lung")
saveCCCRadiomicClass(CCCLungNGTDM, FeatursLungNGTDM, "NGTDM", "Lung")
saveCCCRadiomicClass(CCCLungIQwavelet, FeatursLungIQwavelet, "IQwavelet", "Lung")




















### Graph de variation du COV en fonction du poids ###
ax = sns.boxplot(x='Groupe', y="COV", data=df)
figure= ax.get_figure()
figure.savefig("//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/Figures/t-test_Groupe_COV", dpi=400)

ax = sns.boxplot(x='Groupe', y="COVnorm", data=df)
figure= ax.get_figure()
figure.savefig("//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/Figures/t-test_Groupe_COVnorm", dpi=400)


##################################################
###				 Anayse t-test   		       ###
##################################################

df = pad.read_excel('//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/AnalyseImcIA.xlsx', sheet_name='ACD', usecols="A:P", nrows=113, header=0)
print(stats.ttest_rel(df['COVorignine'], df['COV-IA']))
print(stats.ttest_rel(df['COVOrigni-norm'], df['COV-IA-norm']))


#################################################
### Etude corrélation entre facteurs et poids ###
#################################################

df = pad.read_excel('//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/AnalyseImcIA.xlsx', sheet_name='ACD', usecols="A:P", nrows=113, header=0)

### Corrélation de person ###
print(pearsonr(df['COVorignine'], df["poids"]))
print(pearsonr(df['COV-IA'], df["poids"]))
print(pearsonr(df['COVOrigni-norm'], df["poids"]))
print(pearsonr(df['COV-IA-norm'], df["poids"]))

### Test d'autres corrélations ###
from scipy.stats import spearmanr
print(spearmanr(df['COVorignine'], df["poids"]))
print(spearmanr(df['COV-IA'], df["poids"]))

### Test d'autres corrélations ###
from scipy.stats import kendalltau
print(kendalltau(df['COVorignine'], df["poids"]))
print(kendalltau(df['COV-IA'], df["poids"]))

ax = sns.regplot(x='poids', y='COVorignine', data=df)
figure= ax.get_figure()
figure.savefig("//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/Figures/corr-COVorignine-poids", dpi=400)

ax = sns.regplot(x='poids', y='COV-IA', data=df)
figure= ax.get_figure()
figure.savefig("//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/Figures/corr-COV-IA-poids", dpi=400)

ax = sns.regplot(x='poids', y='COVOrigni-norm', data=df)
figure= ax.get_figure()
figure.savefig("//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/Figures/COVOrigni-norm-poids", dpi=400)

ax = sns.regplot(x='poids', y='COV-IA-norm', data=df)
figure= ax.get_figure()
figure.savefig("//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/Figures/COV-IA-norm-poids", dpi=400)ax = sns.regplot(x='poids', y='COV-IA', data=df)

### Combinaison des deux ###
sns.regplot(x='poids', y='COVorignine', data=df)
sns.regplot(x='poids', y='COV-IA', data=df)

sns.regplot(x='poids', y='COVOrigni-norm', data=df)
sns.regplot(x='poids', y='COV-IA-norm', data=df)



### Test de différence significative entre les corrélations ###

"""
Functions for calculating the statistical significant differences between two dependent or independent correlation
coefficients.
The Fisher and Steiger method is adopted from the R package http://personality-project.org/r/html/paired.r.html
and is described in detail in the book 'Statistical Methods for Psychology'
The Zou method is adopted from http://seriousstats.wordpress.com/2012/02/05/comparing-correlations/
Credit goes to the authors of above mentioned packages!
Author: Philipp Singer (www.philippsinger.info)
"""

from __future__ import division

__author__ = 'psinger'

import numpy as np
from scipy.stats import t, norm
from math import atanh, pow
from numpy import tanh

def rz_ci(r, n, conf_level = 0.95):
	zr_se = pow(1/(n - 3), .5)
	moe = norm.ppf(1 - (1 - conf_level)/float(2)) * zr_se
	zu = atanh(r) + moe
	zl = atanh(r) - moe
	return tanh((zl, zu))

def rho_rxy_rxz(rxy, rxz, ryz):
	num = (ryz-1/2.*rxy*rxz)*(1-pow(rxy,2)-pow(rxz,2)-pow(ryz,2))+pow(ryz,3)
	den = (1 - pow(rxy,2)) * (1 - pow(rxz,2))
	return num/float(den)

def dependent_corr(xy, xz, yz, n, twotailed=True, conf_level=0.95, method='steiger'):
	"""
	Calculates the statistic significance between two dependent correlation coefficients
	@param xy: correlation coefficient between x and y
	@param xz: correlation coefficient between x and z
	@param yz: correlation coefficient between y and z
	@param n: number of elements in x, y and z
	@param twotailed: whether to calculate a one or two tailed test, only works for 'steiger' method
	@param conf_level: confidence level, only works for 'zou' method
	@param method: defines the method uses, 'steiger' or 'zou'
	@return: t and p-val
	"""
	if method == 'steiger':
		d = xy - xz
		determin = 1 - xy * xy - xz * xz - yz * yz + 2 * xy * xz * yz
		av = (xy + xz)/2
		cube = (1 - yz) * (1 - yz) * (1 - yz)

		t2 = d * np.sqrt((n - 1) * (1 + yz)/(((2 * (n - 1)/(n - 3)) * determin + av * av * cube)))
		p = 1 - t.cdf(abs(t2), n - 3)

		if twotailed:
			p *= 2

		return t2, p
	elif method == 'zou':
		L1 = rz_ci(xy, n, conf_level=conf_level)[0]
		U1 = rz_ci(xy, n, conf_level=conf_level)[1]
		L2 = rz_ci(xz, n, conf_level=conf_level)[0]
		U2 = rz_ci(xz, n, conf_level=conf_level)[1]
		rho_r12_r13 = rho_rxy_rxz(xy, xz, yz)
		lower = xy - xz - pow((pow((xy - L1), 2) + pow((U2 - xz), 2) - 2 * rho_r12_r13 * (xy - L1) * (U2 - xz)), 0.5)
		upper = xy - xz + pow((pow((U1 - xy), 2) + pow((xz - L2), 2) - 2 * rho_r12_r13 * (U1 - xy) * (xz - L2)), 0.5)
		return lower, upper
	else:
		raise Exception('Wrong method!')

def independent_corr(xy, ab, n, n2 = None, twotailed=True, conf_level=0.95, method='fisher'):
	"""
	Calculates the statistic significance between two independent correlation coefficients
	@param xy: correlation coefficient between x and y
	@param xz: correlation coefficient between a and b
	@param n: number of elements in xy
	@param n2: number of elements in ab (if distinct from n)
	@param twotailed: whether to calculate a one or two tailed test, only works for 'fisher' method
	@param conf_level: confidence level, only works for 'zou' method
	@param method: defines the method uses, 'fisher' or 'zou'
	@return: z and p-val
	"""

	if method == 'fisher':
		xy_z = 0.5 * np.log((1 + xy)/(1 - xy))
		ab_z = 0.5 * np.log((1 + ab)/(1 - ab))
		if n2 is None:
			n2 = n

		se_diff_r = np.sqrt(1/(n - 3) + 1/(n2 - 3))
		diff = xy_z - ab_z
		z = abs(diff / se_diff_r)
		p = (1 - norm.cdf(z))
		if twotailed:
			p *= 2

		return z, p
	elif method == 'zou':
		L1 = rz_ci(xy, n, conf_level=conf_level)[0]
		U1 = rz_ci(xy, n, conf_level=conf_level)[1]
		L2 = rz_ci(ab, n2, conf_level=conf_level)[0]
		U2 = rz_ci(ab, n2, conf_level=conf_level)[1]
		lower = xy - ab - pow((pow((xy - L1), 2) + pow((U2 - ab), 2)), 0.5)
		upper = xy - ab + pow((pow((U1 - xy), 2) + pow((ab - L2), 2)), 0.5)
		return lower, upper
	else:
		raise Exception('Wrong method!')

print(dependent_corr(.37, .24, .93, 113, method='steiger'))
print(independent_corr(0.37 , 0.24, 113, 113, method='fisher'))
