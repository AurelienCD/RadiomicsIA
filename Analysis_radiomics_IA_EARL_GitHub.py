# coding: utf-8

# Statistical analyses of IA and EARL algorithms impact on radiomics in PET imaging
# Author : Aurélien Corroyer-Dulmont
# Version : 20 January 2021



from scipy import stats
import pandas as pad
import matplotlib.pyplot as plt
import seaborn as sns
import codecs
import math
import statsmodels.formula.api as smf
import statsmodels.api as sm
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import wilcoxon
import statistics

#custom_palette = [sns.xkcd_rgb["medium green"], sns.xkcd_rgb["pale red"], "blue", "orange", "blue","yellow", "purple"]
custom_paletteGeneral = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["pale red"], sns.xkcd_rgb["medium green"], "orange", "blue","yellow", "purple"]
sns.set_palette(custom_paletteGeneral)



########################################
###        Path and dataframe        ###
########################################

Path = 'xxx/'

df = pad.read_excel("xxx/AnalyseRadiomic-IA.xlsx", header=0)
df = pad.read_excel("xxx/AnalyseRadiomic-EARL.xlsx", header=0)
dfSUVMaxVolume = pad.read_excel("xxx/Analyse_SUVMax_Volume.xlsx", header=0)


### To get the values for each location ###
df_liver = df[df[" labels "] == "Liver"]
df_Lung = df[df[" labels "] == "Lung"]
df_BloodPool = df[df[" labels "] == "BloodPool"]
df_Muscle = df[df[" labels "] == "Muscle"]
df_Lesion = df[df[" labels "] == "Lesion"]







																			#####################################
																			### Analysis of Algorithm effect  ###
																			#####################################



### Normality distribution test ###

def normalityDistributionTestShapiro(df):

	### avec test Shapiro
	y=3 
	for i in range(104): 
		NormOrigin = stats.shapiro(df[df.columns[y]][:int(len(df)/2)])
		NormIA = stats.shapiro(df[df.columns[y]][int(len(df)/2):])
		if NormOrigin[1]> 0.05:
			print("La variable " + str(df.columns[y] + " pour les données Origin, ne suit pas une loi normale (p = " + str(NormOrigin[1])) + ")")
		if NormIA[1]> 0.05:
			print("La variable " + str(df.columns[y] + " pour les données IA, ne suit pas une loi normale (p = " + str(NormIA[1])) + ")")
		y+=1
	### avec test Shapiro

normalityDistributionTestShapiro(df)
normalityDistributionTestShapiro(df_liver)
normalityDistributionTestShapiro(df_Lung)
normalityDistributionTestShapiro(df_Muscle)
normalityDistributionTestShapiro(df_Lesion)
normalityDistributionTestShapiro(df_BloodPool)


def normalityDistributionTestAgostino(df):

	### avec test D’Agostino
	y=3 
	for i in range(104): 
		NormOrigin = stats.normaltest(df[df.columns[y]][:int(len(df)/2)])
		NormIA = stats.normaltest(df[df.columns[y]][int(len(df)/2):])
		if NormOrigin[1]> 0.05:
			print("La variable " + str(df.columns[y] + " pour les données Origin, ne suit pas une loi normale (p = " + str(NormOrigin[1])) + ")")
		if NormIA[1]> 0.05:
			print("La variable " + str(df.columns[y] + " pour les données IA, ne suit pas une loi normale (p = " + str(NormIA[1])) + ")")
		y+=1
	### avec test D’Agostino

normalityDistributionTestAgostino(df)
normalityDistributionTestAgostino(df_liver)
normalityDistributionTestAgostino(df_Lung)
normalityDistributionTestAgostino(df_Muscle)
normalityDistributionTestAgostino(df_Lesion)
normalityDistributionTestAgostino(df_BloodPool)

## The distribution is not following a normal distribution, wilcoxon analysis was chosen

### End of Normality distribution test ###


### Wilcoxon analysis ###

def wilcoxonTestSigni(df):
	y=3 
	for i in range(104): 
		wilcoxonResults = wilcoxon(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])
		if wilcoxonResults[1]< 0.05:
			print("La variable " + str(df.columns[y] + " est différent significativement (p = " + str(wilcoxonResults[1])) + ")")
		y+=1

def wilcoxonTestNoSigni(df):
	y=3 
	for i in range(104): 
		wilcoxonResults = wilcoxon(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])
		if wilcoxonResults[1]> 0.05:
			print("La variable " + str(df.columns[y] + " n'est pas différent significativement (p = " + str(wilcoxonResults[1])) + ")")
		y+=1

wilcoxonTestSigni(df_Lesion)
wilcoxonTestNoSigni(df_Lesion)


def savewilcoxonTestAnalysis(df, Loc):
	
	### significant
	y=3 
	for i in range(104): 
		results = wilcoxon(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])
		ax = sns.boxplot(x='Groupe', y=df[df.columns[y]], data=df, showfliers = False) 
		if results[1] < 0.05:
			figure = ax.get_figure()
			figure.savefig(Path + str(Loc) +"/Wilcoxon_test__" + str(df.columns[y]), dpi=400)
			print("Variable " + str(df.columns[y] + " p = " + str(round(results[1],4)))) 
		figure.clear()
		y+=1

	### non significant
	y=3 
	for i in range(104): 
		results = wilcoxon(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])  
		ax = sns.boxplot(x='Groupe', y=df[df.columns[y]], data=df, showfliers = False) 
		if results[1] > 0.05:
			figure = ax.get_figure()
			figure.savefig(Path + str(Loc) +"/PAS SIGNIFICATIF/Wilcoxon___" + str(df.columns[y]), dpi=400)
			print("Variable " + str(df.columns[y] + " p = " + str(round(results[1],4))))
		figure.clear()
		y+=1

savewilcoxonTestAnalysis(df_Lesion, "lesion")

### End of Wilcoxon analysis ###


### To get Median, Min and Max of the difference between original and IA radiomics in lesion ###

def enonceMedianMinAndMax(df):
	
	medianDif = []
	minDif = []
	maxDif = []
	y=3 
	for i in range(104): 
		dfOrigin = df[df.columns[y]][:int(len(df)/2)]
		dfIA = df[df.columns[y]][int(len(df)/2):]
		dfOrigin.reset_index(drop=True, inplace=True)
		dfIA.reset_index(drop=True, inplace=True)
		dfDif = dfIA - dfOrigin 
		medianDif.append(statistics.median(dfDif))
		minDif.append(min(dfDif))
		maxDif.append(max(dfDif))
		y+=1

	y=3
	for i in range(104):
		print(str(round(medianDif[i],2)) + " [" + str(round(minDif[i],2)) + " ; " + str(round(maxDif[i],2)) + "]")
		y+=1

enonceMedianMinAndMax(df_Lesion)

### End of To get Median, Min and Max of the difference between original and IA radiomics in lesion ###







									#############################################################################################
									### 	Correlation study between features values before and after algorithm application  ###
									### 	which radiomics features are stable after algorithm application        		   	  ###
									#############################################################################################

									#####################################################################################
									### Etude corrélation entre facteurs avant et après application de l'aglo d'IA    ###
									### 	est-ce que le paramètre de radiomic et stable après IA        		   	  ###
									#####################################################################################


### Pearson correlation analysis ###

def saveCorrelationFeatures(df, Loc):
	y=3 
	for i in range(104): 
		results = pearsonr(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])
		if results[1] < 0.05:
			ax = sns.regplot(x=df[df.columns[y]][:int(len(df)/2)], y=df[df.columns[y]][int(len(df)/2):], data=df)
			ax.set(xlabel="Original " + str(df.columns[y]), ylabel="IA " + str(df.columns[y]))
			ax.text(0.05, 0.9, "R²=" + str(round(results[0]*results[0],3)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
			ax.axis('equal')
			figure = ax.get_figure()
			figure.savefig(Path + "Corrélation features2/Figures_stabilité_radiomics_" + str(Loc) +"/Effet_algo__" + str(df.columns[y]), dpi=400)
			print("Pour la variable " + str(df.columns[y] + " la valeur de p est inférieur à 0.05, et R² est égal à : " + str(round(results[0]*results[0],3))))
		figure.clear()
		y+=1

saveCorrelationFeatures(df, "toutes_loc")
saveCorrelationFeatures(df_liver, "liver")
saveCorrelationFeatures(df_Lung, "lung")
saveCorrelationFeatures(df_Muscle, "muscle")
saveCorrelationFeatures(df_Lesion, "lesion")
saveCorrelationFeatures(df_BloodPool, "bloodpool")


### To save on the same graph the two correlation curves ###

df1 = pad.read_excel("xxx/AnalyseRadiomic-IA.xlsx", header=0)
df2 = pad.read_excel("xxx/AnalyseRadiomic-EARL.xlsx", header=0)
df_lesion_IA = df1[df1[" labels "] == "Lesion"]
df_lesion_EARL = df2[df2[" labels "] == "Lesion"]


def saveBothCorrelationFeatures(df_lesion_IA, df_lesion_EARL):
	y=3 
	for i in range(104): 
		ax = sns.regplot(x=df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], y=df_lesion_IA[df_lesion_IA.columns[y]][int(len(df_lesion_IA)/2):], data=df_lesion_IA)
		ax = sns.regplot(x=df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], y=df_lesion_EARL[df_lesion_EARL.columns[y]][int(len(df_lesion_EARL)/2):], data=df_lesion_EARL)
		ax.set(xlabel="Original " + str(df_lesion_IA.columns[y]), ylabel="Post-processing " + str(df_lesion_IA.columns[y]))
		figure = ax.get_figure()
		figure.savefig(Path + str(df.columns[y]), dpi=400)
		figure.clear()
		y+=1

saveBothCorrelationFeatures(df_lesion_IA, df_lesion_EARL)


y=6 #SUVmean
resultsIA = pearsonr(df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], df_liver_IA[df_liver_IA.columns[y]][int(len(df_liver_IA)/2):])
resultsEARL = pearsonr(df_liver_EARL[df_liver_EARL.columns[y]][:int(len(df_liver_EARL)/2)], df_liver_EARL[df_liver_EARL.columns[y]][int(len(df_liver_EARL)/2):])

ax = sns.regplot(x=df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], y=df_liver_IA[df_liver_IA.columns[y]][int(len(df_liver_IA)/2):], data=df_liver_IA)
ax = sns.regplot(x=df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], y=df_liver_EARL[df_liver_EARL.columns[y]][int(len(df_liver_EARL)/2):], data=df_liver_EARL)
ax.text(0.5, 0.35, "IA  (R²=" + str(round(resultsIA[0]*resultsIA[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["windows blue"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(0.5, 0.27, "EARL  (R²=" + str(round(resultsEARL[0]*resultsEARL[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["pale red"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.set(xlabel="Original" + str(df_liver_IA.columns[y]), ylabel="Post-processing" + str(df_liver_IA.columns[y]))
plt.xlim(1.5, 4)
plt.ylim(1.5, 4)
plt.show()

resultsIA = pearsonr(df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], df_lesion_IA[df_lesion_IA.columns[y]][int(len(df_lesion_IA)/2):])
resultsEARL = pearsonr(df_lesion_EARL[df_lesion_EARL.columns[y]][:int(len(df_lesion_EARL)/2)], df_lesion_EARL[df_lesion_EARL.columns[y]][int(len(df_lesion_EARL)/2):])

ax = sns.regplot(x=df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], y=df_lesion_IA[df_lesion_IA.columns[y]][int(len(df_lesion_IA)/2):], data=df_lesion_IA)
ax = sns.regplot(x=df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], y=df_lesion_EARL[df_lesion_EARL.columns[y]][int(len(df_lesion_EARL)/2):], data=df_lesion_EARL)
ax.text(0.5, 0.35, "IA  (R²=" + str(round(resultsIA[0]*resultsIA[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["windows blue"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(0.5, 0.27, "EARL  (R²=" + str(round(resultsEARL[0]*resultsEARL[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["pale red"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(0.7, 0.85, "y=x", weight="bold", fontsize=10, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.set(xlabel="Original" + str(df_lesion_IA.columns[y]), ylabel="Post-processing" + str(df_lesion_IA.columns[y]))
plt.xlim(1, 16)
plt.ylim(1, 16)
plt.plot(df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], "k", label="y=x" )
plt.show()
figure = ax.get_figure()
figure.savefig(Path, dpi=400)


y=10 #COV
resultsIA = pearsonr(df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], df_liver_IA[df_liver_IA.columns[y]][int(len(df_liver_IA)/2):])
resultsEARL = pearsonr(df_liver_EARL[df_liver_EARL.columns[y]][:int(len(df_liver_EARL)/2)], df_liver_EARL[df_liver_EARL.columns[y]][int(len(df_liver_EARL)/2):])

ax = sns.regplot(x=df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], y=df_liver_IA[df_liver_IA.columns[y]][int(len(df_liver_IA)/2):], data=df_liver_IA)
ax = sns.regplot(x=df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], y=df_liver_EARL[df_liver_EARL.columns[y]][int(len(df_liver_EARL)/2):], data=df_liver_EARL)
ax.text(0.25, 0.55, "IA  (R²=" + str(round(resultsIA[0]*resultsIA[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["windows blue"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(0.25, 0.48, "EARL  (R²=" + str(round(resultsEARL[0]*resultsEARL[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["pale red"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.set(xlabel="Original" + str(df_liver_IA.columns[y]), ylabel="Post-processing" + str(df_liver_IA.columns[y]))
plt.xlim(0.02, 0.18)
plt.ylim(0.02, 0.18)
plt.show()

ax = sns.regplot(x=df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], y=df_liver_IA[df_liver_IA.columns[y]][int(len(df_liver_IA)/2):], data=df_liver_IA)
ax = sns.regplot(x=df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], y=df_liver_EARL[df_liver_EARL.columns[y]][int(len(df_liver_EARL)/2):], data=df_liver_EARL)
ax.text(0.5, 0.30, "IA  (R²=" + str(round(resultsIA[0]*resultsIA[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["windows blue"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(0.5, 0.23, "EARL  (R²=" + str(round(resultsEARL[0]*resultsEARL[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["pale red"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(0.7, 0.85, "y=x", weight="bold", fontsize=10, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.set(xlabel="Original" + str(df_liver_IA.columns[y]), ylabel="Post-processing" + str(df_liver_IA.columns[y]))
plt.plot(df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], df_liver_IA[df_liver_IA.columns[y]][:int(len(df_liver_IA)/2)], "k", label="y=x" )
plt.xlim(0.08, 0.18)
plt.ylim(0.02, 0.10)
plt.show()


resultsIA = pearsonr(df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], df_lesion_IA[df_lesion_IA.columns[y]][int(len(df_lesion_IA)/2):])
resultsEARL = pearsonr(df_lesion_EARL[df_lesion_EARL.columns[y]][:int(len(df_lesion_EARL)/2)], df_lesion_EARL[df_lesion_EARL.columns[y]][int(len(df_lesion_EARL)/2):])

ax = sns.regplot(x=df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], y=df_lesion_IA[df_lesion_IA.columns[y]][int(len(df_lesion_IA)/2):], data=df_lesion_IA)
ax = sns.regplot(x=df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], y=df_lesion_EARL[df_lesion_EARL.columns[y]][int(len(df_lesion_EARL)/2):], data=df_lesion_EARL)
ax.text(0.1, 0.85, "IA  (R²=" + str(round(resultsIA[0]*resultsIA[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["windows blue"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(0.1, 0.78, "EARL  (R²=" + str(round(resultsEARL[0]*resultsEARL[0],3)) +")", weight="bold", fontsize=12, color=sns.xkcd_rgb["pale red"], horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(0.7, 0.85, "y=x", weight="bold", fontsize=10, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.set(xlabel="Original" + str(df_lesion_IA.columns[y]), ylabel="Post-processing" + str(df_lesion_IA.columns[y]))
plt.xlim(0.1, 0.8)
plt.ylim(0.1, 0.8)
plt.plot(df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], df_lesion_IA[df_lesion_IA.columns[y]][:int(len(df_lesion_IA)/2)], "k", label="y=x" )
plt.show()
figure = ax.get_figure()
figure.savefig(Path, dpi=400)


### End of Pearson correlation analysis ###




###########################################################################################
### 		 														   		   	        ###
### 		 	 Analysis of the Concordance Correlation Coefficient (CCC) 			    ###
### 																   		   	        ###
###########################################################################################

### CCC formula, cf https://nirpyresearch.com/concordance-correlation-coefficient/
def ccc(x,y):
	sxy = np.sum((x - x.mean())*(y - y.mean()))/x.shape[0]      
	rhoc = 2*sxy / (np.var(x) + np.var(y) + (x.mean() - y.mean())**2)      
	return rhoc  


### Storage of the CCC values in list, for the different locations ###

### All loc ###
FeatursAllLocList = []
CCCAllLocValue = []

y=3 
for i in range(104): 
	result = ccc(df[df.columns[y]][:int(len(df)/2)].to_numpy(), df[df.columns[y]][int(len(df)/2):].to_numpy())
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
### All loc ###


### Lesion ###
FeatursLesion = []
CCCLesion = []

y=3 
for i in range(104): 
	result = ccc(df_Lesion[df_Lesion.columns[y]][:int(len(df_Lesion)/2)].to_numpy(), df_Lesion[df_Lesion.columns[y]][int(len(df_Lesion)/2):].to_numpy())
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
for i in range(104): 
	result = ccc(df_liver[df_liver.columns[y]][:int(len(df_liver)/2)].to_numpy(), df_liver[df_liver.columns[y]][int(len(df_liver)/2):].to_numpy())
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
for i in range(104): 
	result = ccc(df_Lung[df_Lung.columns[y]][:int(len(df_Lung)/2)].to_numpy(), df_Lung[df_Lung.columns[y]][int(len(df_Lung)/2):].to_numpy())
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
for i in range(104): 
	result = ccc(df_Muscle[df_Muscle.columns[y]][:int(len(df_Muscle)/2)].to_numpy(), df_Muscle[df_Muscle.columns[y]][int(len(df_Muscle)/2):].to_numpy())  ### pour test non app = stats.ttest_ind
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
for i in range(104): 
	result = ccc(df_BloodPool[df_BloodPool.columns[y]][:int(len(df_BloodPool)/2)].to_numpy(), df_BloodPool[df_BloodPool.columns[y]][int(len(df_BloodPool)/2):].to_numpy())  ### pour test non app = stats.ttest_ind
	FeatursBloodPool.append(str(df_BloodPool.columns[y]))
	CCCBloodPool.append(round(result,3))
	y+=1

plt.figure(figsize=(5,25))
plt.axvline(0.85, 0,1, linewidth=3, color='b')
ax = sns.barplot(x=CCCBloodPool, y=FeatursBloodPool)
ax.set(xlabel="CCC ")
plt.show()
### BloodPool ###



### Save CCC ###
def saveCCCRadiomic(CCC, Featurs):
	df_CustomPalette = pad.DataFrame(CCC, columns = ['CCC'])
	df_CustomPalette['Featurs'] = Featurs
	custom_palette = {}
	for q in set(df_CustomPalette.Featurs):
		val = df_CustomPalette[df_CustomPalette.Featurs == q].CCC
		if val.values < 0.85:
			custom_palette[q] = sns.xkcd_rgb["pale red"]
		else:
			custom_palette[q] = sns.xkcd_rgb["windows blue"]
	plt.figure(figsize=(5,42))
	plt.axvline(0.85, 0,1, linewidth=3, color='b')
	ax = sns.barplot(x=CCC, y=Featurs, palette=custom_palette)
	ax.set(xlabel="CCC ")
	plt.show()

saveCCCRadiomic(CCCLesion, FeatursLesion)
saveCCCRadiomic(CCCLung, FeatursLung)
saveCCCRadiomic(CCCLiver, FeatursLiver)




### To get the CCC with separate class of radiomics ###

### Lesion ###
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
### Lesion ###

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

## To save the radiomics per classes and location
def saveCCCRadiomicClass(CCC, Featurs, Class, Loc):
	df_CustomPalette = pad.DataFrame(CCC, columns = ['CCC'])
	df_CustomPalette['Featurs'] = Featurs
	custom_palette = {}
	for q in set(df_CustomPalette.Featurs):
		val = df_CustomPalette[df_CustomPalette.Featurs == q].CCC
		if val.values < 0.85:
			custom_palette[q] = sns.xkcd_rgb["pale red"]
		else:
			custom_palette[q] = sns.xkcd_rgb["windows blue"]
	plt.figure(figsize=(5,0.36*len(Featurs)))
	plt.axvline(0.85, 0,1, linewidth=3, color='b')
	ax = sns.barplot(x=CCC, y=Featurs, palette=custom_palette)
	ax.set(xlabel="CCC ")
	ax.set_title(Class)
	figure = ax.get_figure()
	figure.savefig(Path + "CCC/"+ str(Loc) + "/CCC_" + str(Loc) + str(Class), dpi=400, bbox_inches='tight')
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

### end of the correlation analyses




### To answer to a reviewer's comment :

### Correlation SUV(Max) to Lesion Volume

results = pearsonr(dfSUVMaxVolume["Max(SUV)"], dfSUVMaxVolume["original_shape_MeshVolume"])
print(results[1]) ## = 2.3033648726612872e-05

ax = sns.regplot(x=dfSUVMaxVolume["volume(ml)"], y=dfSUVMaxVolume["Max(SUV)"], data=dfSUVMaxVolume)
ax.set(xlabel="Tumour volume (ml)", ylabel="SUV Max")
ax.text(0.05, 0.9, "R²=" + str(round(results[0]*results[0],3)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(0.05, 0.8, "p=" + str(round(results[1],5)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
figure = ax.get_figure()
figure.savefig(Path + "Corrélation SUVMax_to_Volume", dpi=400)

## there is a significant correlation (with a very low R²=0.18) but only 3 lesions where smaller than 100 pixels.