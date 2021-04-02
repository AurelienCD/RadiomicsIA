#coding: utf-8


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

#custom_palette = [sns.xkcd_rgb["medium green"], sns.xkcd_rgb["pale red"], "blue", "orange", "blue","yellow", "purple"]
custom_paletteGeneral = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["pale red"], sns.xkcd_rgb["medium green"], "orange", "blue","yellow", "purple"]
sns.set_palette(custom_paletteGeneral)



########################################
###        Path and dataframe        ###
########################################

Path = '//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/test/'
df = pad.read_excel("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/AnalyseRadiomic-IA.xlsx", header=0)
df = pad.read_excel("//s-grp/grp/RADIOPHY/NOUVELLE ARBORESENCE/Imagerie/Projets de RECHERCHE/Médecine Nucléaire/AI TEP SubtleMedical_Article_Front_Oncol/AnalyseRadiomic-EARL.xlsx", header=0)



### To get the values for each location ###
df_liver = df[df[" labels "] == "Liver"]
df_Lung = df[df[" labels "] == "Lung"]
df_BloodPool = df[df[" labels "] == "BloodPool"]
df_Muscle = df[df[" labels "] == "Muscle"]
df_Lesion = df[df[" labels "] == "Lesion"]


																			#####################################
																			### Analysis of Algorithm effect  ###
																			#####################################


def saveTtestAnalysis(df, Loc):
	
	### significant
	y=3 
	for i in range(118): 
		results = stats.ttest_rel(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])  ### pour test non app = stats.ttest_ind
		ax = sns.boxplot(x='Groupe', y=df[df.columns[y]], data=df, showfliers = False) 
		if results[1] < 0.05:
			figure = ax.get_figure()
			figure.savefig(Path + "ttest/" + str(Loc) +"/t-test__" + str(df.columns[y]), dpi=400)
			print("Pour la variable " + str(df.columns[y] + " la valeur de p concernant la différence entre les deux images est de " + str(results[1]))) 
		figure.clear()
		y+=1

	### non significant
	y=3 
	for i in range(118): 
		results = stats.ttest_rel(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])  ### pour test non app = stats.ttest_ind
		ax = sns.boxplot(x='Groupe', y=df[df.columns[y]], data=df, showfliers = False) 
		if results[1] > 0.05:
			figure = ax.get_figure()
			figure.savefig(Path + "ttest/" + str(Loc) +"/PAS SIGNIFICATIF/t-test__" + str(df.columns[y]), dpi=400)
			print("Pour la variable " + str(df.columns[y] + " la valeur de p est supérieur à 0.05, et est égale à : " + str(results[1])))
		figure.clear()
		y+=1

saveTtestAnalysis(df, "toutes_loc")
saveTtestAnalysis(df_liver, "liver")
saveTtestAnalysis(df_Lung, "lung")
saveTtestAnalysis(df_Muscle, "muscle")
saveTtestAnalysis(df_Lesion, "lesion")
saveTtestAnalysis(df_BloodPool, "bloodpool")




									#############################################################################################
									### 	Correlation study between features values before and after algorithm application  ###
									### 	which radiomics features are stable after algorithm application        		   	  ###
									#############################################################################################

									#####################################################################################
									### Etude corrélation entre facteurs avant et après application de l'aglo d'IA    ###
									### 	est-ce que le paramètre de radiomic et stable après IA        		   	  ###
									#####################################################################################


def saveCorrelationFeatures(df, Loc):
	y=3 
	for i in range(118): 
		results = pearsonr(df[df.columns[y]][:int(len(df)/2)], df[df.columns[y]][int(len(df)/2):])
		if results[1] < 0.05:
			ax = sns.regplot(x=df[df.columns[y]][:int(len(df)/2)], y=df[df.columns[y]][int(len(df)/2):], data=df)
			ax.set(xlabel="Origin " + str(df.columns[y]), ylabel="IA " + str(df.columns[y]))
			ax.text(0.05, 0.9, "R²=" + str(round(results[0]*results[0],3)), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
			figure = ax.get_figure()
			figure.savefig(Path + "Corrélation features/Figures_stabilité_radiomics_" + str(Loc) +"/Effet_algo__" + str(df.columns[y]), dpi=400)
			print("Pour la variable " + str(df.columns[y] + " la valeur de p est inférieur à 0.05, et R² est égal à : " + str(round(results[0]*results[0],3))))
		figure.clear()
		y+=1

saveCorrelationFeatures(df, "toutes_loc")
saveCorrelationFeatures(df_liver, "liver")
saveCorrelationFeatures(df_Lung, "lung")
saveCorrelationFeatures(df_Muscle, "muscle")
saveCorrelationFeatures(df_Lesion, "lesion")
saveCorrelationFeatures(df_BloodPool, "bloodpool")




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
	result = ccc(df_Lesion[df_Lesion.columns[y]][:int(len(df_Lesion)/2)].to_numpy(), df_Lesion[df_Lesion.columns[y]][int(len(df_Lesion)/2):].to_numpy())
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
	result = ccc(df_liver[df_liver.columns[y]][:int(len(df_liver)/2)].to_numpy(), df_liver[df_liver.columns[y]][int(len(df_liver)/2):].to_numpy())
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
	result = ccc(df_Lung[df_Lung.columns[y]][:int(len(df_Lung)/2)].to_numpy(), df_Lung[df_Lung.columns[y]][int(len(df_Lung)/2):].to_numpy())
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







def normalityDistributionTestAgostino(df):

	### avec test D’Agostino
	y=3 
	for i in range(118): 
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