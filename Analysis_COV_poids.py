#coding: utf-8


from scipy import stats
import pandas as pad
import matplotlib.pyplot as plt
import seaborn as sns
import codecs
import statsmodels.formula.api as smf
import statsmodels.api as sm

from scipy.stats import pearsonr

## POUR CHANGER COULEUR sns.set_palette("bright")
#custom_palette = [sns.xkcd_rgb["medium green"], sns.xkcd_rgb["pale red"], "blue", "orange", "blue","yellow", "purple"]
custom_palette = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["pale red"], sns.xkcd_rgb["medium green"], "orange", "blue","yellow", "purple"]
sns.set_palette(custom_palette)

##########################
### Analyse effet Algo ###
##########################


##################################################
### Représentation des différences deux à deux ###
##################################################

df = pad.read_excel('//s-grp/grp/RADIOPHY/Personnel/Aurélien Corroyer-Dulmont/Recherche/AI TEP SubtleMedical/article front oncol noise IA/AnalyseImcIA.xlsx', sheet_name='ACD', usecols="A:D", nrows=226, header=0)

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


########################################
### Etude corrélation entre facteurs ###
########################################

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
