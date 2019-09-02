import numpy as np
import pandas as pd
import networkx as nx
import rpy2.robjects as robjects
import os, pickle, itertools
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import sparse

from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from networkx.algorithms import bipartite as bi
import numba

# Suppress R script warnings
import warnings
from rpy2.rinterface import RRuntimeWarning

warnings.filterwarnings("ignore", category=RRuntimeWarning)
pandas2ri.activate()

## @package Bipartite Co-expression Network
# Basic class
class BiCoxNet(object):

    ## Constructor
    # @param mrna: mRNA expression profile: pandas DataFrame (genes x samples)
    # @param serna: seRNA expression profile: pandas DataFrame (se x samples)
    # @param celltype: string of celltype, e.g. 'HES3_GFP_ESC'
    # @param stage_id: int of stage_id
    def __init__(self, mrna, serna, celltype, stage_id):
        try:
            self.mrna_locus = mrna['locus']
            mrna = mrna.drop(['locus'], axis=1)
        except:
            self.mrna_locus = pd.Series(mrna.index, name='locus')
        self.serna_locus = serna['locus']
        serna = serna.drop(['locus'], axis=1)
        self.celltype = celltype
        self.stage_id = stage_id

        self.mrna, self.serna = mrna, serna

    ## Load pre-computed correlation matrix
    # @param filepath: csv file of correlation matrix
    def load_corr(self, filepath):
        self.corr = pd.read_csv(filepath, index_col=0)
        return self

    ## Calculate mRNA-to-seRNA correlation matrix
    #  Filter correlation which p_value is larger than specified threshold
    # @param method: 'bicor', 'pearson(default)', 'spearman'
    # @param p_value: threshold to filter correlation, if 0, then no test
    # Return: mRNA-to-seRNA absolute correlation matrix
    def calculate_corr(self, method='pearson', p_value=0.05, save=False, abs_flag=True):
        self.method = method
        r = robjects.r
        r.source("./network/binet_utils.R", chdir=True)
        
        with localconverter(robjects.default_converter + pandas2ri.converter):
            data1 = robjects.conversion.py2rpy(self.mrna)
            data2 = robjects.conversion.py2rpy(self.serna)
        #data1 = pandas2ri.py2ri(self.mrna)
        #data2 = pandas2ri.py2ri(self.serna)
        with localconverter(robjects.default_converter + pandas2ri.converter):
            if p_value == 0:
                self.corr = robjects.conversion.rpy2py(r.calculate_cor(data1, data2, method=method, abs_flag=abs_flag))
            else:
                self.corr = robjects.conversion.rpy2py(r.calculate_cor_test(data1, data2, p_value, method=method))
        #if p_value == 0:
        #    self.corr = pandas2ri.ri2py(r.calculate_cor(data1, data2, method=method, abs_flag=abs_flag))
        #else:
        #    self.corr = pandas2ri.ri2py(r.calculate_cor_test(data1, data2, p_value, method=method))
        self.corr.index = self.mrna_locus
        self.corr.columns = self.serna_locus

        if save:
            self.corr.to_csv('net_corr_mat{}.csv'.format(self.stage_id), index=True)
        return self

    ## Half-cumulative distribution method applied to correlation matrix(mRNA x seRNA)
    def half_cumulative(self):
        
        return self 

    ## Implement weighted bipartite network projection, weight is corr
    # @param onto: compute which projection onto mrna or serna 
    def projection(self, onto='mrna'):
        
        @numba.jit
        def _calculate(data):
            su0, su1 = data.sum(axis=0), data.sum(axis=1)
            su0[su0 == 0] = 1
            su1[su1 == 0] = 1
            norm = data / su0
            return np.dot(norm, data.T) / su1
        
        if onto == 'mrna':
            s = _calculate(self.corr.values)
        elif onto == 'serna':
            s = _calculate(self.corr.values.T)
        else:
            print("Usage: [object].projection(onto='mrna / serna')")
        return s 
    #==================================================================
    # WGCNA functions
    #==================================================================    
    ## Choose best beta for adjacency function: |cor|^beta
    # @param data: correlation matrix of mRNA to seRNA
    #
    # Criterion: 1. scale-free R^2 close to 1
    #            2. Largest mean connectivity
    #
    def choose_power(self, method='pearson'):
        r = robjects.r
        r.source("./network/binet_utils.R", chdir=True)
        
        with localconverter(robjects.default_converter + pandas2ri.converter):
            data1 = robjects.conversion.py2rpy(self.mrna)
            data2 = robjects.conversion.py2rpy(self.serna)
        #data1 = pandas2ri.py2ri(self.mrna)
        #data2 = pandas2ri.py2ri(self.serna)
        try:
            power_table = r.choose_power(data1, data2, self.celltype, self.stage_id, method=self.method)
        except AttributeError:
            power_table = r.choose_power(data1, data2, self.celltype, self.stage_id, method=method)

        return power_table  
        
    ## Check scale-free property
    # @param beta: The beta parameter you want to check if scale-free property holds
    #
    def check_scale_free(self, beta):
        r = robjects.r
        r.source("./network/binet_utils.R", chdir=True)
        
        with localconverter(robjects.default_converter + pandas2ri.converter):
            data = robjects.conversion.py2rpy(self.corr)
            
        #data = pandas2ri.py2ri(self.corr)
        r.check_scalefree(data, beta, self.celltype, self.stage_id)

    #=======================================================================================================
    # Plotting functions
    #======================================================================================================
    
    ## Plot correlation heatmap
    # @param beta: paramether of adjacency function |cor|^beta
    # @param num_mrna: tuple specify the index of mrna loci to plot(default:50), 'all' indicates plot all loci.
    # @param num_serna: tuple specify the index of serna loci to plot(default:50)
    def plot_heatmap(self, beta, num_mrna=(0,50), num_serna=(0,50), filepath=None):
        f, ax = plt.subplots(1, 1, figsize=(12, 8))
        if num_mrna != 'all' and num_serna != 'all':
            data = (self.corr**beta).iloc[num_mrna[0]:num_mrna[1], num_serna[0]:num_serna[1]]
        elif num_mrna == 'all' and num_serna != 'all':
            data = (self.corr**beta).iloc[:,num_serna[0]:num_serna[1]]
        elif num_serna == 'all' and num_mrna != 'all':
            data = (self.corr**beta).iloc[num_mrna[0]:num_mrna[1],:]
        else:
            data = self.corr**beta

        sns.heatmap(data, vmin=0.0, vmax=1.0) 
        ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)
        ax.set_yticklabels(ax.yaxis.get_majorticklabels(), rotation=0)
        plt.title('{} corr^{} state{} genes{}-{} to seRNA {}'.format(self.method, beta, self.stage_id, num_mrna[0], num_mrna[1], num_serna))
        plt.ylabel('Genes of {} to {}'.format(num_mrna[0], num_mrna[1]))
        plt.xlabel('Super Enhancers')
        f.tight_layout() 
        
        if filepath == None:
            if not os.path.exists('../graph/'+self.celltype):
                os.makedirs('../graph/'+self.celltype)

            plt.savefig(os.path.join('../graph', self.celltype, '{} corr^{} state{} genes{}-{} to seRNA {}'.format(self.method, beta, self.stage_id, num_mrna[0], num_mrna[1], num_serna)))
            plt.close() 
        else:
            plt.savefig(filepath)
            plt.close()

    ## Plot mRNA to seRNA correlation count distribution
    def plot_top_mRNA_hist(self):
        top = []
        label = []
        for i in range(0, 600, 100):
            top.append(self.corr.iloc[i, :]) 
            label.append('Top{} mRNA'.format(i))    
        top = tuple(top)
        n, bins, patches = plt.hist(top, label=label)
        plt.title('state{} top mRNA: count of correlation to seRNA'.format(self.stage_id))
        plt.xlabel('{} correlation'.format(self.method))
        plt.ylabel('Counts')
        plt.legend(loc='best')
        plt.show()
        plt.close() 
