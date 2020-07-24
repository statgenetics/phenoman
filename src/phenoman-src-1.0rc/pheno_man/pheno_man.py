#!/usr/bin/env python3
#
# Copyright (C) 2012-2013 Biao Li (biaol@bcm.edu) and Gao Wang (gaow@bcm.edu)
# License: GNU General Public License (http://www.gnu.org/licenses/)

import os, sys, random, copy, math
from str4r import runR, Str4R
from collections import OrderedDict, Counter
from itertools import izip, chain, product
import re
## global variables
#DATA_DIR = os.path.dirname(sys.modules['pheno_man'].__file__)
#DUP_AA = [x.rstrip().split('\t') for x in open(os.path.join(DATA_DIR, 'AAdups.txt')).readlines()[1:]]
#DUP_EA = [x.rstrip().split('\t') for x in open(os.path.join(DATA_DIR, 'EAdups.txt')).readlines()[1:]]
#DUP_LIST = DUP_AA + DUP_EA

HIST_FOO = '''
gghues <- function(n) {
        hues = seq(15, 375, length=n+1)
        hcl(h=hues, l=65, c=100)[1:n]
}
gghue <- function(n) { gghues(n)[n] }
gghist <- function(dat, stat_foo, vlines=NULL, xname='x') {
        average <- round(mean(dat$x, na.rm=T),4)
        stdev <- round(sd(dat$x, na.rm=T),4)
        #med <- round(median(dat$x, na.rm=T), 4)
        r1 <- round(min(dat$x, na.rm=T),4)
        r2 <- round(max(dat$x, na.rm=T),4)
        suppressMessages(library(ggplot2))
        # convert dat obj from numeric to data.frame
        myplot <- ggplot(dat) +
                        aes_string(x = 'x') +
                        geom_histogram(aes_string(y = '..density..', fill = '..density..'), color = 'white') +
                        scale_fill_gradient('bin mass', low = 'darkolivegreen3', high = colors()[552]) +
                        geom_line(stat = 'density', size = 0.5, linetype = 2, color = 'grey50') +
                        stat_function(fun = stat_foo, color = 'red') +
                        geom_rug(color = 'grey80') +
                        scale_x_continuous(name = paste('\\n', xname)) +
                        scale_y_continuous(name = 'Density\\n') +
                        labs(title = paste('Histogram & fitted vs. theoretical density curves for', xname, '\\n', 'mean = ', toString(average), '; ', 'stdev = ', toString(stdev), '; ', 'range = ', '[', toString(r1), ',', toString(r2), ']', '\\n')) +
                        theme_bw()
        if (!is.null(vlines)) {
        myplot <- myplot + geom_vline(xintercept = vlines, color = '#9ACD32', linetype = 2)
        }
        return(myplot)
}
'''

# R func to perform model selection
SUMMARY_FOO = 'print(summary(dat))'
LM_FOO_both = '''
#options(width=10000)
dat <- na.omit(dat)
delim <- ','
names(dat)[1] <- "y"
if (length(unique(dat[,1])) == 2) {
        md <- glm(y ~ ., data=dat, family='binomial')
} else {
        md <- lm(y ~ ., data=dat)
}
mselect <- step(md, trace=1, direction = 'both')
mresult <- summary(md)
write('INPUT MODEL:', stdout())
write.table(mresult$coefficients, stdout(), sep=delim)
f <- as.numeric(mresult$fstatistic)
if (is.na(f[1])) {
        pval <- 1 - pchisq(mresult$null.deviance-mresult$deviance, mresult$df.null-mresult$df.residual)
} else {
        pval <- 1 - pf(f[1], f[2], f[3])
}
write(paste('p-value:', delim, pval, '\n', sep=''), stdout())

mresult <- summary(mselect)
write('OPTIMIZED MODEL:', stdout())
write.table(mresult$coefficients, stdout(), sep=delim)
f <- as.numeric(mresult$fstatistic)
if (is.na(f[1])) {
        pval <- 1 - pchisq(mresult$null.deviance-mresult$deviance, mresult$df.null-mresult$df.residual)
} else {
        pval <- 1 - pf(f[1], f[2], f[3])
}
write(paste('p-value:', delim, pval, '\n', sep=''), stdout()
'''

LM_FOO_forward = '''
#options(width=10000)
dat <- na.omit(dat)
delim <- ','
names(dat)[1] <- "y"
if (length(unique(dat[,1])) == 2) {
        md <- glm(y ~ ., data=dat, family='binomial')
} else {
        md <- lm(y ~ ., data=dat)
}
mselect <- step(md, trace=1, direction = 'forward')
mresult <- summary(md)
write('INPUT MODEL:', stdout())
write.table(mresult$coefficients, stdout(), sep=delim)
f <- as.numeric(mresult$fstatistic)
if (is.na(f[1])) {
        pval <- 1 - pchisq(mresult$null.deviance-mresult$deviance, mresult$df.null-mresult$df.residual)
} else {
        pval <- 1 - pf(f[1], f[2], f[3])
}
write(paste('p-value:', delim, pval, '\n', sep=''), stdout())

mresult <- summary(mselect)
write('OPTIMIZED MODEL:', stdout())
write.table(mresult$coefficients, stdout(), sep=delim)
f <- as.numeric(mresult$fstatistic)
if (is.na(f[1])) {
        pval <- 1 - pchisq(mresult$null.deviance-mresult$deviance, mresult$df.null-mresult$df.residual)
} else {
        pval <- 1 - pf(f[1], f[2], f[3])
}
write(paste('p-value:', delim, pval, '\n', sep=''), stdout()
'''

LM_FOO_backward = '''
#options(width=10000)
dat <- na.omit(dat)
delim <- ','
names(dat)[1] <- "y"
if (length(unique(dat[,1])) == 2) {
        md <- glm(y ~ ., data=dat, family='binomial')
} else {
        md <- lm(y ~ ., data=dat)
}
mselect <- step(md, trace=1, direction = 'backward')
mresult <- summary(md)
write('INPUT MODEL:', stdout())
write.table(mresult$coefficients, stdout(), sep=delim)
f <- as.numeric(mresult$fstatistic)
if (is.na(f[1])) {
        pval <- 1 - pchisq(mresult$null.deviance-mresult$deviance, mresult$df.null-mresult$df.residual)
} else {
        pval <- 1 - pf(f[1], f[2], f[3])
}
write(paste('p-value:', delim, pval, '\n', sep=''), stdout())

mresult <- summary(mselect)
write('OPTIMIZED MODEL:', stdout())
write.table(mresult$coefficients, stdout(), sep=delim)
f <- as.numeric(mresult$fstatistic)
if (is.na(f[1])) {
        pval <- 1 - pchisq(mresult$null.deviance-mresult$deviance, mresult$df.null-mresult$df.residual)
} else {
        pval <- 1 - pf(f[1], f[2], f[3])
}
write(paste('p-value:', delim, pval, '\n', sep=''), stdout()
'''


# R func to calculate residuals
RESID_FOO = '''
#options(width=10000)
dat <- na.omit(dat)
delim <- ','
names(dat)[1] <- "y"
# residuals only on quantitative traits (lm) not (glm)
md <- lm(y ~ ., data=dat)
residuals = md$residuals
#fitted.values = md$fitted.values
write(paste('residuals:', '\n', sep=''), stdout())
write.table(residuals, stdout(), sep=delim)
#write(paste('fitted.values:', '\n', sep=''), stdout())
#write.table(fitted.values, stdout(), sep=delim)
'''


# utility functions
def typeOfValues(vals):
    '''Figure out type of values and return the type'''
    vals = ['nan' if x.lower() in ('na', 'null', 'none', '') else x for x in vals]
    if len(vals) == 0:
        return 'str'
    try:
        map(int, vals)
        return 'int'
    except:
        try:
            map(float, vals)
            return 'float'
        except:
            return 'str'


def typeOfValues_selectModule(vals):
    types = []
    for val in vals:
        try:
            tmp = int(val)
            types.append(int)
        except:
            try:
                tmp = float(val)
                types.append(float)
            except:
                types.append(str)
    return types
            

def tconvert(x):
    try:
        return int(x)
    except:
        try:
            return float(x)
        except:
            try:
                return eval(x)
            except:
                return x

# missingness filter
def mfilter(frequency):
   try:
       value = float(frequency)
       if not ((value >= 0 and value <= 1) or value == -1):
           msg = "{0} is not valid input. Valid input should fall in range [0, 1] or equal to -1".format(frequency)
           raise ValueError(msg)
   except Exception as e:
       raise argparse.ArgumentTypeError(e)
   return value


def freq(frequency):
    try:
        value = float(frequency)
        if not (value >= 0 and value <= 1):
            msg = "{0} is not valid input. Valid input should fall in range [0, 1]".format(frequency)
            raise ValueError(msg)
    except Exception as e:
        raise argparse.ArgumentTypeError(e)
    return value

def input2dict(inlines, recodeNulls = True):
    '''
    convert input phenotypes to a dictionary
    recodeNulls -- if True, recode all 'NA', 'Null', 'None' (.lower()) into 'nan', otherwise keep original coding for nulls. 
    '''
    # determine header
    if not inlines:
        raise ValueError('No input data')
    headers = inlines[0].rstrip().split('\t')
    if len(headers) == 0:
        raise ValueError('Empty header line. Nothing will be imported.')
    if headers[0] == 'sample_name':
        if len(headers) == 1:
            raise ValueError('No phenotype to be imported')
    else:
        raise ValueError('The input phenotype data must be tab delimited, and starts with a header line '
                                        'with the first column sample_name.')
    pdict = OrderedDict()
    body = [x.rstrip().split('\t') for x in inlines[1:]]
    for k, item in zip(headers, zip(*body)):
        if recodeNulls:
            item = ['nan' if x.lower() in ('na', 'null', 'none', '') else x for x in item]
            item_type = typeOfValues(item)
            pdict[k] = list(map(eval(item_type), item))
        else: # special case for using select module
            types = typeOfValues_selectModule(item)
            pdict[k] = [i(j) for i,j in zip(types, item)]
    return pdict

def pairwise(iterable):
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)

def rsummary(pdict, main = None, covariates = None, option='summary', addResiduals=False, direction='both'):
    '''summary statistics of phenotype dictionary, or model selection or adding residuals, via R'''
    if 'sample_name' in pdict.keys():
        del(pdict['sample_name'])
    # move main phenotype to front
    if main:
        pdict.move_to_end(main, last = False)
    if covariates:
        for idx, k in enumerate(pdict.keys()):
            if k not in covariates and idx > 0:
                del(pdict[k])
    # recode missing values
    for k in pdict.keys():
        pdict[k] = [None if ismissing(x) else x for x in pdict[k]]

    # Run R summary
    rstr = '''
    dat <- as.data.frame({});
    '''.format(Str4R(pdict))
    if option == 'summary':
        rstr += SUMMARY_FOO
    else: # model selection
        if not addResiduals:
            if direction == 'both':
                rstr += LM_FOO_both
            elif direction == 'forward':
                rstr += LM_FOO_forward
            else:
                rstr += LM_FOO_backward
        else: # add residuals
            # return error if 'main' is a case-control trait
            if len(set(pdict[main])) <= 2:
                sys.stderr.write("\nERROR: Cannot calculate residuals for a case-control trait '%s'.\n\n" % main)
                sys.exit(0)
            rstr += RESID_FOO
    out, err = runR(rstr)
    if err:
        sys.stderr.write(err)
        raise ValueError("See **** R ERROR **** above!")
    return out

def rhist(data, phenotype, output = 'histogram.pdf', vlines = None):
    '''draw histogram using ggplot2'''
    if vlines:
        vlines = Str4R(vlines)
    else:
        vlines = 'NULL'
    rdat = Str4R([None if ismissing(x) else x for x in data[phenotype]])
    if output.endswith('.pdf'):
        rstr = '''
                    dat <- as.data.frame({})
                    colnames(dat) <- c('x')
                    stat_foo <- function(x) dnorm(x, mean(dat$x, na.rm=T), sd(dat$x, na.rm=T))
                    vlines <- {}
                    p <- gghist(dat, stat_foo, vlines=vlines, xname={})
                    pdf({},12,9)
                    print(p)
                    dev.off()
                    '''.format(rdat, vlines, repr(phenotype), repr(output))
    elif output.endswith('.png'):
        rstr = '''
                    dat <- as.data.frame({})
                    colnames(dat) <- c('x')
                    stat_foo <- function(x) dnorm(x, mean(dat$x, na.rm=T), sd(dat$x, na.rm=T))
                    vlines <- {}
                    p <- gghist(dat, stat_foo, vlines=vlines, xname={})
                    png({},960,640)
                    print(p)
                    dev.off()
                    '''.format(rdat, vlines, repr(phenotype), repr(output))
    else:
        raise ValueError("rhist can only output figure file in pdf or png format for now, check output file name")
    #
    sys.stdout.write('Generating histogram ...\n')
    out, err = runR(HIST_FOO + rstr)
    sys.stdout.write('Done!\n')
    if err:
        sys.stderr.write(err)
    return

def ismissing(x):
    return (x != x or x == 'nan')

def isOutlier(x, criteria):
    return (x < criteria[0] or x > criteria[1])

def recode_categorical(d, k, refgrp=None):
    if k not in d.keys():
        raise ValueError("Field name '{}' not found".format(k))
    if refgrp and refgrp not in d[k]:
        raise ValueError("Reference group name '{}' not found in field {}".format(refgrp, k))
    mapping = OrderedDict()
    i = 65
    if refgrp:
        mapping[refgrp] = 'A'
        i = 66
    for item in d[k]:
        if item not in mapping.keys():
            mapping[item] = str(chr(i))
            i += 1
    for idx, item in enumerate(d[k]):
        d[k][idx] = mapping[item]
    return d, mapping

def qnorm( p, mean = 0.0, sd = 1.0):
    """
    Modified from the author's original perl code (original comments follow below)
    by dfield@yahoo-inc.com.  May 3, 2004.
    Lower tail quantile for standard normal distribution function.
    This function returns an approximation of the inverse cumulative
    standard normal distribution function.  I.e., given P, it returns
    an approximation to the X satisfying P = Pr{Z <= X} where Z is a
    random variable from the standard normal distribution.
    The algorithm uses a minimax approximation by rational functions
    and the result has a relative error whose absolute value is less
    than 1.15e-9.
    Author:  Peter John Acklam
    Time-stamp:  2000-07-19 18:26:14
    E-mail: pjacklam@online.no
    WWW URL: http://home.online.no/~pjacklam
    """
    if p <= 0 or p >= 1 or p != p:
            # The original perl code exits here, we'll throw an exception instead
            raise ValueError( "Invalid input '%f'. Value must be in open interval (0,1)" % p )
    # Coefficients in rational approximations.
    a = (-3.969683028665376e+01,  2.209460984245205e+02, \
             -2.759285104469687e+02,  1.383577518672690e+02, \
             -3.066479806614716e+01,  2.506628277459239e+00)
    b = (-5.447609879822406e+01,  1.615858368580409e+02, \
             -1.556989798598866e+02,  6.680131188771972e+01, \
             -1.328068155288572e+01 )
    c = (-7.784894002430293e-03, -3.223964580411365e-01, \
             -2.400758277161838e+00, -2.549732539343734e+00, \
              4.374664141464968e+00,  2.938163982698783e+00)
    d = ( 7.784695709041462e-03,  3.224671290700398e-01, \
              2.445134137142996e+00,  3.754408661907416e+00)
    # Define break-points.
    plow  = 0.02425
    phigh = 1 - plow
    # Rational approximation for lower region:
    if p < plow:
            q  = math.sqrt(-2*math.log(p))
            z = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    # Rational approximation for upper region:
    elif phigh < p:
            q  = math.sqrt(-2*math.log(1-p))
            z = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    # Rational approximation for central region:
    else:
            q = p - 0.5
            r = q*q
            z= (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
                            (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    # transform to non-standard:
    return mean + z * sd # !@#$% sorry, just discovered Sep. 9, 2011

def calMean(xx):
    return float(sum(xx))/len(xx)
    
def calStd(xx):
    mean = calMean(xx)
    return (sum([(i-mean)**2 for i in xx])/len(xx))**0.5
    
def addResiduals(dataset, field, residuals, outFileName):
    '''
    add residuals from regression model as a new column 'phenotype_residuals' in 'dataset'
    and dump the dataset to 'outFileName'
    dataset - the dataset which is piped into the regression model
    field - phenotype name
    residuals - list of residuals from regression
    outFileName - output file name of the new dataset containing the column of 'residuals'
    '''
    with open(outFileName, 'w') as outFile:
        dataset[field+'_residuals'] = copy.deepcopy(dataset[field])
        allFields = list(dataset.keys())
        for idx, value in enumerate(residuals):
            dataset[field+'_residuals'][idx] = value
        if 'sample_name' in allFields:
            allFields.remove('sample_name')
        outString = 'sample_name\t'+'\t'.join(allFields)+'\n'
        #
        for idx, name in enumerate(dataset['sample_name']):
            outString += (name+'\t')
            for j, field in enumerate(allFields):
                outString += (str(dataset[field][idx])+'\t') if j < len(allFields)-1 else (str(dataset[field][idx]))
            outString += '\n'
        #
        outFile.write(outString)
    outFile.close()
    sys.stderr.write("\nAdd residuals to '%s' and write new dataset to %s\n\n" % (field, outFileName))
    return

def isEqual(x, y):
    '''
    return True if x is equal to y, particularly allow 'NaN'.lower() to be equal to float('NaN') and 'InF'.lower() to be equal to float('inf')
    '''
    if str(x).lower() in ['nan', 'inf']:
        return str(x).lower() == str(y).lower()
    else:
        return x == y
    
        

class PhenData:
    def __init__(self, fname, recodeNullsInFile=True):
        '''Phenotype data object and member functions to modify the data'''
        with open(fname, 'r') as f:
            inlines = f.readlines()
        self.data = input2dict(inlines, recodeNullsInFile)
    
        #print self.data['adiponectin']
    
    def fsummarize(self):
        '''summarize fields properties'''
        # valid entries count per field
        field_vcounts = []
        # missing entries count per field
        field_mcounts = []
        for field in self.data.keys():
            vcount, mcount = self.countElements(self.data[field])
            field_vcounts.append(vcount)
            field_mcounts.append(mcount)
        return self.data.keys(), field_vcounts, field_mcounts

    def countElements(self, vector):
        '''count valid/missing elements in input vector'''
        mcount = sum([1 for x in vector if ismissing(x)])
        return len(vector) - mcount, mcount

    def findDuplicates(self, fileNameList):
        '''detect duplicates in samples based on a list of files
        with each having two columns/lists,  
        (for each of elements in the list, the first entry is to
        be removed, the 2nd entry is to be remained).
        return a list of sample names to remove'''
        # construct the list of pairs of duplicated individuals
        listDups = [[x.rstrip().split('\t') for x in open(fi).readlines()] for fi in fileNameList]
        dupList = list(chain(*listDups))
        #
        if fileNameList == []:
            dupList = DUP_LIST        
        #
        sampleNames = self.data['sample_name']
        duplicateNames = []
        #
        for dupInfo in dupList:
            remove, keep = dupInfo[0], dupInfo[1]
            # first see if both of them are in the sample, otherwise do nothing
            if remove in sampleNames and keep in sampleNames:
            # if 'remove' is already an element of duplicateNames do nothing
                if remove in duplicateNames:
                    continue
                else:
                    duplicateNames.append(remove)
            else:
                continue
        #
        return duplicateNames

    def findMissing(self, field):
        '''detect missing samples by a field.
        return a list of sample names of having missing values
        on this particular field'''
        nameList = [name for name, value in zip(self.data['sample_name'], self.data[field]) if ismissing(value)]
        return nameList

    def findOutliers(self, field, criteria):
        '''detect outliers by a field based on criteria, and
        return a dict with keys as sample names of outliers (to remove)
        with values as the outlier values of those individuals.
        Note that 'criteria' is a list of 2 values [min, max]
        '''
        if criteria:
            return {name:value for name, value in zip(self.data['sample_name'], self.data[field]) if isOutlier(value, criteria)}
        else:
            return {}

    def removeSamples(self, names):
        '''delete samples based on input sample names
        Will modify self.data directly
        Note that:
        names can be either a list of names (e.g. returned values
        from both func findMissing and func findDuplicates)
        or a dict obj (e.g. returned value from func findOutliers)
        '''
        #
        if type(names) is dict:
            names = names.keys()
        elif type(names) in [list, tuple]:
            pass
        else:
            raise ValueError("arg names need to be a list/tuple of sample names or a dictionary obj with its keys as sample names")
        #
        for name in names:
            try: # if name is not in 'sample_name', do nothing
                idx = self.data['sample_name'].index(name)
                for field in self.data.keys():
                    del self.data[field][idx]
            except:
                print('Warning: sample %s not found in database\n' % name)
                continue
        #
        return

    def log(self, field):
        '''apply log transformation to self.data[field]
        Will add another field named by 'log_field' and modify self.data[log_field]'''
        # Note that: need to remove all missing values before performing log transformation
        self.data["log_"+field] = copy.deepcopy(self.data[field])
        for idx, value in enumerate(self.data["log_"+field]):
            self.data["log_"+field][idx] = math.log(value)
        #
        return
    
    def log10(self, field):
        '''
        apply log10 transformation to self.data[field]. It will add another field named by "log10_field" and modify self.data[log10_field]
        '''
        self.data["log10_"+field] = copy.deepcopy(self.data[field])
        for idx, value in enumerate(self.data["log10_"+field]):
            self.data["log10_"+field][idx] = math.log(value, 10)
        #
        return
    
    def findPercentile(self, field, percentiles):
        '''find values based on input percentiles for self.data[field] and return it
        'percentiles' need to be a list of two values [lower_bound, upper_bound]'''
        # Note that: need to remove all duplicates and missing values before identifying outliers
        #
        fieldInfo = self.data[field]
        minimum, maximum = min(fieldInfo), max(fieldInfo)
        #
        return [(maximum-minimum) * percentiles[0] + minimum, (maximum-minimum) * percentiles[1] + minimum]

    def fillMissing(self, verbose = True, threshold = 1):
        '''for each field in self.data, will fill up missing values
        by average value of the field if len(list(set(self.data[field]))) > 2
        else fill it with min(listset(self.data[field])). Will modify self.data directly.
        Missing value can be checked by condition "ismissing(x)". If verbosity is true,
        write the summary data information in the same manner (but probably different information)
        as "phenoman show fields" via "sys.stdout.write()"
        'threshold' - allow up to 'threshold' (in percentage) missingness
        '''
        allPhenoNames = self.data.keys()[1:]
        propMissingDict = {}
        deletedPhenosDict = {}
        #
        for phenoName in allPhenoNames:
            phenoInfo = self.data[phenoName]
            # if not a numeric trait
            if type(phenoInfo[0]) not in [float, int]:
                self.data[phenoName] = ['NA' if x == 'nan' else x for x in phenoInfo]
                continue
            #
            tmpValues = []
            tmpMissingIdx = []
            for idx, pheno in enumerate(phenoInfo):
                if ismissing(pheno):
                    tmpMissingIdx.append(idx)
                else:
                    tmpValues.append(pheno)
            #
            # check if there are too many missing values
            propMiss = len(tmpMissingIdx)/float(len(phenoInfo))
            if propMiss >= threshold:
                # if proportion of missingness >= threshold
                # save its header info and proportion of missingness into 'deletedPhenosDict' and remove the covariate col with high missingness from self.data
                deletedPhenosDict[phenoName] = propMiss
                del self.data[phenoName]
                #
            else: # fill missing values and save phenotype info and proportion of missingness into 'propMissingDict'
                propMissingDict[phenoName] = [propMiss]
                if len(set(tmpValues)) > 2:
                    avePheno = round(sum(tmpValues)/float(len(tmpValues)),3)
                    #
                elif len(set(tmpValues)) > 0:
                    avePheno = min(tmpValues)
                    #
                else:
                    continue
                for idx in tmpMissingIdx:
                    self.data[phenoName][idx] = avePheno
                # count inds with low valus and with high values
                lowValueCount, highValueCount = 0, 0
                for ind in self.data[phenoName]:
                    if ind <= avePheno:
                        lowValueCount += 1
                    else:
                        highValueCount += 1
                #
                propMissingDict[phenoName].append(lowValueCount)
                propMissingDict[phenoName].append(highValueCount)
        # {end: for}
        if verbose == False:
            return
        #
        outString = 'FIELD\tMISSING_RATE\t#BELOW_AVG\t#ABOVE_AVG\n'
        #
        for key, value in zip(propMissingDict.keys(), propMissingDict.values()):
            outString += (key+'\t')
            outString += ('%.4f' % (value[0]) +'\t')
            outString += (str(value[1]) +'\t')
            outString += (str(value[2]) +'\n')
        #
        warningMessage = '\nWARNING: %s covariate(s) dropped (--filter_missing %.2f):\n\n' % (len(deletedPhenosDict), threshold)
        #
        for key, value in zip(deletedPhenosDict.keys(), deletedPhenosDict.values()):
            warningMessage += '{:>30}\t{}\n'.format((key+'\t'), '%.2f%% missing' % (value * 100))
        #
        if len(deletedPhenosDict.keys()) == 0:
            return outString, ''
        else:
            return outString, warningMessage


    def output(self, fields, required=[]):
        ''' output selected fields of data. If exclude = True
        then the selected fields will be excluded from output
        and the rest of the fields are returned instead.
        return to a string ready to print to screen'''
        #
        fields = required + fields
        warningMessage = ''
        #
        del_fields = []
        for idx, field in enumerate(fields):
            if field not in self.data.keys():
                warningMessage += 'WARNING: %s not found in cleaned data\n' % field
                del_fields.append(field)
        ##
        fields = [x for x in fields if x not in del_fields ]
        outString = 'sample_name\t'+'\t'.join(fields)+'\n'
        #
        for idx, name in enumerate(self.data['sample_name']):
            outString += (name+'\t')
            for j, field in enumerate(fields):
                outString += (str(self.data[field][idx])+'\t') if j < len(fields)-1 else (str(self.data[field][idx]))
            #
            outString += '\n'
        ##
        return outString, warningMessage

    
    def dumpDataset(self, fileName, fields=[]):
        '''
        dump the entire dataset (or selected fields) to an external file
        If 'fields != []', dump selected fields only, otherwise dump all fields
        contained in the dataset.
        '''
        with open(fileName, 'w') as outFile:
            if len(fields) == 0:
                allFields = list(self.data.keys())
            else:
                allFields = fields
            if 'sample_name' in allFields:
                allFields.remove('sample_name')
            out, warn = self.output(allFields)
            outFile.write(out)
        outFile.close()
        sys.stderr.write('\nWrite intermediate dataset to %s\n' % fileName)
        if warn:
            sys.stderr.write("\n"+warn+"\n")
        return
        

    def select(self, fields, criteria, keepSelected=True, outFileName="selectedDataset.txt", sampleSize=0, removeThese=[], keepThese=[], toCases=False, toControls=False, keepOriginal=False, byNames=[]):
        '''
        select individuals that have given values on specific fields.
        Note that length of fields must be equal to length of values, where
        only individuals that have all specified fields following given criteria
        will be selected out.
        params:
            fields - list of field names based on which to select individuals/samples
            criteria - list of criteria with the ith element corresponds to the selection criteria of the ith element in 'fields'
            keepSelected - True (default) to keep those that are satisfying all criteria as selected samples; False to have them removed and return those that are not satisfying all criteria as samples.
            outFileName - output file name to save the selected individuals
            sampleSize - required sample size for selected individuals, if it is less than remained individuals, return a warning message and remove more remained individuals to meet the sample size requirement; otherwise restore some removed one to fill in the remained group to meet the sample size requirement; if it is equal to 0 (default) do nothing.
            removeThese - list of files that contains sample names to be removed under any circumstance
            keepThese - list of files that contains sample names to be kept under all circumstances
            toCases - False (default). If True, set fields[0] in selected individuals to 1 (1 for case).
            toControls - False (default). If True, set fields[0] in selected individuals to 0 (0 for control).
            keepOriginal -False (default). If True, keep original values of fields[0] by adding a new column named by 'fields[0]_original' to the dataset.
            byNames - list of files that contains sample names, which will be selected while ignoring all other conditions
        '''
        ## select individuals by sample names if list(byNames) != []
        listNames = [[x.rstrip() for x in open(fi).readlines()] for fi in byNames]
        listKeep = set(list(chain(*listNames)))
        if len(listKeep) > 0:
            listRemove = [ind for ind in self.data['sample_name'] if ind not in listKeep]
            self.removeSamples(listRemove)
            #
            sys.stdout.write("\nThe following individuals (N = %d) have been selected:\n" % (len(listKeep)))
            for idx, ind in enumerate(listKeep):
                if idx == 5:
                    sys.stdout.write("...\n(%d sample names omitted)" % (len(listKeep)-5))
                    break
                sys.stdout.write(ind+'\n')
            sys.stdout.write('\n')
            sys.stdout.write("\nNumber of individuals selected: %d\n" % len(listKeep))
            ## output the entire remaining dataset
            self.dumpDataset(outFileName)
            return
        ##
        
        flag = 0
        ##
        if len(fields) != len(criteria):
            sys.stderr.write("\nERROR: Number of traits must be equal to number of values to use phenoman select command\n\n")
            flag = 1
        ##
        ## construct lists of sample names from lists of file names 'removeThese' and 'keepThese'
        listRemove = [[x.rstrip() for x in open(fi).readlines()] for fi in removeThese]
        listKeep = [[x.rstrip() for x in open(fi).readlines()] for fi in keepThese]
        removeTheseInds = list(chain(*listRemove))
        keepTheseInds = list(chain(*listKeep))
        ## return error if there is any individual appearing to be in both lists
        tmp = [i[0] for i in product(removeTheseInds, keepTheseInds) if len(set(i)) == 1]
        if len(tmp) > 0:
            sys.stderr.write("\nERROR: the following sample names cannot appear in both 'removethese' and 'keepthese' lists of files: %s\n\n" % ', '.join(tmp))
            flag = 1
        ##
        ## return error if required sample size 'sampleSize' is smaller than len(keepTheseInds)
        if sampleSize is not 0 and sampleSize < len(keepTheseInds):
            sys.stderr.write("\nERROR: required sample size, --samplesize (%d), must be greater than or equal to number of specified individuals to be kept, --keepthese (%d)\n\n" % (sampleSize, len(keepTheseInds)))
            flag = 1
        ##
        ## return error if required sample size 'sampleSize' + len(removeTheseInds) is greater than total number of individuals
        if sampleSize is not 0 and sampleSize + len(removeTheseInds) > len(self.data['sample_name']):
            sys.stderr.write("\nERROR: required sample size, --samplesize (%d), + number of specified individuals to be removed, --removethese (%d), must be less than or equal to the entire sample size (%d)\n\n" % (sampleSize, len(removeTheseInds), len(self.data['sample_name'])))
            flag = 1
        ##
        ## return error if both 'toCases' and 'toControls' are True
        if toCases is True and toControls is True:
            sys.stderr.write("\nERROR: choose to set remained individuals to be either cases (--tobecases) or controls (--tobecontrols) on trait1, cannot specify both\n\n")
            flag = 1
        ##
        if flag:
            sys.exit(0)
        # select/remove samples according to criteria
        criteriaType = [] # 1 for single value, 2 for a range of numbers, 3 for a list of strings to choose from
        for idx, ind in enumerate(criteria):
            tmp = ind.split('|')
            if len(tmp) > 1:
                criteriaType.append(3)
                criteria[idx] = tmp
            else:
                tmp = ind.split(',')
                if len(tmp) == 2:
                    criteriaType.append(2)
                    criteria[idx] = [float(i) for i in tmp]
                else:
                    criteriaType.append(1)
                    try:
                        criteria[idx] = int(ind)
                    except:
                        try:
                            criteria[idx] = float(ind)
                        except:
                            continue
        #
        ## find samples that are not following criteria
        indsToRemove = []
        for idx, name in enumerate(self.data['sample_name']):
            for field, criType, cri in zip(fields, criteriaType, criteria):
                if criType == 1:
                    if not isEqual(self.data[field][idx], cri):
                    #if self.data[field][idx] != cri:
                        indsToRemove.append(name)
                        break
                elif criType == 2:
                    if not (type(self.data[field][idx]) in [float, int]) or self.data[field][idx] < cri[0] or self.data[field][idx] > cri[1]:
                        indsToRemove.append(name)
                        break
                elif criType == 3:
                    types = typeOfValues_selectModule(cri)
                    cri = [i(j) for i,j in zip(types, cri)]
                    if self.data[field][idx] not in cri:
                        indsToRemove.append(name)
                        break
                else:
                    continue
        ## remove those samples if 'keepSelected' is True, otherwise keep those samples and remove all the others
        if keepSelected:
            # keep and/or remove user-specified samples
            indsToRemove = self._keepOrRemoveGivenInds(indsToRemove, removeTheseInds, keepTheseInds)
            self._removeInds(indsToRemove=indsToRemove, sampleSize=sampleSize, removeTheseInds=removeTheseInds, keepTheseInds=keepTheseInds)
        else:
            tmpIndsToRemove = [ind for ind in self.data['sample_name'] if ind not in indsToRemove]
            # keep and/or remove user-specified samples
            indsToRemove = self._keepOrRemoveGivenInds(tmpIndsToRemove, removeTheseInds, keepTheseInds)
            self._removeInds(indsToRemove=indsToRemove, sampleSize=sampleSize, removeTheseInds=removeTheseInds, keepTheseInds=keepTheseInds)
        ##
        ## set fields[0] in remained inds to 1 (cases) or 0 (controls)
        if toCases:
            self.setCoding(coding=1, phenotype=fields[0], keepOriginal=keepOriginal)
        if toControls:
            self.setCoding(coding=0, phenotype=fields[0], keepOriginal=keepOriginal)
        ##
        ## output the entire remaining dataset
        self.dumpDataset(outFileName)
        return
        
        
    def setCoding(self, coding, phenotype, keepOriginal):
        '''
        Set 'phenotype' in remained individuals to be 'coding'.
        If 'keepOriginal' is True, keep original values of 'phenotype' by
        adding a new column named by 'phenotype_original' to the dataset.
        '''
        sampleNames = self.data['sample_name']
        oriPheno = copy.deepcopy(self.data[phenotype])
        self.data[phenotype] = [coding] * len(sampleNames)
        #
        if keepOriginal:
            self.data[phenotype+'_original'] = oriPheno
        #
        return
    
    
    def _keepOrRemoveGivenInds(self, indsToRemove, removeTheseInds, keepTheseInds):
        '''
        add samples that are in list 'removeTheseInds' to list 'indsToRemove' and take those that are in list 'keepTheseInds' out of list 'indsToRemove'
        '''
        allIndNames = self.data['sample_name']
        # remove specified individuals
        tmp = [ind for ind in removeTheseInds if (ind in allIndNames and ind not in indsToRemove)]
        indsToRemove.extend(tmp)
        # keep specified individuals
        tmp = [ind for ind in keepTheseInds if (ind in allIndNames and ind in indsToRemove)]
        [indsToRemove.remove(ind) for ind in tmp]
        #
        return indsToRemove
    
    
    def _removeInds(self, indsToRemove, sampleSize, removeTheseInds, keepTheseInds):
        '''
        subfunction of primary function 'self.select', see param list of func 'select' for more details
        '''
        # check sampleSize
        if sampleSize > len(self.data['sample_name']):
            sys.stderr.write("\nERROR: required sample size after selection (%d) cannot be greater than the sample size of the original input dataset (%d)\n" % (sampleSize, len(self.data['sample_name'])))
            sys.exit(0)
        #
        remainedSampleSize = len(self.data['sample_name']) - len(indsToRemove)
        # compare remainedSampleSize with required 'sampleSize'
        if sampleSize in [0, remainedSampleSize]:
            # choose all if 'sampleSize == 0'
            pass
        elif sampleSize < remainedSampleSize:
            sys.stderr.write("\nWARNING: selected/remained number of indivduals (%d) is greater than the required sample size (%d), some individuals are randomly chosen from remained ones to be further removed to meet the sample size requirement. If you need to keep all remained individuals please rerun the last command WITHOUT specifying --samplesize argument\n\n" % (remainedSampleSize, sampleSize))
            numDiff = remainedSampleSize - sampleSize
            remainedInds = [ind for ind in self.data['sample_name'] if ind not in indsToRemove]
            # keep those that are in 'keepTheseInds'
            [remainedInds.remove(ind) for ind in keepTheseInds if ind in remainedInds]
            indsToFurtherRemove = random.sample(remainedInds, numDiff)
            indsToRemove += indsToFurtherRemove
            #
            sys.stderr.write("WARNING: these individuals (N = %d) are randomly chosen to add in the removal list:\n" % (len(indsToFurtherRemove)))
            for idx, ind in enumerate(indsToFurtherRemove):
                if idx == 5:
                    sys.stderr.write("...\n(%d sample names omitted)" % (len(indsToFurtherRemove)-5))
                    break
                sys.stderr.write(ind+'\n')
            sys.stderr.write('\n')
            #
        else: # sampleSize > remainedSampleSize, need to restore some samples from removed ones to remained ones
            numDiff = sampleSize - remainedSampleSize
            # avoid restoring those that are in 'removeTheseInds'
            indsToRemoveLess = copy.deepcopy(indsToRemove)
            [indsToRemoveLess.remove(ind) for ind in removeTheseInds if ind in indsToRemoveLess]
            restoredInds = random.sample(indsToRemoveLess, numDiff)
            sys.stderr.write("\nWARNING: The following individuals ought to be removed according to the criteria but are randomly chosen from removed ones to be kept in the remained dataset in order to meet the sample size requirement (%d). If you still want them to be removed, please rerun the last command WITHOUT specifying --samplesize argument\n\n" % sampleSize)
            sys.stderr.write("WARNING: these individuals (N = %d) are restored:\n" % (len(restoredInds)))
            for idx, ind in enumerate(restoredInds):
                if idx == 5:
                    sys.stderr.write("...\n(%d sample names omitted)" % (len(restoredInds)-5))
                    break
                sys.stderr.write(ind+'\n')
            sys.stderr.write('\n')
            #
            [indsToRemove.remove(ind) for ind in restoredInds]
        ##
        self.removeSamples(indsToRemove)
        #
        sys.stdout.write("\nThe following individuals (N = %d) have been removed:\n" % (len(indsToRemove)))
        for idx, ind in enumerate(indsToRemove):
            if idx == 5:
                sys.stdout.write("...\n(%d sample names omitted)" % (len(indsToRemove)-5))
                break
            sys.stdout.write(ind+'\n')
        sys.stdout.write('\n')
        sys.stdout.write("\nNumber of individuals removed: %d\n" % len(indsToRemove))
        sys.stdout.write("Number of individuals remained/selected: %d\n" % len(self.data['sample_name']))   
        
    
    def massage(self, field, lower, upper, outFileName, percentile):
        '''
        Reset field values of individuals that are out of specified bounds into boundary values.
            field - trait of interest
            lower - lower limit
            upper - upper limit
            percentile - boolean parameter, if or not 'lower and upper' are values in percentiles
            outFileName - output file name to save the intermediate data set
        '''
        # remove individuals that are missing values on the trait 'field'
        pheno_missing = self.findMissing(field)
        sys.stdout.write('There are %d individuals that are missing values on phenotype %s and are removed\n\n' % (len(pheno_missing), str(field)))
        self.removeSamples(pheno_missing)
        # find the values on the lower and/or upper limits and replace those that are out of limits with limit values
        if percentile:
            if lower is None:
                lower = 0
            if upper is None:
                upper = 1
            limitValues = self.findPercentile(field, [lower, upper])
        else:
            if lower is None:
                lower = min(self.data[field])
            if upper is None:
                upper = max(self.data[field])
            limitValues = [lower, upper]
        #
        print ("sample_name\t%s\twinsorized_%s" % (field, field))
        numInds = 0
        for idx in range(len(self.data[field])):
            if self.data[field][idx] > limitValues[1]:
                print ("%s\t%.3f\t%.3f" % (self.data['sample_name'][idx], self.data[field][idx], limitValues[1]))
                self.data[field][idx] = limitValues[1]
                numInds += 1
            elif self.data[field][idx] < limitValues[0]:
                print ("%s\t%.3f\t%.3f" % (self.data['sample_name'][idx], self.data[field][idx], limitValues[0]))
                self.data[field][idx] = limitValues[0]
                numInds += 1
            else:
                continue
        sys.stdout.write('There are %d individuals that have been winsorized on phenotype %s\n\n' % (numInds, field))
        ##
        # output the entire data set
        allFields = list(self.data.keys())
        allFields.remove('sample_name')
        outFile = open(outFileName, 'w')
        out, warn = self.output(allFields)
        outFile.write(out)
        outFile.close()
        if warn:
            sys.stderr.write(warn + '\n')
        #
        return
        #pass
    
    def scale(self, field):
        '''
        Apply scaling to self.data[field] (x-mean) and add an additional field named
        by 'field_scaled' and modify self.data[field_scaled]
        '''
        # note: need to remove all missing values before performing scaling
        self.data[field+"_scaled"] = copy.deepcopy(self.data[field])
        mean = calMean(self.data[field])
        for idx, value in enumerate(self.data[field+"_scaled"]):
            self.data[field+"_scaled"][idx] = value - mean
        sys.stdout.write('Add a new column %s_scaled to the dataset\n\n' % field)
        return
        
    def standardize(self, field):
        '''
        Apply standardization to self.data[field] (x-mean)/std and add an
        additional field named by 'field_standardized' and modify
        self.data[field_standardized]
        '''
        # note: need to remove all missing values before performing standardization
        self.data[field+"_standardized"] = copy.deepcopy(self.data[field])
        mean = calMean(self.data[field])
        std = calStd(self.data[field])
        for idx, value in enumerate(self.data[field+"_standardized"]):
            self.data[field+"_standardized"][idx] = (value - mean) / std
        sys.stdout.write('Add a new column %s_standardized to the dataset\n\n' % field)
        return
        
    def normalize(self, field):
        '''
        Apply normalization to self.data[field], (x-x_min)/(x_max-x_min), and add an
        additional field named by 'field_normalized' and modify
        self.data[field_normalized]
        '''
        # note: need to remove all missing values before performing normalization
        self.data[field+"_normalized"] = copy.deepcopy(self.data[field])
        minimum = min(self.data[field])
        maximum = max(self.data[field])
        for idx, value in enumerate(self.data[field+"_normalized"]):
            self.data[field+"_normalized"][idx] = (value - minimum) / (maximum - minimum)
        sys.stdout.write('Add a new column %s_normalized to the dataset\n\n' % field)
        return
    
    def qnormalize(self, field, bounds):
        '''
        Apply Gaussian quantile normalization to self.data[field], and add an
        additional field named by 'field_qnormalized' and modify
        self.data[field_qnormalized]
        '''
        lp, hp = bounds[0], bounds[1]
        if not (lp < hp and 0 < lp < 1 and 0 < hp < 1):
            sys.stderr.write("\nERROR: arg --qnormalize should be specified by a range of probabilities, e.g. 0.025 0.975\n")
        # note: need to remove all missing values before performing Gaussian quantile normalization
        self.data[field+"_qnormalized"] = copy.deepcopy(self.data[field])
        tmp = copy.deepcopy(self.data[field])
        tmp.sort()
        inc = (hp-lp)/(len(tmp)-1)
        for idx, value in enumerate(self.data[field+"_qnormalized"]):
            incIndex = tmp.index(value)
            prob = lp + inc * incIndex
            self.data[field+"_qnormalized"][idx] = qnorm(prob)
        sys.stdout.write('Add a new column %s_qnormalized to the dataset\n\n' % field)
        return
    
    
    #######
    
 

from subprocess import PIPE, Popen

def set_count(words):
    '''return value: uniq elements in list and their respective counts'''
    c = Counter(words)
    r1, r2 = list(zip(*sorted(c.most_common())))
    return list(r1), list(r2)

class CombinedDummyCoder:
    '''Combine & dummy code covariates for ESP analysis'''
    def __init__(self, param, fn):
        self.fn = fn
        self.sample_col = param.sample_col
        self.display_name = param.header
        self.mp = param.merger_rule
        try:
            self.cols = list(map(int, param.cols))
        except:
            raise ValueError("column ids should be integers")
        self.order_by = param.order_by
        self.v = param.verbose
        self.as_letter = param.as_letter
        self.use_r = param.use_r
        self.style = param.style
        if self.use_r and self.style == 'triangle':
            sys.stderr.write("WARNING: R will only offer 'diagonal' style coding. Disable --use_r if you want to use triangle style\n")
        self.dat, self.guest_dat = self.read(self.fn, self.sample_col, self.cols, self.order_by)
        self.merger_rule = self.merge()
        self.mapping, self.cdat = self.convert(self.dat)

    def run(self):
        if self.use_r:
            rstr = self.codeR(self.cdat)
            dat = self.R(rstr)
        else:
            dat = self.code(self.cdat, self.mapping, self.style)
        out = [[y.replace('"', '') for y in x.split(',')] for x in dat.strip().split('\n')]
        coding, counts = set_count(['\t'.join(x[1:]) for x in out])
        self.codings = {x[0]:'N = {:<8}'.format(y) + x[1:] for x, y in zip(coding, counts)}
        if self.v:
            sys.stderr.write("\nCombined dummy codes:\n")
            self.writecoding(self.mapping, self.codings)
            if self.merger_rule:
                sys.stderr.write("\nMerged groups:\n")
                self.writemerger(self.merger_rule)
        # revert names to raw strings
        if not self.as_letter:
            rmapping = {}
            for k in self.mapping.keys():
                rmapping[self.mapping[k]] = '+'.join(k)
            for idx, item in enumerate(out):
                out[idx][1] = rmapping[out[idx][1]]
        #
        head = ['sample_name', '{}'.format(self.display_name)]
        for i in range(1, len(out[0]) - 1):
            head.append('{}_'.format(self.display_name) + str(i))
        output = '\t'.join(head) + '\t' + '\t'.join(self.guest_dat['sample_name'])
        if len(out) + 1 != len(self.guest_dat):
            raise ValueError("BUG: covariate data not matching dummy codes")
        for item in out:
            output += '\n' + '\t'.join(item) + '\t' + '\t'.join(self.guest_dat[item[0]])
            
        return output

    def read(self, fn, scol, cols, oby):
        '''read input data and return a list of data like [[sample_name, (col1, col2 ...)] ...]
        where the tuple is sorted by the col specified'''
        if scol <= 0:
            raise ValueError("sample name column not valid")
        if not cols:
            raise ValueError("columns for covariates not specified")
        try:
            dat = [tuple([x.strip().split('\t')[scol-1]] + [y for idx, y in enumerate(x.strip().split('\t')) if idx+1 in cols]) for x in open(fn).readlines()]
            guest_dat = {x.strip().split('\t')[scol-1]:[y for idx, y in enumerate(x.strip().split('\t')) if not idx+1 == scol] for x in open(fn).readlines()}
        except Exception as e:
            raise ValueError("failed to read input data: {}".format(e))
        #
        if dat[0][0] != 'sample_name':
            raise ValueError('The first line of input phenotype file must be header line with sample_name column named "sample_name"')
        dat = dat[1:]
        if oby:
            if oby not in cols:
                raise ValueError("wrong order_by option")
            idx = cols.index(oby)
            dat = sorted(dat, key=lambda d: d[idx+1])
        return [[x[0], tuple(x[1:])] for x in dat], guest_dat

    def merge(self):
        '''merge combined dummy variable groups'''
        if not self.mp:
            return {}
        mdict = OrderedDict()
        for item in self.mp:
            try:
                lhs, rhs = item.split('=')
                if not (rhs.startswith('(') or rhs.endswith(')')):
                    raise ValueError('() not match!')
                else:
                    rhs = rhs[1:-1].split(',')
                lhs = lhs.strip()
                rhs = [tuple(x.strip().split('+')) for x in rhs]
                for r in rhs:
                    if len(r) != len(self.cols):
                        raise ValueError('Invalid -m/--merger_rule. Please make sure "?? + ?? + ..." matches the --cols input.')
                    mdict[r] = tuple([lhs])
            except Exception as e:
                raise ValueError('Invalid input: "{}".\nExample: --merger_rule "a=(rr + jj, tt + uu, pp + mm)" "b=(xx + ii, yy + zz)"'.format(' '.join(self.mp)))
        seen = []
        for i in range(len(self.dat)):
            if self.dat[i][1] in list(mdict.keys()):
                seen.append(self.dat[i][1])
                self.dat[i][1] = mdict[self.dat[i][1]]
        for item in list(mdict.keys()):
            if not item in seen:
                raise ValueError('Invalid input: group combination "{}" not found in input data'.format(' + '.join(item)))
        return mdict

    def convert(self, dat):
        '''return a list of list of strings ready be formatted for R'''
        d = OrderedDict()
        i = 65
        # find unique coding
        for item in dat:
            if item[1] not in d.keys():
                d[item[1]] = str(chr(i))
                i += 1
        dat = [[x[0], d[x[1]]] for x in dat]
        return d, dat

    def writecoding(self, mapping, codings):
        for item in mapping.keys():
            sys.stderr.write('{:<30}-->{:>5}\t{:>8}\n'.format('+'.join(item), mapping[item], codings[mapping[item]].replace('\t', ' ')))
        sys.stderr.write('\n')

    def writemerger(self, merger):
        for item in merger.keys():
            sys.stderr.write('{:<40}>>>{:>15}\n'.format('+'.join(item), ' '.join(merger[item])))
        sys.stderr.write('\n')

    def code(self, dat, mapping, style):
        '''dummy coding by diagnal or triangle style'''
        codings = {}
        l = len(mapping) - 1
        i = -1
        for item in mapping.values():
            codings[item] = ['0'] * l
            if i >= 0:
                if style == 'triangle': codings[item][:(i+1)] = ['1'] * (i+1)
                else: codings[item][i] = '1'
            i += 1
        return '\n'.join([','.join(x + codings[x[1]]) for x in dat])

    def codeR(self, dat):
        '''return a string ready to pipe to R'''
        rstr = 'colnames(dat) <- c("names", "groups")\n'
        rstr += 'mat <- model.matrix(lm(runif(length(dat$groups))~dat$groups))\n'
        rstr += 'write.table(cbind(as.character(dat$names), as.character(dat$groups), apply(as.matrix(mat[,-1]), 2, as.character)), stdout(), row.names=F, col.names=F, sep=",")'
        if self.v and self.use_r:
            sys.stdout.write('== R implementation ==\n' + rstr+'\n======================\n\n')
        rstr = 'dat <- as.data.frame(cbind(' + ', '.join(['c(' + ', '.join(list(map(repr, x))) + ')' for x in list(zip(*dat))]) + '))\n' + rstr
        return rstr

    def R(self, rstr):
        ''' Run R command and parse the output '''
        tc = Popen(["R", '--subordinate', '--vanilla'], stdin = PIPE, stdout = PIPE, stderr = PIPE)
        out, error = tc.communicate(rstr.encode(sys.getdefaultencoding()))
        if (tc.returncode):
            raise ValueError(" (exception captured from R) \n{0}".format(error.decode(sys.getdefaultencoding())))
        else:
            if error:
                sys.stderr.write("WARNING message from R: \n{}\n".format(error.decode(sys.getdefaultencoding())))
        return out.decode(sys.getdefaultencoding())


#if __name__ == '__main__':
#   DUP_AA = [x.rstrip().split('\t') for x in open('AAdups.txt').readlines()[1:]]
#   DUP_EA = [x.rstrip().split('\t') for x in open('EAdups.txt').readlines()[1:]]
#   DUP_LIST = DUP_AA + DUP_EA
   ##
#          #obj = PhenData('AA.txt')
#          #
#          #sys.exit(1)
#          #
#          ### test 1: if # samples in AA.txt + # samples in EA.txt =? # samples in ALL.txt
#          #objAA = PhenData('AA.txt')
#          #objEA = PhenData('EA.txt')
#          #objAll = PhenData('ALL.txt')
#          #print len(objAA.data['sample_name'])+len(objEA.data['sample_name']), len(objAll.data['sample_name'])
#          ### test 2: test PhenData.fillMissing func
#          #obj = PhenData('AA.txt')
#          #propMissing, deletedPheno = obj.fillMissing(verbose=None, threshold=1)
#          #
#          ## test 3: test PhenData.output
#          obj = PhenData('AA.txt')
#          # <1> remove missingness for the primary trait
#          missingOnes = obj.findMissing(field = 'esp_dbp_baseline')
#          obj.removeSamples(names = missingOnes)
#          output1 = obj.output(fields=['esp_dbp_baseline', 'esp_stroke_famhx_baseline', 'esp_waist_baseline', 'esp_lpa'])
#          # <2> remove duplicates
#          dups = obj.findDuplicates()
#          obj.removeSamples(names = dups)
#          output2 = obj.output(fields=['esp_dbp_baseline', 'esp_stroke_famhx_baseline', 'esp_waist_baseline', 'esp_lpa'])
#          # <3> remove outliers
#          outliers = obj.findOutliers(field = 'esp_dbp_baseline', criteria = obj.findPercentile(field='esp_dbp_baseline', percentiles=[.1, .9]))
#          obj.removeSamples(names = outliers)
#          output3 = obj.output(fields=['esp_dbp_baseline', 'esp_stroke_famhx_baseline', 'esp_waist_baseline', 'esp_lpa'])
#          # <4> fill missingness for covariates
#          obj.fillMissing(verbose=None, threshold=.8)
#          output4 = obj.output(fields=['esp_dbp_baseline', 'esp_stroke_famhx_baseline', 'esp_waist_baseline', 'esp_lpa'])
