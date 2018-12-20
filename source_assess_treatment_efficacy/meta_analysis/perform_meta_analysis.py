# -*- coding: utf-8 -*-

"""
.. module:: perform_meta_analysis
    :synopsis: module performing a meta-analysis
 
.. moduleauthor:: Aurore Bussalb <aurore.bussalb@mensiatech.com>
"""

import numpy as np
import scipy.stats as scp
import pandas as pd
import warnings
import matplotlib.pyplot as plt


def _effect_size_ppc(n_treatment, n_control, mean_post_test_treatment, mean_pre_test_treatment, mean_pre_test_control, mean_post_test_control,
                    std_pre_test_treatment, std_pre_test_control):   
    """Computes the pre post control effect size (Scott B. Morris (2008), also called the effect size between "Estimating Effect Sizes From Pretest-Posttest Control Group Designs 
    and under a random effects model", Organizational Research Methods (Equation 8)).
    
    Parameters
    ----------
    n_treatment: int
        Number of patients included in the treatment group.

    n_control: int
        Number of patients included in the control group.

    mean_post_test_treatment: float
        Mean score after the treatment.

    mean_pre_test_treatment: float
        Mean score before the treatment.

    mean_pre_test_control: float
        Mean score before the treatment in the control group.

    mean_post_test_control: float
        Mean score after the treatment in the control group.
           
    std_pre_test_treatment: float
        Standard deviation of the mean score before the treatment.

    std_post_test_treatment: float
        Standard deviation of the mean score after the treatment.

    Returns
    -------
    effect_size: float
        Value estimating the efficacy of the treatment.
        If it's negative, the result is in favor of the treatment.
        
    """     

    S_within = np.sqrt(((n_treatment - 1)*std_pre_test_treatment**2 + (n_control - 1)*std_pre_test_control**2)/
                            (n_treatment + n_control - 2))
    d = ((mean_post_test_treatment - mean_pre_test_treatment) - (mean_post_test_control - mean_pre_test_control))/S_within
    
    # Correction factor for small sample size. This correction factor is close to 1 unless the degree of freedom is very
    # small (<10), see Borenstein, Introdution to meta-analysis, 2009.
    if (n_treatment + n_control - 2) < 10:  
        warnings.warn('Since the sample size is too small, a correction factor is applied to the effect size')
        correction_factor = 1 - (3/(4*(n_treatment + n_control - 2) - 1))
        effect_size = d*correction_factor 
    else:
        effect_size = d
    
    return effect_size


def _standard_error_effect_size(n_treatment, n_control, effect_size, pre_post_correlation):    
    """Scott B. Morris (2008) "Estimating Effect Sizes From Pretest-Posttest Control Group Designs and under 
    a random effects model", Organizational Research Methods (Equation 25).
    
    Parameters
    ----------
    n_treatment: int
        Number of patients included in the treatment group.

    n_control: int
    Number of patients included in the control group.

    effect_size: float
        Value estimating the efficacy of the treatment.
        If it's negative, the result is in favor of the treatment.

    pre_post_correlation: float
        Pearson correlation of the pre-test and post-test values (i.e the pooled within-groups Pearson correlation.
     
    Returns
    -------
    standard_error_ES: float 
        Standard error of the effect size.

    variance_ES: float
        Variance of the effect size.     
        
    """
    
    # Correction factor for small sample size. This correction factor is close to 1 unless the degree of freedom is very
    # small (<10), see Borenstein, Introdution to meta-analysis, 2009. 
    if (n_treatment + n_control - 2) < 10: 
        correction_factor = 1 - (3/(4*(n_treatment + n_control - 2) - 1))
        warnings.warn('Since the sample size is too small, a correction factor is applied to the variance of the effect size')
    else:
        correction_factor = 1
    
    # Variance 
    variance_ES = (2*(correction_factor**2)*(1 - pre_post_correlation)*((n_treatment + n_control)/
                               (n_treatment*n_control))*((n_treatment + n_control - 2)/(n_treatment + n_control - 4))*
                               (1 + ((effect_size**2)/(2*(1 - pre_post_correlation)*((n_treatment + n_control)/
                               (n_treatment*n_control))))) - effect_size**2)
    
    # Standard Error
    standard_error_ES = np.sqrt(variance_ES) 
    
    return standard_error_ES

   
def run_meta_analysis(df, scale_to_reverse=[], pre_post_correlation=0.5):
    """Performs a meta analysis with the formulae described in Scott B. Morris (2008) "Estimating Effect Sizes From Pretest-
    Posttest Control Group Designs and under a random effects model", *Organizational Research Methods* and in Borenstein (2009)
    *Introduction to meta-analysis*. These formulae are the same as the ones used in Cortese et al., 2016. 

    A negative effect size favours the treatment. 
  
    Parameters
    ----------
    df: pandas.DataFrame
        Parents, teachers or clinicians ratings required to perform the meta-analysis.
        This dataframe corresponds to one of those obtained with the ``import_csv_for_meta_analysis`` module.
        If you want to run the meta-analysis on parent assessments enter ``df_values_parents``, to run it on teacher assessments
        enter ``df_values_teachers``, and to run on clinicians assessments, run ``df_values_clinicians``
        Each row corresponds to a study, the disease symptoms are assessed by parents, teachers, or clinicians.
        Columns are: mean_post_test_treatment, mean_post_test_control, mean_pre_test_treatment, mean_pre_test_control, n_treatment, 
        n_control, std_post_test_treatment, std_post_test_control, std_pre_test_treatment, std_pre_test_control, raters for each study.  

    scale_to_reverse: list of str, optional
        List of strings listing the clinical scales having a positive correlation with symptoms of the disease; 
        i.e increasing when a patient gets better.
    
    pre_post_correlation: float, default = 0.5
        Pearson correlation of the pre-test and post-test values (i.e the pooled within-groups Pearson correlation). Set to 0.5 by
        default (see Cuijpers et al., 2016 and Balk et al., 2012 "Empirical Assessment of Within-Arm Correlation Imputation in Trials 
        of Continuous Outcomes").                  

    Returns
    -------
    df_results_per_study: pandas.DataFrame 
        Results per study.
        Rows of the dataframe correspond to the studies, columns correspond to the effect size of the study, its standard 
        error, its 95% confidence interval, and the weight of the study.

    df_results: pandas.DataFrame
        Global results.
        It contains the summary effect, its 95% confidence interval, its variance, its standard error, its p-value, 
        the between studies variance (Tau²), the heterogeneity (I²), its p-value, and the Chi2 value.

    Notes
    -----
        Effect sizes computed for each study correspond to the effect sizes between subjects. Thus, the studies included in the meta-analysis 
        must be controlled and provide pre and post scores for treatment and control groups.
        
    """
   
    # Creation of the dataframe for total results
    index = ['Results']
    df_results = pd.DataFrame(index=index)
    
    
    # Compute the effect size    
    df['effect_size'] = df[
                ['n_treatment', 'n_control', 'mean_post_test_treatment', 
                 'mean_pre_test_treatment', 'mean_pre_test_control', 
                 'mean_post_test_control', 'std_pre_test_treatment', 'std_pre_test_control']
                          ].apply(lambda row:_effect_size_ppc(**row), axis=1)
    
    
    # Compute the standard error of the effect size
    df['standard_error_ES'] = df[
                        ['n_treatment', 'n_control', 'effect_size']
                                ].apply(lambda row:_standard_error_effect_size(row['n_treatment'], row['n_control'],
                                                         row['effect_size'], pre_post_correlation), axis=1)

     
    # Check if all the scales measure the desease severity the same way (high score = more symptomps) and homogenize
    for scale_name in scale_to_reverse:
        df['effect_size'][ df['score_name']==scale_name ] *= -1
    
    
    # All the following equations come from M. Borenstein and L. Hedges (2009) Introduction to Meta-Analysis
    
    # 95% Confidence interval (Equations 8.3 and 8.4)    
    df['confidence_interval_of_the_ES'] = df[
                                ['effect_size', 'standard_error_ES']].apply(lambda row: (
                                 row['effect_size'] - 1.96*row['standard_error_ES'],
                                 row['effect_size'] + 1.96*row['standard_error_ES']), axis=1)
    
    
    # Compute the inverse of the variance = weight under a fixed effect model (Equation 11.2)
    df['weight_fixed_model'] = 1/(df['standard_error_ES']**2)

    
    # Computation of Tau²: between studies variance
    
    ## Compute degrees of freedom (Equation 12.4)
    degrees_of_freedom = len(df.index) - 1 
    
    ## Compute Q (Equation 12.3)  
    Q = (df['weight_fixed_model']*df['effect_size']**2).sum() - ((df['weight_fixed_model']*df['effect_size']).sum())**2/df['weight_fixed_model'].sum()
    df_results['Chi2'] = Q
        
    ## P value of the heterogeneity
    # To know if heterogeneity is statistically significant, we can use Q and degrees of freedom
    # Null hypothesis: all studies share a common effect size
    # Under the null hypothesis, Q will follow a central chi-squared distribution
    df_results['p-value Heterogeneity'] = 1 - scp.chi2.cdf(Q, degrees_of_freedom)

    ## Compute C (Equation 12.5)
    C = df['weight_fixed_model'].sum() - ((df['weight_fixed_model']**2).sum()/df['weight_fixed_model'].sum())
    
    ## Tau² (Equation 12.2)
    # When Tau2 is negative, we put it at zero (this negative value is due to sampling issues, 
    # when the observed dispersion is less than we would expect by chance, see Borenstein)
    Tau2 = (Q - degrees_of_freedom)/C    
    if Tau2 < 0:
        Tau2 = 0
    df_results['Tau2'] = Tau2 

        
    # Compute the weight of each study under a random effects model
    ## Compute the weights (Equation 12.6)
    df['weight'] = 1/(df['standard_error_ES']**2 + Tau2)
    ## In percentage
    df['percentage_weight'] = (df['weight']*100)/df['weight'].sum()
    

    # Summary effect (Equation 12.7)
    df_results['Summary Effect'] = (df['effect_size']*df['percentage_weight']).sum()/df['percentage_weight'].sum()

    # Variance and SE of the summary effect (Equations 12.8 and 12.9)
    df_results['Variance Summary Effect'] = 1/df['weight'].sum()
    df_results['Standard Error Summary Effect'] = np.sqrt(df_results['Variance Summary Effect'])

    # 95% Confidence interval (Equations 12.10 and 12.11)
    df_results['95% Confidence Interval of the Summary Effect'] = df_results[
                                         ['Summary Effect', 'Standard Error Summary Effect']].apply(lambda row: (
                                           row['Summary Effect'] - 1.96*row['Standard Error Summary Effect'], 
                                           row['Summary Effect'] + 1.96*row['Standard Error Summary Effect']), axis=1)
   
    # P value for the summary effect (Equations 12.12 and 12.14)
    # Null hypothesis: control group and treatment group have no different effect
    z = df_results['Summary Effect']/df_results['Standard Error Summary Effect']
    df_results['p-value'] = 2*(1 - scp.norm.cdf(abs(z)))
       
        
    # Heterogeneity (Equation 16.9)
    I2 = (((Q - degrees_of_freedom))/Q)*100
    if I2 < 0:
        I2 = 0
    df_results['Heterogeneity'] = I2 
    
    
    # Creation of the dataframe with results by studies
    df_results_per_study = pd.DataFrame({'Year': df['year'],
                                        'Effect size': df['effect_size'],
                                        'Standard Error of the ES': df['standard_error_ES'],
                                        '95% Confidence interval of the ES': df['confidence_interval_of_the_ES'],
                                        'Weight': df['percentage_weight']},
                                         index=df.index)

    return df_results_per_study, df_results, df['effect_size']

if __name__ == '__main__':
    meta_analysis('values_total_meta_analysis.csv', 'Parents') 


def forest_plot(df_results_per_study, df_results):
    """Creates a forest plot.
    
    Parameters
    ----------
    df_results_per_study: pandas.DataFrame
        Results per study.
        Dataframe obtained after performing the meta-analysis with ``run_meta_analysis``.
        Rows of the dataframe correspond to the studies, columns correspond to the effect size of the study, its standard 
        error, its 95% confidence interval, and the weight of the study.

    df_results: pandas.DataFrame
        Global results.
        It contains the summary effect, its 95% confidence interval, its variance, its standard error, its p-value, 
        the between studies variance (Tau²), the heterogeneity (I²), its p-value, and the Chi2 value.
        
    Returns
    -------
    forest_plot: matplotlib.figure
        Graphical representation of the meta-analysis' results.
        Representation of the effect size and its 95% confidence interval for each study.
        
    """

    # Sort data so that studies with bigger effect size are in the top of the forest plot
    df_results_per_study = df_results_per_study.sort_values(df_results_per_study.columns[1], ascending=[True])
    
    # Conversion to lists for the plotting
    ES = df_results_per_study['Effect size'].tolist()
    weight = df_results_per_study['Weight'].tolist()
    names = df_results_per_study.index.tolist()
    names = [i[0] for i in names]

    # Preparing for the plotting
    ## Confidence Interval
    lower_limit = []
    upper_limit = []
    for confidence_interval in df_results_per_study['95% Confidence interval of the ES']:
        lower_limit.append(confidence_interval[0])
        upper_limit.append(confidence_interval[1])
    lower_limit_summary = df_results['95% Confidence Interval of the Summary Effect'][0][0]
    upper_limit_summary = df_results['95% Confidence Interval of the Summary Effect'][0][1]
    
    # Add the confidence interval of the summary effect to others
    lower_limit.extend([lower_limit_summary])
    lower_limit.reverse() # the summary effect must be at the bottom
    upper_limit.extend([upper_limit_summary]) 
    upper_limit.reverse()
    
    # Add the summary effect at the other effects size
    names.append('Summary Effect')
    names.reverse()
    ES.extend(df_results['Summary Effect'])
    ES.reverse()  
    
    # Make the effect size representation more visible (squares are bigger)
    weight = [i * 5 for i in weight]

    # Graphic
    y = np.array(range(1,len(names)+1))
    forest_plot = plt.figure()
    plt.yticks(y, names)
    # Vertical line in zero
    plt.axvline(0, color = 'k')   
    # Plot Confidence Interval
    for i in range(0,len(names)): 
        plt.plot([lower_limit[i], upper_limit[i]], [y[i],y[i]], color = 'g')
    # Plot effect sizes
    plt.scatter(ES[1:len(names)], y[1:len(names)], s=weight[1:len(names)],
                marker = 's', color = 'b')
    plt.scatter(ES[0], y[0], s=100, marker = 'D', color = 'b')
    plt.xlabel('Effect size')
    plt.title('Standard Mean Difference, 95% Confidence Interval', fontweight = "bold")
    
    return forest_plot 

