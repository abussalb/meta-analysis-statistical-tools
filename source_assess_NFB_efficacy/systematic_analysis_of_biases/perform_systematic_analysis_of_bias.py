# -*- coding: utf-8 -*-

"""
.. module:: perform_systematic_analysis_of_bias
    :synopsis: module performing a SAOB
 
.. moduleauthor:: Aurore Bussalb <aurore.bussalb@mensiatech.com>
"""

import warnings
import pandas as pd
import numpy as np
import random
from sklearn.linear_model import LassoCV, LassoLarsIC
from sklearn.model_selection import LeaveOneOut
import statsmodels.api as sm
from sklearn import tree
import graphviz


def effect_size_within_subjects(mean_post_test_NFB, mean_pre_test_NFB, std_post_test_NFB, std_pre_test_NFB):
    """Computes effects sizes inside a group, this effect size reflects the evolution inside a group between pre and post test. The
    fomula used comes from Cohen J. (1988), *Statistical Power Analysis for the Behavioral Sciences*. 
    
    Parameters
    ----------
    n_NFB: int
        Number of patients included in the Neurofeedback (NFB) group.
    
    mean_post_test_NFB: float
        Mean score after the treatment (NFB here).
    
    mean_pre_test_NFB: float
        Mean score before the treatment (NFB here).
               
    std_pre_test_NFB: float
        Standard deviation of the mean score before the treatment (Neurofeedback here).

    std_post_test_NFB: float
        Standard deviation of the mean score before the treatment (Control here).                         

    Returns
    -------
    effect_size: float
        Value estimating the efficacy of NFB.
        If it's negative, the result is in favor of NFB.

    Notes
    -----
        Effect sizes computed for each study correspond to the effect sizes within subjects. Thus, the effect size is 
        computed from the pre and post test values of the treatment group. So here a control group is not neccessary.
    
    """

    Effect_size_NFB = (mean_post_test_NFB - mean_pre_test_NFB)/np.sqrt((std_pre_test_NFB**2 + std_post_test_NFB**2)/2)
    
    return Effect_size_NFB


def preprocess_factors(df):
    """Preprocesses factors before running the SAOB.

    Factors with too many missing values and with too many identical observations will be removed.
    Besides, values can be standardized. The categorial variables will be codded as dummies.

    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe containing all observations in rows, factors and also values to compute the effect size within subjects in columns. 
        It is obtained after the import of the csv file containing all data by ``import_csv_for_factors``.

    Returns
    -------
    X: pandas.DataFrame
        Preprocessed dataframe containing all observations in rows and factors in columns.
        Factors with too many missing values and with too many identical observations have been removed.
        Categorical variables are coded in dummies.
        Values here are standardized. 

    X_non_standardized: pandas.DataFrame
        Preprocessed dataframe containing all observations in rows and factors in columns.
        Factors with too many missing values and with too many identical observations have been removed.
        Categorical variables are coded in dummies.
        Values here are not standardized. 

    """

    # Remove factors with too few observations    
    df_number_of_nans = df.isnull().sum()
    columns_to_remove_nans = df_number_of_nans[(df_number_of_nans > round(len(df)*20/100) + 1)]
    df = df.drop(columns_to_remove_nans.index.values, axis=1)
    
    # Dataframe containing only factors
    X = df.drop(['mean_post_test_NFB', 'mean_pre_test_NFB','n_NFB', 'raters', 'score_name', 'std_post_test_NFB',
                 'std_pre_test_NFB', 'effect_size_NFB'], axis=1)

    # Turn into dummy variables the categorical variables 
    categorical_factors = list(set(X.columns) - set(X._get_numeric_data().columns))
    X = pd.get_dummies(X, columns=categorical_factors, drop_first=True)

    # Remove categorical factors too homogeneous  
    bool_cols = [col for col in X if X[col].dropna().value_counts().index.isin([0,1]).all()]
    X_categorical = X[bool_cols]
    X_categorical_count = X_categorical.apply(pd.value_counts)
    columns_to_remove_no = X_categorical_count.iloc[0][(X_categorical_count.iloc[0] > round(len(df)*80/100) + 1)]
    columns_to_remove_yes = X_categorical_count.iloc[1][(X_categorical_count.iloc[1] > round(len(df)*80/100) + 1)]
    columns_to_remove = columns_to_remove_no.index.tolist() + columns_to_remove_yes.index.tolist()
    X = X.drop(columns_to_remove, axis=1)  
        
    # Put -1 instead of NaNs in continuous factors
    X.fillna(value=-1, inplace=True)  

    # Standardization of the independent variables for WLS and Lasso
    X_non_standardized = X
    X = (X - X.mean())/X.std()

    return X, X_non_standardized 


def weighted_linear_regression(df, X, y):
    """Performs Weighted Least Squares.

    Dependent variable = effect size within subjects; independent variable = factors. 
    P independent variables and n observations.
    model: WXB = Wy; W (nxn) diagonal matrix of weights, y (nx1) column vector
    of dependent variables, X (nxP) matrix of independent variables, B (Px1) 
    column vector of coefficients.

    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe containing all observations in rows, factors in columns and also values to compute the effect size within subjects. 
        It is obtained after the import of the csv file containing all data by ``import_csv_for_factors``.

    X: pandas.DataFrame
        Preprocessed dataframe containing all observations in rows and factors in columns (the independent variables). 
        Factors with too many missing values and with too many identical observations have been removed. Besides, values have been standardized. 
        Categorical variables are coded in dummies.
        This dataframe is obtained thanks to the ``preprocess_factors function``.

    y: pandas.Series
        Effect size within subjects computed for each observation (the dependent variable).

    Returns
    -------
    summary: statsmodels.iolib.summary.Summary
        Summary of the WLS regression.
        In particular values of coefficients, associated pvalues, F-statistics, Prob(Omnibus), skew and kurtosis.

    """

    # Find the number of scales per study
    df['number_of_scales'] = df.index.value_counts()

    # Compute the weight of each ES: number_of_NFB_patients_in_the_study/number_of_scales_in_the_study
    df['weight'] = df['n_NFB']/df['number_of_scales']
    W = np.diag(df['weight'])

    # Get rank of the moment matrix and its condition number: it has to be full rank,
    # to have eigen values > 0  and a high condition number to be invertible
    rank_X = np.linalg.matrix_rank(X)
    X_transpose = X.transpose()
    W_transpose = W.transpose()
    rank_W = np.linalg.matrix_rank(W)
    moment_matrix = np.linalg.multi_dot([X_transpose,W_transpose,W,X])
    moment_matrix_rank = np.linalg.matrix_rank(moment_matrix)
    eigen_values_moment_matrix = np.linalg.eigvals(moment_matrix)
    condition_number = np.linalg.cond(moment_matrix)
    
    if moment_matrix.shape[0] == moment_matrix_rank:
        print('Moment matrix is invertible')
    else: 
        warnings.warn('Moment matrix is not invertible, be carefull while interpreting the results')
    
    # WLS
    weighted_regression = sm.WLS(y, sm.add_constant(X), weights=df['weight'])
    weighted_regression_fit = weighted_regression.fit()
    summary = weighted_regression_fit.summary()

    return summary
    
    
def ordinary_linear_regression(X, y):
    """Performs Ordinary Least Squares.

    Dependent variable = effect size within subjects; independent variable = factors. 
    P independent variables and n observations.
    model: XB = y; y (nx1) column vector of dependent variables, X (nxP) matrix
    of independent variables, B (Px1) column vector of coefficients.

    Parameters
    ----------
    X: pandas.DataFrame
        Preprocessed dataframe containing all observations in rows and factors in columns (the independent variables). 
        Factors with too many missing values and with too many identical 
        observations have been removed. Besides, values have been standardized. 
        Categorical variables are coded in dummies.
        This dataframe is obtained thanks to the ``preprocess_factors`` function.

    y: pandas.Series
        Effect size within subjects computed for each observation (the dependent variable).

    Returns
    -------
    summary_ols: statsmodels.iolib.summary.Summary
        Summary of the OLS regression.
        In particular values of coefficients, associated pvalues, F-statistics, Prob(Omnibus), skew and kurtosis.

    """

 	# Get rank of the moment matrix and its condition number: it has to be full rank,
    # to have eigen values > 0  and a high condition number to be invertible
    rank_X = np.linalg.matrix_rank(X)
    X_transpose = X.transpose()
    moment_matrix = np.dot(X_transpose,X)
    moment_matrix_rank = np.linalg.matrix_rank(moment_matrix)
    eigen_values_moment_matrix = np.linalg.eigvals(moment_matrix)
    condition_number = np.linalg.cond(moment_matrix)
    
    if moment_matrix.shape[0] == moment_matrix_rank:
        print('Moment matrix is invertible')
    else:
        warnings.warn('Moment matrix is not invertible, be carefull while interpreting the results')

    # Run the OLS 
    regression = sm.OLS(y, sm.add_constant(X))
    regression_fit = regression.fit()
    summary_ols = regression_fit.summary()     
     
    return summary_ols


def regularization_lassocv(X, y):
    """Performs Lasso linear model with iterative fitting along a regularization path.
    The best model is selected by cross-validation.

    Dependent variable = effect size within subjects; independent variable = factors. 
    P independent variables and n observations.

    Parameters
    ----------
    X: pandas.DataFrame
        Preprocessed dataframe containing all observations in rows and factors in columns (the independent variables). 
        Factors with too many missing values and with too many identical observations have been removed. Besides, values have been standardized. 
        Categorical variables are coded in dummies.
        This dataframe is obtained thanks to the ``preprocess_factors`` function.

    y: pandas.Series
        Effect size within subjects computed for each observation (the dependent variable).

    Returns
    -------
    coeff: pandas.DataFrame
        Results of the Lasso.
        Column with coefficients obtained after regularization and the names of the associated factors.

    mse_test: numpy.ndarray
        Mean square error for the test set on each fold, varying alpha.
        Shape (n_alphas, n_folds). 

    alphas: numpy.ndarray
        The grid of alphas used for fitting.
        Shape (n_alphas).

    alpha: float
        The amount of penalization chosen by cross validation.

    """

    # Cross validation (leave one out) to choose the tuning parameter alpha 
    loo = LeaveOneOut()
    lassocv = LassoCV(alphas = None, cv = loo) # leave one out method, internal cross validation: it performs cv on the 
                                               # data it receives
    
    # Fit all the data with the alpha found 
    lassocv.fit(X, y)
    alpha = lassocv.alpha_
    alphas = lassocv.alphas_
    mse_test = lassocv.mse_path_

    # Coefficients that don't explain the model are now reduced to exactly zero
    coeff = pd.DataFrame({'Factors':X.columns, 'Coefficients':lassocv.coef_})

    return coeff, mse_test, alphas, alpha


def regularization_lassoAIC(X, y):
    """Performs Lasso model fit with Lars using AIC for model selection.

    Dependent variable = effect size within subjects; independent variable = factors. 
    P independent variables and n observations.
    model: XB = y; y (nx1) column vector of dependent variables, X (nxP) matrix
    of independent variables, B (Px1) column vector of coefficients.

    Parameters
    ----------
    X: pandas.DataFrame
        Preprocessed dataframe containing all observations in rows and factors in columns (the independent variables). 
        Factors with too many missing values and with too many identical 
        observations have been removed. Besides, values have been standardized. 
        Categorical variables are coded in dummies.
        This dataframe is obtained thanks to the ``preprocess_factors`` function.

    y: pandas.Series
        Effect size within subjects computed for each observation (the dependent variable).

    Returns
    -------
    coeff_aic: pandas.DataFrame
        Results of the Lasso.
        Column with coefficients obtained after regularization and the names of the associated factors.

    """
    
    model = LassoLarsIC(criterion='aic') 
    model.fit(X, y) 
    coeff_aic = pd.DataFrame({'Factors':X.columns, 'Coefficients':model.coef_})

    return coeff_aic


def decision_tree(X_non_standardized, y): 
    """Computes a Decision tree.

    Non linear and hierarchical method.

    Parameters
    ----------
    X: pandas.DataFrame
        Preprocessed dataframe containing all observations in rows and factors in columns (the independent variables). 
        Factors with too many missing values and with too many identical 
        observations have been removed. Besides, values have been standardized. 
        Categorical variables are coded in dummies.
        This dataframe is obtained thanks to the ``preprocess_factors`` function.

    y: pandas.Series
        Effect size within subjects computed for each observation (the dependent variable).

    Returns
    -------
    decision_tree: pdf
        Decision tree obtained.

    """  

    # Decision tree (criterion: mean square error)
    clf = tree.DecisionTreeRegressor(criterion='mse', min_samples_leaf=6)
    clf.fit(X_non_standardized, y)
    score_decision_tree = clf.score(X_non_standardized, y)
    print('R² decision tree', score_decision_tree)
   
    # Visualization
    dot_data = tree.export_graphviz(clf, feature_names=X_non_standardized.columns, out_file=None)
    graph = graphviz.Source(dot_data)
    graph.render('decision_tree', view=True)