# -*- coding: utf-8 -*-

"""
.. module:: import_csv_for_meta_analysis
    :synopsis: module importing data for meta-analysis
 
.. moduleauthor:: Aurore Bussalb <aurore.bussalb@mensiatech.com>
"""

import warnings
import csv
import pandas as pd


def _common_read(csv_file, raters):
    """Reads data from a csv file containing parents and sometimes teachers' severity assessments of ADHD in children. This csv file 
    contains all the data required to perform a meta-analysis.     
    
    Args:
        csv_file (str): name or localisation of the csv file that contains all the values of RCT with a pretest posttest design
        required to perform a meta analysis.
            The csv file must have a specific form: 
            - nine columns: Author, Year, Group, Score Name, Number of patients, Raters, Time, Mean, Std;
            - each study has 4 lines if only parent assessments are available, 8 if teacher assessments are provided too: the Author 
            name (repeated 4 or 8 times), Year (repeated 4 or 8 times), Group (NFB or control), Score Name, Raters (Parents or Teachers), Time 
            (pre and post), Mean and Std;
            - for each author, the 4 lines of the Group column have to be filled as follows: NFB, NFB, Control, Control (pattern repeated one more
            time if teachers assessment is available);
            - for each author, the 4 lines of the Raters column have to be filled as follows: Parents, Parents, Parents, Parents and if teachers' assessment
            is available 4 more lines are added (Teachers, Teachers, Teachers, Teachers);
            - for each author, the 4 lines of the Time column have to be filled as follows: pre, post, pre, post (pattern repeated one more
            time if teachers assessment is available);
            - Mean and Std correspond to the clinical score extracted from studies.
               
        raters (str): 'Teachers' or 'Parents'.
            Person assessing ADHD symptoms.
                           
    Returns:
        df_values (pandas.DataFrame): dataframe used to perform the meta-analysis.
            Each row corresponds to a study
            columns correspond to mean_post_test_NFB, mean_post_test_control, mean_pre_test_NFB, mean_pre_test_control, n_NFB, 
            n_control, std_post_test_NFB, std_post_test_control, std_pre_test_NFB, std_pre_test_control, raters for each study. 
        
    """
    
    # Remove NaNs
    remove_nans =  lambda iterable: [x for x in iterable if not pd.isnull(x)]
    
    # Import a csv file as a dataframe   
    df = pd.read_csv(csv_file)

    # Name of studies with an evaluation 
    name_studies = []
    score_name = []
    indices = df[(df['Group'] == "NFB") & (df['Time'] == "pre") & (df['Raters'] == raters)].index 
    for item in indices:
        if not pd.isnull(df.loc[item,'Mean']) and not pd.isnull(df.loc[item,'Std']):
            name_studies.append(df.loc[item,'Author'])
            score_name.append(df.loc[item,'Score Name'])    
    
    
    # Extract NFB values
    NFB_indices_pre = df[(df['Group'] == "NFB") & (df['Time'] == "pre") & (df['Raters'] == raters)].index
    n_NFB = df.loc[NFB_indices_pre, 'Number of patients']
    year = df.loc[NFB_indices_pre, 'Year']
    mean_pre_test_NFB = df.loc[NFB_indices_pre, 'Mean']
    std_pre_test_NFB = df.loc[NFB_indices_pre, 'Std']
    
    NFB_indices_post = df[(df['Group'] == "NFB") & (df['Time'] == "post") & (df['Raters'] == raters)].index
    mean_post_test_NFB = df.loc[NFB_indices_post, 'Mean']
    std_post_test_NFB = df.loc[NFB_indices_post, 'Std']
    
    ## Delete NaN values
    n_NFB = remove_nans(n_NFB)
    year = remove_nans(year)
    mean_pre_test_NFB = remove_nans(mean_pre_test_NFB)
    mean_post_test_NFB = remove_nans(mean_post_test_NFB)
    std_pre_test_NFB = remove_nans(std_pre_test_NFB)
    std_post_test_NFB = remove_nans(std_post_test_NFB)
    
    
    # Extract Control values
    control_indices_pre = df[(df['Group'] == "Control") & (df['Time'] == "pre") & (df['Raters'] == raters)].index
    n_control = df.loc[control_indices_pre, 'Number of patients']
    year_control = df.loc[control_indices_pre, 'Year']
    mean_pre_test_control = df.loc[control_indices_pre, 'Mean']
    std_pre_test_control = df.loc[control_indices_pre, 'Std']
    
    control_indices_post = df[(df['Group'] == "Control") & (df['Time'] == "post") & (df['Raters'] == raters)].index
    mean_post_test_control = df.loc[control_indices_post, 'Mean']
    std_post_test_control = df.loc[control_indices_post, 'Std']
    
    ## Delete NaN values
    n_control = remove_nans(n_control)
    mean_pre_test_control = remove_nans(mean_pre_test_control)
    mean_post_test_control = remove_nans(mean_post_test_control)
    std_pre_test_control = remove_nans(std_pre_test_control)
    std_post_test_control = remove_nans(std_post_test_control)

   
    # Creation of the data frame containing the results
    df_values = pd.DataFrame({'year': year,
                            'n_control': n_control,
                            'mean_pre_test_control': mean_pre_test_control,
                            'mean_post_test_control': mean_post_test_control,
                            'std_pre_test_control': std_pre_test_control,
                            'std_post_test_control': std_post_test_control,
                            'n_NFB': n_NFB,
                            'score_name': score_name,
                            'mean_pre_test_NFB': mean_pre_test_NFB,
                            'mean_post_test_NFB': mean_post_test_NFB,
                            'std_pre_test_NFB': std_pre_test_NFB,
                            'std_post_test_NFB': std_post_test_NFB,
                            'raters' : raters},
                            index=[name_studies])
    
    return df_values


def import_csv(csv_file, raters=''): 
    """Imports data from a csv file containing parents and sometimes teachers' severity assessments of ADHD in children. This csv file 
    contains all the data required to perform a meta-analysis. It is possible to import parents' ratings or teachers' only but also both.
    
    Args:
        csv_file (str): name or localisation of the csv file that contains all the values of RCT with a pretest posttest design required to perform a meta analysis.
            The csv file must have a specific form: 
             - nine columns: Author, Year, Group, Score Name, Number of patients, Raters, Time, Mean, Std;
             - each study has 4 lines if only parent assessments are available, 8 if teacher assessments are provided too: the Author 
               name (repeated 4 or 8 times), Year (repeated 4 or 8 times), Group (NFB or control), Score Name, Raters (Parents or Teachers), Time 
               (pre and post), Mean and Std;
             - for each author, the 4 lines of the Group column have to be filled as follows: NFB, NFB, Control, Control (pattern repeated one more
               time if teachers assessment is available);
             - for each author, the 4 lines of the Raters column have to be filled as follows: Parents, Parents, Parents, Parents and if teachers' assessment
               is available 4 more lines are added (Teachers, Teachers, Teachers, Teachers);
             - for each author, the 4 lines of the Time column have to be filled as follows: pre, post, pre, post (pattern repeated one more
               time if teachers assessment is available);
             - Mean and Std correspond to the clinical score extracted from studies.
               
        raters (str): optional 'Teachers' or 'Parents'.
            Person assessing ADHD symptoms, if no raters are precised then all values will be returned. 
                           
    Returns:
        df_values_parents (pandas.DataFrame): parents's ratings required to perform the meta-analysis.
            It will be returned if ``raters = 'Parents'``.
            Each row corresponds to a study, ADHD symptoms are assessed by parents
            columns correspond to mean_post_test_NFB, mean_post_test_control, mean_pre_test_NFB, mean_pre_test_control, n_NFB, 
            n_control, std_post_test_NFB, std_post_test_control, std_pre_test_NFB, std_pre_test_control, raters for each study.  
    
        df_values_teachers (pandas.DataFrame): teachers's ratings required to perform the meta-analysis.
            It will be returned if ``raters = 'Teachers'``.
            Each row corresponds to a study, ADHD symptoms are assessed by teachers
            columns correspond to mean_post_test_NFB, mean_post_test_control, mean_pre_test_NFB, mean_pre_test_control, n_NFB, 
            n_control, std_post_test_NFB, std_post_test_control, std_pre_test_NFB, std_pre_test_control for each study. 

        .. note:: both dataframes will be returned if no rater is precised.
        
    """
    
    # Parents
    if raters == 'Parents':
        df_values_parents = _common_read(csv_file, raters='Parents')
        return df_values_parents
        
    # Teachers
    elif raters == 'Teachers':
        df_values_teachers = _common_read(csv_file, raters='Teachers')
        return df_values_teachers 
        
    # All
    elif raters == '':
        df_values_parents = _common_read(csv_file, raters='Parents')
        df_values_teachers = _common_read(csv_file, raters='Teachers')
        
        return df_values_parents, df_values_teachers 
    
    # Error
    else:
        warnings.warn('Raters are either teachers or parents')
              
if __name__ == '__main__':
    import_csv_for_meta_analysis('values_inattention_meta_analysis.csv')   



