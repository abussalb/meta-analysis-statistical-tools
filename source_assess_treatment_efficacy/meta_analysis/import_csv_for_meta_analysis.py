# -*- coding: utf-8 -*-

"""
.. module:: import_csv_for_meta_analysis
    :synopsis: module importing data for meta-analysis
 
.. moduleauthor:: Aurore Bussalb <aurore.bussalb@mensiatech.com>
"""

import warnings
import pandas as pd


def _common_read(csv_file, raters):
    """Reads data from a csv file containing parents and sometimes teachers and clinicians severity assessments of a disease. This csv file 
    contains all the data required to perform a meta-analysis.
   
    Parameters
    ----------
    csv_file: str
        Name or localisation of the csv file that contains all the values of clinical trials required to perform the SAOB.
        The csv file must have a specific form: 
        - eight first columns: Author, Year, Score Name, Number of patients, Raters (Parents, Teachers, or Clinicians), Time (pre and post), Mean, Std;
        - other columns correspond to factors to analyze;
        - each study has at least 2 lines: Author name, Year, Score Name, Number of patients, Raters (Parents, Teachers, or Clinicians) repeated twice, 
          Time (pre and post), Mean (at pre-test and post-test), Std (at pre-test and post-test). If Teachers assessments are available, 
          two lines must be added (same in the case of clinicians);
        - sometimes, for a study, several clinical scales are available, they are all entered in the csv file;
        - for each author, the 2 lines of the Raters column have to be filled as follows: Parents, Parents, if teachers' assessment is available
          2 more lines are added (Teachers, Teachers), same in the case of clinicians (Clinicians, Clinicians);
        - for each author, the 2 lines of the Time column have to be filled as follows: pre, post (pattern repeated one more time if teachers or
          clinicians assessment is available);
        - Mean and Std correspond to the clinical score extracted from studies at pre-test and post-test.
               
    raters: str, 'Teachers', 'Parents' or 'Clinicians'
        Person assessing the disease symptoms.
                           
    Returns
    -------
    df_values: pandas.DataFrame
        Dataframe used to perform the SAOB.
        Each row corresponds to a study rated by a specific rater on a specific scale.
        Columns are: mean_post_test_treatment, mean_post_test_control, mean_pre_test_treatment, mean_pre_test_control, n_treatment.
        n_control, std_post_test_treatment, std_post_test_control, std_pre_test_treatment, std_pre_test_control, raters for each study.    

    """
    
    # Import a csv file as a dataframe   
    df = pd.read_csv(csv_file)

    # Name of studies with an evaluation and its associated score
    indices_name_studies = df.loc[ (df['Group'] == "Treatment")
                & (df['Time'] == "pre")
                & (df['Raters'] == raters),
                ['Mean', 'Std']
                ].isnull().sum(axis=1).index
    name_studies = df.loc[indices_name_studies, 'Author']
    score_name = df.loc[indices_name_studies, 'Score Name']
    

    # Extract treatment values
    treatment_indices_pre = df[
                    (df['Group'] == "Treatment") 
                  & (df['Time'] == "pre") 
                  & (df['Raters'] == raters)
                  ].index
    n_treatment = df.loc[treatment_indices_pre, 'Number of patients']
    year = df.loc[treatment_indices_pre, 'Year']
    mean_pre_test_treatment = df.loc[treatment_indices_pre, 'Mean']
    std_pre_test_treatment = df.loc[treatment_indices_pre, 'Std']
    
    treatment_indices_post = df[
                    (df['Group'] == "Treatment") 
                  & (df['Time'] == "post") 
                  & (df['Raters'] == raters)
                  ].index
    mean_post_test_treatment = df.loc[treatment_indices_post, 'Mean']
    std_post_test_treatment = df.loc[treatment_indices_post, 'Std']
    
    
    # Extract Control values
    control_indices_pre = df[
                    (df['Group'] == "Control") 
                  & (df['Time'] == "pre") 
                  & (df['Raters'] == raters)
                  ].index
    n_control = df.loc[control_indices_pre, 'Number of patients']
    year_control = df.loc[control_indices_pre, 'Year']
    mean_pre_test_control = df.loc[control_indices_pre, 'Mean']
    std_pre_test_control = df.loc[control_indices_pre, 'Std']
    
    control_indices_post = df[
                    (df['Group'] == "Control") 
                  & (df['Time'] == "post") 
                  & (df['Raters'] == raters)
                  ].index
    mean_post_test_control = df.loc[control_indices_post, 'Mean']
    std_post_test_control = df.loc[control_indices_post, 'Std']

   
    # Creation of the data frame containing the results
    df_values = pd.DataFrame({'year': year.tolist(),
                            'n_control': n_control.tolist(),
                            'mean_pre_test_control': mean_pre_test_control.tolist(),
                            'mean_post_test_control': mean_post_test_control.tolist(),
                            'std_pre_test_control': std_pre_test_control.tolist(),
                            'std_post_test_control': std_post_test_control.tolist(),
                            'n_treatment': n_treatment.tolist(),
                            'score_name': score_name.tolist(),
                            'mean_pre_test_treatment': mean_pre_test_treatment.tolist(),
                            'mean_post_test_treatment': mean_post_test_treatment.tolist(),
                            'std_pre_test_treatment': std_pre_test_treatment.tolist(),
                            'std_post_test_treatment': std_post_test_treatment.tolist()},
                            index=[name_studies])

    df_values = df_values[ df_values.isna().sum(axis=1)==0 ]
    
    return df_values


def import_csv(csv_file, raters=''): 
    """Imports data from a csv file containing parents and sometimes teachers and clinicians severity assessments of a disease. This csv file 
    contains all the data required to perform a meta-analysis. It is possible to import parents' ratings, teachers' or clinicians' only but also the three
    of them.
    
    Parameters
    ----------
    csv_file: str
        Name or localisation of the csv file that contains all the values of RCT with a pretest posttest design required to perform a meta analysis.
        The csv file must have a specific form: 

        - nine columns: Author, Year, Group, Score Name, Number of patients, Raters, Time, Mean, Std;
        - each study has 4 lines if only parent assessments are available, 8 if teacher assessments are provided too, and 12 if clinicians' 
          are available: the Author name (repeated 4, 8, or 12 times),
          Year (repeated 4 or 8 times), Group (Treatment or control), Score Name, Raters (Parents, Teachers, or Clinicians), Time (pre and post), 
          Mean and Std;
        - for each author, the 4 lines of the Group column have to be filled as follows: Treatment, Treatment, Control, Control 
          (pattern repeated one more time if teachers or clinicians assessment is available);
        - for each author, the 4 lines of the Raters column have to be filled as follows: Parents, Parents, Parents, Parents and if teachers' assessment 
          is available 4 more lines are added (Teachers, Teachers, Teachers, Teachers), same for clinicians;
        - for each author, the 4 lines of the Time column have to be filled as follows: pre, post, pre, post (pattern repeated one moretime if teachers assessment is available);
        - Mean and Std correspond to the clinical score extracted from studies.
           
    raters: str, optional 'Teachers', 'Parents' or 'Clinicians'
        Person assessing a disease symptoms, if no raters are precised then all values will be returned. 
                           
    Returns
    -------
    df_values_parents: pandas.DataFrame
        Parents' ratings required to perform the meta-analysis.
        It will be returned if ``raters = 'Parents'``.
        Each row corresponds to a study, the disease symptoms are assessed by parents.
        Columns are: mean_post_test_treatment, mean_post_test_control, mean_pre_test_treatment, mean_pre_test_control, n_treatment.
        n_control, std_post_test_treatment, std_post_test_control, std_pre_test_treatment, std_pre_test_control, raters for each study.  

    df_values_teachers: pandas.DataFrame
        Teachers' ratings required to perform the meta-analysis.
        It will be returned if ``raters = 'Teachers'``.
        Each row corresponds to a study, the disease symptoms are assessed by teachers.
        Columns are: mean_post_test_treatment, mean_post_test_control, mean_pre_test_treatment, mean_pre_test_control, n_treatment.
        n_control, std_post_test_treatment, std_post_test_control, std_pre_test_treatment, std_pre_test_control for each study. 

    df_values_clinicians: pandas.DataFrame
        Clinicians' ratings required to perform the meta-analysis.
        It will be returned if ``raters = 'Clinicians'``.
        Each row corresponds to a study, the disease symptoms are assessed by clinicians.
        Columns are: mean_post_test_treatment, mean_post_test_control, mean_pre_test_treatment, mean_pre_test_control, n_treatment. 
        n_control, std_post_test_treatment, std_post_test_control, std_pre_test_treatment, std_pre_test_control for each study. 

    Notes
    -----   
        The three dataframes will be returned if no rater is precised.
        
    """
    
    # Parents
    if raters == 'Parents':
        df_values_parents = _common_read(csv_file, raters='Parents')
        return df_values_parents
        
    # Teachers
    elif raters == 'Teachers':
        df_values_teachers = _common_read(csv_file, raters='Teachers')
        return df_values_teachers 

    # Clinicians
    elif raters == 'Clinicians':
        df_values_clinicians = _common_read(csv_file, raters='Clinicians')
        return df_values_clinicians 
        
    # All
    elif raters == '':
        df_values_parents = _common_read(csv_file, raters='Parents')
        df_values_teachers = _common_read(csv_file, raters='Teachers')
        df_values_clinicians = _common_read(csv_file, raters='Clinicians')
        
        return df_values_parents, df_values_teachers, df_values_clinicians 
    
    # Error
    else:
        warnings.warn('Raters are either teachers, parents, or clinicians')
              
if __name__ == '__main__':
    import_csv_for_meta_analysis('values_inattention_meta_analysis.csv')   



