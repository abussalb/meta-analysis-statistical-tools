# -*- coding: utf-8 -*-

"""
.. module:: import_csv_for_factors
    :synopsis: module importing data for SAOB

.. moduleauthor:: Aurore Bussalb <aurore.bussalb@mensiatech.com>
"""

import pandas as pd


def _common_read(csv_file, raters):
    """Reads data from a csv file containing parents and sometimes teachers' severity assessments of ADHD in children and factor values. 
    This csv file contains all the data required to perform the SAOB. 
   
    Parameters
    ----------
    csv_file: str
        Name or localisation of the csv file that contains all the values of clinical trials required to perform teh SAOB.
        The csv file must have a specific form: 
        - eight first columns: Author, Year, Score Name, Number of patients, Raters (Parents or Teachers), Time (pre and post), Mean, Std;
        - other columns correspond to factors to analyze;
        - each study has at least 2 lines: Author name, Year, Score Name, Number of patients, Raters (Parents or Teachers) repeated twice, 
        Time (pre and post), Mean (at pre-test and post-test), Std (at pre-test and post-test). If Teachers assessments are available,
        two lines must be added;
        - sometimes, for a study, several clinical scales are available, they are all entered in the csv file;
        - for each author, the 2 lines of the Raters column have to be filled as follows: Parents, Parents, if teachers' assessment
        is available 2 more lines are added (Teachers, Teachers);
        - for each author, the 2 lines of the Time column have to be filled as follows: pre, post (pattern repeated one more
        time if teachers assessment is available);
        - Mean and Std correspond to the clinical score extracted from studies at pre-test and post-test.
               
    raters: str, 'Teachers' or 'Parents'
        Person assessing ADHD symptoms.
                           
    Returns
    -------
    df_values: pandas.DataFrame
        Dataframe used to perform the SAOB.
        Each row corresponds to a study rated by a specific rater on a specific scale.
        Columns correspond to mean_post_test_NFB, mean_post_test_control, mean_pre_test_NFB, mean_pre_test_control, n_NFB, 
        n_control, std_post_test_NFB, std_post_test_control, std_pre_test_NFB, std_pre_test_control, raters for each study.    

    """
    
    # Import a csv file as a dataframe   
    df = pd.read_csv(csv_file)

    # Name of studies with an evaluation 
    name_studies = []
    score_name = []
    indices = df[(df['Time'] == "pre") & (df['Raters'] == raters)].index 
    for item in indices:
        if not pd.isnull(df.loc[item,'Mean']) and not pd.isnull(df.loc[item,'Std']):
            name_studies.append(df.loc[item,'Author'])
            score_name.append(df.loc[item,'Score Name'])    
    
    # Extract NFB values 
    NFB_indices_pre = df[(df['Time'] == "pre") & (df['Raters'] == raters)].index
    n_NFB = df.loc[NFB_indices_pre, 'Number of patients']
    mean_pre_test_NFB = df.loc[NFB_indices_pre, 'Mean']
    std_pre_test_NFB = df.loc[NFB_indices_pre, 'Std']
    
    NFB_indices_post = df[(df['Time'] == "post") & (df['Raters'] == raters)].index
    mean_post_test_NFB = df.loc[NFB_indices_post, 'Mean']
    std_post_test_NFB = df.loc[NFB_indices_post, 'Std']
    
    # Extract factors
    pblind = df.loc[NFB_indices_pre, 'Probably Blind']
    number_of_sessions = df.loc[NFB_indices_pre, 'Number of sessions']
    SMR = df.loc[NFB_indices_pre, 'SMR']
    theta_up = df.loc[NFB_indices_pre, 'Theta up']
    beta_up_central = df.loc[NFB_indices_pre, 'Beta up central']
    theta_down = df.loc[NFB_indices_pre, 'Theta down']
    beta_up_frontal = df.loc[NFB_indices_pre, 'Beta up frontal']
    SCP = df.loc[NFB_indices_pre, 'SCP']
    on_drugs = df.loc[NFB_indices_pre, 'On drugs during NFB/assessments']
    age_min = df.loc[NFB_indices_pre, 'Age min']
    age_max = df.loc[NFB_indices_pre, 'Age max']
    randomization = df.loc[NFB_indices_pre, 'Randomization']
    IRB = df.loc[NFB_indices_pre, 'Institutional Review Board']
    transfer_phase = df.loc[NFB_indices_pre, 'Transfer phase']
    transfer_card = df.loc[NFB_indices_pre, 'Transfer card']
    EOG_correction = df.loc[NFB_indices_pre, 'EOG correction']
    artifact_correction_based_on_amplitude = df.loc[NFB_indices_pre, 'Artifact correction based on amplitude']
    thresholding = df.loc[NFB_indices_pre, 'Thresholding']
    session_pace = df.loc[NFB_indices_pre, 'Session pace (per week)']
    session_length = df.loc[NFB_indices_pre, 'Session length (min)']
    treatment_length = df.loc[NFB_indices_pre, 'Treatment length (weeks)']
    more_than_one_active_electrode = df.loc[NFB_indices_pre, '>1 active electrode']
    EEG_quality = df.loc[NFB_indices_pre, 'EEG quality']
    control_group = df.loc[NFB_indices_pre, 'Control group']
   
    # Creation of the data frame containing the results
    df_values = pd.DataFrame({'n_NFB': n_NFB.tolist(),
                            'score_name': score_name,
                            'mean_pre_test_NFB': mean_pre_test_NFB.tolist(),
                            'mean_post_test_NFB': mean_post_test_NFB.tolist(),
                            'std_pre_test_NFB': std_pre_test_NFB.tolist(),
                            'std_post_test_NFB': std_post_test_NFB.tolist(),
                            'raters' : raters,
                            'pblind' : pblind.tolist(),
                            'number_of_sessions' : number_of_sessions.tolist(),
                            'SMR' : SMR.tolist(),
                            'theta_up' : theta_up.tolist(),
                            'theta_down' : theta_down.tolist(),
                            'beta_up_central': beta_up_central.tolist(),
                            'beta_up_frontal': beta_up_frontal.tolist(),
                            'SCP': SCP.tolist(), 
                            'on_drugs': on_drugs.tolist(),
                            'age_min': age_min.tolist(),
                            'age_max': age_max.tolist(),
                            'randomization': randomization.tolist(),
                            'IRB': IRB.tolist(),
                            'transfer_phase': transfer_phase.tolist(),
                            'transfer_card': transfer_card.tolist(),
                            'EOG_correction': EOG_correction.tolist(),
                            'artifact_correction_based_on_amplitude': artifact_correction_based_on_amplitude.tolist(), 
                            'thresholding': thresholding.tolist(),
                            'session_pace': session_pace.tolist(),
                            'session_length': session_length.tolist(),
                            'treatment_length': treatment_length.tolist(),
                            'more_than_one_active_electrode': more_than_one_active_electrode.tolist(),
                            'EEG_quality': EEG_quality.tolist(),
                            'control_group': control_group.tolist()},
                             index=[name_studies])
    
    return df_values


def import_csv(csv_file):
    """Imports data from a csv file containing parents and sometimes teachers' severity assessments of ADHD in children and factor values. 
    This csv file contains all the data required to perform the SAOB. 
    
    Parameters
    ----------
    csv_file: str
        Name or localisation of the csv file that contains all the values of clinical trials required to perform teh SAOB.
        The csv file must have a specific form: 

        - eight first columns: Author, Year, Score Name, Number of patients, Raters (Parents or Teachers), Time (pre and post), Mean, Std;
        - other columns correspond to factors to analyze;
        - each study has at least 2 lines: Author name, Year, Score Name, Number of patients, Raters (Parents or Teachers) repeated twice, 
          Time (pre and post), Mean (at pre-test and post-test), Std (at pre-test and post-test). If Teachers assessments are available, 
          two lines must be added;
        - sometimes, for a study, several clinical scales are available, they are all entered in the csv file;
        - for each author, the 2 lines of the Raters column have to be filled as follows: Parents, Parents, if teachers' assessment is available
          2 more lines are added (Teachers, Teachers);
        - for each author, the 2 lines of the Time column have to be filled as follows: pre, post (pattern repeated one more time if teachers 
          assessment is available);
        - Mean and Std correspond to the clinical score extracted from studies at pre-test and post-test.

    Returns
    -------
    df_values_parents: pandas.DataFrame
        Parents's ratings required to perform the SAOB.
        Each row corresponds to a study and a specific clinical scale, ADHD symptoms are assessed by parents,
        Columns correspond to mean_post_test_NFB, mean_post_test_control, mean_pre_test_NFB, mean_pre_test_control, n_NFB, 
        n_control, std_post_test_NFB, std_post_test_control, std_pre_test_NFB, std_pre_test_control, raters and factors. 

    df_values_teachers: pandas.DataFrame
        Teachers's ratings required to perform the SAOB.
        Each row correspond to a study and a specific clinical scale,, ADHD symptoms are assessed by teachers,
        Columns correspond to mean_post_test_NFB, mean_post_test_control, mean_pre_test_NFB, mean_pre_test_control, n_NFB, 
        n_control, std_post_test_NFB, std_post_test_control, std_pre_test_NFB, std_pre_test_control, raters and factors. 
        
    """
        
    # Import parents ans teachers values
    df_values_parents = _common_read(csv_file, raters='Parents')
    df_values_teachers = _common_read(csv_file, raters='Teachers')
        
    return df_values_parents, df_values_teachers 
        
if __name__ == '__main__':
    import_csv_for_factors('values_total_meta_analysis_all_factors.csv') 
    

