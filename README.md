# meta-analysis-statistical-tools

This repository describes and offers a module for running traditional meta-analysis of randomized controlled trials (RCT) and a systematic analysis of biases (SAOB), an approach meant to identify the individual contribution of methodological and technical biases to the intervention efficacy (see *Bussalb et al., 2019*). 

# Description

The repository has the following structure:
* source-assess-treatment-efficacy
  * meta_analysis: package to perform meta-analysis and plots
  * systematic_analysis_of_biaises: package to perform systematic analysis of biaises
* documentation
  * package-documentation: package documentation in html
  * article-html: original article in html 
  * article-latex-source: tex source of original article
* example
  * meta-analysis: notebook that uses the meta_analysis package
  * systematic_analysis_of_biaises: notebook that uses the systematic_analysis_of_biases package 

Please feel free to re-use, suggest improvement, and contribute. 
Please cite *Bussalb et al., 2019*.

# Installation 

## Obtain source code

Obtain a local copy of this repository:

```git clone git@https://github.com/AuroreBussalb/meta-analysis-statistical-tools.git``` 

## Install requirements

Install the requirements noted in the ```requirements.txt``` file in the virtual environment ```treatment-efficacy-env```:

```conda env create -f environment.yml``` (if graphviz has not been installed, run ```pip install graphviz``` after activating the environment) 

or create an environment:

```conda create --name treatment-efficacy-env python=3```

activate the environment: 

```activate treatment-efficacy-env```

and install the requirements:

```pip install -r requirements.txt```

# Setup

1. Activate the environment:

```activate treatment-efficacy-env```

2. Run from a command line from the root directory (where the setup.py file is):

```pip install .```

3. You can uninstall this package by doing:

```pip uninstall source_assess_treatment_efficacy```

