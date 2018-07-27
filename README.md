# meta-analysis-statistical-tools

This repository describes and offers a module for running traditional meta-analysis of randomized controlled trials (RCT) and an systematic analysis of biases (SAOB) an approach meant to identify the individual contribution of methodological and technical biases to the intervention efficacy. 

# Description

The repository has the following structure:
* source-assess-NFB-efficacy
  * meta_analysis: package to perform meta-analysis and plots
  * systematic_analysis_of_biaises: package to perform systematic analysis of biaises
* tests: unitary tests and validation compatible with pytest
* documentation
  * package-documentation: html package documentation
  * article-html: html of the original article 
  * article-latex-source: tex source of original article
* example
  * meta-analysis: notebook that uses the meta_analysis package
  * systematic_analysis_of_biaises: notebook that uses the systematic_analysis_of_biases package 

Please feel free to re-use, suggest improvement, and contribute. 
Please cite *Bussalb et al. (complete ref here)*.

# Installation 

## Obtain source code

Obtain a local copy of this repository:

```git clone git@https://github.com/AuroreBussalb/meta-analysis-statistical-tools.git``` 

## Install requirements

Install the requirements noted in the ```requirements.txt``` file in the virtual environment assess-NFB-efficacy-env:

```conda env create -f environment.yml``` 

or create an environment:

```conda create --name NFB-efficacy-env python=3```

activate the environment: 

```activate NFB-efficacy-env```

and install the requirements:

```pip install -r requirements.txt```

# Setup

1. Activate the environment:

```activate NFB-efficacy-env```

2. Run from a command line from the root directory (where the setup.py file is):

```pip install .```

or 

```python setup.py install```

3. You can uninstall this package by doing:

```pip uninstall source_assess_NFB_efficacy```

