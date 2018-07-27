#!/usr/bin/env python

from setuptools import setup

setup(name='source_assess_NFB_efficacy',
      description='Perform meta-analysis and systematic analysis of biases',
      url='https://github.com/AuroreBussalb/meta-analysis-statistical-tools',
      author='Aurore Bussalb',
      author_email='aurore.bussalb@mensiatech.com',
      packages=['source_assess_NFB_efficacy/meta_analysis', 'source_assess_NFB_efficacy/systematic_analysis_of_biases'],
      zip_safe=False)