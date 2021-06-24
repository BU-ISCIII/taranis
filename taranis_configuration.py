#!/usr/bin/env python3

import os

###  Settings configuration  for logging
taranis_dir = os.path.dirname(os.path.realpath(__file__))
###LOGGING_CONFIGURATION = taranis_dir
###LOGGING_CONFIGURATION = taranis_dir + '/logging_config.ini'
LOGGING_CONFIGURATION = os.path.join(taranis_dir,'logging_config.ini')
print(LOGGING_CONFIGURATION)

## Settings configuration for create schema functionality

