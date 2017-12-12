# this script produces a version of the catalog of papers without the keywords used for searching each theory,
# in lower case, missing symbols and grammatical markers.

import re
import nltk

from nltk.tokenize import word_tokenize

text = open('catalog.txt', 'r').read()
new_text = text.lower()  # lower case
new_text = re.sub(r'\([Cc]\).+\.', '', new_text)  # it removes copyright notices.
new_text = re.sub(r'[^\w\s]s', '', new_text)  # removes possessives
new_text = re.sub(r'[^\w\s]', '', new_text)  # removes punctuantion and symbols

new_text = re.sub(r'actr', '', new_text)  # controls act-r cognition
new_text = re.sub(r'bayesian', '', new_text)  # controls bayesian
new_text = re.sub(r'computational', '', new_text)  # controls computational
new_text = re.sub(r'connectionism', '', new_text)  # controls bayesian
new_text = re.sub(r'distributed', '', new_text)  # controls distributed
new_text = re.sub(r'dynamical', '', new_text)  # controls dynamical systems
new_text = re.sub(r'ecological', '', new_text)  # controls ecological
new_text = re.sub(r'embodied', '', new_text)  # controls embodied
new_text = re.sub(r'enactive', '', new_text)  # controls enactive
new_text = re.sub(r'parallel\sdistributed\sprocessing', '', new_text)  # controls PDP.
new_text = re.sub(r'PDP', '', new_text)  # also controls PDP.

with open('clean_catalog.txt', 'w') as outputdb:
    outputdb.write(new_text)
