#this script takes the clean database produced by rep_clean_database and produces a two different versions of the database: a stemmed version and a lemmatized version, as rep_lemmatized_catalog.txt and rep_stemmed_catalog.txt#


import nltk
import copy
from nltk.tokenize import word_tokenize
from nltk.stem.porter import PorterStemmer
from nltk.stem import WordNetLemmatizer
porter = PorterStemmer()
lemmat = WordNetLemmatizer()

def flatten_everything(list_of_lists):
    ###transforms a list of lists of strings into one tab and new line delimited string for writing, and returns it###
    flat_list = ['\t'.join(article) for article in list_of_lists]
    new_string = '\n'.join(flat_list)
    return new_string

def write_cat(flat_catalog, file_name):
###takes string.txt as file name and a catalog that is just a long string and turns it into a .txt file###
    with open(file_name, 'w') as outputcat:
        outputcat.write(flat_catalog)

def stem_and_lemmat(catalog, stem_catalog, lem_catalog):
###takes original catalog (list of lists) and returns two new catalogs (list of lists): one stemmed and one lemmatized, as objects stem_cat and lem_cat. then writes two new catalogs in .txt###
    n = 0
    abstract_list = [article[6] for article in catalog]
    for abstract in abstract_list:
        new_abstract = word_tokenize(abstract) #tokenize
        new_abstract = nltk.pos_tag(new_abstract) #tag
        new_abstract = [word_tag[0] for word_tag in new_abstract if word_tag[1] not in tag_dictionary and len(word_tag[0]) > 2] #filters
        stem_catalog[n][6] = ' '.join([porter.stem(word) for word in new_abstract])
        lem_catalog[n][6] = ' '.join([lemmat.lemmatize(word) for word in new_abstract])
        n = n + 1
    stem_catalog = flatten_everything(stem_catalog)
    lem_catalog = flatten_everything(lem_catalog)
    write_cat(stem_catalog, 'rep_stemmed_catalog.txt')
    write_cat(lem_catalog, 'rep_lemmatized_catalog.txt')

tag_dictionary = {'CC':0, 'CD':1, 'DT':2, 'EX':3, 'IN':4, 'LS':5, 'MD':6, 'PDT':7, 'POS':8, 'PRP':9, 'PRP$':10, 'RB':11, 'RBS':12, 'RBR':13, 'RP':14, 'SYM':15, 'TO':16, 'UH':17, 'WDT':18, 'WP':19, 'WP$':20, 'WRB':21} #removes everything but adjectives,

whole_text = open('rep_clean_catalog.txt', 'r').read().split('\n')
whole_text = [article.split('\t') for article in whole_text]
whole_text = whole_text[:-1] #deletes the last item of the list, artifact of the catalog generator putting an empty line at the end of the document.
stem_cat = copy.deepcopy(whole_text)
lem_cat = copy.deepcopy(whole_text)

stem_and_lemmat(whole_text, stem_cat, lem_cat)