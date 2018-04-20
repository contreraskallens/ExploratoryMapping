# this script creates a document by (unique) term of the catalog specified as database (default:
# 'lemmatized_catalog.txt'). Writes it down as 'document_by_term.txt'#

import re


def flatten(mylist):  # flatten one level deep list
    new_list = []
    for item in mylist:
        n = 0
        while n < len(item):
            new_list.append(item[n])
            n = n + 1
    return new_list


database = open('lemmatized_catalog.txt', 'r').read().split('\n')  # open lemmatized catalog

abstract_dictionary = {int(lista[0]): lista[6] for lista in
                       [line.split('\t') for line in database]}  # makes a dictionary with the data of the catalog

term_list = [term for term in (sorted(set(flatten([abstract.split(' ') for abstract in abstract_dictionary.values()]))))
             if len(term) > 3]  # makes a list of the unique terms in the abstracts longer than 3 characters

word_dictionary = {key: abstract.split(' ') for key, abstract in
                   abstract_dictionary.items()}  # makes a dictionary of the words that each abstract has.

dbt = {}  # generates the table of values with documents as rows and terms as columns.

for i in range(len(abstract_dictionary.values())):
    value_list = []  # the values will go here
    value_list.extend([0] * len(term_list))
    dbt[i + 1] = value_list  # pads list with 0s
    words = word_dictionary[i + 1]  # words of the abstract
    for word in words:
        if len(word) > 3:
            term_index = term_list.index(word)
            dbt[i + 1][term_index] += 1

with open('document_by_term_lemmatized.txt', 'w') as output_file:
    for term in term_list:
        output_file.writelines('\t' + term)
    output_file.writelines('\n')
    for key, matrix in dbt.items():
        output_file.writelines(str(key))
        for value in matrix:
            output_file.writelines('\t' + str(value))
        output_file.writelines('\n')