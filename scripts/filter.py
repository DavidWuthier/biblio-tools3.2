#! /usr/bin/env python

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2017
"""

# usage: filter.py -i DIR -o DIR [-v]
# generic filter function to adapt according to your needs

import os
import time
import argparse

## ##################################################
## ##################################################
## ##################################################
def filter_corpus(in_dir,out_dir,verbose):

  ## INITIALIZATION
  t1=time.time()
  # files to filter (will not be taken into account if they don't exist)
  filenames=["articles", "addresses", "authors", "countries", "institutions", "journals", "keywords", "references", "subjects", "subjects2", "AUaddresses", "cities", "fundingtext", "abstracts"] 

  ## ON WHICH ITEMS DO YOU WANT TO FILTER? 
  # 'AU':'author','CU':'country/territory','I':'institution','DT':'document type','LA':'language','Y':'publication year','S':'subject category','S2':'subject subcategory' (scopus only!!),'J':'publication source','K':'keyword','AK':'authors\' keyword','TK':'title word','R':'reference','RJ':'reference source'
  FILTER_ON=["Y","J"]

  # combine all filters by "union" (one of the conditions needs to be checked)  or "intersection" (all conditions needs to be respected)
  COMBINE="intersection"

  ## FILTER ACCORDING TO WHAT? - Beware, these are case sensitive!
  ymin=2010;ymax=2020;                                 
  keywords=["network", "bibliographic"]
  subjects=["Physics", "Computer Sciences"]
  subsubjects=["Education"]
  pubsources=["NATURE", "SCIENCE", "PHYS REV E"]
  authors=["Grauwin S"]
  countries=["France", "Germany", "UK", "Italy", "Spain", "Belgium", "Netherlands"]
  institutions=["CNRS"]
  references=["Grauwin S, 2009, P NATL ACAD SCI USA, 106"]
  refsources=["P NATL ACAD SCI USA"]
  languages=["English","French"]
  doctypes=['Article', 'Review', 'Proceedings Paper', 'Article', 'Article in Press', 'Review', 'Conference Paper', 'Conference Review']



  ## ##################################################
  ## SELECT ID OF PUBLICATIONS TO KEEP
  KEEPID=dict()

  with open(os.path.join(in_dir, "articles.dat") ,"r") as file: data_lines=file.read().split("\n")[:-1]
  nbpub=len(data_lines)  
  if ("Y" in FILTER_ON): 
    print("...selecting publications with selected publication years")
    KEEPID["Y"]=dict([(l.split("\t")[0],'') for l in data_lines if (int(l.split("\t")[2])<=ymax and int(l.split("\t")[2])>=ymin)])
  if ("J" in FILTER_ON):
    print("...selecting publications with selected sources")
    KEEPID["J"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[3] in pubsources)])
  if ("DT" in FILTER_ON):
    print("...selecting publications with selected doctypes")
    KEEPID["DT"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[8] in doctypes)])
  if ("LA" in FILTER_ON):
    print("...selecting publications with selected languages")
    KEEPID["LA"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[9] in languages)])
  del data_lines

  if ("S" in FILTER_ON):
    print("...selecting publications with selected subjects")
    with open(os.path.join(in_dir, "subjects.dat") ,"r", encoding='utf8') as file: data_lines=file.read().split("\n")[:-1]  
    KEEPID["S"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[1] in subjects)])
    del data_lines

  if ("S2" in FILTER_ON):
    print("...selecting publications with selected subsubjects")
    with open(os.path.join(in_dir, "subjects2.dat") ,"r", encoding='utf8') as file: data_lines=file.read().split("\n")[:-1]  
    KEEPID["S2"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[1] in subsubjects)])
    del data_lines

  if ("AU" in FILTER_ON):
    print("...selecting publications with selected authors")
    with open(os.path.join(in_dir, "authors.dat") ,"r", encoding='utf8') as file: data_lines=file.read().split("\n")[:-1]  
    KEEPID["AU"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[2] in authors)])
    del data_lines

  if ("CU" in FILTER_ON):
    print("...selecting publications with countries")
    with open(os.path.join(in_dir, "countries.dat") ,"r", encoding='utf8') as file: data_lines=file.read().split("\n")[:-1]  
    KEEPID["CU"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[2] in countries)])
    del data_lines

  if ("I" in FILTER_ON):
    print("...selecting publications with selected institutions")
    with open(os.path.join(in_dir, "institutions.dat") ,"r", encoding='utf8') as file: data_lines=file.read().split("\n")[:-1]  
    KEEPID["I"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[2] in institutions)])
    del data_lines

  if ("K" in FILTER_ON or "TK" in FILTER_ON or "AK" in FILTER_ON):
    with open(os.path.join(in_dir, "keywords.dat") ,"r", encoding='utf8') as file: data_lines=file.read().split("\n")[:-1]  
    if ("K" in FILTER_ON):
      print("...selecting publications with selected keywords")
      KEEPID["K"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[1]=="IK" and l.split("\t")[2] in keywords)])
    if ("TK" in FILTER_ON):
      print("...selecting publications with selected title words")
      KEEPID["TK"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[1]=="TK" and l.split("\t")[2] in keywords)])
    if ("AK" in FILTER_ON):
      print("...selecting publications with selected authors' jeywords")
      KEEPID["AK"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[1]=="AK" and l.split("\t")[2] in keywords)])
    del data_lines

  if ("R" in FILTER_ON or "RJ" in FILTER_ON):
    with open(os.path.join(in_dir, "references.dat") ,"r", encoding='utf8') as file: data_lines=file.read().split("\n")[:-1]  
    if ("R" in FILTER_ON):
      print("...selecting publications with selected references")
      KEEPID["R"]=dict([(l.split("\t")[0],'') for l in data_lines if ( ', '.join([l.split("\t")[k] for k in range(1,6)]).replace(',0,0','').replace(', 0, 0','').replace(', 0','') in references)])
    if ("RJ" in FILTER_ON):
      print("...selecting publications with selected ref sources")
      KEEPID["RJ"]=dict([(l.split("\t")[0],'') for l in data_lines if (l.split("\t")[3] in refsources)])
    del data_lines

  ## combine the filtering conditions by union / intersection
  TOKEEP={}
  if COMBINE=="intersection":
    for pubid in KEEPID[FILTER_ON[0]]:
      cond=True
      for k in range(1,len(FILTER_ON)):
        if pubid not in KEEPID[FILTER_ON[k]]: cond=False
      if cond: TOKEEP[pubid]=''
  if COMBINE=="union":
    for x in FILTER_ON:
      for pubid in KEEPID[x]:TOKEEP[pubid]=''

  ## how much?
  print (".. %d publications selected out of %d (%.2f%%)" % (len(TOKEEP), nbpub, len(TOKEEP)*100.0/nbpub) )

  ## ##################################################
  ## PROCEED
  for fff in filenames:
    if (os.path.exists(os.path.join(in_dir, fff+".dat"))):
      print ("filtering %s" % fff)
      with open(os.path.join(in_dir, fff+".dat") ,"r", encoding='utf8') as file: 
        data_lines=file.read().split("\n")[:-1] 
      with open(os.path.join(out_dir, fff+".dat") ,"w", encoding='utf8') as outfile: 
        for l in data_lines:
          if l.split("\t")[0] in TOKEEP: outfile.write("%s\n" % l)
        del data_lines
  #
  with open(os.path.join(in_dir, "database.dat") ,"r", encoding='utf8') as file: 
      foo=file.read()
  with open(os.path.join(out_dir, "database.dat") ,"w", encoding='utf8') as outfile: 
      outfile.write("%s\n" % foo)
  del foo

  ## ##################################################
  ## END  
  print ('END - total time needed: %ds' % (time.time()-t1))
  return

## ##################################################
## ##################################################
## ##################################################



def main():
# usage: describe_corpus.py [-h] [--version] -i DIR [-v]
#
# optional arguments:
#   -h, --help            show this help message and exit
#   --version             show program's version number and exit
#   -i DIR, --input_dir DIR input directory name
#   -r 
  # Parse line options.
  # Try to have always the same input options
  parser = argparse.ArgumentParser(description = 'parser')

  parser.add_argument('--version', action='version', version='%(prog)s 1.1')
  
  parser.add_argument("-i", "--input_dir", nargs=1, required=True,
          action = "store", dest="in_dir",
          help="input directory name",
          metavar='DIR')

  parser.add_argument("-o", "--output_dir", nargs=1,
          action = "store", dest="out_dir",
          help="output directory name",
          metavar='DIR')

  parser.add_argument("-v", "--verbose",
          action = "store_true", dest="verbose",
          default = False,
          help="verbose mode [default %(default)s]")

  #Analysis of input parameters
  args = parser.parse_args()
  
  if (not os.path.exists(args.in_dir[0])):
      print ("Error: Input directory does not exist: ", args.in_dir[0])
      exit()

  if not os.path.exists(args.out_dir[0]): 
    os.makedirs(args.out_dir[0])
    print ("Output directory %s, which did not previously exist, was created" % args.out_dir[0])


  ##      

  filter_corpus(args.in_dir[0],args.out_dir[0],args.verbose)

  return

## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################
