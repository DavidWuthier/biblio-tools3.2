#! /usr/bin/python -v
# -*- coding: utf-8 -*-

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2017
"""

# usage: parser.py -i DIR [-o DIR] [-e]
# 

import os
import numpy
import math
import time
import argparse
import Utils.Utils as Utils
from fnmatch import fnmatch

## ##################################################
# a dictionnary of common_words used when extracting title words
common_words = ['a', 'able', 'about', 'across', 'after', 'all', 'almost', 'also', 'am', 'among', 'an', 'and', 'any', 'are', 'as', 'at', 'be', 'because', 'been', 'between', 'but', 'by', 'can', 'cannot', 'could', 'dear', 'did', 'do', 'does', 'during','either', 'else', 'ever', 'every', 'for', 'from', 'get', 'got', 'had', 'has', 'have', 'he', 'her', 'hers', 'him', 'his', 'how', 'however', 'i', 'if', 'in', 'into', 'is', 'it', 'its', 'just', 'least', 'let', 'like', 'likely', 'may', 'me', 'might', 'most', 'must', 'my', 'neither', 'no', 'nor', 'not', 'of', 'off', 'often', 'on', 'only', 'or', 'other', 'our', 'own', 'rather', 'said', 'say', 'says', 'she', 'should', 'since', 'so', 'some', 'than', 'that', 'the', 'their', 'them', 'then', 'there', 'these', 'they', 'this', 'tis', 'to', 'too', 'twas', 'us', 'wants', 'was', 'we', 'were', 'what', 'when', 'where', 'which', 'while', 'who', 'whom', 'why', 'will', 'with', 'would', 'yet', 'you', 'your']
french_words = ['un','une','au','aux','entre','a','le','la','les','du','de','des','mais','ou','et','dans','avec','en','sur','sous','avant','apres','vers','par','pendant','depuis','pour','chez','est','ont']
common_words = common_words + french_words;
punctuation = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '.', '/', ':', ';', '<', '=', '>', '?', '@', '[', '\\', ']', '^', '_', '`', '{', '|', '}', '~', ' - ']
# I kept the '-' out of this list
## ##################################################

## ##################################################
## ##################################################
def biblio_parser(in_dir,out_dir,database,expert):

  ## INITIALIZATION
  t1=time.time()

  #.. detecting raw files
  if database == 'wos': pattern = "*.txt"
  if database == 'scopus': pattern = "*.csv"
  print ("..Analysing files %s/%s" % (in_dir,pattern) )
  srclst=[]
  for path, subdirs, files in os.walk(in_dir):
    for name in files:
        if fnmatch(name, pattern):
            srclst.append( os.path.join(path, name))
  print ("....%d '%s' files detected" % (len(srclst),pattern))


  #.. prep empty parsed files
  outfilenames={"A":"articles","AU":"authors","K":"keywords","S":"subjects","R":"references","CU":"countries","AD":"addresses","I":"institutions"}
  if expert:
    #output all info that might be useful in advanced projects 
    outfilenames.update({ "AB":"abstracts" ,"CI":"cities"})
    if database=="wos": outfilenames.update({"FT":"fundingtext", "RP":"AUaddresses"})
  if database=="scopus": outfilenames.update({"S2":"subjects2"})
  outf=dict()
  for elt in outfilenames: outf[elt] = open(os.path.join(out_dir, outfilenames[elt]+".dat"),'w', encoding='utf-8')
  
  #.. keep origin of data
  dst0  = os.path.join(out_dir, "database.dat") 
  with open(dst0,'w') as ff:
    if database == 'wos': ff.write('Web Of Science')
    if database == 'scopus': ff.write('Scopus')

  # SCOPUS JOURNAL CATEGORIES
  #... scopus export data do not contain any subject category. Here we upload an official scopus list from february 2015 listing the SUBJCAT each journal correspond to (multiple category per journal possible)
  if (database =='scopus'):
    print ("..upload Scopus publication sources' categories from auxiliary file")
    journalCATS=dict()
    issnCATS=dict()
    code_cat=dict()
    appearsXtimes=[]
    script_dir = os.path.dirname(__file__)
    #.. read cat_catcode file
    filename = os.path.join(script_dir, 'Utils/scopus_cat_codes.txt')
    with open(filename,'r', encoding='utf8') as fd:
      for line in fd.readlines():
        line = line.strip('\n')
        foo=line.split('\t')
        code_cat[foo[1]]=foo[0]
    #.. read journal_catcode file
    filename = os.path.join(script_dir, 'Utils/scopus_journals_issn_cat.txt')
    with open(filename,'r', encoding='utf8') as fd:
      for line in fd.readlines():
        line = line.strip('\n')
        foo=line.split('\t')
        jrn=foo[0]
        issn=foo[1]
        bar=foo[2].replace(' ','').split(';')
        if (jrn not in journalCATS and len(foo[2])>1): journalCATS[jrn]=[[],[]]
        if (issn not in issnCATS and len(issn)>1 and len(foo[2])>1): issnCATS[issn]=[[],[]]
        else: appearsXtimes.append(jrn)
        for i in range(len(bar)):
          if len(bar[i])>0: 
            cat=code_cat[bar[i]]
            journalCATS[jrn][0].append(cat) 
            if len(issn)>1: issnCATS[issn][0].append(cat)  
            cat=code_cat[bar[i][0:2]+'00'].replace('(all)','').replace('General ','')
            journalCATS[jrn][1].append(cat) 
            if len(issn)>1: issnCATS[issn][1].append(cat) 
    for jrn in journalCATS: journalCATS[jrn]=[list(set(journalCATS[jrn][0])),list(set(journalCATS[jrn][1]))]
    for issn in issnCATS: issnCATS[issn]=[list(set(issnCATS[issn][0])),list(set(issnCATS[issn][1]))]
    del code_cat;     
    # create dict of journals / sources not in that category file        
    journalNOT=dict()

    #print("..%d sources appear more than once in the scopus's source/category file" % len(list(set(appearsXtimes))))



  ## #############################################################
  ## TREAT DATA
  #.. some parameters to count errors / filtered publis
  kompt_publis = 0
  kompt_refs = 0
  kompt_corrupt_refs = 0
  kompt_refs_with_DOI = 0
  #.. "id" will be the new id for articles in the dataset, "unique_ids" tracks the WOS or SCOPUS unique id that we use to remove duplicates
  id = int(-1)
  UNIQUE_IDS = dict()

  # parse and clean metadata about each article
  for src in srclst:
      pl = Utils.ArticleList()
      pl.read_file(src,database)
      print ("..processing %d publications in file %s" % (len(pl.articles), src)) 
      kompt_publis+=len(pl.articles)
      if (len(pl.articles) > 0):
          for article in pl.articles:
            if article.UT not in UNIQUE_IDS:
              UNIQUE_IDS[article.UT] = ''
              id = id + 1
              #article 
              foo = article.AU.split('; ')
              firstAU = foo[0].replace(',','')
              if (article.J9=='' and article.SO==''): article.J9='[]';
              if (article.J9==''): article.J9=article.SO;
              if (article.DT==''): article.DT='[]';
              if (database=='scopus'): article.J9=article.SO;
              outf["A"].write("%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (id,firstAU,article.PY,article.J9,article.VL,article.BP,article.DI,article.PT,article.DT,article.LA,article.TC,article.TI,article.PD,article.UT))
              #journals
              if "J" in outf: outf["J"].write("%d\t%s\t%s\n" % (id, article.SO,article.SN))
              #authors
              if(article.AU != ""): 
                  foo = article.AU.split('; ')
                  for i in range(len(foo)):
                      #... check the upper or lower format of letter to have the author name in format "Grauwin S"
                      foobar=normauthor(foo[i])
                      outf["AU"].write("%d\t%d\t%s\n" % (id,i,foobar))
              #keywords
              if(article.DE != ""):
                  #.. author keywords
                  foo = article.DE.split('; ')
                  for f in foo:
                      outf["K"].write("%d\tAK\t%s\n" % (id,f.upper()))
              if(article.ID != ""):
                  #.. WOS or SCOPUS keywords 
                  foo = article.ID.split('; ')
                  for f in foo:
                      outf["K"].write("%d\tIK\t%s\n" % (id,f.upper()))
              if(article.TI != ""):
                  #.. title words (excluding the common words listed on the top of this file)
                  foo = article.TI
                  #... remove ponctuations !"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~
                  for p in punctuation: foo = foo.replace(p,'')
                  foo = foo.split(' ')
                  for f in foo:
                    bar = f.lower()
                    if bar not in common_words and len(bar)>0:
                      outf["K"].write("%d\tTK\t%s\n" % (id, bar.upper()))
              #subjects
              if database=='wos':
                if(article.WC != ""):
                  foo = article.WC.split('; ')
                  for cat in foo: outf["S"].write("%d\t%s\n" % (id,cat))
              if database=='scopus':
                if article.SO in journalCATS:
                  foo=journalCATS[article.SO][1]
                  for cat in foo: outf["S"].write("%d\t%s\n" % (id,cat))
                  foo=journalCATS[article.SO][0]
                  for cat in foo: outf["S2"].write("%d\t%s\n" % (id,cat))
                elif article.SN in issnCATS:
                  foo=issnCATS[article.SN][1]
                  for cat in foo: outf["S"].write("%d\t%s\n" % (id,cat))
                  foo=issnCATS[article.SN][0]
                  for cat in foo: outf["S2"].write("%d\t%s\n" % (id,cat))
                else:
                  if (article.SO, article.SN) not in journalNOT: journalNOT[(article.SO, article.SN)]=0
                  journalNOT[(article.SO, article.SN)]+=1;
                if ( (article.SO in journalCATS and 'Education' in journalCATS[article.SO][0]) or (article.SN in issnCATS and 'Education' in issnCATS[article.SN][0]) or 'education' in article.SO.lower() or 'learning' in article.SO.lower() ): outf["S"].write("%d\tEducation*\n" % (id))
              # #references
              if(article.CR != ""):
                foo = article.CR.split('; ')
                for i in range(len(foo)):
                  if len(foo[i])>3:
                    ref=Utils.Ref()
                    ref.parse_ref(foo[i],database)
                    kompt_refs += 1 
                    if (ref.year > 0): 
                      if (ref.DOI !=""): kompt_refs_with_DOI+=1
                      outf["R"].write("%d\t%s\t%d\t%s\t%s\t%s\t%s\n" % (id,ref.firstAU,ref.year,ref.journal,ref.volume,ref.page,ref.DOI))
                    if (ref.year == 0):
                      #print (foo[i])
                      kompt_corrupt_refs += 1  
              #countries / cities / institutions
              if(article.C1 != ""):
                  adresse = article.C1
                  aux1 = adresse.find('[')
                  aux2 = adresse.find(']')
                  while (aux1 < aux2) and aux1 > -1:
                      aux = adresse[aux1:min(aux2+2,len(adresse))]
                      adresse = adresse.replace(aux,'')
                      aux1 = adresse.find('[')
                      aux2 = adresse.find(']')
                  foo = adresse.split('; ')
                  for i in range(len(foo)):
                      #... "institutions" will keep everything listed between commas in the adresses, except from cities and countries
                      foo[i] = foo[i].replace(', ',',')
                      bar = foo[i].split(',') 
                      ll = len(bar)
                      for j in range(ll - 2):
                        outf["I"].write("%d\t%d\t%s\n" % (id,i,bar[j]))
                      #... country
                      country = bar[ll-1].replace('.','').replace(';','').upper()
                      lll=len(bar[ll-1])
                      #....  if pb, indicate country X
                      if lll<2: country='X'
                      #.... put all USA states together under the label 'USA'
                      usa_states=['AL','AK','AZ','AR','CA','NC','SC','CO','CT','ND','SD','DE','FL','GA','HI','ID','IL','IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MS','MO','MT','NE','NV','NH','NJ','NM','NY','OH','OK','PA','RI','TN','TX','UT','VT','VA','WV','WA','WI','WY','DC'];
                      usa_states2=[f+' ' for f in usa_states];
                      if (country[lll-3:lll] == 'USA' or country in usa_states or country[0:3] in usa_states2): country = 'USA'
                      #.... put England, Wales, North Ireland, Scotland in UK
                      if country in ['ENGLAND','WALES','NORTH IRELAND','SCOTLAND','UNITED KINGDOM']: country='UK'
                      if country not in ['USA','UK']: 
                        country=" ".join(w.capitalize() for w in country.lower().split())
                      if (database =='scopus' and country == 'USA'):country='United States'
                      if (database =='scopus' and country == 'UK'):country='United Kingdom'
                      outf["CU"].write("%d\t%d\t%s\n" % (id,i,country))
                      #... cities
                      if "CI" in outf:
                        if country=="Australia":city=bar[ll-3]
                        else: city=bar[ll-2]
                        for nb in range(10): city=city.replace(str(nb),'');
                        city=" ".join(w.capitalize() for w in city.split(' ') if w!='F-')
                        city=" ".join(w for w in city.split(' ') if '-' not in w)
                        if len(city) >0:
                          outf["CI"].write("%d\t%d\t%s\n" % (id,i,city))
                      #... addresses
                      add=""
                      for j in range(ll - 1): 
                        add+="%s, " % (bar[j])
                      add +=country
                      outf["AD"].write("%d\t%d\t%s\n" % (id,i,add))

              # abstract (the 2 conditions refer to wos and scopus way to indicate the absence of abstract)
              if "AB" in outf:
                if(article.AB != "") and (article.AB != "[No abstract available]"):
                  outf["AB"].write("%d\t%s\n" % (id,article.AB))

              # funding text
              if "FT" in outf:
                if(article.FX != ""):
                  outf["FT"].write("%d\t%s\n" % (id,article.FX))

              # RP
              if "RP" in outf:
                RPfoo=article.C1.replace('; [',';\t[').split('\t')
                for elt in RPfoo:
                  aux1 = elt.find('[')
                  aux2 = elt.find(']')
                  if (aux1 < aux2 and aux1>-1):
                    auxA = elt[aux1:min(aux2+2,len(elt))]
                    AD=elt.replace(auxA,'').replace(';','')
                    auxA=auxA.replace('[','').replace(']','').split('; ')
                    for a in auxA:
                      outf["RP"].write("%d\t%s\t%s\n"% (id, normauthorB(a), AD))  

      #.. delete the stuff dealing from the raw file just treated from memory        
      del pl

  ## END
  #... how many papers?
  print (".......")
  print ("..%d parsed publications in total (filtered out of %d publications)" % (id + 1, kompt_publis) )
  print ("..%d refs in total, %d with a DOI" % (kompt_refs, kompt_refs_with_DOI))
  #... error in refs?
  if kompt_refs > 0: print ("..%d inadequate refs out of %d (%f%%) have been rejected by this parsing process (no publication year: unpublished, in press, arxiv ...) " % (kompt_corrupt_refs, kompt_refs, (100.0 * kompt_corrupt_refs) / kompt_refs))
  else: print ('..no references found in your dataset. Check whether you extracted your data properly!')
  #... journal with unknown cat?
  if database =='scopus':
    print ("..%d publis from %d distinct sources with unknown subject category" % (sum([journalNOT[jr] for jr in journalNOT]), len(journalNOT)))
    if len(journalNOT)>0:
      foo=len([jr for jr in journalNOT if journalNOT[jr]>=5])
      print ("..%d of these sources with more than 5 papers in the corpus" % (foo))
      filename = 'scopus_sources_without_categories.dat'
      print ("...all these publication sources are listed in %s" % filename)
      with open(filename,'w') as out:
        out.write('Publication source\tISSN\tNpapers\n')
        fb=list(journalNOT.items())
        fb.sort(key=lambda e:-e[1])
        for elt in fb:
          try: 
            out.write("%s\t%s\t%d\n" % (elt[0][0],elt[0][1],elt[1]))
          except:
            out.write("[encoding problem]\t%d\n" % (elt[1]))
  #... close output files
  for elt in outfilenames: outf[elt].close()


  #... end
  t2=time.time()
  print ('..time needed: %ds' % (t2-t1))
  return

## ##################################################
## ##################################################

def normauthor(foo):
    foo = foo.replace(',','')
    aux1 = foo.rfind(' ')
    aux2 = len(foo)
    foobar = foo.lower().capitalize()
    if aux1 > 0: 
        s1 = foobar[aux1:aux2]
        s2 = s1.upper() 
        foobar = foobar.replace(s1,s2)
    aux = foobar.find('-')
    if aux > 0: 
        bar1 = foobar[aux:aux+2]
        bar2 = '-' + foobar[aux+1].upper()
        foobar = foobar.replace(bar1,bar2)
    aux = foobar.find(' ')
    if aux > 0 and (aux!=len(foobar)-1): 
        bar1 = foobar[aux:aux+2]
        bar2 = ' ' + foobar[aux+1].upper()
        foobar = foobar.replace(bar1,bar2)
    return foobar

def normauthorB(foo):
  bar=foo.lower().split(', ')
  foobar=' '.join([e.capitalize() for e in bar[0].split(' ') if len(e)>0])
  foobar='-'.join([f[0].upper()+f[1:] for f in foobar.split('-') if len(f)>0])
  if len(bar)>1: foobar += ' ' + ''.join([e[0].upper() for e in bar[1].replace('-',' ').split(' ') if len(e)>0])

  return foobar

## ##################################################
## ##################################################
## ##################################################

def main():
# usage: parser.py [-h] [--version] -i DIR -o DIR [-v]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   --version             show program's version number and exit
#   -o DIR, --output_dir DIR
#                         output directory name
#   -i DIR, --input_dir DIR
#                         input directory name
  # Parse line options.
  # Try to have always the same input options
  parser = argparse.ArgumentParser(description = 'parser')

  parser.add_argument('--version', action='version', version='%(prog)s 1.1')
  
  parser.add_argument("-i", "--input_dir", nargs=1, required=True,
          action = "store", dest="in_dir",
          help="input directory name",
          metavar='DIR')
          
  parser.add_argument("-o", "--output_dir", nargs=1, required=False,
          action = "store", dest="out_dir",
          help="output directory name",
          default = "Desktop/",
          metavar='DIR')

  parser.add_argument("-d", "--database",
          action = "store", dest="database",
          default = 'wos',
          help="database [default %(default)s]",
          metavar='string')  

  parser.add_argument("-e", "--expert",
          action = "store_true", dest="expert",
          default = False,
          help="expert mode [default %(default)s]")

  #Analysis of input parameters
  args = parser.parse_args()
  
  if (not os.path.exists(args.in_dir[0])):
    print ("Error: Input directory does not exist: ", args.in_dir[0])
    exit()


  if args.out_dir == 'Desktop/':
    args.out_dir = args.in_dir

  if (not os.path.exists(args.out_dir[0])):
    args.out_dir[0]=args.in_dir[0]
    print ("Error: Output directory does not exist: ", args.out_dir[0])
    exit()

  if args.database not in ['wos','scopus']:
    print ("Error: database must be either 'wos' or 'scopus'")
    exit()

  ##      

  biblio_parser(args.in_dir[0],args.out_dir[0],args.database,args.expert)

  return


    
## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################

