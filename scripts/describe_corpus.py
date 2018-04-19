#! /usr/bin/env python

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2017
"""

# usage: describe_corpus.py -i DIR [-v]
# 

import os
import time
import math
import numpy
import argparse
import itertools
from collections import Counter


## ##################################################
## ##################################################
## ##################################################

def describe_corpus(in_dir, out_dir, verbose):

  ## INITIALIZATION
  t1=time.time()
  if verbose: print ("INITIALIZATION")
  # empty dicts for the "top 100 coocurrence network that will be shown in BiblioMaps"
  myNODES=[]
  myLINKS=[]
  with open(os.path.join(in_dir, "articles.dat") ,"r", encoding='utf-8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1]  
  # total nb_art  
  nb_art=len(data_lines)
  print ("... %d publis" % nb_art)
  # keep track of the database
  with open(os.path.join(in_dir, "database.dat"),'r') as ff: thedatabase=ff.read().strip('\n')
  # file that will record all proba distribs  
  outfile_itemuse=os.path.join(out_dir,'DISTRIBS_itemuse.json')
  with open(outfile_itemuse,'w') as out: out.write('{\n\t"database":"%s",\n\t"N":%d' % (thedatabase, nb_art))


  ## ##################################################
  ## AUX FUNCTIONS
  def do_analysis(stuff,filename,nick):
    # HOW MUCH ITEMS BY PUBLI?
    stuff.sort(key=lambda e:e[0])
    bah=[(my_id, len(set(list(item)))) for my_id, item in itertools.groupby(stuff,key=lambda e:e[0])]
    bah.sort(key=lambda x:x[1])
    itemdistrib=[(nb,len(list(v))) for nb,v in itertools.groupby(bah,key=lambda e:e[1])]
    with open(outfile_itemuse,'a') as out: out.write(',\n\t"q%s":[%s, %s]'% ( nick,[elt[0] for elt in itemdistrib],[elt[1] for elt in itemdistrib]))

    # HOW MUCH TIME AN ITEM APPEARS?
    stuff.sort(key=lambda e:e[1])
    # duplicates are removed: we count in how many distinct publis an item appears 
    foo=[(xx,len(set(list(art_id)))) for xx,art_id in itertools.groupby(stuff,key=lambda e:e[1])]
    foo.sort(key=lambda x:-x[1])
    # generate co-occurrence network
    keep=[foo[i] for i in range(min(100,len(foo)))]
    generate_cooc(stuff,keep,nick, in_dir);
    #
    nb_with=len(list(set([e[0] for e in stuff])))
    if nb_with < nb_art:
      foo+=[('none available', nb_art-nb_with)]
      foo.sort(key=lambda x:-x[1])
    #
    outfile = os.path.join(out_dir, filename)
    # output the top used items
    mythr=len(foo);kompt=1;
    while(mythr>10000):
      mythr=len([elt for elt in foo if elt[1]>kompt])
      kompt+=1

    with open(outfile,'w',encoding='utf-8-sig') as out:
      out.write("item,count,f\n")
      for i in range(mythr):
        (x,v)=foo[i]
        out.write('"%s",%d,%.2f\n' % (x,v,100.0*v/nb_art))
    # output the total distrib
    distrib=[(nb,len(list(v))) for nb,v in itertools.groupby(foo,key=lambda e:e[1])]
    distrib.sort()
    xx=[elt[0] for elt in distrib]
    yy=[elt[1] for elt in distrib]
    yy=list(sum(yy)-numpy.cumsum([0]+yy))
    yy=yy[0:-1]
    with open(outfile_itemuse,'a') as out: out.write(',\n\t"p%s":[%s, %s]'% (nick,xx,yy))

    # free memory
    del stuff, xx, yy, distrib, keep, foo

    return

  ## ################

  def generate_cooc(stuff,keep,nick,in_dir):
    if nick in ["AU","CU","S","S2","K","R","RJ","I","AK","TK"]:
      nodeID={};aux={};comm={};
      # nodes
      for i in range(len(keep)):
        nodeID[keep[i][0]]=i
        myNODES.append('{"type":"%s","name":%d,"item":"%s","size":%d}' % (nick,i,keep[i][0],keep[i][1]))
      # links
      for elt in stuff:
        if elt[1] in nodeID: 
          if elt[0] not in aux: aux[elt[0]]=[]
          aux[elt[0]].append(elt[1])
      for ee in aux:
        for itema in aux[ee]:
          for itemb in [ff for ff in aux[ee] if ff > itema]:
             if (itema, itemb) not in comm: comm[(itema, itemb)]=0
             comm[itema,itemb]+=1
      for foo in comm:
        myLINKS.append('{"type":"%s","source":%d,"target":%d,"Ncooc":%d}' % (nick,nodeID[foo[0]],nodeID[foo[1]],comm[foo]))

      # free memory
      del aux, nodeID, comm
    return

  ## ################

  def treat_item(item, nick, which):
    if verbose: print ("... dealing with %s" % (item) )
    with open(os.path.join(in_dir, item+".dat") ,"r",encoding='utf-8') as file:
      # dat file have one trailing blank line at end of file
      data_lines=file.read().split("\n")[:-1] 
    aux = [(l.split("\t")[0], l.split("\t")[which]) for l in data_lines]
    do_analysis(aux,'freq_'+item+'.dat',nick)  
    # free memory
    del data_lines 
    # end
    return

  ## ##################################################
  ## TREAT DATA 
  if verbose: print ("... dealing with years + journals + doctypes + languages" )
  # data_lines was already filtered to keep only publi in the correct time-window
  #.. Y
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines]
  do_analysis(aux,'freq_years.dat','Y')
  #.. J
  aux = [(l.split("\t")[0], l.split("\t")[3]) for l in data_lines]
  do_analysis(aux,'freq_journals.dat','J')
  #.. D
  aux = [(l.split("\t")[0], l.split("\t")[8]) for l in data_lines]
  do_analysis(aux,'freq_doctypes.dat','DT') 
  #.. L
  aux = [(l.split("\t")[0], l.split("\t")[9]) for l in data_lines]
  do_analysis(aux,'freq_languages.dat','LA')    
  #.. free memory
  del data_lines
  
  treat_item('authors', 'AU', 2)
  treat_item('subjects', 'S', 1)
  treat_item('institutions', 'I', 2)
  treat_item('countries', 'CU', 2)
  #treat_item('cities', 'CI', 2) 
  with open(os.path.join(in_dir, "database.dat"),'r') as ff:
    if ff.read()=='Scopus': treat_item('subjects2', 'S2', 1)

  if verbose: print ("... dealing with keywords + title words") 
  with open(os.path.join(in_dir, "keywords.dat") ,"r",encoding='utf-8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1]
  #.. IK
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines if l.split("\t")[1]=='IK' ]
  do_analysis(aux,'freq_keywords.dat','K')
  #.. AK
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines if l.split("\t")[1]=='AK' ]
  do_analysis(aux,'freq_authorskeywords.dat','AK')  
  #.. TK
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines if l.split("\t")[1]=='TK' ]
  do_analysis(aux,'freq_titlewords.dat','TK')  
  #.. free memory
  del data_lines 


  if verbose: print ("... dealing with references + references journals") 
  with open(os.path.join(in_dir, "references.dat") ,"r",encoding='utf-8', errors='ignore') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1]
  #.. R
  aux = [(l.split("\t")[0],", ".join(l.split("\t")[1:5]).replace(',0,0','').replace(', 0, 0','').replace(', 0','')) for l in data_lines]
  do_analysis(aux,'freq_references.dat','R')  
  #.. RJ
  aux = [(l.split("\t")[0], l.split("\t")[3]) for l in data_lines]
  do_analysis(aux,'freq_refjournals.dat','RJ')  
  #.. free memory
  del data_lines

  ## ##################################################
  ## END  
  with open(os.path.join(out_dir,'coocnetworks.json'), 'w') as out:
    out.write('{\n"nodes":[\n\t%s\n],\n"links":[\n\t%s]\n}' % (',\n\t'.join([elt for elt in myNODES]).replace('\\',''), ',\n\t'.join([elt for elt in myLINKS]) ))
  with open(outfile_itemuse,'a') as out: out.write('\n}')
  if (verbose): print ('Time needed: %ds' % (time.time()-t1))
  return

## ##################################################
## ##################################################
## ##################################################

def main():
# usage: describe_corpus.py [-h] [--version] -i DIR [-v] [-H INT] [-O]
#
# optional arguments:
#   -h, --help            show this help message and exit
#   --version             show program's version number and exit
#   -i DIR, --input_dir DIR input directory name 
  # Parse line options.
  # Try to have always the same input options
  parser = argparse.ArgumentParser(description = 'parser')

  parser.add_argument('--version', action='version', version='%(prog)s 1.1')
  
  parser.add_argument("-i", "--input_dir", nargs=1, required=True,
          action = "store", dest="in_dir",
          help="input directory name",
          metavar='DIR')
          
  parser.add_argument("-v", "--verbose",
          action = "store_true", dest="verbose",
          default = False,
          help="verbose mode [default %(default)s]")

  # Analysis of input parameters
  args = parser.parse_args()
  if (not os.path.exists(args.in_dir[0])):
      print ("Error, Input directory does not exist: ", args.in_dir[0])
      exit()

  # Prep out folder
  out_dir=os.path.join(args.in_dir[0], 'freqs')
  if not os.path.exists(out_dir): os.makedirs(out_dir)

  # Proceed with main function
  describe_corpus(args.in_dir[0],out_dir,args.verbose)

  return

## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################
