#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2012
   All rights reserved.
   BSD license.
   .... If you are using these scripts, please cite our "Scientometrics" paper:
   .... S Grauwin, P Jensen, Mapping Scientific Institutions. Scientometrics 89(3), 943-954 (2011)
"""

import os
import sys
import glob
import numpy
import argparse
import itertools
import math
import json
import Utils.Utils as Utils

## ##################################################
## ########################### Comm tables
## ##################################################

def comm_tables(in_dir,parti,thr,verbose):
  ## this function is written in order to account for overlapping or not overlapping clusters
  partition=parti.copy()
  # transform partition values into lists if they are not (thanks to this, the following will work whether articles belong to one or several clusters)
  if type(list(partition.values())[0]) is not list:
    for elt in partition: partition[elt]=[partition[elt]]

  # list clusters and their sizes: below we are only interested by articles within communities of size > thr
  comm_size = dict();
  list_nodes = dict();
  list_comm = [] ## list of communities of size > thr
  allcom=[]; [allcom.extend(y) for y in partition.values()];
  for com in set(allcom):
    list_nodes[com] = dict()
    aa=list(set([node for node in partition.keys() if com in partition[node]]))
    for a in aa: list_nodes[com][a]='';
    comm_size[com] = len(list_nodes[com])
    if comm_size[com] > thr: list_comm.append( com )

  ## 
  list_id = dict(); ## list of id -- articles within the BC network and within a community of size > thr
  for id_art in partition:
    for com in partition[id_art]:
      if comm_size[com] > thr:
        if id_art not in list_id: list_id[id_art] = []
        list_id[id_art].append(com)

  ## AUX FUNCTIONS
  def extract_refstuff(mylist, mysize):
    mylist.sort()
    # h index
    h=1
    if (len(mylist)==0): 
      h=0
    else: 
      while mylist[-h] > h: h+=1
    # how many ref are used 10% / 5% / 2% times?
    x10=len([a for a in mylist if a > 0.1*mysize]);
    x5=len([a for a in mylist if a > 0.05*mysize]);
    x2=len([a for a in mylist if a > 0.02*mysize]);

    return [h,x10,x5,x2]

  def extractMostFreq(aux,howmuch,quant):
    #.. keep only stuff about papers in BC network
    #aux.sort(key=lambda e:e[0])
    #aux=[e for e in aux ]; #if int(e[0]) in list_id];
    #.. count the % of publis with info existing/available
    fooaux=dict((elt[0],'') for elt in aux)  
    NOitemin={}
    for pub in list_id:
      if str(pub) not in fooaux: NOitemin[str(pub)]=''
    #.. total occurrences of a given item
    aux.sort(key=lambda e:e[1])
    freqTOT=[(xx,len(set(list(art_id)))) for xx,art_id in itertools.groupby(aux,key=lambda e:e[1])]
    freqTOT=dict(freqTOT)
    #.. number of papers with items
    NN=len(set([elt[0] for elt in aux]))
    #.. occurrences of a given item in each cluster
    stuffX=dict(); available_for=dict();
    if quant: quantR=dict()
    foo=dict();
    for com in list_comm: foo[com]=[];
    aux=[e for e in aux if int(e[0]) in list_id];
    for e in aux:
      for com in list_id[int(e[0])]: foo[com].append(e)
    for com in list_comm:
      foo[com].sort(key=lambda e:e[1])
      L=[(xx,len(set(list(art_id)))) for xx,art_id in itertools.groupby(foo[com],key=lambda e:e[1])]
      L.sort(key=lambda e:-e[1])
      # ...
      stuffX[com] = dict();
      cs = comm_size[com];
      # count the % of pub with available items
      available_for[com]=len(list(set([e[0] for e in foo[com]])))*100.0/cs;
      #
      if quant:
        LL = [elt[1] for elt in L]; 
        quantR[com]=extract_refstuff(LL, cs); 
      for i in range(min(howmuch,len(L))):
        item = L[i][0];
        f = L[i][1] * 1.0 / cs; 
        p = freqTOT[item] * 1.0 / NN; 
        if (p < 1 and L[i][0]!='none available'): sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffX[com][i] = [item.replace('&','\\&'), f*100, sigma]
    if quant:
      return stuffX, available_for, quantR  
    else: 
      return stuffX, available_for

  ## TREAT DATA - EXTRACT ITEMS
  BIGSTUFF={}; AVAIL={};

  # KEYWORDS & TITLE WORDS
  if verbose: print ("....most frequent keywords / title words")
  with open(os.path.join(in_dir, "keywords.dat"), "r", encoding='utf8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1] 
  auxa = [(l.split("\t")[0], l.split("\t")[1], l.split("\t")[2]) for l in data_lines]
  aux = [(a[0], a[2]) for a in auxa if a[1]=='IK']
  (BIGSTUFF['K'], AVAIL['K']) = extractMostFreq(aux,20,0)
  aux = [(a[0], a[2]) for a in auxa if a[1]=='TK']
  (BIGSTUFF['TK'], AVAIL['TK']) = extractMostFreq(aux,20,0)
  #.. free memory
  del data_lines, aux

  # SUBJECTS
  if verbose: print ("....most frequent subjects")
  with open(os.path.join(in_dir, "subjects.dat"), "r", encoding='utf8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1] 
  aux = [(l.split("\t")[0], l.split("\t")[1]) for l in data_lines]
  (BIGSTUFF['S'], AVAIL['S']) = extractMostFreq(aux,20,0)
  #.. free memory
  del data_lines 

  with open(os.path.join(in_dir, "database.dat"),'r') as ff:
    if ff.read()=='Scopus': 
      with open(os.path.join(in_dir, "subjects2.dat"), "r", encoding='utf8') as file:
        # dat file have one trailing blank line at end of file
        data_lines=file.read().split("\n")[:-1] 
      aux = [(l.split("\t")[0], l.split("\t")[1]) for l in data_lines]
      (BIGSTUFF['S2'], AVAIL['S2']) = extractMostFreq(aux,20,0)
      #.. free memory
      del data_lines 

  # JOURNALS / YEARS
  if verbose: print ("....most frequent journals")
  with open(os.path.join(in_dir, "articles.dat"), "r", encoding='utf8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1] 
  aux = [(l.split("\t")[0], l.split("\t")[3]) for l in data_lines]
  (BIGSTUFF['J'], AVAIL['J']) = extractMostFreq(aux,20,0)
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines]
  (BIGSTUFF['Y'], AVAIL['Y']) = extractMostFreq(aux,20,0)  
  #.. free memory
  del data_lines   

  # AUTHORS
  if verbose: print ("....most frequent authors")
  with open(os.path.join(in_dir, "authors.dat"), "r", encoding='utf8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1] 
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines]
  (BIGSTUFF['A'], AVAIL['A']) = extractMostFreq(aux,20,0)
  #.. free memory
  del data_lines  

  # INSTITUTIONS
  optionI='D'
  if optionI=='D':
    if verbose: print ("....most frequent institutions")
    with open(os.path.join(in_dir, "institutions.dat"), "r", encoding='utf8') as file:
      # dat file have one trailing blank line at end of file
      data_lines=file.read().split("\n")[:-1] 
    aux = [(l.split("\t")[0], l.split("\t")[2].upper()) for l in data_lines]
    (BIGSTUFF['I'], AVAIL['I']) = extractMostFreq(aux,20,0)
    #.. free memory
    del data_lines    
  if optionI=='Dclean':
    if verbose: print ("....most frequent institutions")
    with open(os.path.join(in_dir, "institutionsC.dat"), "r", encoding='utf8') as file:
      # dat file have one trailing blank line at end of file
      data_lines=file.read().split("\n")[:-1] 
    aux = [(l.split("\t")[0], l.split("\t")[2].upper()) for l in data_lines]
    (BIGSTUFF['I'], AVAIL['I']) = extractMostFreq(aux,20,0)
    #.. free memory
    del data_lines  
  if optionI=='fI':
    if verbose: print ("....most frequent institutions")
    with open(os.path.join(in_dir, "addresses.dat") ,"r") as file:
      # dat file have one trailing blank line at end of file
      data_lines=file.read().split("\n")[:-1] 
    aux = [(l.split("\t")[0], l.split("\t")[2].split(",")[0].upper()) for l in data_lines]
    (BIGSTUFF['I'], AVAIL['I']) = extractMostFreq(aux,20,0)
    #.. free memory
    del data_lines  
  

  # COUNTRIES
  if verbose: print ("....most frequent countries")
  with open(os.path.join(in_dir, "countries.dat"), "r", encoding='utf8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1] 
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines]
  (BIGSTUFF['C'], AVAIL['C']) = extractMostFreq(aux,20,0)
  #.. free memory
  del data_lines  


  def myref(spl):
    if spl[4]=='0':
      ref=", ".join(spl[1:4])
    else: 
      ref=", ".join(spl[1:6])
    return ref

  # REFS / REF JOURNALS
  if verbose: print ("....most frequent refs / ref journals")
  with open(os.path.join(in_dir, "references.dat"), "r", encoding='utf8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1] 
  aux = [(l.split("\t")[0],myref(l.split("\t"))) for l in data_lines]
  (BIGSTUFF['R'], AVAIL['R'], quantR) = extractMostFreq(aux,20,1)
  #
  aux = [(l.split("\t")[0], l.split("\t")[3]) for l in data_lines]
  (BIGSTUFF['RJ'], AVAIL['RJ']) = extractMostFreq(aux,20,0)  
  #.. free memory
  del data_lines  

  ## END
  return (quantR, BIGSTUFF, AVAIL)


## ##################################################
## ########################### Comm groups (compute % of subjects / countries in clusters)
## ##################################################
def comm_groups(in_dir,parti,thr,verbose):
  ## 
  ## all the codes below are written to account for overlapping clusters!
  ##
  partition=parti.copy()
  # transform partition values into lists if they are not (this function need to work in case when article belong to one or several clusters)
  if type(list(partition.values())[0]) is not list:
    for elt in partition: partition[elt]=[partition[elt]]

  # Communities sizes - we are mostly interested by articles within communities of size > thr
  comm_size = dict();
  list_nodes = dict();
  list_comm = [] ## list of communities of size > thr
  allcom=[]; [allcom.extend(y) for y in partition.values()];
  for com in set(allcom):
    list_nodes[com] = dict()
    aa=list(set([node for node in partition.keys() if com in partition[node]]))
    for a in aa: list_nodes[com][a]='';
    comm_size[com] = len(list_nodes[com])
    if comm_size[com] > thr: list_comm.append( com )

  ## 
  list_id = dict(); ## list of id -- articles within the BC network and within a community of size > thr
  for id_art in partition:
    for com in partition[id_art]:
      if comm_size[com] > thr:
        if id_art not in list_id: list_id[id_art] = []
        list_id[id_art].append(com)

  def extractS(aux):
    #.. keep only stuff about papers in BC network
    aux.sort(key=lambda e:e[0])
    aux=[e for e in aux if int(e[0]) in list_id];
    #.. total occurrences of a given item
    aux.sort(key=lambda e:e[1])
    freqTOT=[(xx,len(set(list(art_id)))) for xx,art_id in itertools.groupby(aux,key=lambda e:e[1])]
    #.. sort list of keys in alphabetical order 
    freqTOT=dict(freqTOT)
    listX=list(freqTOT.keys())
    listX.sort()
    #.. number of papers with items
    NN=len(set([elt[0] for elt in aux]))
    #.. occurrences of a given item in each cluster
    stuffX=dict();
    foo=dict();
    for com in list_comm: foo[com]=[];
    aux=[e for e in aux if int(e[0]) in list_id];
    for e in aux:
      for com in list_id[int(e[0])]: foo[com].append(e)
    for com in list_comm:
      foo[com].sort(key=lambda e:e[1])
      L=[(xx,len(set(list(art_id)))) for xx,art_id in itertools.groupby(foo[com],key=lambda e:e[1])]
      L=dict(L);
      # ...
      stuffX[com] = [];
      for kk in listX:
        if kk in L: stuffX[com].append(L[kk])
        else: stuffX[com].append(0) 

    allX=[]
    for x in listX: allX.append([x,1.0*freqTOT[x]/NN ]);    
    return (allX, stuffX)

  ## TREAT DATA - EXTRACT ITEMS

  # SUBJECTS
  if verbose: print ("....subjects")
  with open(os.path.join(in_dir, "subjects.dat") ,"r") as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1] 
  aux = [(l.split("\t")[0], l.split("\t")[1]) for l in data_lines]

  """
  # for special analysis / scopus
  with open(os.path.join(in_dir, "subsubjects.dat") ,"r") as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1] 
  aux = aux + [(l.split("\t")[0], l.split("\t")[1]) for l in data_lines]
  """

  (listS, stuffS) = extractS(aux)
  #.. free memory
  del data_lines 

  ## END
  return (listS, stuffS)

## ##################################################
## ########################### my_write (to report tables in .tex document)
## ##################################################
def my_write(f_out,txttype,my_stuff,com,nbline):
    if com in my_stuff:
      f_out.write("{\\bf %s }& f(\%%) & $\sigma$\\\\\n\hline\n" % ( txttype ))
      for i in range(min(nbline,len(my_stuff[com]))):
        f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( my_stuff[com][i][0].replace('[','$[$').replace(']','$]$'), my_stuff[com][i][1], my_stuff[com][i][2]) )         
      for i in range(max(0,nbline-len(my_stuff[com]))):
        f_out.write(" &  & \\\\\n")
      f_out.write("\hline\n")    
    else:
      if txttype != "Subject":
        f_out.write("%s & f(\%%) & $\sigma$\\\\\n\hline\n" % ( txttype ))
        for i in range(nbline): f_out.write(" &  & \\\\\n")


## ##################################################
## ########################### most representative / cited articles / authors
## ##################################################
def comm_AR(in_dir,parti,thr,dgcl_id,verbose):
  ## INI
  LL=20 # length of output lists
  CR_papers=dict()
  CR_authors=dict()
  limtitle = 120
  
  partition=parti.copy()
  # transform partition values into lists if they are not (this function need to work in case when article belong to one or several clusters)
  if type(list(partition.values())[0]) is not list:
    for elt in partition: partition[elt]=[partition[elt]]

  if dgcl_id=='NaN':
    dgcl_id=dict()
    for elt in partition: dgcl_id[elt]=1;

  #.. cluster sizes
  cluster_size = dict();
  allcom=[]; [allcom.extend(y) for y in partition.values()];
  for com in set(allcom):
    list_nodes = [node for node in partition.keys() if com in partition[node]]
    cluster_size[com] = len(list_nodes)
  #..
  stuff_papers = dict()
  stuff_authors = dict()
  for com in cluster_size: 
    if cluster_size[com]>thr: 
      stuff_papers[com]=[]
      stuff_authors[com]=dict()

  ## INPUT DATA
  # all authors
  my_auth=dict()
  src2 = os.path.join(in_dir, "authors.dat") 
  pl = Utils.Author()
  pl.read_file(src2)
  for l in pl.authors:
    if l.id not in my_auth: my_auth[l.id]=l.author
    else: my_auth[l.id]+= ', ' + l.author

  # abstract
  """
  my_abs=dict()
  src8  = os.path.join(in_dir, "abstracts.dat") 
  pl = Utils.Abstract()
  pl.read_file(src8)
  for l in pl.abstracts:
    if l.id not in my_abs: 
      my_abs[l.id]=l.abstract
  """

  # article
  src1  = os.path.join(in_dir, "articles.dat") 
  pl = Utils.Article()
  pl.read_file(src1)  
  for l in pl.articles:
    if l.id in partition:
      ## dealt with cases when article belong to several com
      for com in partition[l.id]:
        if com in stuff_papers: 
          footitle=l.title.replace('&','\&')
          if len(footitle)>limtitle: 
            aux = footitle[0:limtitle].rfind(' ')
            footitle=footitle[0:aux] + "..."
          if l.id in my_auth: 
            authors=my_auth[l.id]
            for auth in authors.split(', '):
              if auth not in stuff_authors[com]:stuff_authors[com][auth]=[0,0,0]
              stuff_authors[com][auth][0]+=int(l.times_cited)
              stuff_authors[com][auth][1]+=1
              stuff_authors[com][auth][2]+=dgcl_id[l.id]
          else: authors= ''
          stuff_papers[com].append([l.firstAU, l.year, footitle, l.journal.replace('&','\&'), l.volume, l.doctype, int(l.times_cited), dgcl_id[l.id], authors ])  

  # TREAT DATA 
  #.. define average degree of an author's papers in each cluster
  #for com in stuff_authors: 
  #  for auth in stuff_authors[com]: stuff_authors[com][auth][2]*=1.0/stuff_authors[com][auth][1]

  #.. prep dict
  for com in stuff_papers: 
    CR_papers[com]=dict()
    CR_authors[com]=dict()
    for KK in ['MC','MC_K90','MC_K95','MR','MR_TC90','MR_TC95', 'MR_TCsup5']: CR_papers[com][KK]=dict()
    for KK in ['MC','MC_K50','MC_K80','MR','MR_TC50','MR_TC80']: CR_authors[com][KK]=dict()

  #.. compute most cited , most cited with d>avg(d), most representative, most representative with TC > avg(TC) for papers and authors in clusters 
  for com in stuff_papers:
    # PAPERS
    foo = stuff_papers[com]
    """
    auxTC=[elt[6] for elt in foo]
    auxTC.sort()
    avgTC=sum(auxTC)/len(foo)
    TC90=auxTC[int(round(len(foo)*0.9))-1]
    TC95=auxTC[int(round(len(foo)*0.95))-1]
    auxK=[elt[7] for elt in foo]
    auxK.sort()
    avgK=sum(auxK)/len(foo)
    K90=auxK[int(round(len(foo)*0.9))-1]
    K95=auxK[int(round(len(foo)*0.95))-1]
    """

    # most cited papers
    foo.sort(key=lambda e:-e[6])
    for k in range(min(LL,len(foo))):
      CR_papers[com]['MC'][k]=foo[k] 
    """
    # most cited papers with ...
    foof=[elt for elt in foo if elt[7]>K90]
    foof.sort(key=lambda e:-e[6])
    CR_papers[com]['MC_K90']['p']=K90
    for k in range(min(LL,len(foof))):
      CR_papers[com]['MC_K90'][k]=foof[k] 
    #  
    foof=[elt for elt in foo if elt[7]>K95]
    foof.sort(key=lambda e:-e[6])
    CR_papers[com]['MC_K95']['p']=K95
    for k in range(min(LL,len(foof))):
      CR_papers[com]['MC_K95'][k]=foof[k]  
    """     
    # most representative papers (in terms of degree)
    foo.sort(key=lambda e:-e[7])
    for k in range(min(LL,len(foo))):
      CR_papers[com]['MR'][k]=foo[k] 
    """
    # most representative papers cited at least 5 times
    foof=[elt for elt in foo if elt[6]>5]
    foof.sort(key=lambda e:-e[7])
    for k in range(min(LL,len(foof))):
      CR_papers[com]['MR_TCsup5'][k]=foof[k]       
    # most representative papers (in terms of degree) with ...
    foof=[elt for elt in foo if elt[6]>TC90]
    foof.sort(key=lambda e:-e[7])
    CR_papers[com]['MR_TC90']['p']=TC90
    for k in range(min(LL,len(foof))):
      CR_papers[com]['MR_TC90'][k]=foof[k] 
    #
    foof=[elt for elt in foo if elt[6]>TC95]
    foof.sort(key=lambda e:-e[7])
    CR_papers[com]['MR_TC95']['p']=TC95
    for k in range(min(LL,len(foof))):
      CR_papers[com]['MR_TC95'][k]=foof[k] 
    """
      

    # AUTHORS
    foo = list(stuff_authors[com].items())
    """
    auxTC=[elt[1][0] for elt in foo]
    auxTC.sort()
    avgTC=sum(auxTC)/len(foo)
    TC50=auxTC[int(round((len(foo)-1)*0.5))]
    TC80=auxTC[int(round((len(foo)-1)*0.8))]
    auxK=[elt[1][2] for elt in foo]
    auxK.sort()
    avgK=sum(auxK)/len(foo)
    K50=auxK[int(round((len(foo)-1)*0.5))]
    K80=auxK[int(round((len(foo)-1)*0.8))]
    """

    # most cited authors:
    foo.sort(key=lambda e:-e[1][0])
    for k in range(min(LL,len(foo))):
      CR_authors[com]['MC'][k]=foo[k]
    """
    # most cited authors with ...:
    foof=[elt for elt in foo if elt[1][2]>K50]
    foof.sort(key=lambda e:-e[1][0])
    CR_authors[com]['MC_K50']['p']=K50
    for k in range(min(LL,len(foof))):
      CR_authors[com]['MC_K50'][k]=foof[k]
    #
    foof=[elt for elt in foo if elt[1][2]>K80]
    foof.sort(key=lambda e:-e[1][0])
    CR_authors[com]['MC_K80']['p']=K80
    for k in range(min(LL,len(foof))):
      CR_authors[com]['MC_K80'][k]=foof[k]   
    """   
    # most representative authors (in terms of degree)
    foo.sort(key=lambda e:-e[1][2])
    for k in range(min(LL,len(foo))):
      CR_authors[com]['MR'][k]=foo[k]
    """
    # most representative authors with ...:
    foof=[elt for elt in foo if elt[1][0]>TC50]
    foof.sort(key=lambda e:-e[1][2])
    CR_authors[com]['MR_TC50']['p']=TC50
    for k in range(min(LL,len(foof))):
      CR_authors[com]['MR_TC50'][k]=foof[k]
    # 
    foof=[elt for elt in foo if elt[1][0]>TC80]
    foof.sort(key=lambda e:-e[1][2])
    CR_authors[com]['MR_TC80']['p']=TC80
    for k in range(min(LL,len(foof))):
      CR_authors[com]['MR_TC80'][k]=foof[k]
    """

  return (CR_papers, CR_authors)


## ##################################################

def my_writeMCP(f_out,my_stuff,nbline):
    for i in range(min(nbline,len(my_stuff)-1)):
        artref= "%s, %d, %s(%s)" % (my_stuff[i][0].replace('[','$[$').replace(']','$]$'), my_stuff[i][1], my_stuff[i][3], my_stuff[i][4])
        f_out.write("%d & %d & %s & %s\\\\\n" % ( my_stuff[i][7], my_stuff[i][6], artref, my_stuff[i][2]) ) 
    f_out.write("\hline\n")      

def my_writeMCA(f_out,my_stuff,nbline):
    for i in range(min(nbline,len(my_stuff)-1)):
        f_out.write("%s & %d & %d & %d \\\\\n" % ( my_stuff[i][0].replace('[','$[$').replace(']','$]$'), my_stuff[i][1][1], my_stuff[i][1][0], my_stuff[i][1][2]  ) ) 
    f_out.write("\hline\n")  


## ##################################################
## ########################### Direction / Glue
## ##################################################

def direction_and_glue(G,in_dir,out_dir,partition,bcthr,thr,verbose,txt):

  #. initialize
  list_nodes = dict();
  comm_size = dict();
  for com in set(list(partition.values())) :
    list_nodes[com] = [nodes for nodes in partition.keys() if partition[nodes] == com]
    comm_size[com] = len(list_nodes[com])
  #    
  list_comm = [] ## list of communities of size > thr
  list_id = dict(); ## list of id -- articles within the BC network and within a community of size > thr
  for com in comm_size: 
    if comm_size[com] > thr: list_comm.append( com )
  for id_art in partition:
    com = partition[id_art]
    if comm_size[com] > thr: list_id[id_art] = ''   

  #. cluster ref
  if verbose: print ("....prep clusters refs table")
  freqR  = dict(); ## records the freq of each ref in each community of size > thr
  norm = dict(); ## tot of previous doc for normalization
  for com in list_comm: 
    freqR[com] = dict()
    norm[com] = 0

  with open(os.path.join(in_dir, "references.dat") ,"r") as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1] 
  aux = [(l.split("\t")[0],", ".join(l.split("\t")[1:]).replace(',0,0','')) for l in data_lines]
  del data_lines;
  for elt in aux:
      reference = elt[1];
      if elt[0] in list_id: 
        for com in partition[elt[0]]:
          if reference not in freqR[com]: freqR[com][reference] = 0
          freqR[com][reference] += 1
  del aux; 
  #      
  for com in list_comm:
    for ref in freqR[com]:
      norm[com] += freqR[com][ref]*freqR[com][ref]
    norm[com] = math.sqrt(norm[com]) 
  #
  name = "%s_refs_thr%d.csv" % (txt,bcthr)
  dst = os.path.join(out_dir, name)      
  f_r = open(dst,'w')
  f_r.write('com\tref\tn\n')      
  for com in list_comm:
    aux=freqR[com].items()
    aux.sort(key=lambda e:-e[1]) 
    for kpt in range(min(len(aux),50)):
      f_r.write('%d\t%s\t%d\n' % (com,aux[kpt][0],aux[kpt][1]))      
  f_r.close() 

  #. direction of links btw clusters
  if verbose: print ("....prep coupling refs table")
  name = "%s_couplingrefs_thr%d.csv" % (txt,bcthr)
  dst = os.path.join(out_dir, name)
  f_cr = open(dst,'w')
  f_cr.write('com1\tcom2\tref\tn12\tn1\tn2\n')
  direction = dict()
  for com1 in list_nodes:
    for com2 in list_nodes:
      size1 = len(list_nodes[com1]); size2 = len(list_nodes[com2]);
      if size1 > thr and size2 > thr and com1 > com2:
        normsh=0; cos12=0; cos1=0; cos2=0;
        shared_ref=[ref for ref in freqR[com2] if ref in freqR[com1]]
        aux=dict()
        for ref in shared_ref:
          cos12 += freqR[com1][ref]*freqR[com2][ref]
          cos1 += freqR[com1][ref]*freqR[com1][ref]*freqR[com2][ref]
          cos2 += freqR[com2][ref]*freqR[com1][ref]*freqR[com2][ref]
          normsh += freqR[com1][ref]*freqR[com2][ref]*freqR[com1][ref]*freqR[com2][ref]
          aux[ref] = freqR[com1][ref]*freqR[com2][ref]
        normsh=math.sqrt(normsh) 
        if normsh >0:
          cos12 *= 1.0/(norm[com1]*norm[com2])
          cos1 *= 1.0/(norm[com1]*normsh)
          cos2 *= 1.0/(norm[com2]*normsh)
        direction[(com1,com2)]=[cos1,cos2,cos12]
        #.. top 50 coupling ref
        coupling=aux.items()
        coupling.sort(key=lambda e:-e[1])
        if cos1<cos2:
          for kpt in range(min(len(coupling),50)):
            ref=coupling[kpt][0]
            f_cr.write('%d\t%d\t%s\t%d\t%d\t%d\n' % (com1, com2, ref, aux[ref], freqR[com1][ref], freqR[com2][ref]))
        else:
          for kpt in range(min(len(coupling),50)):
            ref=coupling[kpt][0]
            f_cr.write('%d\t%d\t%s\t%d\t%d\t%d\n' % (com2, com1, ref, aux[ref], freqR[com2][ref], freqR[com1][ref]))    
  f_cr.close()      

  """
  #. output links with local glue > 0
  list_ref=dict()
  # list_ref['Gardner H, 1983, FRAMES MIND THEORY M']=''
  if len(list_ref) >0:
    if verbose: print ("....compute local glue")
    for r in list_ref: print r
    name = "%s_links_dir_glue.txt" % (txt)
    dst = os.path.join(out_dir, name)
    f_gephi = open(dst,'w')
    f_gephi.write("edgedef>node1 VARCHAR,node2 VARCHAR,num_links DOUBLE,weight DOUBLE,logweight DOUBLE,cos1 DOUBLE,cos2 DOUBLE")
    for kk in range(len(list_ref)):
      f_gephi.write(',locglue%s DOUBLE' % kk)
    f_gephi.write('\n')  
    for com1 in list_nodes:
      for com2 in list_nodes:
        size1 = len(list_nodes[com1]); size2 = len(list_nodes[com2]);
        if size1 > thr and size2 > thr and com1 > com2:
          W = 0; N = 0;
          for id1 in list_nodes[com1]:
            for id2 in list_nodes[com2]:
              if id2 in G.edge[id1]: 
                N += 1
                W += G.edge[id1][id2]['weight']
          W *= 1000.0 / (size1 * size2)
          cos1=direction[(com1,com2)][0]
          cos2=direction[(com1,com2)][1]
          if W > 0.000001:
            if cos1 < cos2:
              f_gephi.write("%d,%d,%d,%1.9f,%1.2f,%.5f,%1.5f" % (com1, com2, N, W, 6 + math.log(W)/math.log(10),cos1,cos2))
            else:
              f_gephi.write("%d,%d,%d,%1.9f,%1.2f,%.5f,%1.5f" % (com2, com1, N, W, 6 + math.log(W)/math.log(10),cos2,cos1))  
            #
            for r in list_ref:
              z=0
              if r in freqR[com1] and r in freqR[com2]: z=1 
              f_gephi.write(',%d' % z)
            f_gephi.write('\n')  
    ## ... end
    f_gephi.close() 
  """
  
  #. compute global glue
  if verbose: print ("....compute global glue")
  glue = dict()
  for com1 in list_nodes:
    for com2 in list_nodes:
      size1 = len(list_nodes[com1]); size2 = len(list_nodes[com2]);
      if size1 > thr and size2 > thr and com1 > com2:
        #.. num_links
        N = 0;
        for id1 in list_nodes[com1]:
          for id2 in list_nodes[com2]:
            if id2 in G.edge[id1]: 
              N += 1
        #..add all ref in glue
        if N>0:
          for ref in freqR[com1]:
            if ref in freqR[com2]:
              if ref not in glue: glue[ref]=0
              glue[ref] += freqR[com1][ref]*freqR[com2][ref]*1.0/N

  #. normalize factor
  ZZ=sum(list(glue.values()))
  #. export the 1000 more gluing refs   
  stuff = dict()
  L = glue.items()
  L.sort(key=lambda e:-e[1])
  for i in range(min(1000,len(L))):
    r=L[i][0]
    g=L[i][1]*100.0/ZZ
    stuff[i]=[r,g]

  #. output
  fooname="%s_glue_thr%d.txt" % (txt,bcthr)
  filename = os.path.join(out_dir, fooname)
  f_out = open(filename,'w')
  if verbose: print ("....output most gluing references")
  for i in range(len(stuff)):
    f_out.write('%s\t%.3f\n' % (stuff[i][0],stuff[i][1]))
  f_out.close() 

  ## END
  return direction
       

## ##################################################
## ##################################################


def main():
# usage: BC.py [-h] [--version] -i DIR -o DIR -p FILE [-v]
#
# optional arguments:
#   -i DIR, --input_dir DIR input directory name
#   -o DIR, --output_dir DIR input directory name
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
          metavar='DIR')
          
  parser.add_argument("-p", "--partition",
          action = "store", dest="partition")

  parser.add_argument("-thr", "--thr", nargs=1, type=int,
          action = "store", dest="thr", default=[10],
          help="",
          metavar='thr')

  parser.add_argument("-v", "--verbose",
          action = "store_true", dest="verbose",
          default = False,
          help="verbose mode [default %(default)s]")

  #Analysis of input parameters
  args = parser.parse_args()
  
  if (not os.path.exists(args.in_dir[0])):
      print ("Error: Input directory does not exist: ", args.in_dir[0])
      exit()

  if (not os.path.exists(args.partition[0])):
      print ("Error: Input partition does not exist: ", args.partition[0])
      exit()

  if (not os.path.exists(args.out_dir[0])):
      print ("Error: Output directory does not exist: ", args.out_dir[0])
      exit()

  ##

  partition = dict();      

  f_in = open(args.partition,'r')
  for line in f_in.readlines():
    foo = line.stript().split('\t')
    if len(foo) == 2:
      partition[int(foo[0])] = int(foo[1])
  f_in.close()
  
  ##
  comm_tables(args.in_dir[0],partition,thr,args.verbose)

  return


    
## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################
