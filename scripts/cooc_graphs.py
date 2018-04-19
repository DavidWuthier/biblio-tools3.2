#! /usr/bin/env python

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2017
"""

# usage: cooc_graph.py -i DIR [-v]
# 

import os
import time
import math
import argparse
import itertools

## ##########################################################################
## SOME VARIABLES
# default threshold value: "K":50 means that the keyword co-occurrence network will features keywords appearing in at least 50 publications of the corpus
# if you don't want a certain type of item, put a riducously high threshold, eg 1000000
COOC_THR={"AU":50,"CU":1,"I":50,"S":1,"S2":1,"K":50,"AK":50,"TK":50,"R":30,"RJ":200,"Y":1,"J":20,"DT":1,"LA":1}
## ##########################################################################


## ##################################################
## ##################################################
## ##################################################

def cooc_graph(in_dir, ini_suff, timeWND, prephet, verbose):

  THR=COOC_THR
  # for heteregeneous graph
  global HET_SUFF, HET_TABLE;
  HET_TABLE={};
  HET_SUFF='';
  # default color for gephi display
  color_nodes = {'Y': '255,255,0', 'J': '150,0,150', 'AU': '20,50,255', 'K': '255,0,255', 'AK': '255,0,255', 'TK': '205,0,205', 'S': '50,0,150', 'S2': '50,0,150', 'R': '255,0,0', 'RJ': '255,97,0', 'I': '0,255,0', 'CU': '0,255,255', 'LA': '0,180,0', 'DT': '0,180,0'}


  ## INITIALIZATION
  t1=time.time()
  if verbose: print ("INITIALIZATION")
  with open(os.path.join(in_dir, "articles.dat") ,"r", encoding='utf-8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1]  
  # filter if necessary
  if timeWND!=[]: 
    data_lines=[l for l in data_lines if (int(l.split("\t")[2])>=timeWND[0] and int(l.split("\t")[2])<=timeWND[1] )]
    tokeep=dict()
    for l in data_lines: tokeep[l.split("\t")[0]]=''
    if len(tokeep)==0: return
  # total nb_art  
  nb_art=len(data_lines)
  print ("... %d publis" % nb_art)

  # prep out folder
  outname='gdffiles'
  out_dir=os.path.join(in_dir,outname)
  if not os.path.exists(out_dir): os.makedirs(out_dir)


  ## ##################################################
  ## AUX FUNCTIONS
  def generate_cooc(stuff,nick):
    global HET_TABLE, HET_SUFF
    # HOW MUCH TIME AN ITEM APPEARS?
    stuff.sort(key=lambda e:e[1])
    # duplicates are removed: we count in how many distinct publis an item appears 
    foo=[(xx,len(set(list(art_id)))) for xx,art_id in itertools.groupby(stuff,key=lambda e:e[1])]
    foo.sort(key=lambda x:-x[1])
    # keep only the items appearing more than THR
    keep=[foo[i] for i in range(len(foo)) if foo[i][1]>=THR[nick] ]

    if (len(keep)>0):
      if prephet: HET_SUFF+='_'+nick+str(THR[nick])
      outfilename=os.path.join(out_dir,'cooc_'+nick+ini_suff)
      if timeWND != []: outfilename+='_'+str(timeWND[0])+'_'+str(timeWND[1])
      outfilename+='_thr'+str(THR[nick])+'.gdf'
      with open(outfilename, 'w') as f_gephi:
        nodeID={};aux={};comm={};
        # nodes
        f_gephi.write("nodedef>name VARCHAR,label VARCHAR,type VARCHAR,width DOUBLE,height DOUBLE,size DOUBLE,color VARCHAR\n")
        for i in range(len(keep)):
          nodeID[keep[i][0]]=i
          size=keep[i][1]
          f_gephi.write("%d,'%s',%s,%f,%f,%d,'%s'\n" % (i, keep[i][0].replace('&','-'), nick, math.sqrt(size),math.sqrt(size), size, color_nodes[nick]))
        # links
        f_gephi.write("edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE,nb_cooc DOUBLE\n")
        for elt in stuff:
          if elt[1] in nodeID: 
            if elt[0] not in aux: aux[elt[0]]=[]
            aux[elt[0]].append(nodeID[elt[1]])
            if prephet: 
              if (nick,elt[1]) not in HET_TABLE: HET_TABLE[(nick, elt[1])]=[]
              HET_TABLE[(nick, elt[1])].append( elt[0] ) 
        for ee in aux:
          for itema in aux[ee]:
            for itemb in [ff for ff in aux[ee] if ff > itema]:
               if (itema, itemb) not in comm: comm[(itema, itemb)]=0
               comm[itema,itemb]+=1
        for (a,b) in comm:
          f_gephi.write("%d,%d,%.6f,%d\n" % (a, b, comm[(a,b)]/math.sqrt(keep[a][1]*keep[b][1]), comm[(a,b)]))
      # no cooccurence between years or publication sources: remove these files  
      if nick in ['J','Y']: os.remove(outfilename)

      # free memory
      del aux, stuff, keep, nodeID, comm

    return

  ## ################
  def output_hetgraph():
    global HET_TABLE, HET_SUFF
    print("... Prep het graph")
    outfilename=os.path.join(out_dir,'cooc_heterogeneous'+ini_suff)
    if timeWND != []: outfilename+='_'+str(timeWND[0])+'_'+str(timeWND[1])
    outfilename+=HET_SUFF+'.gdf'
    with open(outfilename, 'w') as f_gephi:
      nodeID={};aux={};comm={};size={};
      # nodes
      print(".... nodes")
      f_gephi.write("nodedef>name VARCHAR,label VARCHAR,type VARCHAR,width DOUBLE,height DOUBLE,size DOUBLE,color VARCHAR\n")  
      i=0
      for (t,x) in HET_TABLE:
        nodeID[(t,x)]=i
        size[i]=len(HET_TABLE[(t,x)])
        f_gephi.write("%d,'%s',%s,%f,%f,%d,'%s'\n" % (i, x.replace('&','-'), t, math.sqrt(size[i]),math.sqrt(size[i]), size[i], color_nodes[t]))
        i+=1
      # links
      e = len(HET_TABLE) * len(HET_TABLE) / 2; ee = 0; p=5; 
      print(".... links")
      f_gephi.write("edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE,nb_cooc DOUBLE\n")
      for elt in HET_TABLE:
        for pub in HET_TABLE[elt]:
          if pub not in aux: aux[pub]=[]
          aux[pub].append(nodeID[elt])
      for pub in aux:
        for itema in aux[pub]:
          for itemb in [ff for ff in aux[pub] if ff > itema]:
             if (itema, itemb) not in comm: comm[(itema, itemb)]=0
             comm[itema,itemb]+=1
      for (a,b) in comm:
        f_gephi.write("%d,%d,%.6f,%d\n" % (a, b, comm[(a,b)]/math.sqrt(size[a]*size[b]), comm[(a,b)]))     

    # free memory
    del aux, comm, nodeID, HET_TABLE

    return


  ## ################

  def treat_item(item, nick, which):
    if verbose: print ("... dealing with %s" % (item) )
    with open(os.path.join(in_dir, item+".dat") ,"r",encoding='utf-8') as file:
      # dat file have one trailing blank line at end of file
      data_lines=file.read().split("\n")[:-1] 
    if (timeWND == []): aux = [(l.split("\t")[0], l.split("\t")[which]) for l in data_lines]
    else: aux = [(l.split("\t")[0], l.split("\t")[which]) for l in data_lines if l.split("\t")[0] in tokeep]
    generate_cooc(aux,nick)  
    # free memory
    del data_lines 
    # end
    return

  ## ##################################################
  ## TREAT DATA 
  if prephet:    
    if verbose: print ("... dealing with years + journals + doctypes + languages" )
    # data_lines was already filtered to keep only publi in the correct time-window
    #.. Y
    aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines]
    generate_cooc(aux,'Y')
    #.. J
    aux = [(l.split("\t")[0], l.split("\t")[3]) for l in data_lines]
    generate_cooc(aux,'J')

  #.. D
  aux = [(l.split("\t")[0], l.split("\t")[8]) for l in data_lines]
  generate_cooc(aux,'DT') 
  #.. L
  aux = [(l.split("\t")[0], l.split("\t")[9]) for l in data_lines]
  generate_cooc(aux,'LA')    
  #.. free memory
  del data_lines
  
  treat_item('authors', 'AU', 2)
  treat_item('subjects', 'S', 1)
  treat_item('institutions', 'I', 2)
  treat_item('countries', 'CU', 2)
  with open(os.path.join(in_dir, "database.dat"),'r') as ff:
    if ff.read()=='Scopus': treat_item('subjects2', 'S2', 1)

  if verbose: print ("... dealing with keywords + title words") 
  with open(os.path.join(in_dir, "keywords.dat") ,"r",encoding='utf-8') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1]
  if timeWND!=[]: data_lines=[l for l in data_lines if (l.split("\t")[0] in tokeep) ]
  #.. IK
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines if l.split("\t")[1]=='IK' ]
  generate_cooc(aux,'K')
  #.. AK
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines if l.split("\t")[1]=='AK' ]
  generate_cooc(aux,'AK')  
  #.. TK
  aux = [(l.split("\t")[0], l.split("\t")[2]) for l in data_lines if l.split("\t")[1]=='TK' ]
  generate_cooc(aux,'TK')  
  #.. free memory
  del data_lines 


  if verbose: print ("... dealing with references + references journals") 
  with open(os.path.join(in_dir, "references.dat") ,"r",encoding='utf-8', errors='ignore') as file:
    # dat file have one trailing blank line at end of file
    data_lines=file.read().split("\n")[:-1]
  if timeWND!=[]: data_lines=[l for l in data_lines if (l.split("\t")[0] in tokeep) ] 
  #.. R
  aux = [(l.split("\t")[0],", ".join(l.split("\t")[1:5]).replace(',0,0','').replace(', 0, 0','').replace(', 0','')) for l in data_lines]
  generate_cooc(aux,'R')  
  #.. RJ
  aux = [(l.split("\t")[0], l.split("\t")[3]) for l in data_lines]
  generate_cooc(aux,'RJ')  
  #.. free memory
  del data_lines

  ## ##################################################
  ## OUTPUT HETEROGENEOUS GRAPH IF ASKED
  if prephet: output_hetgraph()


  ## ##################################################
  ## END    
  print ('.. Time needed: %ds' % (time.time()-t1))
  return

## ##################################################
## ##################################################
## ##################################################


def main():
# usage: cooc_graph.py [-h] [--version] -i DIR [-v] [-H INT] [-O]
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

  parser.add_argument("-hg", "--prephet",
          action = "store_true", dest="prephet",
          default = False,
          help="prepare heterogeneous graph? [default %(default)s]")
          
  parser.add_argument("-v", "--verbose",
          action = "store_true", dest="verbose",
          default = False,
          help="verbose mode [default %(default)s]")

  #Analysis of input parameters
  args = parser.parse_args()
  if (not os.path.exists(args.in_dir[0])):
      print ("Error, Input directory does not exist: ", args.in_dir[0])
      exit()

  cooc_graph(args.in_dir[0], '', [], args.prephet, args.verbose)


  return

## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################
