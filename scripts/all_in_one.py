#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2017
"""

# usage: allinone.py -i DIR [-o DIR] [-d scopus] [-e]
# -i dir indicates the folder where your wos / scopus data is
# -o dir indicates the main folder where all the files will be created


import os
import time
import numpy
import math
import argparse
import itertools
from collections import Counter
import json
import Utils.Utils as Utils
import Utils.BCUtils as BCUtils
import Utils.community as community
import networkx as nx
from fnmatch import fnmatch


from biblio_parser import biblio_parser
from describe_corpus import describe_corpus
from biblio_coupling import *
from cooc_graphs import cooc_graph

## ON WHICH TIME PERIOD DO YOU WISH TO MAKE ANALYSIS? 
# Leave timeWND=[] for all years in corpus, create an array of ymin-ymax otherwise, eg 
# if you want successive windows of 5 years between 1990 and 2019 (ie 90-94, 95-99, 00-04, etc), write sth like
# timeWND=[[yy,yy+4] for yy in range(1990,2019,5)]
# if you want windows of 10 years between 1990 and 2019 with a 5 years overlap (ie 90-99, 95-04, 00-09, etc), write sth like
# timeWND=[[yy,yy+9] for yy in range(1990,2015,5)]
# in general, if you want windows of N years between YMIN and YMAX with a X years overlap, edit the following:
# timeWND=[[yy,yy+N-1] for yy in range(YMIN,YMAX,N-X)]
timeWND=[]   
  

## WHICH ANALYSIS SHOULD WE PERFORM? (remove those you don't want from the list) 
# 'PARSE' will do the first parsing of the data (you NEED to do it at once)
# 'CD' stands for corpus description on each time period
# 'BC' stands for bibliographic coupling on each time period
# 'CO' stands for co-occurrence graph on each time period 
# 'MATCH' will perform a SIMPLE matching analysis between the BC clusters of different time periods. You may change the parameter THETA: we keep only links between clusters of diff periods whose Kessler's cosines are above this threshold. Default: 0.05
#
whichanalysis=["PARSE", "CD", "BC", "CO", "MATCH"]
THETA=0.05

## ##################################################
## ##################################################

def allinone(dataraw_dir, main_dir, database, expert):
  ## INITIALIZATION
  tzero=time.time()  
  
  ## 1st step: parse & describe the whole corpus
  if "PARSE" in whichanalysis:
    print("PARSE THE WHOLE CORPUS")
    biblio_parser(dataraw_dir,main_dir,database,expert)
    print("DESCRIBE THE WHOLE CORPUS")
    out_dir=os.path.join(main_dir, 'freqs')
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    describe_corpus(main_dir, out_dir, False)
  

  ## 2nd step: do the required analysis on whole corpus if no time window was given
  if timeWND==[]: 
    if "BC" in whichanalysis: 
      print ("BC ANALYSIS ON WHOLE CORPUS")
      BC_network(main_dir,main_dir,'',False,False,False)
    if "CO" in whichanalysis:
      print ("COOC ANALYSIS ON WHOLE CORPUS")
      cooc_graph(main_dir,'', [], True, False)


  ## 3rd step: filter, describe and analyse the corpus in the different time windows if any
  if timeWND!=[]: print("studied periods:", timeWND)
  if (len(timeWND)>0 and not os.path.exists(os.path.join(main_dir, 'dataPeriods')) ): os.makedirs(os.path.join(main_dir, 'dataPeriods'))
  if ("CD" in whichanalysis or "BC" in whichanalysis or "CO" in whichanalysis):
    for period in timeWND:
      print ("DEALING WITH TIME PERIOD %d-%d" % (period[0],period[1]))
      suff='_'+str(period[0])+'_'+str(period[1])
      # filter data
      print("..FILTER")
      datadir=os.path.join(main_dir, 'dataPeriods/data'+suff)
      if not os.path.exists(datadir): 
        os.makedirs(datadir)
        PUBYfilter(main_dir, datadir, period)
      # describe data
      if "CD" in whichanalysis: 
        print("..DESCRIBE")
        out_dir=os.path.join(main_dir, 'freqs'+suff)
        if not os.path.exists(out_dir): os.makedirs(out_dir)
        describe_corpus(datadir, out_dir, False)
      # bc analysis?
      if "BC" in whichanalysis:
        print("..BC ANALYSIS")
        BC_network(datadir,main_dir,suff,False,False,False)
      if "CO" in whichanalysis:
        print ("COOC ANALYSIS")
        cooc_graph(datadir,suff,[],True,False)

  ## 4th step: matching BC clusters from different time periods
  # Right now, IT'S MADE IN A "QUICK AND DIRTY" WAY: compute the bc links between clusters of different time periods (and NOT the links between clusters of the same period), keep only those above a given threshold (this way only the strongest links will remain), define the predecessor / successor of each cluster as the one in the period before / after with which the links is higher, build streams or "meta-clusters" as chain of clusters which are successor/predessesors to one another.
  # ==> this is a SIMPLE attempt at building an history, not a definitive one!  
  # ==> the meta-clusters IDs will overwrite the "group" data in the json files used for the colors of the clusters in the BiblioMaps viz.  
  if "MATCH" in whichanalysis:
    print("..MATCHING CLUSTERS OF DIFF. PERIODS")
    # upload (top) partitions limited to the list of clusters kept in the viz
    theparts=dict()
    for period in timeWND:
      suff='_'+str(period[0])+'_'+str(period[1])
      datadir=os.path.join(main_dir, 'dataPeriods/data'+suff)
      theparts[suff]=extract_part(datadir, main_dir, suff)
    # compute the BC network between these clusters
    (GG,LABEL)=prep_bcnetw(main_dir,theparts,THETA)
    # detect meta-partition
    metapart=detect_metapart(GG, LABEL) 
    # overwrite the groups in BCclusters_xxx.json files
    for period in timeWND:
      suff='_'+str(period[0])+'_'+str(period[1])
      newgroup=dict()
      for elt in [ee for ee in LABEL if ee.split('Z')[1]==suff]:
        newgroup[int(elt.split('Z')[0])]=metapart[LABEL[elt]]
      cfile=os.path.join(main_dir, 'jsonfiles/BCclusters'+suff+part_suffixBC+'.json')
      updategroups(cfile,newgroup)


  #... LOG
  with open(os.path.join(main_dir,'AA_log.txt'),'a') as out:
    mye='';
    if (expert): mye='-e';
    out.write("\n\n%s GMT\nall_in_one.py -i %s -o %s -d %s %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), dataraw_dir, main_dir, database, mye))
    if "CD" in whichanalysis and timeWND!=[]:
      out.write("TO EXPLORE THE CORPUSES OF ALL PERIODS IN BIBLIOMAPS:\nCopy / paste these lines on appropriate place (should be line 116) in the 'corpus.html' file + switch the 'howmany' variable (should be on line 100) from 'one' to 'several'\n")
      out.write("\ttimeWND=%s\n" % timeWND)
    if "BC" in whichanalysis and timeWND!=[]:
      out.write("TO EXPLORE THE BC CLUSTERS OF ALL PERIODS IN BIBLIOMAPS:\n Copy / paste these lines on appropriate place (should be line 205) in the 'corpus.html' file + switch the 'howmany' variable (should be on line 190) from 'one' to 'several'\n")
      out.write("\ttimeWND=%s\n" % timeWND)
      out.write("\tfiles_suffix='%s'\n" % part_suffixBC)


  ## END
  print ('END of allinone - TOTAL TIME NEEDED: %ds' % (time.time()-tzero))

  return

## ##########################################################################################################
## ##########################################################################################################

def PUBYfilter(in_dir,out_dir,per):
  ## INITIALIZATION
  ymin=per[0]
  ymax=per[1]
  # files to filter (will not be taken into account if they don't exist)
  filenames=["articles", "addresses", "authors", "countries", "institutions", "journals", "keywords", "references", "subjects", "subjects2", "AUaddresses", "cities", "fundingtext", "abstracts"] 

  ## ##################################################
  ## SELECT ID OF PUBLICATIONS TO KEEP
  with open(os.path.join(in_dir, "articles.dat") ,"r") as file: data_lines=file.read().split("\n")[:-1]
  nbpub=len(data_lines)  
  TOKEEP=dict([(l.split("\t")[0],'') for l in data_lines if (int(l.split("\t")[2])<=ymax and int(l.split("\t")[2])>=ymin)])
  del data_lines
  print ("..%d publications selected out of %d (%.2f%%)" % (len(TOKEEP), nbpub, len(TOKEEP)*100.0/nbpub) )

  ## ##################################################
  ## PROCEED
  for fff in filenames:
    if (os.path.exists(os.path.join(in_dir, fff+".dat"))):
      print (".... filtering %s" % fff)
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
  return

## ##################################################
## ##################################################

def extract_part(datadir, maindir, suff):
  # upload list of cluster to keep
  ffin = open(os.path.join(maindir, "jsonfiles/BCdefaultVAR"+suff+part_suffixBC+".json"),'r')
  defdata = json.loads(ffin.read())
  tokeep=dict([(elt[0],'') for elt in defdata["TOPpositions"]])
  #
  #... upload previously computed partition
  fooname="partitions%s.txt" % (suff+part_suffixBC)
  ffin = open(os.path.join(datadir, fooname),'r')
  foo = ffin.read()
  lines = foo.split('\n')
  louvain_partition = json.loads(lines[0])
  parttokeep={}
  aux=str(len(louvain_partition)-2)
  for key in louvain_partition[aux]:
    if louvain_partition[aux][key] in tokeep:
      parttokeep[int(key)]=louvain_partition[aux][key]
  # end
  return parttokeep

## ##################################################
## ##################################################

def prep_bcnetw(main_dir,theparts,THETA):
  print ("..Create the 'Bibliographic Coupling' weight table")
  nR = dict() # store the number of ref of each publication
  ref_table = dict() # store the id of articles using a given ref
  BC_table = dict() # store the number of common refs between pairs of publications

  print ("....loading refs table")
  with open(os.path.join(main_dir, "references.dat") ,"r", encoding='utf8') as file: 
    data_lines=file.read().split("\n")[:-1] 
  for l in data_lines:
    foo =', '.join([l.split("\t")[k] for k in range(1,6)]).replace(',0,0','').replace(', 0, 0','').replace(', 0','')
    pubid=int(l.split("\t")[0])
    pubyear=int(l.split("\t")[2])
    if foo in ref_table: ref_table[foo].append( pubid )
    else: ref_table[foo] = [pubid]
    if pubid not in nR: nR[pubid]=0
    nR[pubid] += 1
  del data_lines

  print ("....detecting common references")
  for foo in ref_table:
    if len(ref_table[foo]) >= 2:
      for i in ref_table[foo]:
          for j in ref_table[foo]:
              if (i<j):
                  if i not in BC_table: BC_table[i] = dict()
                  if j not in BC_table[i]: BC_table[i][j] = 0
                  BC_table[i][j] += 1 
  del ref_table

  ## PREP NETWORK 
  print ("....define graph in networkx format")
  cluster_links=dict()
  LABEL=dict(); kompt=0;

  for partA in theparts:
    for i in theparts[partA]:
      foo=str(theparts[partA][i])+'Z'+partA
      if foo not in LABEL:
        LABEL[foo]=kompt; kompt+=1;
      if i in BC_table:
        for partB in [part for part in theparts if part!=partA]:
          for j in [jj for jj in BC_table[i] if jj in theparts[partB]]:
            w_ij = (1.0 * BC_table[i][j]) / math.sqrt(nR[i] * nR[j])
            foo=(str(theparts[partA][i])+'Z'+partA, str(theparts[partB][j])+'Z'+partB)
            if foo not in cluster_links: 
              cluster_links[foo]=[0,0]; cluster_links[(foo[1],foo[0])]=[0,0]
            cluster_links[foo][0]+=w_ij
            cluster_links[foo][1]+=1
      # in case of overlapping periods, we need to add potential self-loops from a publication to itself, EVEN IF i not in BC_table
      for partB in [part for part in theparts if part!=partA]:
        if i in theparts[partB]:
          foo=(str(theparts[partA][i])+'Z'+partA, str(theparts[partB][i])+'Z'+partB)
          if foo not in cluster_links: 
            cluster_links[foo]=[0,0]; cluster_links[(foo[1],foo[0])]=[0,0]
          cluster_links[foo][0]+=1
          cluster_links[foo][1]+=1
  del(BC_table)

  ## keep only links above THETA
  GG=nx.Graph()

  for (a,b) in cluster_links:
    if (a<=b):
      W=cluster_links[(a,b)][0]+cluster_links[(b,a)][0]
      N=cluster_links[(a,b)][1]+cluster_links[(b,a)][1]
      if W/N>=THETA:
        GG.add_edge(LABEL[a], LABEL[b], weight=W/N)

  return (GG, LABEL)

## ##################################################
## ##################################################

def detect_metapart(GG, LABEL):
  # louvain partition ==> problem: may group together clusters on same period with same predecessor!
  # foodendogram = community.generate_dendogram(GG, part_init=None)
  # metapart = community.partition_at_level(foodendogram,len(foodendogram)-1)
  # mod=community.modularity(metapart, GG)
  # print ("....top-clusters partition in %d meta-clusters, Q=%.4f" % ( len(set(metapart.values())), mod))
 
  # list all clusters by periods 
  allsuff=['_'+str(period[0])+'_'+str(period[1]) for period in timeWND]
  byperiod=dict([(suff,[]) for suff in allsuff])
  for lab in LABEL: byperiod[lab.split('Z')[1]].append(LABEL[lab])  
  for per in allsuff:
    print(".... %s: %d top clusters" % (per, len(byperiod[per]) ))

  # aux function
  def findmax(mylinks, suff):
    lookin=[(com, mylinks[com]["weight"]) for com in mylinks if com in byperiod[suff]]
    lookin.sort(key=lambda e:-e[1])
    if len(lookin)>0: return lookin[0][0]
    else: return -1

  # find successors / predecessors
  pred=dict();pred[-1]=-1;pred2=dict();pred2[-1]=-1;
  succ=dict();succ[-1]=-1;succ2=dict();succ2[-1]=-1;
  TT=len(allsuff)
  for per in range(TT):
    for com in byperiod[allsuff[per]]:
      pred[com]=-1; succ[com]=-1;
      pred2[com]=-1; succ2[com]=-1;
      if com in GG.edge:
        # predecessors
        if per>0: pred[com]=findmax(GG.edge[com], allsuff[per-1])
        if per>1: pred2[com]=findmax(GG.edge[com], allsuff[per-2])
        # successors
        if per<TT-1: succ[com]=findmax(GG.edge[com], allsuff[per+1])
        if per<TT-2: succ2[com]=findmax(GG.edge[com], allsuff[per+2])
      #print(com, pred[com], succ[com])

  """
  foo=[(xx,len(set(list(aux)))) for xx,aux in itertools.groupby(list(pred.items()),key=lambda e:e[1])]
  print ("%d 2-splits" % len([elt for elt in foo if elt[1]==2]))
  print ("%d 3p-splits" % len([elt for elt in foo if elt[1]>2]))

  foo=[(xx,len(set(list(aux)))) for xx,aux in itertools.groupby(list(succ.items()),key=lambda e:e[1])]
  print ("%d 2-merges" % len([elt for elt in foo if elt[1]==2]))
  print ("%d 3p-merges" % len([elt for elt in foo if elt[1]>2]))
  """

  # metapart = streams of succ / pred
  def buildstream(com):
    stream=[com]
    c=com;
    while (c in succ and pred[succ[c]]==c): 
      stream.append(succ[c]); c=succ[c]; 
    c=com;
    while (c in pred and succ[pred[c]]==c):
      stream.append(pred[c]); c=pred[c];
    return stream
  ##
  metapart={}
  metakompt=0;
  streams=[];
  statstream=[0 for kk in range(TT)]
  for com in LABEL:
    if LABEL[com] not in metapart:
      stream=buildstream(LABEL[com])
      streams.append(stream)
      for cc in stream: metapart[cc]=metakompt
      metakompt+=1
      statstream[len(stream)-1]+=1

  for kk in range(TT):
    if statstream[kk]>0: print(".... %d streams spread over %d periods" % (statstream[kk],kk+1))

  return metapart

## ##################################################
## ##################################################

def updategroups(myfile, newgroup):

  with open(myfile, 'r') as data_file: data = json.load(data_file)
  old_new=dict()  
  for n in data['nodes']: old_new[n['group']]=newgroup[n['id_top']]  
  #
  with open(myfile, 'r') as data_file: stuff = data_file.read()
  stuff=stuff.replace('"group":','"groupX":')
  #  
  for g in old_new:
    aa='"groupX":%d,' % g
    bb='"group":%d,' % old_new[g]
    stuff=stuff.replace(aa,bb)
  # 
  with open(myfile, 'w') as f_out: f_out.write('%s' % stuff)

  return

## ##################################################
## ##################################################
## ##################################################

def main():
# usage: parser.py [-h] [--version] -i DIR [-o DIR] [-v] [-e]
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
          default = "blah",
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

  """
  if (not os.path.exists(args.out_dir[0])):
    args.out_dir[0]=args.in_dir[0]
    print ("Error: Output directory does not exist: ", args.out_dir[0])
    exit()
  """

  if args.out_dir == 'blah':
    args.out_dir = args.in_dir

  if args.database not in ['wos','scopus']:
    print ("Error: database must be either 'wos' or 'scopus'")
    exit()

  ##      

  allinone(args.in_dir[0],args.out_dir[0],args.database,args.expert)

  return


    
## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################