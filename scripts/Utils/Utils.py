#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2012
   All rights reserved.
   BSD license.
"""

import os
import sys
import glob
import numpy
import argparse
import unidecode

def remove_accents(input_str):
    return unidecode.unidecode(input_str)

## ##################################################
## ##################################################
## ##################################################

class Biblio_line:
    def __init__(self):
        self.PT = "" ## Publication Type (J=Journal; B=Book; S=Series)
        self.AU = "" ## Authors
        self.BA = "" ## Book Authors
        self.BE = "" ## Book Editor
        self.GP = "" ## Book Group Authors
        self.AF = "" ## Author Full Name
        self.BF = "" ##         
        self.CA = "" ## Group Authors
        self.TI = "" ## Document Title
        self.SO = "" ## Publication Name
        self.SE = "" ## Book Series Title
        self.BS = "" ## Book Series Subtitle
        self.LA = "" ## Language
        self.DT = "" ## Document Type
        self.CT = "" ## Conference Title 
        self.CY = "" ## Conference Date 
        self.CL = "" ## Conference Location 
        self.SP = "" ## Conference Sponsors 
        self.HO = "" ## Conference Host
        self.DE = "" ## Author Keywords
        self.ID = "" ## Keywords Plus
        self.AB = "" ## Abstract
        self.C1 = "" ## Author Address
        self.RP = "" ## Reprint Address
        self.EM = "" ## E-mail Address
        self.RI = "" ##
        self.OI = "" ##
        self.FU = "" ## Funding Agency and Grant Number
        self.FX = "" ## Funding Text
        self.CR = "" ## Cited References
        self.NR = "" ## Cited Reference Count
        self.TC = "" ## Times Cited
        self.Z9 = "" ## 
        self.PU = "" ## Publisher
        self.PI = "" ## Publisher City
        self.PA = "" ## Publisher Address
        self.SN = "" ## ISSN
        self.EI = "" ##
        self.BN = "" ## ISBN
        self.J9 = "" ## 29-Character Source Abbreviation
        self.JI = "" ## ISO Source Abbreviation
        self.PD = "" ## Publication Date
        self.PY = 0 ## Year Published
        self.VL = "" ## Volume
        self.IS = "" ## Issue
        self.PN = "" ## Part Number
        self.SU = "" ## Supplement
        self.SI = "" ## Special Issue
        self.BP = "" ## Beginning Page
        self.EP = "" ## Ending Page
        self.AR = "" ## Article Number
        self.DI = "" ## Digital Object Identifier (DOI)
        self.D2 = "" ## 
        self.PG = "" ## Page Count
        self.WC = "" ## Web of Science Category
        self.SC = "" ## Subject Category
        self.GA = "" ## Document Delivery Number
        self.UT = "" ## Unique Article Identifier

    def parse_line(self, line, defCols, numCols, database):
        """
        parse a line of the WoS or Scopus output file  
        """

        if database =='wos':
            s = line.split("\t");
            if len(s)<=numCols+1 and len(s)> numCols-2:
                if(s[defCols['PT']]=='J'): self.PT = 'Journal' ## Publication Type (J=Journal; B=Book; S=Series)
                if(s[defCols['PT']]=='B'): self.PT = 'Book' 
                if(s[defCols['PT']]=='S'): self.PT = 'Series' 
                self.AU = s[defCols['AU']] ## Authors
                self.TI = s[defCols['TI']] ## Document Title
                self.SO = s[defCols['SO']] ## Publication Name
                self.SN = s[defCols['SN']] ## Source ISSN
                self.DT = s[defCols['DT']] ## Document Type
                self.DE = s[defCols['DE']] ## Author Keywords
                self.ID = s[defCols['ID']] ## Keywords Plus
                self.C1 = s[defCols['C1']] ## Author Address
                self.RP = s[defCols['RP']] ## Reprint Address
                self.CR = s[defCols['CR']] ## Cited References
                self.TC = s[defCols['TC']] ## Times Cited
                self.DI = s[defCols['DI']] ## Doi
                self.J9 = s[defCols['J9']] ## 29-Character Source Abbreviation
                if s[defCols['PY']].isdigit(): self.PY = int(s[defCols['PY']])
                else:               self.PY = 0  ## Year Published
                self.PD = s[defCols['PD']] ## Publication date
                self.VL = s[defCols['VL']] ## Volume
                self.IS = s[defCols['IS']] ## Issue
                self.BP = s[defCols['BP']] ## Beginning Page
                self.WC = s[defCols['WC']] ## Web of Science Category
                self.UT = s[defCols['UT']] ## Unique Article Identifier
                self.FX = s[defCols['FX']] ## Funding text
                self.AB = s[defCols['AB']] ## Abstract
                self.LA = s[defCols['LA']] ## Language
                #self.line = line

        elif database == 'scopus':
            aux=0
            s = []
            line=line.replace('""',"'")
            while aux < len(line) and len(s)<numCols-1:
              # we are at 'aux'. The next column can be delimited by "" or not
              if line[aux]=='"':
                aux2=line.find('",',aux)
                s.append(line[aux+1:aux2])
                aux=aux2+2
              else:
                aux2=line.find(',',aux)
                s.append(line[aux:aux2])
                aux=aux2+1
            s.append(line[aux:])

            if len(s)==numCols:
                self.AU = remove_accents(s[defCols['Authors']].replace(',',';').replace('.','')) ## Authors
                self.TI = s[defCols['Title']] ## Document Title
                self.SO = s[defCols['Source title']] ## Publication Name
                self.DT = s[defCols['Document Type']] ## Document Type
                self.DE = s[defCols['Author Keywords']] ## Author Keywords
                self.ID = s[defCols['Index Keywords']] ## Index Keywords
                self.C1 = remove_accents(s[defCols['Affiliations']]) ## Affiliations
                self.RP = remove_accents(s[defCols['Authors with affiliations']])  ## 
                self.CR = remove_accents(s[defCols['References']]) ## Cited References
                self.TC = s[defCols['Cited by']] ## Times Cited
                self.DI = s[defCols['DOI']] ## Doi
                self.J9 = s[defCols['Abbreviated Source Title']] ##  Source Abbreviation
                if s[defCols['Year']].isdigit(): self.PY = int(s[defCols['Year']])
                else:               self.PY = 0  ## Year Published
                self.VL = s[defCols['Volume']] ## Volume
                self.IS = s[defCols['Issue']] ## Issue
                self.BP = s[defCols['Page start']] ## Beginning Page
                self.UT = s[defCols['EID']] ## Unique Article Identifier
                self.AB = s[defCols['Abstract']] ## Abstract
                self.LA = s[defCols['Language of Original Document']] ## Language
                self.SN = s[defCols['ISSN']] ## Source ISSN


## ##################################################

def defColumns(line,database):

  # initialize
  defCols = dict();
  #
  if database == 'wos':
    Cols = ['PT', 'AU', 'TI', 'SO', 'SN', 'DT', 'DE', 'ID', 'AB', 'C1', 'RP', 'CR', 'TC', 'J9', 'PD', 'PY', 'LA', 'VL', 'IS', 'BP', 'WC', 'FX', 'DI', 'UT']
    foo = line.replace('\ufeff','').replace('\xef\xbb\xbf','').split('\t')
  elif database == 'scopus':
    Cols = ['Authors', 'Title', 'Source title', 'Document Type', 'Author Keywords', 'Index Keywords', 'Abstract', 'Affiliations', 'Authors with affiliations', 'References', 'Cited by', 'Abbreviated Source Title', 'Year', 'Volume', 'Issue', 'Page start', 'Language of Original Document', 'DOI', 'EID', 'ISSN']
    foo = line.replace('\xef\xbb\xbf','').replace('\ufeff','').split(',')

  # match columns number in "line"
  for i in range(len(foo)):
    if foo[i] in Cols: 
      defCols[foo[i]] = i
  numCols = len(foo)
  return (defCols, numCols)

## ##################################################

## ####ff=ope##############################################
## ##################################################
## ##################################################

class ArticleList:

    def __init__(self):
        self.articles      = []      # articles list
 
    def read_file(self,filename,database):

        articles_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename,'rU',encoding='utf-8') 
                #fd = open(filename,'rU',encoding="ISO-8859-1") 
            else:
                fd = sys.stdin
            # read
            aux = 0
            for line in fd.readlines():
                line = line.strip('\n') # removes \n
                if (line != ""):
                  if (aux == 1): # do not take 1st line into account! 
                    wline = Biblio_line()
                    wline.parse_line(line, defCols, numCols, database)
                    articles_list.append( wline )
                  if (aux == 0): # define columns thanks to 1st line
                    (defCols, numCols) = defColumns( line, database )
                    aux = 1
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.articles   = articles_list

## ##################################################
## ##################################################
## ##################################################

class Article:

    def __init__(self):
        self.id          = 0
        self.firstAU     = ""
        self.pubdate     = ""       
        self.year        = 0
        self.journal     = ""
        self.volume      = ""
        self.page        = ""
        self.doi         = ""
        self.pubtype     = ""
        self.doctype     = ""
        self.language    = ""
        self.times_cited = 0
        self.title       = ""
        self.uniqueID    = ""
 
        self.articles      = []      # liste d'articles


    def read_file(self,filename):
        """
        Lecture des articles
        """
        articles_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename, encoding='utf8')
            else:
                fd = sys.stdin
            # read
            aux = 0
            for line in fd.readlines():
                line = line.strip() # removes \n
                if (line != ""):
                    s = line.split("\t")
                    aline = Article()
                    aline.id = int(s[0])
                    if(len(s)>1): aline.firstAU = s[1]
                    if(len(s)>2): aline.year = int(s[2]) 
                    if(len(s)>3): aline.journal = s[3] 
                    if(len(s)>4): aline.volume = s[4] 
                    if(len(s)>5): aline.page = s[5] 
                    if(len(s)>6): aline.doi  = s[6]
                    if(len(s)>7): aline.pubtype = s[7]
                    if(len(s)>8): aline.doctype = s[8]
                    if(len(s)>9): aline.language = s[9]
                    if(len(s)>10): 
                        try: aline.times_cited = int('0'+s[10])
                        except: aline.times_cited = 0
                    if(len(s)>11): aline.title = s[11]
                    if(len(s)>12): aline.pubdate = s[12]
                    if(len(s)>13): aline.uniqueID = s[13]

                    if len(aline.journal) == 0:
                      aline.journal='%s - %s %d, %s ' % (aline.pubtype, aline.firstAU, aline.year, aline.title)
   
                    articles_list.append( aline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.articles   = articles_list

## ##################################################
## ##################################################
## ##################################################

class Author:

    def __init__(self):
        self.id     = 0
        self.rank   = 0       
        self.author = ""
 
        self.authors  = []      # liste 


    def read_file(self,filename):
        """
        Lecture des articles
        """
        alines_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename)
            else:
                fd = sys.stdin
            # read
            for line in fd.readlines():
                line = line.strip() # removes \n
                if (line != ""):
                    s = line.split("\t")
                    aline = Author()
                    aline.id = int(s[0])
                    aline.rank = int(s[1])
                    aline.author = s[2]  
                    alines_list.append( aline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.authors   = alines_list

## ##################################################
## ##################################################
## ##################################################

class Country:

    def __init__(self):
        self.id     = 0
        self.rank   = 0       
        self.country = ""
 
        self.countries  = []      # liste 


    def read_file(self,filename):
        """
        Lecture des articles
        """
        clines_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename)
            else:
                fd = sys.stdin
            # read
            for line in fd.readlines():
                line = line.strip() # removes \n
                if (line != ""):
                    s = line.split("\t")
                    cline = Country()
                    cline.id = int(s[0])
                    cline.rank = int(s[1])
                    cline.country = s[2].lower().capitalize()  
   
                    clines_list.append( cline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.countries   = clines_list

## ##################################################
## ##################################################
## ##################################################

class City:

    def __init__(self):
        self.id     = 0
        self.rank   = 0       
        self.city = ""
 
        self.cities  = []      # liste 


    def read_file(self,filename):
        """
        Lecture des articles
        """
        clines_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename)
            else:
                fd = sys.stdin
            # read
            for line in fd.readlines():
                line = line.strip() # removes \n
                if (line != ""):
                    s = line.split("\t")
                    cline = City()
                    cline.id = int(s[0])
                    cline.rank = int(s[1])
                    cline.city = s[2]  
   
                    clines_list.append( cline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.cities   = clines_list

## ##################################################
## ##################################################
## ##################################################

class Institution:

    def __init__(self):
        self.id     = 0
        self.rank   = 0       
        self.institution = ""
 
        self.institutions  = []      # liste 


    def read_file(self,filename):
        """
        Lecture des articles
        """
        ilines_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename)
            else:
                fd = sys.stdin
            # read
            for line in fd.readlines():
                line = line.strip() # removes \n
                if (line != ""):
                    s = line.split("\t")
                    iline = Institution()
                    if len(s)==3:
                      iline.id = int(s[0])
                      iline.rank = int(s[1])
                      iline.institution = s[2].upper()  
   
                      ilines_list.append( iline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.institutions   = ilines_list

## ##################################################
## ##################################################
## ##################################################

class Keyword:

    def __init__(self):
        self.id      = 0
        self.ktype   = ""       
        self.keyword = ""
 
        self.keywords  = []      # liste 


    def read_file(self,filename):
        """
        Lecture des articles
        """
        klines_list = []
        try:
            # open
            if filename != 'st.lower().capitalize()din':
                fd = open(filename)
            else:
                fd = sys.stdin
            # read
            for line in fd.readlines():
                line = line.strip() # removes \n
                if (line != ""):
                    s = line.split("\t")
                    kline = Keyword()
                    kline.id = int(s[0])
                    kline.ktype = s[1]
                    if len(s)>2:
                      kline.keyword = s[2].upper()  

                      klines_list.append( kline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.keywords   = klines_list

## ##################################################
## ##################################################
## ##################################################

class Ref:

    def __init__(self):
        self.id      = 0
        self.firstAU = ""       
        self.year    = 0
        self.journal = ""
        self.volume  = "0"
        self.page    = "0"
        self.DOI     = ""
        self.refs      = []      # liste de refs

    def parse_ref(self, ref, database):
        """
        parse a ref of the WoS or Scopus format  
        """
        if database == 'wos':
            s = ref.split(', ')
            if(len(s)>0):
                s[0]=s[0].replace('. ','').replace('.','')
                aux=s[0].split(' ')
                aux[-1]=aux[-1].replace('-','')
                foo =' '.join([elt.capitalize() for elt in aux[0:len(aux)-1]] + [aux[len(aux)-1]])
                foo = '-'.join([elt[0].upper()+elt[1:] for elt in foo.split('-') if len(elt)>0])
                self.firstAU = foo
            fooy=0
            if(len(s)>1): 
                if s[1].isdigit(): self.year = int(s[1])
                else:              self.year = 0; fooy=1
            if (foo.isdigit() and int(foo)>1500 and int(foo)<2100) :
            	self.firstAU = "[Anonymous]"
            	self.year = int(foo)
            	fooy=1
            if(len(s)>2-fooy): self.journal = s[2-fooy]
            for x in range(3-fooy,len(s)):
                if(s[x][0]=='V'): self.volume  = s[x].replace('V','')
                elif(s[x][0]=='P'): self.page  = s[x].replace('P','')
                elif(s[x][0]=='p'): self.page  = s[x].replace('p','')
                elif(s[x][0:3]=="DOI"): self.DOI = s[x].replace('DOI ','')

        if database == 'scopus':
            aux1=ref.find('(')
            aux2=ref.find(')',aux1)
            aux3=ref.find(',',aux2+2)
            if (aux1 >-1 and aux3==-1):aux3=len(ref)-1
            # sometimes there are strange "authors" / 
            while (aux3==aux2+1 and ref.find('(',aux3) > -1) or ((aux2 > -1) and (aux1 > -1) and ref[aux1+1:aux2].isdigit()==False ):   
              #print '%d %d %d' % (aux1,aux2,aux3)
              aux1=ref.find('(',aux2)
              aux2=ref.find(')',aux1)
              aux3=ref.find(',',aux2+2)
              if (aux1 >-1 and aux3==-1):aux3=len(ref)-1
            if (aux1 > -1) and ref[aux1+1:aux2].isdigit(): self.year=int(ref[aux1+1:aux2])
            if (aux2 > -1) and (aux3 > -1): self.journal=ref[aux2+2:aux3].upper()
            if (aux2 > -1) and (aux3 == -1): self.journal=ref[aux2+2:].upper()
            aux0=ref.find(',',ref.find(',')+1)
            if aux0 < aux1 and aux0 > -1:
              self.firstAU=ref[0:aux0].replace(',','').replace('.','')
            #
            aux = ref.find(', ,')
            # if this is present, ref is in 'book format': no volume / page; else it is in 'article format'
            if (aux == -1) and (aux3 > -1) and (aux3 < len(ref)-1):
                aux4=ref.find(',',aux3+1)
                if aux4 > -1:
                    foo=ref[aux3+2:aux4]
                    if ('pp' not in foo) and (foo.split(' ')[0].isdigit()): 
                        self.volume=foo.split(' ')[0]
                    if ('pp' in foo) and (foo.replace('pp. ','').replace('-',' ').split(' ')[0].isdigit()):
                        self.page=foo.replace('pp. ','').replace('-',' ').split(' ')[0]
                aux5=ref.find(',',aux4+1)
                if aux5 == -1: aux5=len(ref)-1 
                if aux4 > -1:
                    foo=ref[aux4+2:aux5]
                    if ('pp' in foo) and (foo.replace('pp. ','').replace('-',' ').split(' ')[0].isdigit()):
                        self.page=foo.replace('pp. ','').replace('-',' ').split(' ')[0]

    def read_file(self,filename):
        """
        Lecture des refs
        """
        refs_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename, encoding='utf8')
            else:
                fd = sys.stdin
            # read
            for line in fd.readlines():
                line = line.strip('\n') # removes \n
                if (line != ""):
                    s = line.split("\t")
                    refline = Ref()
                    refline.id = int(s[0])
                    refline.firstAU = s[1]
                    refline.year = int(s[2]) 
                    refline.journal = s[3] 
                    refline.volume = s[4] 
                    refline.page = s[5] 
   
                    refs_list.append( refline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.refs   = refs_list

## ##################################################
## ##################################################
## ##################################################

class Subject:

    def __init__(self):
        self.id      = 0
        self.subject = ""       
 
        self.subjects  = []      # liste 


    def read_file(self,filename):
        """
        Lecture des articles
        """
        slines_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename)
            else:
                fd = sys.stdin
            # read
            for line in fd.readlines():
                line = line.strip() # removes \n
                if (line != ""):
                    s = line.split("\t")
                    sline = Subject()
                    sline.id = int(s[0])
                    sline.subject = s[1] 
   
                    slines_list.append( sline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.subjects   = slines_list

## ##################################################
## ##################################################
## ##################################################

class Abstract:

    def __init__(self):
        self.id      = 0
        self.abstract = ""       
 
        self.abstracts  = []      # liste 


    def read_file(self,filename):
        #"
        #Lecture des abstract
        #"
        alines_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename)
            else:
                fd = sys.stdin
            # read
            for line in fd.readlines():
                line = line.strip() # removes \n
                if (line != ""):
                    s = line.split("\t")
                    aline = Abstract()
                    aline.id = int(s[0])
                    aline.abstract = s[1] 
   
                    alines_list.append( aline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.abstracts   = alines_list

## ##################################################
## ##################################################
## ##################################################

class Labo:

    def __init__(self):
        self.id   = 0      
        self.labo = ""
 
        self.labos  = []      # liste 


    def read_file(self,filename):
        """
        Lecture des labos
        """
        llines_list = []
        try:
            # open
            if filename != 'stdin':
                fd = open(filename)
            else:
                fd = sys.stdin
            # read
            for line in fd.readlines():
                line = line.strip() # removes \n
                if (line != ""):
                    s = line.split("\t")
                    lline = Labo()
                    if len(s)==2:
                      lline.id = int(s[0])
                      lline.labo = s[1]
   
                      llines_list.append( lline )
            # close  
            if filename != 'stdin':
                fd.close()
        except IOError:
            print ("file does not exist")

        self.labos   = llines_list

## ##################################################
## ##################################################
## ##################################################
    
## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################

