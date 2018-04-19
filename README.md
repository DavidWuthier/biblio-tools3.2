
# BiblioTools

This is a fork of the BiblioTools/BiblioMaps project that can be found at http://www.sebastian-grauwin.com/bibliomaps.

It's meant to be run on Linux. It includes an example project with a make file to simplify the workflow.

## Dependencies

For Python:

```
pip install numpy
pip install argparse
pip install json
pip install itertools
pip install unidecode
pip install collections
pip install subprocesses
pip install networkx
```

For Python 3:

```
sudo apt-get install python3-numpy
sudo apt-get install python3-unidecode
sudo apt-get install python3-networkx
```

Apache2 server for BiblioMaps:

```
sudo apt-get install apache2
```

## Workflow

Add your .csv files in the "raw" folder, then

```
make
sudo make install
```

Open Firefox to the address "localhost", select "bibliomaps" and that's it!

## Official description

   BiblioTools is a set of python scripts that can be used to analyse bibligraphic data. essentially, the scripts take Web of Science or Scopus bibliographic data files as input and prepare formatted output files that can be used either in (the graph vizualisation tool) Gephi or in (the web interative visualition platform) BiblioMaps.

   BiblioMaps is a set of html / css / javascript files that can be used to deploy an interactive visualisation platform on your (local or remote) website.

   BiblioTools and BiblioMap will allow you to
   ** perform standard statistical (frequency) analysis of the bibliographic items (keywords, subjects, journals, authors, references, etc) of the publications within your corpus and visualize the results in a web interface
   ** prepare standard co-occurrence networks (co-authors, co-citations, heterogeneous networks, etc) to be visualized in Gephi or BiblioMaps
   ** prepare a Bibliographic Coupling analysis: construction of the network, detection of hierarchical topical clusters, interactive visualisation of the resulting clusters.

   For more details, see the tutorials, examples and references on my website: http://www.sebastian-grauwin.com/bibliomaps

   -----------------------------------------------------------------------------

   Author : Sebastian Grauwin (seb.grauwin@gmail.com)
   Version: November 2017

   -----------------------------------------------------------------------------

   BiblioTools is a freely available software: you can redistribute it and/or modify it as you want, but feedbacks or freelance jobs proposals are always appreciated! In particular, if you find a bug, please send me a bug report including if necessary the input files. Note that the python scripts are expecting a friendly use and do not make much verifications about the arguments you may provide.

   BiblioTools codes are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY, express or implied, including warranties of merchantability and fitness for a particular purpose.
   
   -----------------------------------------------------------------------------
