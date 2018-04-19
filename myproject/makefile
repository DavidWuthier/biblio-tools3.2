DB=scopus

all: description coupling cooc

parsing:
	mkdir -p processed
	python3 ../scripts/biblio_parser.py -i ./raw -o ./processed -d $(DB)

description: parsing
	python3 ../scripts/describe_corpus.py -i ./processed -v

coupling: parsing
	python3 ../scripts/biblio_coupling.py -i ./processed -v

cooc: parsing
	python3 ../scripts/cooc_graphs.py -i ./processed -v -hg

clean:
	rm -f ./scopus_sources_without_categories.dat
	rm -rf ./processed

install:
	rm -f /var/www/html/index.html
	cp -r ../html/bibliomaps /var/www/html
	mkdir -p /var/www/html/bibliomaps/data
	cp -r ./processed/freqs /var/www/html/bibliomaps/data
	cp ./processed/jsonfiles/* /var/www/html/bibliomaps/data

uninstall:
	rm -rf /var/www/html/bibliomaps
