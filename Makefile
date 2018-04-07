setup: build-all 

build-all: build-cripsr-sites build-offtarget

build-cripsr-sites:
	cd crispr_sites && $(MAKE)

build-offtarget:
	cd offtarget && go build -o offtarget

# Re-enable after merge into master, which contains ash
#install-ash:
#	cd ash && $(MAKE)

.PHONY: test

test:
	cd crispr_sites && make tests

