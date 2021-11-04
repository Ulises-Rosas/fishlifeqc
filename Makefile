SHELL := /bin/bash

BIN_DIR = ext_bin

NAME    = `python3 setup.py --name`
VERSION = `python3 setup.py --version`
FULL_NAME  = $(NAME)-$(VERSION)
COMPRESSED = $(FULL_NAME)'.tar.gz'

all: dist_dir update_bins upload

dist_dir:
	python3 setup.py sdist bdist_wheel

update_bins:
	tar -xvzf dist/$(COMPRESSED)
	rm -rf $(FULL_NAME)/$(NAME).egg-info
	rm -rf $(FULL_NAME)/ext_bin
	cp -rf ext_bin $(FULL_NAME)
	tar -cvzf $(COMPRESSED) $(FULL_NAME)
	rm -rf dist
	rm -rf build
	rm -rf $(NAME).egg-info


upload:
	twine upload $(COMPRESSED)
	rm $(COMPRESSED)


