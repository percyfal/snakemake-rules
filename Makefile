.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help
define BROWSER_PYSCRIPT
import os, webbrowser, sys
try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts


clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

lint: ## check style with flake8
	flake8 snakemake_rules tests

test: ## run tests quickly with the default Python
	py.test


test-all: ## run tests on every Python version with tox
	tox

coverage: ## check code coverage quickly with the default Python
	coverage run --source snakemake_rules -m pytest

		coverage report -m
		coverage html
		$(BROWSER) htmlcov/index.html

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/snakemake_rules.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ snakemake_rules
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

release: clean ## package and upload a release
	python setup.py sdist upload
	python setup.py bdist_wheel upload

conda: clean ## package and upload a conda release
	conda build conda
	anaconda upload $(shell conda build conda --output)

id:=$(shell docker images | grep -v REPOSITORY | head -1 | awk '{print $$3}')
DOCKERUSER:=percyfal
DOCKERTAG:=snakemake-rules
VERSION:=$(shell python setup.py --version)

docker: clean ## Make docker image
	docker build -t ${DOCKERTAG} .

dockerupload: clean ## Make docker image and upload
	docker tag ${id} ${DOCKERUSER}/${DOCKERTAG}:latest
	docker login -u ${DOCKERUSER}
	docker push ${DOCKERUSER}/${DOCKERTAG}

dockerrelease: clean ## Make docker image and upload, tagging with latest version
	docker build -t ${DOCKERTAG} .
	docker tag ${id} ${DOCKERUSER}/${DOCKERTAG}:latest
	docker tag ${id} ${DOCKERUSER}/${DOCKERTAG}:$(VERSION)
	docker login -u ${DOCKERUSER}
	docker push ${DOCKERUSER}/${DOCKERTAG}

dist: clean ## builds source and wheel package
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean ## install the package to the active Python's site-packages
	python setup.py install
