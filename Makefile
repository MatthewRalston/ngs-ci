docs:
	sphinx-apidoc -f -o docs ngsci
	cd docs && make html && cd ..

.PHONY: all docs



