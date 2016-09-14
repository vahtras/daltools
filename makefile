test:
	python -m nose -x
debug:
	python -m nose -x --pdb
coverage:
	python -m nose --with-coverage --cover-package daltools --cover-html
