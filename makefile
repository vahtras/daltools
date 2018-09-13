test:
	python -m pytest -x 2>&1 | tee errors.err
debug:
	python -m pytest -x --pdb
coverage:
	python -m pytest --cov=daltools --cov-report="html"
