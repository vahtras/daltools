test:
	python -m pytest -x 2>&1 | tee errors.err
debug:
	python -m pytest -x --pdb
coverage:
	python -m pytest tests --cov=daltools --cov-report=html --cov-report=term-missing
