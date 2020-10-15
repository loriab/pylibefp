.DEFAULT_GOAL := all
isort = isort -rc pylibefp
black = black pylibefp
autoflake = autoflake -ir --remove-all-unused-imports --ignore-init-module-imports --remove-unused-variables qcelemental

.PHONY: format
format:
#	$(autoflake)
	$(isort)
	$(black)

.PHONY: lint
lint:
	$(isort) --check-only
	$(black) --check

.PHONY: mypy
mypy:
	mypy pylibefp

.PHONY: test
test:
	pytest -v --cov=pylibefp/

.PHONY: clean
clean:
	rm -rf `find . -name __pycache__`
	rm -f `find . -type f -name '*.py[co]' `
	rm -f `find . -type f -name '*~' `
	rm -f `find . -type f -name '.*~' `
	rm -rf .cache
	rm -rf .pytest_cache
	rm -rf .mypy_cache
	rm -rf htmlcov
	rm -rf *.egg-info
	rm -f .coverage
	rm -f .coverage.*
	rm -rf build
	rm -rf dist
	rm -f qcelemental/*.c qcelemental/*.so
