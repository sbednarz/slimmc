SLIMMC ?= slimmc
PYTHON ?= python3
.PHONY: all clean
all:
	PYTHONPATH=../.. $(PYTHON) run_regressions.py --slimmc $(SLIMMC)
clean:
	find . -type d -name results -prune -exec rm -rf {} +
