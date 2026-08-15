# Root workflow for the slimmc family.
# Component-specific commands remain in homo/, copo/, and cli/.

PYTHON ?= python
NIM ?= nim
MPLCONFIGDIR ?= /tmp/pyslimmc-mpl
MPLBACKEND ?= Agg
PYTHON_RUN = env MPLBACKEND="$(MPLBACKEND)" MPLCONFIGDIR="$(MPLCONFIGDIR)" $(PYTHON)
PYTHON_CHECK = $(PYTHON) scripts/check_python_env.py

.PHONY: help info python-info check-versions check-release-config check-docs check-makefiles build build-engines build-all debug test-fast test test-devel test-full test-release test-pyslimmc test-pyslimmc-opt test-integration test-engines test-cli test-resolved-model \
	test-phase-a test-phase-b test-phase-c test-phase-d test-phase-e \
	test-validation test-depropagation test-terminal-microstructure test-homo-copo-equivalence \
	clean clean-generated

help:
	@echo "Slimmc build and test targets"
	@echo "============================="
	@echo ""
	@echo "Build:"
	@echo "  make build"
	@echo "      Build slimmc and slimmc-summary."
	@echo ""
	@echo "  make build-engines"
	@echo "      Build standalone homo and copo developer engines."
	@echo ""
	@echo "  make build-all"
	@echo "      Build user binaries and standalone developer engines."
	@echo ""
	@echo "  make debug"
	@echo "      Build debug variants."
	@echo ""
	@echo "  make info"
	@echo "      Show configured tools and component versions."
	@echo ""
	@echo "  make python-info"
	@echo "      Show the selected Python interpreter and dependency locations."
	@echo ""
	@echo "Tests:"
	@echo "  make test-fast"
	@echo "      Run fast pyslimmc tests only."
	@echo ""
	@echo "  make test"
	@echo "      Run unit tests plus 92 real CLI/engine/pyslimmc/opt integrations."
	@echo ""
	@echo "  make test-devel"
	@echo "      Run the complete development regression: phases A-E, technical"
	@echo "      validation, depropagation, terminal/penultimate/microstructure"
	@echo "      validation, and homo-copo equivalence."
	@echo ""
	@echo "  make test-full"
	@echo "      Run standard and complete development regression; never examples."
	@echo ""
	@echo "  make test-release"
	@echo "      Run release metadata, Git provenance, Makefile, and documentation gates."
	@echo ""
	@echo "Development subtests:"
	@echo "  make test-phase-a     Run phase A for homo and copo."
	@echo "  make test-phase-b     Run phase B for homo and copo."
	@echo "  make test-phase-c     Run phase C for homo and copo."
	@echo "  make test-phase-d     Run phase D for homo and copo."
	@echo "  make test-phase-e     Run phase E for homo and copo."
	@echo "  make test-validation  Run technical validation for homo and copo."
	@echo "  make test-depropagation"
	@echo "      Run detailed copo depropagation validation."
	@echo "  make test-terminal-microstructure"
	@echo "      Run copo terminal, penultimate, and microstructure validation."
	@echo "  make test-homo-copo-equivalence"
	@echo "      Compare homo with effective one-monomer copo."
	@echo ""
	@echo "Cleanup:"
	@echo "  make clean"
	@echo "      Remove compiled binaries, caches, and generated test files."
	@echo ""
	@echo "Information:"
	@echo "  make help"
	@echo "      Show this help."

info:
	@echo "PYTHON=$(PYTHON)"
	@(command -v $(PYTHON) >/dev/null 2>&1 && $(PYTHON) --version) || echo "Python not found"
	@echo "NIM=$(NIM)"
	@(command -v $(NIM) >/dev/null 2>&1 && $(NIM) --version | head -n 1) || echo "Nim not found"
	@echo "Slimmc=$$(cat VERSION)"
	@echo "pyslimmc=$$(sed -n 's/^__version__ = \"\([^\"]*\)\"/\1/p' pyslimmc/_version.py)"
	@echo "pyslimmc-opt=$$(sed -n 's/^__version__ = \"\([^\"]*\)\"/\1/p' pyslimmc_opt/__init__.py)"

check-versions:
	$(PYTHON) scripts/check_versions.py

check-release-config:
	$(PYTHON) scripts/check_release_config.py

check-docs:
	$(PYTHON) scripts/check_documentation.py

check-makefiles:
	$(PYTHON) scripts/check_makefile_references.py

python-info:
	$(PYTHON_CHECK) --require numpy matplotlib

build: cli-build summary-build

build-engines: homo-build copo-build

build-all: build build-engines

homo-build:
	$(MAKE) -C homo build NIM="$(NIM)" PYTHON="$(PYTHON)"

copo-build:
	$(MAKE) -C copo build NIM="$(NIM)" PYTHON="$(PYTHON)"

cli-build:
	$(MAKE) -C cli build NIM="$(NIM)" PYTHON="$(PYTHON)"

summary-build:
	$(MAKE) -C cli summary NIM="$(NIM)" PYTHON="$(PYTHON)"

debug:
	$(MAKE) -C homo debug NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(MAKE) -C copo debug NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(MAKE) -C cli debug NIM="$(NIM)" PYTHON="$(PYTHON)"


test-run-id:
	$(NIM) c -r --hints:off common/tests/test_run_id.nim
	$(NIM) c -r --hints:off common/tests/test_model_contract.nim

test-homo-results-v1:
	$(MAKE) -C homo build TARGET=slimmc-stage-h NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(PYTHON) common/tests/check_homo_storage_v1_stage_c.py
	$(PYTHON) common/tests/check_homo_storage_v1_actions.py
	$(PYTHON) common/tests/check_homo_storage_v1_finalization.py
	$(PYTHON) common/tests/check_homo_storage_v1_chains.py
	$(PYTHON) common/tests/check_homo_storage_v1_moments.py
	$(PYTHON) common/tests/check_homo_storage_v1_validation.py


test-resolved-model:
	$(MAKE) -C homo build TARGET=slimmc-resolved-test NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(MAKE) -C copo build NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(PYTHON) common/tests/check_resolved_model_metadata.py \
		--homo-engine homo/slimmc-resolved-test --copo-engine copo/slimmc-copo

test-pyslimmc:
	$(PYTHON_CHECK) --require numpy matplotlib pytest
	$(PYTHON) -m compileall -q pyslimmc
	$(PYTHON_RUN) -m pytest -q pyslimmc/tests

test-pyslimmc-opt:
	$(PYTHON_CHECK) --require numpy scipy sklearn pytest
	$(PYTHON) -m compileall -q pyslimmc_opt
	$(PYTHON_RUN) -m pytest -q pyslimmc_opt/tests

test-integration: build
	$(PYTHON_CHECK) --require numpy scipy sklearn matplotlib pytest
	$(PYTHON_RUN) -m pytest -q tests/integration

test-engines:
	$(MAKE) -C homo test NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(MAKE) -C copo test NIM="$(NIM)" PYTHON="$(PYTHON)"

test-cli:
	$(MAKE) -C cli test NIM="$(NIM)" PYTHON="$(PYTHON)"

# Shared phase entry points. These preserve the component-specific suites while
# giving developers one stable command at repository root.
test-phase-a: cli-build
	$(PYTHON_RUN) homo/tests/validation/phase_a/check_phase_a.py --engine bin/slimmc
	$(PYTHON_RUN) copo/tests/validation/phase_a/check_phase_a.py --engine bin/slimmc


test-phase-b: cli-build
	$(PYTHON_RUN) homo/tests/validation/phase_b/check_phase_b.py --engine bin/slimmc
	$(PYTHON_RUN) copo/tests/validation/phase_b/check_phase_b.py --engine bin/slimmc


test-phase-c: cli-build
	$(PYTHON_RUN) homo/tests/validation/phase_c/check_phase_c.py --engine bin/slimmc
	$(PYTHON_RUN) copo/tests/validation/phase_c/check_phase_c.py --engine bin/slimmc


test-phase-d: cli-build
	$(PYTHON_RUN) homo/tests/validation/phase_d/check_phase_d.py --engine bin/slimmc
	$(PYTHON_RUN) copo/tests/validation/phase_d/check_phase_d.py --engine bin/slimmc


test-phase-e: cli-build
	$(PYTHON_RUN) homo/tests/validation/phase_e/check_phase_e.py --engine bin/slimmc
	$(PYTHON_RUN) copo/tests/validation/phase_e/check_phase_e.py --engine bin/slimmc



# Copolymer-specific detailed suites with no direct homo counterpart.
test-depropagation: cli-build
	$(PYTHON_RUN) copo/tests/validation/depropagation/check_depropagation.py --engine bin/slimmc


test-terminal-microstructure: cli-build
	$(PYTHON_RUN) copo/tests/validation/terminal_microstructure/check_terminal_microstructure.py --engine bin/slimmc


# Cross-engine black-box equivalence. The standalone binaries are required so
# the same paired chemistry is forced through each engine independently.
test-homo-copo-equivalence:
	$(MAKE) -C homo build TARGET=slimmc-equivalence NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(MAKE) -C copo build NIM="$(NIM)" PYTHON="$(PYTHON)"
	@work_dir="$$(mktemp -d)"; \
	trap 'rm -rf "$$work_dir"' EXIT; \
	cp homo/slimmc-equivalence "$$work_dir/slimmc-homo"; \
	cp copo/slimmc-copo "$$work_dir/slimmc-copo"; \
	chmod +x "$$work_dir/slimmc-homo" "$$work_dir/slimmc-copo"; \
	$(PYTHON_RUN) common/tests/equivalence_homo_copo/check_equivalence.py \
		--homo-engine "$$work_dir/slimmc-homo" --copo-engine "$$work_dir/slimmc-copo"

# Fast feedback while editing Python API code.
test-fast: test-pyslimmc

# Standard local test suite.
test: check-versions check-release-config check-docs check-makefiles test-run-id test-pyslimmc test-pyslimmc-opt test-engines test-cli test-resolved-model test-integration
	@echo "slimmc family test: PASS"

# Technical validation for both engines.
test-validation: test-phase-a test-phase-b test-phase-c test-phase-d test-phase-e
	@echo "slimmc family test-validation: PASS"

# Complete development regression.
test-devel:
	$(MAKE) test-validation NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(MAKE) test-depropagation NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(MAKE) test-terminal-microstructure NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(MAKE) test-homo-copo-equivalence NIM="$(NIM)" PYTHON="$(PYTHON)"
	@echo "slimmc family test-devel: PASS"

# Complete local regression.
test-full:
	$(MAKE) test NIM="$(NIM)" PYTHON="$(PYTHON)"
	$(MAKE) test-devel NIM="$(NIM)" PYTHON="$(PYTHON)"
	@echo "slimmc family test-full: PASS"

# Release packaging gate. This intentionally avoids rerunning the complete
# numerical regression already owned by test-full.
test-release: build check-versions check-release-config check-makefiles
	$(PYTHON) scripts/check_documentation.py --engine ./bin/slimmc --require-engine
	./bin/slimmc --version
	$(PYTHON) scripts/check_build_provenance.py bin/slimmc --require-git
	@echo "slimmc family test-release: PASS"

clean:
	$(MAKE) -C homo clean PYTHON="$(PYTHON)"
	$(MAKE) -C copo clean PYTHON="$(PYTHON)"
	$(MAKE) -C cli clean PYTHON="$(PYTHON)"
	rm -f common/tests/test_run_id common/tests/test_run_id.exe
	rm -f common/tests/test_model_contract common/tests/test_model_contract.exe
	rm -f common/tests/test_results_writers common/tests/test_results_writers.exe
	rm -f common/tests/test_storage_manifest common/tests/test_storage_manifest.exe
	rm -f homo/slimmc-resolved-test homo/slimmc-resolved-test.exe
	rm -f homo/slimmc-equivalence homo/slimmc-equivalence.exe
	$(MAKE) clean-generated PYTHON="$(PYTHON)"

clean-generated:
	$(PYTHON) scripts/clean_generated.py
