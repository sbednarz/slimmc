# Standalone copo-engine developer workflow.

NIM ?= nim
PYTHON ?= python
APP ?= slimmc-copo
SRC ?= src/slimmc_copo.nim
NIMFLAGS_COMMON ?= --hints:off

# Build provenance. Values are compiled into the executable; runtime never
# consults the user's current Git repository. Callers may override these.
BUILD_GIT_COMMIT ?= $(shell git rev-parse HEAD 2>/dev/null)
BUILD_GIT_TAG ?= $(shell git describe --tags --exact-match HEAD 2>/dev/null)
BUILD_GIT_DIRTY ?= $(shell if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then test -z "$$(git status --porcelain)" && echo false || echo true; else echo unknown; fi)
BUILD_TIMESTAMP_UTC ?= $(shell date -u +%Y-%m-%dT%H:%M:%SZ)
PROVENANCE_DEFINES = $(if $(BUILD_GIT_COMMIT),-d:BuildGitCommit=$(BUILD_GIT_COMMIT),) $(if $(BUILD_GIT_TAG),-d:BuildGitTag=$(BUILD_GIT_TAG),) -d:BuildGitDirty=$(BUILD_GIT_DIRTY) -d:BuildTimestampUtc=$(BUILD_TIMESTAMP_UTC)
NIM_PATH ?= --path:src
VERSION := $(shell sed -n '1p' ../VERSION)
VERSION_DEFINE = -d:SlimmcVersion=$(VERSION)

.PHONY: all help info build debug test clean

all: build

help:
	@echo "Slimmc copo engine $(VERSION)"
	@echo "  make build  build the standalone release engine"
	@echo "  make debug  build the standalone debug engine"
	@echo "  make test   run copo Nim unit tests"
	@echo "  make info   show configured tools and version"
	@echo "  make clean  remove generated component artifacts"

info:
	@echo "VERSION=$(VERSION)"
	@echo "NIM=$(NIM)"
	@(command -v $(NIM) >/dev/null 2>&1 && $(NIM) --version | head -n 1) || echo "Nim not found"
	@echo "PYTHON=$(PYTHON)"
	@(command -v $(PYTHON) >/dev/null 2>&1 && $(PYTHON) --version) || echo "Python not found"

build:
	$(NIM) c $(NIMFLAGS_COMMON) $(NIM_PATH) $(VERSION_DEFINE) $(PROVENANCE_DEFINES) -d:release --opt:speed -o:$(APP) $(SRC)
	@test -f "$(APP).exe" || chmod +x "$(APP)"

debug:
	$(NIM) c $(NIMFLAGS_COMMON) $(NIM_PATH) $(VERSION_DEFINE) $(PROVENANCE_DEFINES) --debugger:native -o:$(APP)-debug $(SRC)
	@test -f "$(APP)-debug.exe" || chmod +x "$(APP)-debug"

test:
	$(NIM) c -r $(NIMFLAGS_COMMON) $(NIM_PATH) tests/run_all_tests.nim

clean:
	rm -rf nimcache src/nimcache tests/nimcache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -f $(APP) $(APP).exe $(APP)-debug $(APP)-debug.exe
	rm -f tests/run_all_tests tests/run_all_tests.exe
