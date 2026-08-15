# slimmc cli dispatcher (homo/copo)
#
# Single dispatcher binary: links both engines directly in-process (no
# sibling-binary exec at runtime). Architecture and routing contracts are
# documented in docs/development/ARCHITECTURE.md and docs/reference/CLI.md.

NIM ?= nim
TARGET_DIR ?= ../bin
SLIMMC_TARGET ?= $(TARGET_DIR)/slimmc
PYTHON ?= python
NIMFLAGS_COMMON ?= --hints:off

# Build provenance. Values are compiled into the executable; runtime never
# consults the user's current Git repository. Callers may override these.
BUILD_GIT_COMMIT ?= $(shell git rev-parse HEAD 2>/dev/null)
BUILD_GIT_TAG ?= $(shell git describe --tags --exact-match HEAD 2>/dev/null)
BUILD_GIT_DIRTY ?= $(shell if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then test -z "$$(git status --porcelain)" && echo false || echo true; else echo unknown; fi)
BUILD_TIMESTAMP_UTC ?= $(shell date -u +%Y-%m-%dT%H:%M:%SZ)
PROVENANCE_DEFINES = $(if $(BUILD_GIT_COMMIT),-d:BuildGitCommit=$(BUILD_GIT_COMMIT),) $(if $(BUILD_GIT_TAG),-d:BuildGitTag=$(BUILD_GIT_TAG),) -d:BuildGitDirty=$(BUILD_GIT_DIRTY) -d:BuildTimestampUtc=$(BUILD_TIMESTAMP_UTC)
NIMCACHE ?= nimcache
CLI_VERSION ?= $(shell sed -n '1p' ../VERSION)
VERSION_DEFINES = -d:CliVersion=$(CLI_VERSION) -d:SlimmcVersion=$(CLI_VERSION)

.PHONY: all help build summary debug test check-all clean

# Keep a bare `make` equivalent to `make build` without GNU .DEFAULT_GOAL.
all: build

help:
	@echo "slimmc cli dispatcher Makefile targets"
	@echo "  make build       build release executable ($(SLIMMC_TARGET))"
	@echo "  make summary     build standalone slimmc-summary"
	@echo "  make debug       build debug executable"
	@echo "  make test        run dispatch-logic unit tests (no engine build needed)"
	@echo "  make check-all   build + test + end-to-end smoke test"
	@echo "  make clean       remove build artifacts"

build:
	mkdir -p $(TARGET_DIR)
	$(NIM) c $(NIMFLAGS_COMMON) $(VERSION_DEFINES) $(PROVENANCE_DEFINES) --nimcache:$(NIMCACHE)/slimmc --path:../homo/src --path:../copo/src -d:release --opt:speed \
		-o:$(SLIMMC_TARGET) slimmc_cli.nim
	@test -f "$(SLIMMC_TARGET).exe" || chmod +x "$(SLIMMC_TARGET)"

summary:
	mkdir -p $(TARGET_DIR)
	$(NIM) c $(NIMFLAGS_COMMON) -d:AppVersion=$(CLI_VERSION) $(PROVENANCE_DEFINES) --nimcache:$(NIMCACHE)/summary -d:release -o:$(TARGET_DIR)/slimmc-summary slimmc_summary.nim
	@test -f "$(TARGET_DIR)/slimmc-summary.exe" || chmod +x "$(TARGET_DIR)/slimmc-summary"

debug:
	mkdir -p $(TARGET_DIR)
	$(NIM) c $(NIMFLAGS_COMMON) $(VERSION_DEFINES) $(PROVENANCE_DEFINES) --nimcache:$(NIMCACHE)/debug --path:../homo/src --path:../copo/src -d:debug --debugger:native \
		-o:$(SLIMMC_TARGET)-debug slimmc_cli.nim

test:
	$(NIM) c --hints:off -r -o:tests/run_dispatch_tests tests/test_dispatch.nim

check-all: build test
	$(PYTHON) tests/test_dispatch_smoke.py $(SLIMMC_TARGET)

clean:
	rm -rf $(NIMCACHE)
	rm -f $(SLIMMC_TARGET) $(SLIMMC_TARGET).exe $(SLIMMC_TARGET)-debug $(SLIMMC_TARGET)-debug.exe
	rm -f slimmc slimmc.exe slimmc-debug slimmc-debug.exe
	rm -f $(TARGET_DIR)/slimmc-summary $(TARGET_DIR)/slimmc-summary.exe
	rm -f tests/run_dispatch_tests tests/run_dispatch_tests.exe
