from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def text(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def main() -> None:
    homo_parser = text("homo/src/slimmc_parser.nim")
    copo_parser = text("copo/src/copo_parser.nim")
    homo_kmc = text("homo/src/slimmc_kmc.nim")
    copo_kmc = text("copo/src/copo_kmc.nim")

    for source in (homo_parser, copo_parser):
        assert 'of "kmc_volume"' in source
        assert 'of "init_volume"' in source
        assert 'tokens[2] == "repeat"' in source or 'toks[2] == "repeat"' in source
        assert 'tokens[4] == "every"' in source or 'toks[4] == "every"' in source
        assert 'feed action syntax: feed NAME VOLUME_ML' in source

    for source in (homo_kmc, copo_kmc):
        assert 'of eaFeed:' in source
        assert 'newPhysical = m.currentVolumeMl + doseMl' in source
        assert 'newV = oldV * newPhysical / m.currentVolumeMl' in source
        assert 'remaining == 0' in source

    models = list(ROOT.rglob("*.model"))
    assert models
    assert not any("param volume " in p.read_text(encoding="utf-8", errors="ignore") for p in models)

    print("semibatch core source checks: PASS")


if __name__ == "__main__":
    main()
