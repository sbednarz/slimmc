import std/[os, strutils]
import ../storage_manifest

let root = getTempDir() / "slimmc_storage_manifest_test"
if dirExists(root): removeDir(root)
createDir(root)
createDir(root / "state")
createDir(root / ".work")
writeFile(root / "input.model", "param seed 1\n")
writeFile(root / "state" / "t.npy", "abc")
writeFile(root / "run_metadata.json", "ignored")
writeFile(root / "RESULTS_COMPLETE", "ignored")
writeFile(root / ".work" / "temporary", "ignored")

let first = writeStorageManifest(root)
doAssert first.fileCount == 2
doAssert fileExists(root / "checksums.sha256")
doAssert first.text.contains("input.model")
doAssert first.text.contains("state/t.npy")
doAssert not first.text.contains("run_metadata.json")
doAssert not first.text.contains("RESULTS_COMPLETE")
doAssert not first.text.contains(".work/")

let second = writeStorageManifest(root)
doAssert second.hash == first.hash
doAssert second.text == first.text

writeFile(root / "state" / "t.npy", "abcd")
let third = writeStorageManifest(root)
doAssert third.hash != first.hash

removeDir(root)
echo "storage manifest: PASS"
