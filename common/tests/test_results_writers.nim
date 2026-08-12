import std/[json, os, strutils]
import ../[result_json, npy_writer, results_types]

let root = getCurrentDir() / "common" / "tests" / "generated_results_writer"
if dirExists(root): removeDir(root)
createDir(root)

let metadata = %*{
  "run_id": "run_000001",
  "storage": StorageName,
  "storage_format_version": StorageFormatVersion,
  "run_status": $rsRunning,
  "seed": "18446744073709551615"
}
writePrettyJsonAtomic(root / "run_metadata.json", metadata)
let metadataText = readFile(root / "run_metadata.json")
doAssert metadataText.endsWith("\n")
doAssert metadataText.contains("\n  \"run_id\"")
doAssert parseJson(metadataText)["run_status"].getStr == "running"
doAssert not fileExists(root / "run_metadata.json.tmp")
metadata["run_status"] = %"completed"
writePrettyJsonAtomic(root / "run_metadata.json", metadata)
doAssert parseJson(readFile(root / "run_metadata.json"))["run_status"].getStr == "completed"

let jsonlRecords = [
  %*{"record_type": "schema_header", "schema_version": "1.0.0"},
  %*{"record_type": "table", "name": "snapshots"}
]
writeJsonlAtomic(root / "schema.jsonl", jsonlRecords)
let jsonlLines = readFile(root / "schema.jsonl").strip().splitLines()
doAssert jsonlLines.len == 2
for line in jsonlLines:
  doAssert line.len > 0
  discard parseJson(line)

writeNpyUint64(root / "uint64.npy", [0'u64, 1'u64, high(uint64)])
writeNpyUint32(root / "uint32.npy", [0'u32, 1'u32, high(uint32)])
writeNpyFloat64(root / "float64.npy", [0.0, -1.25, 1.0e-18])
writeNpyBool(root / "bool.npy", [false, true, true])
writeNpyUint64(root / "empty_uint64.npy", newSeq[uint64]())
writeNpyFloat64(root / "empty_float64.npy", newSeq[float64]())

for name in ["uint64.npy", "uint32.npy", "float64.npy", "bool.npy", "empty_uint64.npy", "empty_float64.npy"]:
  doAssert fileExists(root / name)
  doAssert not fileExists(root / (name & ".tmp"))

echo root
