## Shared public constants and small types for Slimmc Storage.

const
  StorageName* = "slimmc-storage"
  StorageFormatVersion* = "1.2.0"
  CanonicalByteOrder* = "little"
  CanonicalNpyVersion* = "1.0"
  AvogadroConstantMolInv* = 6.02214076e23

type
  RunStatus* = enum
    rsRunning,
    rsCompleted,
    rsFailed,
    rsInterrupted

proc `$`*(status: RunStatus): string =
  case status
  of rsRunning: "running"
  of rsCompleted: "completed"
  of rsFailed: "failed"
  of rsInterrupted: "interrupted"
