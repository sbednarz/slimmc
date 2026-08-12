## Cross-platform atomic publication of a completed sibling temporary file.

import std/os

when defined(windows):
  import std/widestrs
  import windows/winlean

proc atomicReplaceFile*(tmpPath, finalPath: string) =
  ## Replace `finalPath` atomically with the already closed `tmpPath`.
  ## Both paths must reside on the same filesystem; callers use a sibling temp.
  when defined(windows):
    const MovefileWriteThrough = 0x8'i32
    let ok = moveFileExW(
      newWideCString(tmpPath),
      newWideCString(finalPath),
      DWORD(MOVEFILE_REPLACE_EXISTING or MovefileWriteThrough)
    )
    if ok == 0:
      raiseOSError(osLastError(), "cannot atomically replace " & finalPath)
  else:
    # POSIX rename replaces an existing non-directory destination atomically.
    moveFile(tmpPath, finalPath)
