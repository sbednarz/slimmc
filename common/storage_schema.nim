## Runtime constructors for Slimmc Storage schema records.
## These avoid expanding large `%*{...}` JSON macros in engine modules.

import std/json

proc storageRecord*(fields: openArray[(string, JsonNode)]): JsonNode =
  result = newJObject()
  for (key, value) in fields:
    result[key] = value
