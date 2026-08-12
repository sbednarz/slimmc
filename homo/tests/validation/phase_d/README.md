# Homo chemistry phase D

Phase D completes the remaining homopolymer chemistry coverage without changing
the model language or public API.

- H15: macro initiation
- H16: transfer to a hydrogen donor / chain-transfer agent
- H17: termination by radical capper/inhibitor (`term_x`)
- H18: elementary unimolecular, unlike-bimolecular, and same-species bimolecular reactions
- H19: elementary-reaction efficiency, including productive/nonproductive counters

Each mechanism has fast Nim tests and full `.model -> engine -> storage -> pyslimmc`
integration checks. Run with `make test-phase-d` from `homo/`.
