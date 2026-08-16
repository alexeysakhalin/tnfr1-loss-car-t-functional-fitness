# Independent validation

`recalculate_checkmate_survival.py` independently reconstructs the prespecified
nivolumab-only CheckMate survival models from the checksum-verified Braun
supplementary workbook. It writes aggregate diagnostics and does not stage
patient-level data for version control.

The primary analysis uses 181 RNA-profiled nivolumab-treated tumors from
CM-009, CM-010 and CM-025. Event coding is retained exactly as supplied by the
source workbook.
