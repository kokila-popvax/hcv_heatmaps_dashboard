"""
hcv_data.py
-----------
Pure data-processing logic for the HCV cross-neutralization dashboard.
Depends only on pandas / numpy / re / ast  -> fully unit-testable without Streamlit.

This module mirrors the V15/V16 Colab logic:
  * filter to HCV experiments, drop control groups + antibody/plasma samples
  * parse the Dilution / Avg_Neut_percent_corrected arrays
  * extract % neutralization at a chosen dilution
  * compute log10(IC50), flagging NN (No Neutralization)
  * assign Prime / Boost1 / Boost2 buckets from a hardcoded BUCKET_MAP
  * classify each construct into one of the 10 subgroups via SUBGROUP_RULES

>>> THE TWO BLOCKS YOU MUST REPLACE WITH YOUR AUTHORITATIVE V16 VERSIONS <<<
    1. BUCKET_MAP      (search for: ===== PASTE BUCKET_MAP =====)
    2. SUBGROUP_RULES  (search for: ===== PASTE SUBGROUP_RULES =====)
The starter contents below are reconstructed fragments so the app runs out of
the box; anything not covered falls through to 'Unknown' bucket / 'Uncategorized'
subgroup and the UI reports how many rows that affected.
"""

from __future__ import annotations

import ast
import re
import numpy as np
import pandas as pd


# ============================================================
# CONSTANTS  (same as V15/V16)
# ============================================================

HCV_EXPERIMENTS = {"PVXE141", "PVXE142", "PVXE167A", "PVXE167B", "PVXE174", "PVXE181"}
DROP_GROUPS = {"GNG", "G1"}
DROP_SAMPLE_TYPES = {"antibody", "ab", "mab", "plasma"}

UNCATEGORIZED = "10. Uncategorized"
BUCKETS = ["Prime", "Boost1", "Boost2"]


# ============================================================
# ===== PASTE BUCKET_MAP =====
# Keys are (Experiment, PSVX_number_as_str, Day_as_int) -> bucket name.
# Replace the starter block below with your full V16 BUCKET_MAP.
# ============================================================

# ============================================================
# ===== BUCKET_MAP (from V15 script) =====
# Keys are (Experiment, PSVX_number_as_str, Day_as_int) -> bucket name.
# ============================================================

BUCKET_MAP: dict[tuple[str, str, int], str] = {}


def _add(exp: str, psvx: str, day: int, bucket: str) -> None:
    BUCKET_MAP[(str(exp), str(psvx), int(day))] = bucket


# ------ PVXE141 ------
for day in [21, 28]:
    _add("PVXE141", "179", day, "Prime")
    _add("PVXE141", "180", day, "Prime")
_add("PVXE141", "179", 112, "Boost1")
_add("PVXE141", "180", 112, "Boost1")
_add("PVXE141", "180", 119, "Boost1")
_add("PVXE141", "180", 126, "Boost1")

# ------ PVXE142 ------
_add("PVXE142", "179", 28, "Prime")
_add("PVXE142", "179", 112, "Boost1")
_add("PVXE142", "180", 28, "Prime")
_add("PVXE142", "180", 112, "Boost1")
_add("PVXE142", "180", 119, "Boost1")
_add("PVXE142", "180", 126, "Boost1")

# ------ PVXE167A ------
_add("PVXE167A", "230", 28,  "Prime")
_add("PVXE167A", "242", 28,  "Prime")
_add("PVXE167A", "265", 28,  "Prime")
_add("PVXE167A", "225", 49,  "Boost1")
_add("PVXE167A", "228", 56,  "Boost1")
_add("PVXE167A", "229", 49,  "Boost1")
_add("PVXE167A", "230", 56,  "Boost1")
_add("PVXE167A", "253", 28,  "Prime")
_add("PVXE167A", "253", 56,  "Boost1")
_add("PVXE167A", "253", 98,  "Boost2")
_add("PVXE167A", "260", 56,  "Boost1")
_add("PVXE167A", "265", 49,  "Boost1")
_add("PVXE167A", "235", 98,  "Boost2")
_add("PVXE167A", "236", 98,  "Boost2")
_add("PVXE167A", "246", 98,  "Boost2")
_add("PVXE167A", "250", 91,  "Boost2")
_add("PVXE167A", "253", 91,  "Boost2")
_add("PVXE167A", "254", 28,  "Prime")
_add("PVXE167A", "254", 56,  "Boost1")
_add("PVXE167A", "254", 98,  "Boost2")
_add("PVXE167A", "265", 98,  "Boost2")

# ------ PVXE167B ------
_add("PVXE167B", "242", 28,  "Prime")
_add("PVXE167B", "228", 56,  "Boost1")
_add("PVXE167B", "242", 56,  "Boost1")
_add("PVXE167B", "253", 56,  "Boost1")
_add("PVXE167B", "242", 98,  "Boost2")
_add("PVXE167B", "246", 98,  "Boost2")
_add("PVXE167B", "250", 91,  "Boost2")
_add("PVXE167B", "254", 98,  "Boost2")

# ------ PVXE181 ------
_add("PVXE181", "258", 21,  "Prime")
_add("PVXE181", "262", 21,  "Prime")
_add("PVXE181", "272", 14,  "Prime")
_add("PVXE181", "262", 35,  "Boost1")
_add("PVXE181", "214", 35,  "Boost1")
_add("PVXE181", "259", 35,  "Boost1")
_add("PVXE181", "260", 35,  "Boost1")
_add("PVXE181", "270", 70,  "Boost2")
_add("PVXE181", "273", 70,  "Boost2")
_add("PVXE181", "271", 70,  "Boost2")
_add("PVXE181", "272", 49,  "Boost2")

# ------ PVXE174 ------
_add("PVXE174", "257", 28,  "Prime")
_add("PVXE174", "260", 28,  "Prime")
_add("PVXE174", "254", 56,  "Boost1")
_add("PVXE174", "257", 140, "Boost2")

# ============================================================
# CUSTOM CONSTRUCT LABELS FOR PVXE181 G2, G5, G6
# These groups received cumulative injections. The label for each
# bucket reflects ALL constructs present at that timepoint.
# ============================================================

_G2_P  = "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72"
_G2_B1 = "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154"
_G2_B2 = "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154\n+ 1a138 + 1b34 + UKNP5.2.1 + UKNP1.11.6 + 1b25"

_G5_P  = "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=1 ILP"
_G5_B1 = "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=1 ILP\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154 + HVR HPF Reactivity=1 ILP"
_G5_B2 = "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=1 ILP\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154 + HVR HPF Reactivity=1 ILP\n+ 1a138 + 1b34 + UKNP5.2.1 + UKNP1.11.6 + 1b25 + HVR HPF Reactivity=1 ILP"

_G6_P  = "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=3 ILP"
_G6_B1 = "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=3 ILP\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154 + HVR HPF Reactivity=3 ILP"
_G6_B2 = "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=3 ILP\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154 + HVR HPF Reactivity=3 ILP\n+ 1a138 + 1b34 + UKNP5.2.1 + UKNP1.11.6 + 1b25 + HVR HPF Reactivity=3 ILP"

CUSTOM_LABELS = {
    ("PVXE181", "G2", "Prime"):  _G2_P,
    ("PVXE181", "G2", "Boost1"): _G2_B1,
    ("PVXE181", "G2", "Boost2"): _G2_B2,
    ("PVXE181", "G5", "Prime"):  _G5_P,
    ("PVXE181", "G5", "Boost1"): _G5_B1,
    ("PVXE181", "G5", "Boost2"): _G5_B2,
    ("PVXE181", "G6", "Prime"):  _G6_P,
    ("PVXE181", "G6", "Boost1"): _G6_B1,
    ("PVXE181", "G6", "Boost2"): _G6_B2,
}

# ============================================================
# ===== END BUCKET_MAP =====
# ============================================================

# ============================================================
# SUBGROUP CLASSIFICATION
# Source of truth: "Organized List" sheet of HCV_Constructs_List.xlsx
#
# Strategy:
#   1. EXACT_LOOKUP  — normalised exact match (case-insensitive, whitespace
#                      collapsed).  Covers every named construct in the sheet.
#   2. KEYWORD_RULES — ordered regex fallback for constructs whose sheet name
#                      differs slightly from the IC50 sheet (typos, extra
#                      spaces, short bare-strain names like "UKNP1.16.3").
#   3. Everything else → "9. Misc. Individual Components"
# ============================================================

SG1  = "1. Soluble Envelope Proteins"
SG2  = "2. E1E2 Heterodimers & Combinations"
SG3  = "3. De novo Designed Constructs"
SG4  = "4. Multi-Component Prime-Boost"
SG5  = "5. HVR Immunogens"
SG6  = "6. E1/E2 TMD Variants"
SG7  = "7. Complex HVR Arrays + Adjuvants"
SG8  = "8. Misc. Individual Components"
SG9 = "9. Uncategorized"

SUBGROUP_ORDER = [SG1, SG2, SG3, SG4, SG5, SG6, SG7, SG8, SG9]

# ── Exact-name lookup (all names taken verbatim from the Excel sheet) ────────
# Key is normalised: lower-cased + whitespace collapsed to single space + stripped.
# Value is the subgroup string.

_EXACT_MAP_RAW = {
    # 1. Soluble Envelope Proteins
    "H77C_sE2_Ferritin":                                            SG1,
    "H77C_sE1_Gaussia soluble E1 non-VLP":                         SG1,
    "H77C_sE2_Gaussia: soluble E2 non-VLP":                        SG1,
    "AMS0232_sE2_Gaussia_nonVLP":                                   SG1,
    "H77C_sE1_Gaussia_FibAlpha":                                    SG1,
    "H77C_sE2_Gaussia_FibAlpha":                                    SG1,
    "H77C_sE2_Gaussia_E2p_FibAlpha":                                SG1,
    "H77C_sE2F442NYT_Gaussia_ferritin":                             SG1,
    "H77C_sE2_ΔHVR_HPF":                                           SG1,
    "H77C_sE2_FYQ442NYT_HPF":                                       SG1,
    "H77C_sE1_Gaussia_ferritin_FibAlpha":                           SG1,
    "H77C_sE2_Gaussia_ferritin_FibAlpha":                           SG1,

    # 2. E1E2 Heterodimers & Combinations
    "H77C_sE1E2.LZ_Ferritin":                                       SG2,
    "H77C_sE1E2GS3_Gaussia_ferritin":                               SG2,
    "H77C_sE1_Gaussia soluble E1 non-VLP + H77C_sE2_Gaussia: soluble E2 non-VLP mix": SG2,
    "H77C_sE1E2_Gaussia_ferritin s:soluble":                        SG2,
    "LNP mix of: H77C_sE1_Gaussia_FibAlpha + H77C_sE2_Gaussia_FibAlpha": SG2,
    "H77C_sE1E2_Gaussia_furinP2A_FibAlpha":                         SG2,
    "H77C_sE1.Ferritin||sE2.Ferritin_FurinP2A_FibAlpha":            SG2,
    "LNP mix of: H77C_sE1_Gaussia_ferritin_FibAlpha + H77C_sE2_Gaussia_ferritin_FibAlpha": SG2,
    "mRNA mix of: H77C_sE1_Gaussia_ferritin_FibAlpha + H77C_sE2_Gaussia_ferritin_FibAlpha (2ug dose)": SG2,
    "mRNA mix of: H77C_sE1_Gaussia_ferritin_FibAlpha + H77C_sE2_Gaussia_ferritin_FibAlpha (5ug dose)": SG2,
    "H77C_sE1.LZ.ferritin||E2.LZ.ferritin_Gaussia_furinP2A_FibAlpha": SG2,
    "H77C_sE1.LZ||sE2.LZ.ferritin_Gaussia_furinP2A_FibAlpha":      SG2,
    "H77_sE1_LIF_sE2":                                              SG2,
    "1a154 (H77_w7) E1_TMD||E2 _TMD IVT construct":                 SG6,
    "1a154 E1_TMD||E2_TMD IVT construct":                           SG6,
    "1a154_E1E2_TMD_IVT":                                            SG6,
    "UKNP3.1.2 E1_TMD||E2 _TMD IVT construct":                      SG6,
    "UKNP3.1.2 E1_TMD||E2_TMD IVT construct":                       SG6,
    "UKNP3.1.2_E1E2_TMD_IVT":                                       SG6,
    "UKNP1.11.6 E1_TMD||E2 _TMD IVT construct":                     SG6,
    "UKNP1.11.6 E1_TMD||E2_TMD IVT construct":                      SG6,
    "UKNP1.11.6_E1E2_TMD_IVT":                                      SG6,
    "1b58 E1_TMD||E2 _TMD IVT construct":                           SG6,
    "1b58 E1_TMD||E2_TMD IVT construct":                            SG6,
    "1b58_E1E2_TMD_IVT":                                             SG6,
    "UKNP1.18.1 E1_TMD||E2 _TMD IVT construct":                     SG6,
    "UKNP1.18.1 E1_TMD||E2_TMD IVT construct":                      SG6,
    "UKNP1.18.1_E1E2_TMD_IVT":                                      SG6,
    "1b34 E1_TMD||E2 _TMD IVT construct":                           SG6,
    "1b34 E1_TMD||E2_TMD IVT construct":                            SG6,
    "1b34_E1E2_TMD_IVT":                                             SG6,
    "1a72 E1_TMD||E2 _TMD IVT construct":                            SG6,
    "1a72 E1_TMD||E2_TMD IVT construct":                             SG6,
    "1a72_E1E2_TMD_IVT":                                             SG6,
    "UKNP1.10.1_E1E2_TMD_IVT":                                      SG6,
    "UKNP1.9.1_E1E2_TMD_IVT":                                       SG6,
    "1a123_E1E2_TMD_IVT":                                            SG6,
    "UKNP1.16.3 E1_TMD||E2 _TMD IVT construct":                     SG6,
    "UKNP1.16.3 E1_TMD||E2_TMD IVT construct":                      SG6,
    "UKNP1.16.3_E1E2_TMD_IVT":                                      SG6,
    "1b25 E1_TMD||E2 _TMD IVT construct":                            SG6,
    "1b25 E1_TMD||E2_TMD IVT construct":                             SG6,
    "1b25_E1E2_TMD_IVT":                                             SG6,
    "UKNP4.2.2 E1_TMD||E2 _TMD IVT construct":                      SG6,
    "UKNP4.2.2 E1_TMD||E2_TMD IVT construct":                       SG6,
    "UKNP4.2.2_E1E2_TMD_IVT":                                       SG6,
    "UKNP5.2.1 E1_TMD||E2 _TMD IVT construct":                      SG6,
    "UKNP5.2.1 E1_TMD||E2_TMD IVT construct":                       SG6,
    "UKNP5.2.1_E1E2_TMD_IVT":                                       SG6,
    "1a138_E1E2_TMD_IVT":                                            SG6,
    # Bare strain names (appear in IC50 sheet without "E1_TMD||E2_TMD" suffix)
    "UKNP3.1.2":   SG2,
    "UKNP1.18.1":  SG2,
    "UKNP4.2.2":   SG2,
    "1b58":        SG2,
    "1a72":        SG2,
    "UKNP1.16.3":  SG2,
    "UKNP1.10.1":  SG2,
    "UKNP1.9.1":   SG2,
    "1a123":       SG2,
    "1a154":       SG2,
    "1a138":       SG2,
    "1b34":        SG2,
    "UKNP5.2.1":   SG2,
    "UKNP1.11.6":  SG2,
    "1b25":        SG2,

    # 3. De novo Designed Constructs
    "4n0y_C1_372_dldesign_0_cycle1_af2pred_6_strict_gaussia_ferritin IGH526": SG3,
    "4n0y_C1_72_dldesign_0_cycle1_af2pred_2_strict_gaussia_ferritin IGH526":  SG3,
    "fib_a_UTR_gaussia_HCV_E1_S1_ferritin_HBA_UTR from Corriea et al":        SG3,
    "19B3_GL_6BZV_AS412_89_9_af2out_gaussia_ferritin 19B3GL":                 SG3,
    "19B3_GL_6BZV_AS412_63_8_af2out_gaussia_ferritin 19B3GL":                 SG3,
    "fib_a_UTR_gaussia_HCV_E2_S2_1_ferritin_HBA_UTR from Corriea et al":      SG3,
    "8w0w_C1_219_hcab64_4_rank_6_Gaussia_ferritin":                            SG3,
    "8w0w_C2_213_hcab64_3_rank_1_ferritinVLP":                                 SG3,
    "8w0w_C1_471_hcab64_3_rank_3_ferritinVLP":                                 SG3,
    "PopVax AR4A 1234 – design_1-183_10_ROG_RFD (Run 1) with ferritin":        SG3,
    "PopVax AR4A 1234 – design_2-136_9_ROG_RFD (Run 1) with ferritin":         SG3,
    "PopVax AR4A 1234 – design_2-353_4_ROG_hotspot_RFD (Run 2) with ferritin": SG3,
    "PopVax AR4A 1234 – design_2-43_4_woROG_RFD (Run 1 Additional) with ferritin": SG3,
    "Bruno's E1 scaffolded on HbAg":                     SG3,
    "Bruno's E2 scaffolded on HbAg":                     SG3,
    "PopVax's E2 scaffolded on HPF":                     SG3,

    # 4. Multi-Component Prime-Boost — new labels use \n between rounds (not |)
    "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72":                                          SG4,  # G2 Prime
    "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154": SG4,  # G2 Boost1
    "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154\n+ 1a138 + 1b34 + UKNP5.2.1 + UKNP1.11.6 + 1b25": SG4,  # G2 Boost2
    "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=1 ILP":               SG4,  # G5 Prime
    "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=1 ILP\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154 + HVR HPF Reactivity=1 ILP": SG4,  # G5 Boost1
    "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=1 ILP\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154 + HVR HPF Reactivity=1 ILP\n+ 1a138 + 1b34 + UKNP5.2.1 + UKNP1.11.6 + 1b25 + HVR HPF Reactivity=1 ILP": SG4,  # G5 Boost2
    "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=3 ILP":               SG4,  # G6 Prime
    "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=3 ILP\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154 + HVR HPF Reactivity=3 ILP": SG4,  # G6 Boost1
    "UKNP3.1.2 + UKNP1.18.1 + UKNP4.2.2 + 1b58 + 1a72 + HVR HPF Reactivity=3 ILP\n+ UKNP1.16.3 + UKNP1.10.1 + UKNP1.9.1 + 1a123 + 1a154 + HVR HPF Reactivity=3 ILP\n+ 1a138 + 1b34 + UKNP5.2.1 + UKNP1.11.6 + 1b25 + HVR HPF Reactivity=3 ILP": SG4,  # G6 Boost2

    # 5. HVR Immunogens
    "HVR HPF Reactivity=1 ILP":                          SG5,
    "HVR HPF Reactivty=3 ILP mixture":                   SG5,   # note typo "Reactivty"
    "HVR HPF Reactivity=3 ILP":                          SG5,   # alternate spelling
    "HVR_wILP1_HPF":                                     SG5,
    "HVR_wILP3_HPF LNPmix":                              SG5,
    "HVR Most reactive five No ferritin":                SG5,
    "HVR Most reactive five Ferritin":                   SG5,
    "HVR Most reactive five 8mer Ferritin":              SG5,
    "HVR Most reactive five 8 mer no Ferritin":          SG5,
    "HVR Most reactive Single 8 mer Ferritin":           SG5,
    "HVR Most reactive Single 8mer + PADRE Ferritin":    SG5,
    "HVR Most reactive Single Ferritin":                 SG5,
    "HVR Most reactive Single PADRE Ferritin":           SG5,
    "HVR_wILP1_HPF+1a154_E1E2_TMD_IVT LNP mix":         SG5,
    "HVR_wILP3_HPF +1a154_E1E2_TMD_IVT LNP mix":        SG5,
    "HVR_wILP1_HPF+H77C_sE2_ΔHVR_HPF LNP mix":          SG5,
    "HVR_wILP3_HPF+H77C_sE2_ΔHVR_HPF LNP mix":          SG5,
    "HVR_wILP1_HPF+H77C_sE2_FYQ442NYT_HPF":              SG5,
    "HVR_wILP3_HPF+H77C_sE2_FYQ442NYT_HPF":              SG5,
    "H77mbE1E2_native_codon optimized + HVRs_wILP3max_21aa_HPF_G4Ssep": SG5,
    "H77mbE1E2_native_codon optimized + HVR_wILP3_HPF (a and b)":       SG5,

    # 6. E1/E2 TMD Variants
    "H77mbE2_Codon optimized":                           SG6,
    "H77_sE2-HATMD_Codon Optimized":                     SG6,
    "H77mbE1E2_native_codon optimized":                  SG6,
    "H77_sE1_HATMD-sE2_HATMD_codon optimized":           SG6,
    "H77mbE1E2_native_FYQ442NYT_codon optimized":        SG6,
    "Drew's glycosylated H77 sE2 on Ferritin":           SG6,
    "Mansun's no HVR no stalk H77 sE2 on Ferritin":      SG6,

    # 7. Complex HVR Arrays + Adjuvants
    "HVRs_wILP3max_21aa_HPF_G4Ssep":                     SG7,
    "HVRs_wILP3max_21aa_HATMD_G4Ssep":                   SG7,
    "HVRs_Feld_HATMD_G4Ssep":                            SG7,
    "HVRs_Feld_HPF_G4Ssep":                              SG7,
    "HVRs_wILP3max_1G4S_H77mbE1_ΔΗVR_E2_native":        SG7,
    "HVRs_wILP3max_HATMD_1G4S_H77mbE1_ΔΗVR_E2_native":  SG7,
    "HVRs_wILP3max_1G4S_H77mbE1_ΔΗVR_E2_native+Adju-Phos":       SG7,
    "HVRs_wILP3max_HATMD_1G4S_H77mbE1_ΔΗVR_E2_native+Adju-Phos": SG7,
    "HVRs_Feld_HPF_G4Ssep+C,I,I-FA":                    SG7,
    "HVRs_wILP3max_21aa_HPF_G4Ssep+C,I,I-FA":           SG7,
    "HVRs_Feld_HPF_G4Ssep+Adju-phos":                   SG7,
    "HVRs_Feld_HATMD_G4Ssep+C,I,I-FA":                  SG7,
    "HVRs_wILP3max_21aa_HPF_G4Ssep+Adju-phos":          SG7,
    "HVRs_wILP3max_21aa_HATMD_G4Ssep+C,I,I-FA":         SG7,
    "HVRs_wILP3max_21aa_HATMD_G4Ssep+Adju-Phos":        SG7,
    "GaussiaSP-HVRs_top16_mostreactive_21aa_sep1G4S-2G4S-HPF+C,I,I-FA":  SG7,
    "GaussiaSP-HVRs_top16_mostreactive_21aa_sep1G4S-2G4S-HPF+Adju-phos": SG7,
    "GaussiaSP-HVRs_top16_mostreactive_21aa_sep1G4S-HATMD-C,I,I-FA":     SG7,
    "GaussiaSP-HVRs_top16_mostreactive_21aa_sep1G4S-HATMD+Adju-Phos":    SG7,
}

# Normalise keys: lower + collapse whitespace + remove space-before-underscore
# e.g. "E2 _TMD" → "E2_TMD" so Excel and IC50-sheet spellings match
def _norm(s):
    s = re.sub(r'\s+', ' ', str(s).strip()).lower()
    s = re.sub(r'\s+_', '_', s)   # "E2 _TMD" → "E2_TMD"
    s = re.sub(r'_\s+', '_', s)   # "E2_ TMD" → "E2_TMD" (defensive)
    return s

EXACT_LOOKUP = {_norm(k): v for k, v in _EXACT_MAP_RAW.items()}

# ── Keyword fallback rules (ordered, first match wins) ───────────────────────
# Used only when exact lookup fails.  Handles:
#   • PVXE181 cumulative labels (contain " | ")
#   • Minor spelling variants not listed above
#   • Constructs described differently in the IC50 sheet vs the Excel sheet
KEYWORD_RULES = [
    # 4. Multi-Component Prime-Boost — newline separator is now the marker
    # (pipe "|" removed from labels; "\n" separates rounds in new format)
    (r'\n.*UKNP',                                        SG4),
    # Catch single-round Prime strings (no newline yet)
    (r'(?:Prime|Boost[12]).*UKNP3\.1\.2',               SG4),
    (r'UKNP3\.1\.2.*UKNP1\.18\.1.*UKNP4\.2\.2',         SG4),

    # 7. Complex HVR Arrays — HVRs_ (plural) variants not listed above
    (r'HVRs_.*[Aa]dju',                                  SG7),
    (r'HVRs_.*[Cc],\s*[Ii],\s*[Ii]',                    SG7),
    (r'HVRs_.*H77mb',                                    SG7),
    (r'GaussiaSP.*HVRs',                                 SG7),
    (r'HVRs_(?:wILP|Feld)',                              SG7),

    # 5. HVR Immunogens — HVR (singular) variants
    (r'HVR[\s_](?:HPF|wILP|Most)',                       SG5),
    (r'H77mbE1E2_native.*HVR',                           SG5),

    # 3. De novo Designed
    (r'dldesign|af2pred|af2out|hcab64',                  SG3),
    (r'PopVax\s+AR4A',                                   SG3),
    (r'19B3_GL',                                         SG3),
    (r'fib_a_UTR_gaussia_HCV',                           SG3),

    # 6. E1/E2 TMD Variants
    (r'H77mbE[12]|H77_s[Ee][12][-_]HATMD|H77mbE1E2',    SG6),
    (r"Drew.s glycosylated|Mansun.s no HVR",             SG6),
    (r'FYQ442NYT|ΔHVR|sE2F442NYT',                       SG6),

    # 6. E1/E2 TMD Variants — individual strain IVT constructs (moved from SG2)
    (r'E1_TMD.*E2.*TMD|E1E2.*TMD|_TMD_IVT',             SG6),
    # 2. E1E2 Heterodimers — mixtures and non-TMD heterodimers
    (r'sE1E2|sE1\.LZ|E1\.LZ\.ferritin',                  SG2),
    (r'(?:LNP|mRNA)\s+mix\s+of:',                        SG2),
    (r'H77_sE1_LIF_sE2',                                  SG2),
    # Bare short strain codes
    (r'^(?:UKNP[\d.]+|1[ab]\d+)$',                        SG2),

    # 1. Soluble Envelope Proteins
    (r'sE[12]_Gaussia|sE2_Ferritin|AMS0232_sE2',         SG1),
    (r'sE2_ΔHVR_HPF|sE2_FYQ442NYT_HPF|sE2F442NYT',       SG1),
]

_compiled_kw = [(re.compile(pat, re.IGNORECASE), sg) for pat, sg in KEYWORD_RULES]

def classify_subgroup(construct_desc):
    s = str(construct_desc).strip()
    # 1. Exact match (normalised)
    hit = EXACT_LOOKUP.get(_norm(s))
    if hit:
        return hit
    # 2. Keyword fallback
    for pattern, subgroup in _compiled_kw:
        if pattern.search(s):
            return subgroup
    return SG9   # "9. Misc. Individual Components" — not Uncategorized


# Optional: (Experiment, Group, Bucket) -> cumulative label override (e.g. PVXE181 SG4).
# Replace with your V16 CUSTOM_LABELS if you want the cumulative regimen labels.
# Optional: (Experiment, Group, Bucket) -> cumulative label override (e.g. PVXE181 SG4).
# From V15: PVXE181 G2, G5, G6 received cumulative injections. The label reflects ALL constructs
# present in the animal at that timepoint.
CUSTOM_LABELS: dict[tuple[str, str, str], str] = {
    # G2 — no HVR component
    ("PVXE181", "G2", "Prime"):  _G2_P,
    ("PVXE181", "G2", "Boost1"): _G2_B1,
    ("PVXE181", "G2", "Boost2"): _G2_B2,
    # G5 — HVR HPF Reactivity=1 ILP added at every stage
    ("PVXE181", "G5", "Prime"):  _G5_P,
    ("PVXE181", "G5", "Boost1"): _G5_B1,
    ("PVXE181", "G5", "Boost2"): _G5_B2,
    # G6 — HVR HPF Reactivity=3 ILP added at every stage
    ("PVXE181", "G6", "Prime"):  _G6_P,
    ("PVXE181", "G6", "Boost1"): _G6_B1,
    ("PVXE181", "G6", "Boost2"): _G6_B2,
}
# ============================================================
# ===== END SUBGROUP_RULES =====
# ============================================================


# Optional: PSV strain -> genotype, used only to group/label columns.
# Fill in your 15-strain panel (4x1a, 3x1b, 8 rare) for genotype grouping.
PSV_GENOTYPE: dict[str, str] = {}


# ============================================================
# HELPERS
# ============================================================


def find_col(df: pd.DataFrame, *patterns: str) -> str | None:
    """Return the first column whose name matches any of the regex patterns."""
    for p in patterns:
        for c in df.columns:
            if re.search(p, str(c), re.IGNORECASE):
                return c
    return None


def safe_parse_array(s) -> list[float]:
    """Parse a stringified list like '[10.0, np.float64(30.0), ...]' -> [floats]."""
    if isinstance(s, (list, tuple)):
        try:
            return [float(x) for x in s]
        except Exception:
            return []
    s = str(s).strip()
    if not s or s.lower() in ("nan", "none"):
        return []
    s = re.sub(r"np\.float64\(([^)]+)\)", r"\1", s)
    try:
        return [float(x) for x in ast.literal_eval(s)]
    except Exception:
        return []


def get_neut_at_dilution(dil_str, neut_str, target: float, tol: float = 0.01):
    """% neutralization at `target` dilution, or None if absent/mismatched."""
    dilutions = safe_parse_array(dil_str)
    neuts = safe_parse_array(neut_str)
    if not dilutions or not neuts or len(dilutions) != len(neuts):
        return None
    for d, n in zip(dilutions, neuts):
        if abs(d - target) < tol:
            return n
    return None


def _is_nn(raw_ic50) -> bool:
    """True if the IC50 cell denotes No Neutralization (NN, blank, or 0)."""
    s = str(raw_ic50).strip().lower()
    if s in ("nn", "no neutralization", "no neut"):
        return True
    if s in ("", "nan", "none", "-", "na", "n/a"):
        return True
    try:
        return float(s.replace(",", "")) == 0.0
    except Exception:
        return False


def _to_float(raw) -> float:
    """Parse a numeric cell, tolerating thousands separators (e.g. '2,000.0')."""
    try:
        return float(str(raw).strip().replace(",", ""))
    except Exception:
        return np.nan


# ============================================================
# PREPARE: raw sheet -> tidy standardized frame
# ============================================================

STANDARD_COLS = [
    "Experiment", "Group", "PSV", "Day", "PSVX_num", "Construct_Description",
    "Bucket_Type", "Subgroup", "IC50", "Is_NN", "Log10_IC50",
    "DilutionArray", "NeutArray",
]


def build_construct_lookup(df_const: pd.DataFrame) -> dict[tuple[str, str], str]:
    """
    Build a {(psvx_bare, group): construct_description} lookup from the
    HCV Constructs List sheet (mirrors the notebook's resolve_construct logic).
    """
    lookup: dict[tuple[str, str], str] = {}
    if df_const is None or df_const.empty:
        return lookup
    # filter to HCV rows if a Program column exists
    prog_col = next((c for c in df_const.columns
                     if re.search(r"^program$", c, re.I)), None)
    df_hcv = df_const[df_const[prog_col].astype(str).str.strip().str.upper() == "HCV"].copy() \
             if prog_col else df_const.copy()
    psvx_col = next((c for c in df_hcv.columns
                     if re.search(r"psvx.?id|psvx.?no|^psvx$", c, re.I)), None)
    grp_col  = next((c for c in df_hcv.columns
                     if re.search(r"^group$", c, re.I)), None)
    desc_col = next((c for c in df_hcv.columns
                     if re.search(r"construct.?desc", c, re.I)), None)
    if not all([psvx_col, grp_col, desc_col]):
        return lookup
    for _, row in df_hcv.iterrows():
        psvx_bare = (str(row[psvx_col]).strip()
                     .upper().replace("PSVX-", "").replace("PSVX", "").strip())
        psvx_bare = re.sub(r"\s*\([CRrc]\)\s*", "", psvx_bare).strip()
        grp = str(row[grp_col]).strip()
        desc = str(row[desc_col]).strip()
        if psvx_bare and desc and desc.lower() not in ("nan", "none", ""):
            lookup[(psvx_bare, grp)] = desc
            # also store without group so we can fallback
            lookup.setdefault((psvx_bare, ""), desc)
    return lookup


def prepare_dataframe(df_ic50: pd.DataFrame,
                      construct_lookup: dict | None = None,
                      corrected_ic50: bool = True) -> tuple[pd.DataFrame, dict]:
    """
    Turn the raw IC50 sheet into a tidy frame with STANDARD_COLS.
    construct_lookup: optional {(psvx_bare, group) -> description} from the
    HCV Constructs List sheet (build with build_construct_lookup()).
    Returns (tidy_df, info) where info carries diagnostics + detected column names.
    """
    info: dict = {}
    df = df_ic50.copy()
    df.columns = [str(c).strip() for c in df.columns]

    C_SAMPLE = find_col(df, r"sample.?type")
    C_GROUP = find_col(df, r"^group$", r"\bgroup\b")
    C_PSVX = find_col(df, r"psvx.?no", r"psvx_no", r"^psvx")
    C_PSV = find_col(df, r"^psv$", r"\bpseudovirus\b")
    C_DAY = find_col(df, r"^day$", r"\bday\b")
    C_EXP = find_col(df, r"^experiment$", r"\bexperiment\b", r"pvxe")
    C_DESC = find_col(df, r"construct.?description", r"\bconstruct\b", r"immunogen")
    C_DIL = find_col(df, r"dilution")
    # Pick corrected or raw columns based on caller preference
    if corrected_ic50:
        C_IC50  = find_col(df, r"^ic50_corrected$",            r"ic50.?corrected", r"corrected.?ic50")
        C_NEUT  = find_col(df, r"^avg_neut_percent_corrected$", r"avg.?neut.?percent.?corrected")
    else:
        C_IC50  = find_col(df, r"^ic50$")
        C_NEUT  = find_col(df, r"^avg_neut_percent$")
    # Fall back if preferred column is absent
    if C_IC50 is None:
        C_IC50 = find_col(df, r"ic50")
    if C_NEUT is None:
        C_NEUT = find_col(df, r"avg.?neut.?percent", r"neutrali")

    info["columns"] = dict(sample=C_SAMPLE, group=C_GROUP, psvx=C_PSVX, psv=C_PSV,
                           day=C_DAY, experiment=C_EXP, desc=C_DESC,
                           dilution=C_DIL, neut=C_NEUT, ic50=C_IC50)
    info["n_raw"] = len(df)

    # --- filters (identical to V15/V16) ---
    if C_SAMPLE:
        mask = df[C_SAMPLE].astype(str).str.strip().str.lower().isin(DROP_SAMPLE_TYPES)
        df = df[~mask]
    if C_GROUP:
        df = df[~df[C_GROUP].astype(str).str.strip().isin(DROP_GROUPS)]
    if C_EXP:
        df = df[df[C_EXP].astype(str).str.strip().isin(HCV_EXPERIMENTS)]
    info["n_after_filter"] = len(df)

    out = pd.DataFrame(index=df.index)
    out["Experiment"] = df[C_EXP].astype(str).str.strip() if C_EXP else ""
    out["Group"] = df[C_GROUP].astype(str).str.strip() if C_GROUP else ""
    out["PSV"] = df[C_PSV].astype(str).str.strip() if C_PSV else ""
    out["Day"] = pd.to_numeric(df[C_DAY], errors="coerce") if C_DAY else np.nan
    # Start with C_DESC fallback; will be overridden below if construct_lookup provided
    out["Construct_Description"] = (df[C_DESC].astype(str).str.strip()
                                    if C_DESC else "")

    if C_PSVX:
        out["PSVX_num"] = (df[C_PSVX].astype(str).str.strip()
                           .str.replace(r"(?i)psvx-?", "", regex=True)
                           .str.replace(r"\s*\([CRrc]\)\s*", "", regex=True)
                           .str.strip())
    else:
        out["PSVX_num"] = ""

    # Resolve construct descriptions from the construct lookup sheet (mirrors notebook)
    if construct_lookup:
        def _resolve(r):
            psvx = str(r["PSVX_num"]).strip().upper()
            grp  = str(r["Group"]).strip()
            desc = construct_lookup.get((psvx, grp)) or construct_lookup.get((psvx, ""))
            if desc:
                return desc
            # fall back to whatever is already in C_DESC
            return r["Construct_Description"]
        out["Construct_Description"] = out.apply(_resolve, axis=1)

    # arrays + % neut
    if C_DIL:
        out["DilutionArray"] = df[C_DIL].apply(safe_parse_array)
    else:
        out["DilutionArray"] = [[] for _ in range(len(df))]
    if C_NEUT:
        out["NeutArray"] = df[C_NEUT].apply(safe_parse_array)
    else:
        out["NeutArray"] = [[] for _ in range(len(df))]

    # IC50 + log10 + NN
    if C_IC50:
        out["Is_NN"] = df[C_IC50].apply(_is_nn)
        out["IC50"] = df[C_IC50].apply(_to_float)
    else:
        out["Is_NN"] = False
        out["IC50"] = np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        log10 = np.log10(out["IC50"].where(out["IC50"] > 0))
    out["Log10_IC50"] = log10

    # bucket + subgroup
    def _bucket(r):
        try:
            key = (str(r["Experiment"]), str(r["PSVX_num"]), int(r["Day"]))
        except (ValueError, TypeError):
            return "Unknown"
        return BUCKET_MAP.get(key, "Unknown")

    out["Bucket_Type"] = out.apply(_bucket, axis=1)
    out["Subgroup"] = out["Construct_Description"].apply(classify_subgroup)

    # apply cumulative-label overrides if provided
    if CUSTOM_LABELS:
        def _relabel(r):
            return CUSTOM_LABELS.get(
                (r["Experiment"], r["Group"], r["Bucket_Type"]),
                r["Construct_Description"],
            )
        out["Construct_Description"] = out.apply(_relabel, axis=1)
        # Re-classify subgroup now that construct descriptions may have changed
        out["Subgroup"] = out["Construct_Description"].apply(classify_subgroup)

    info["n_unknown_bucket"] = int((out["Bucket_Type"] == "Unknown").sum())
    info["n_uncategorized"] = int((out["Subgroup"] == UNCATEGORIZED).sum())
    info["all_dilutions"] = sorted({d for arr in out["DilutionArray"] for d in arr})
    info["experiments"] = sorted(out["Experiment"].unique().tolist())
    info["groups"] = sorted(g for g in out["Group"].unique().tolist() if g)
    info["psvs"] = sorted(out["PSV"].unique().tolist())
    info["subgroups_present"] = [s for s in SUBGROUP_ORDER
                                 if s in set(out["Subgroup"].unique())]

    return out.reset_index(drop=True), info


# ============================================================
# COMPUTE VIEW: tidy frame -> per (bucket, construct, PSV) cell table
# ============================================================

def compute_view(df: pd.DataFrame, metric: str,
                 dilution: float | str | None = None) -> pd.DataFrame:
    """
    Aggregate to one row per (Bucket_Type, Construct_Description, PSV).
    metric == 'pct_neut' -> value = mean % neut at `dilution`
                            dilution may be a float or the special string 'all'
                            (in 'all' mode, value = mean across all dilution points)
    metric == 'log10_ic50' -> value = mean log10(IC50) over non-NN replicates
    Output columns: Bucket_Type, Construct_Description, PSV, value, n_tested, all_nn
    """
    d = df[df["Bucket_Type"].isin(BUCKETS)].copy()

    if metric == "pct_neut":
        if dilution is None:
            raise ValueError("dilution is required for metric='pct_neut'")

        if dilution == "all":
            # Explode every (dilution, neut) pair into separate rows then average
            rows = []
            for _, r in d.iterrows():
                dils = safe_parse_array(r["DilutionArray"])
                neuts = safe_parse_array(r["NeutArray"])
                for dv, nv in zip(dils, neuts):
                    rows.append({
                        "Bucket_Type": r["Bucket_Type"],
                        "Construct_Description": r["Construct_Description"],
                        "PSV": r["PSV"],
                        "cell_value": float(nv),
                    })
            if not rows:
                return pd.DataFrame(columns=["Bucket_Type", "Construct_Description",
                                             "PSV", "value", "n_tested", "all_nn"])
            exploded = pd.DataFrame(rows)
            grp = (exploded.groupby(["Bucket_Type", "Construct_Description", "PSV"],
                                    as_index=False)
                           .agg(value=("cell_value", "mean"),
                                n_tested=("cell_value", "size")))
        else:
            d["cell_value"] = [
                get_neut_at_dilution(dil, neut, float(dilution))
                for dil, neut in zip(d["DilutionArray"], d["NeutArray"])
            ]
            d = d[d["cell_value"].notna()]
            grp = (d.groupby(["Bucket_Type", "Construct_Description", "PSV"], as_index=False)
                     .agg(value=("cell_value", "mean"),
                          n_tested=("cell_value", "size")))
        grp["all_nn"] = False
        return grp

    # log10_ic50
    grp = (d.groupby(["Bucket_Type", "Construct_Description", "PSV"])
             .agg(n_tested=("Is_NN", "size"),
                  n_nn=("Is_NN", "sum"),
                  value=("Log10_IC50", "mean"))
             .reset_index())
    grp["all_nn"] = grp["n_nn"] == grp["n_tested"]
    grp.loc[grp["all_nn"], "value"] = np.nan
    return grp[["Bucket_Type", "Construct_Description", "PSV", "value", "n_tested", "all_nn"]]


# ============================================================
# PIVOTS + BREADTH
# ============================================================

def cell_status(value: float, tested: bool, all_nn: bool,
                metric: str, mode: str, threshold: float,
                ge: bool = True) -> str:
    """
    Display status for one cell:
      'not_tested' | 'nn' | 'miss' | 'hit'
    mode == 'threshold' -> binary hit/miss (+ nn / not_tested)
    mode == 'gradient'  -> 'hit' wherever a numeric value exists (+ nn / not_tested)
    """
    if not tested:
        return "not_tested"
    if metric == "log10_ic50" and all_nn:
        return "nn"
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "nn" if metric == "log10_ic50" else "not_tested"
    if metric == "pct_neut" and value <= 0:
        return "nn"
    if mode == "threshold":
        hit = (value >= threshold) if ge else (value > threshold)
        return "hit" if hit else "miss"
    return "hit"


def build_pivots(df_filtered_or_view, view=None, bucket=None, metric='log10_ic50', mode='none', threshold=None, ge=None):
    # Support both build_pivots(f, view, ...) and build_pivots(view, ...)
    if view is None:
        view = df_filtered_or_view
    """Build pivot tables for heatmap visualization from the compute_view output."""
    ge = ge if ge is not None else True
    threshold = threshold or 0.0

    v = view.copy()
    if bucket and bucket != "All":
        v = v[v['Bucket_Type'] == bucket]

    v['status'] = v.apply(
        lambda r: cell_status(
            r['value'], r['n_tested'] > 0, r['all_nn'],
            metric, mode, threshold, ge
        ),
        axis=1,
    )

    value_pivot = v.pivot_table(
        index='Construct_Description', columns='PSV',
        values='value', aggfunc='mean',
    )
    status_pivot = v.pivot_table(
        index='Construct_Description', columns='PSV',
        values='status', aggfunc=lambda x: 'hit' if (x == 'hit').any() else x.iloc[0],
    )
    counts = v.groupby(['Construct_Description', 'PSV'])['n_tested'].sum()

    # Ensure identical shape
    idx = value_pivot.index.union(status_pivot.index)
    cols = value_pivot.columns.union(status_pivot.columns)
    value_pivot = value_pivot.reindex(index=idx, columns=cols)
    status_pivot = status_pivot.reindex(index=idx, columns=cols).fillna('not_tested')

    return value_pivot, status_pivot, counts

def apply_custom_construct_labels(df):
    """Apply custom labels for PVXE181 G2, G5, G6 multi-round constructs"""
    def get_custom_label(row):
        key = (row['Experiment'], row['Group'], row['Bucket_Type'])
        return CUSTOM_LABELS.get(key, row['Construct_Description'])
    
    df['Construct_Description'] = df.apply(get_custom_label, axis=1)
    return df


def get_curve(df: pd.DataFrame, construct: str, psv: str,
              buckets: list[str] | None = None) -> pd.DataFrame:
    """
    Mean % neutralization vs dilution for one construct x PSV (across replicates).
    Returns a DataFrame with columns [dilution, pct_neut].
    """
    d = df[(df["Construct_Description"] == construct) & (df["PSV"] == psv)]
    if buckets:
        d = d[d["Bucket_Type"].isin(buckets)]
    rows = []
    for dil, neut in zip(d["DilutionArray"], d["NeutArray"]):
        if len(dil) == len(neut) and len(dil) > 0:
            rows.extend(zip(dil, neut))
    if not rows:
        return pd.DataFrame(columns=["dilution", "pct_neut"])
    curve = (pd.DataFrame(rows, columns=["dilution", "pct_neut"])
               .groupby("dilution", as_index=False)["pct_neut"].mean()
               .sort_values("dilution"))
    return curve

def parse_ic50(val):
    """Parse IC50 value, handling NN (No Neutralization)"""
    s = str(val).strip().upper()
    
    # Check for NN or empty values
    if s in ('NN', '', 'N/A', 'NA', 'NAN', 'NONE', '-', '0'):
        return 0.0, True  # (log10_value, is_NN_flag)
    
    try:
        v = float(s.replace(',', ''))
        if v > 0:
            return (np.log10(v), False)
        else:
            return (0.0, True)
    except ValueError:
        return (0.0, True)
