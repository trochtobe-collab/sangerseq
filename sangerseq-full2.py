#!/usr/bin/env python3
"""
Sanger Sequencing Pipeline — Full Automation (Steps 1–5)
Replaces: SnapGene (Step 1), BioEdit CAP assembly (Step 2), MEGA alignment (Steps 3–5)

FILENAME CONVENTION:
    {year}_{cohort}_{patient}_{visit}_{amplicon}-{direction}_{runID}.ab1
    Example: 2025_Adult_083_1_A6-FW_6445495.ab1
        Sample key  -> 2025_Adult_083_1   (groups all amplicons for one sample)
        Amplicon    -> A6                  (determines assembly order)
        Direction   -> FW or RV            (RV reads are reverse-complemented)
        Run ID      -> 6445495             (ignored)

FOLDER STRUCTURE:
    input/
        Raw_Data/
            A_strains/    <- RSV A .ab1 files
            B_strains/    <- RSV B .ab1 files
        references/
            ref_RSV_A.fasta
            ref_RSV_B.fasta

OUTPUTS:
    output/
        01_amplicon_reads/     <- Step 1: strand-corrected .fa per amplicon
        02_assembled_contigs/  <- Step 2: assembled contig per sample
        03_analysis_reports/   <- Steps 3–5: per-sample reports + master CSV
        04_extracted_F_protein/<- Extracted F-protein sequences (DNA + AA)

DEPENDENCIES:
    pip install biopython

─────────────────────────────────────────────────────────────────────
DESIGN NOTES
─────────────────────────────────────────────────────────────────────

IUPAC BASE HANDLING
    Ambiguous IUPAC codes (R Y S W K M B D H V) emitted by the base caller
    encode exactly which two bases the caller was uncertain between.
    This information is preserved through the full pipeline:

    Step 1  Each IUPAC base is resolved to its highest-peak candidate using
            raw .ab1 chromatogram data (DATA9–DATA12 channels, PLOC2 index),
            deposited as lowercase (low-confidence).  The full candidate frozenset
            is stored on the SeqRecord (iupac_candidates dict).

    Step 2  Pileup consensus uses cross-strand IUPAC intersection:
            FW = R {A,G}, RV = M {A,C}  →  {A,G} ∩ {A,C} = {A}  → confirmed A.
            A single-element intersection is promoted to uppercase (confirmed).
            Contradictory sets → N (uncovered).  Unresolved multi-element → ambiguous.

    Report  Positions resolved to lowercase (still ambiguous) are tracked
            separately from genuinely uncovered (N/gap) positions.

TRIMMING  (Mott / Trimmomatic SLIDINGWINDOW style — end-only, hard floor 500 bp)
    Algorithm    : Strictly end-only sliding window (equiv. Trimmomatic
                   SLIDINGWINDOW:10:15).  The 5′ end is scanned left→right
                   and the 3′ end right→left; the scan STOPS the moment it
                   finds the first window of TRIM_WINDOW bases whose mean
                   Phred reaches TRIM_MIN_QUAL.  A quality dip anywhere in
                   the interior of the read NEVER causes trimming.
    TRIM_WINDOW  : 10 bp  — standard Sanger window width.
    TRIM_MIN_QUAL: 15     — Q15 ≈ 96.8% per-base accuracy (lenient end-trim).
    Internal mask: bases below Q15 temporarily masked to N for anchor
                   alignment only; original base is still deposited.
    MIN_USABLE_BP: 100 bp — discard read only if fewer bp survive trim.
    MIN_KEEP_BP  : 500 bp — hard floor.  Quality trimming will NEVER reduce
                   a read below 500 bp; boundaries expand outward from the
                   anchor windows if needed, with a warning logged.

SEMI-GLOBAL ALIGNER (Step 3)
    Gap open penalty = -20.  Sanger reads are physically continuous; any gap
    in the sample side of the alignment is almost certainly an artefact of
    low-quality bases, not a real deletion.  The high penalty forces the aligner
    to accept mismatches rather than opening a gap, preventing the cascade of
    false mutations caused by a single-base frame shift.

FRAMESHIFT REPAIR (three-strategy cascade)
    After building the corrected CDS, the pipeline checks for internal stop
    codons and attempts repair using three strategies in order:

    Strategy 1 — Window scan (fast, most targeted):
        Scans ±60 bp upstream of each stop codon and tries DEL/INS at every
        position in that window.  Catches the vast majority of single-indel
        artefacts near the error site.

    Strategy 2 — Low-quality pairs (handles compound frameshifts):
        Takes the 20 lowest-Phred positions in the CDS and tries all pairwise
        combinations of DEL+DEL, INS+INS, DEL+INS, INS+DEL.  Handles +2/-2
        and ±1 compound shifts caused by two nearby artefact bases.

    Strategy 3 — Full CDS single scan (exhaustive fallback):
        Tries DEL and INS at every position in the entire CDS.  Only reached
        if both Strategy 1 and 2 fail.

    If any strategy eliminates all stops, the repair is applied and flagged for
    manual confirmation.  If only a partial reduction is achieved, the best
    partial result is applied and flagged.  If all strategies fail, specific
    chromatogram positions are reported for manual inspection.

PREVIOUSLY FIXED BUGS (preserved for reference)
    FIX 1  corrected_cds uppercased before .translate() — BioPython ignores
           lowercase, emitting X for every codon containing one.
    FIX 2  corrected_chars never includes gap '–' characters.
    FIX 3  ref_filled_dna_pos built from corrected_chars indices, not gapped
           alignment positions.
"""

# =============================================================================
#
#   ██████╗  █████╗ ██████╗  █████╗ ███╗   ███╗███████╗████████╗███████╗██████╗ ███████╗
#   ██╔══██╗██╔══██╗██╔══██╗██╔══██╗████╗ ████║██╔════╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
#   ██████╔╝███████║██████╔╝███████║██╔████╔██║█████╗     ██║   █████╗  ██████╔╝███████╗
#   ██╔═══╝ ██╔══██║██╔══██╗██╔══██║██║╚██╔╝██║██╔══╝     ██║   ██╔══╝  ██╔══██╗╚════██║
#   ██║     ██║  ██║██║  ██║██║  ██║██║ ╚═╝ ██║███████╗   ██║   ███████╗██║  ██║███████║
#   ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
#
#   ╔══════════════════════════════════════════════════════════════════════════╗
#   ║              PIPELINE PARAMETERS — EDIT THIS SECTION ONLY              ║
#   ║     All settings that affect pipeline behaviour live here.             ║
#   ║     You do NOT need to touch anything below to tune the pipeline.      ║
#   ╚══════════════════════════════════════════════════════════════════════════╝
#
# =============================================================================
#
# HOW TO READ THIS SECTION
# ─────────────────────────
#   Each parameter has:
#     • A plain-English explanation of what it does
#     • The effect of making it HIGHER vs LOWER
#     • A recommended value for typical Sanger RSV data
#     • A preset table where relevant
#
#   If you are unsure, just use the recommended values — they have been
#   validated against real Sanger chromatograms from this lab.
#
# =============================================================================
#
#
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  SECTION 1 — READ QUALITY TRIMMING  (Steps 1 & 2)                     │
# │                                                                         │
# │  The pipeline trims low-quality bases from the ENDS of each read       │
# │  before assembly. Interior bases are NEVER trimmed — only the ragged   │
# │  ends produced by the Sanger chemistry are removed.                    │
# │                                                                         │
# │  Method: sliding window scan from each end inward. The scan stops the  │
# │  moment it finds a window of TRIM_WINDOW bases whose average Phred     │
# │  quality score reaches TRIM_MIN_QUAL.                                   │
# └─────────────────────────────────────────────────────────────────────────┘

# ── How wide is the sliding window? ──────────────────────────────────────────
#
#   The window scans this many bases at once when checking end quality.
#   A wider window is more stable but less sensitive to short bad patches.
#   10 is the standard for Sanger sequencing.
#
#   ↑ Higher  → smoother, less sensitive to a single bad base at the end
#   ↓ Lower   → more aggressive trimming of short noisy end patches
#   Recommended: 10
#
TRIM_WINDOW = 10


# ── How good does the end need to be before we stop trimming? ────────────────
#
#   Phred quality score threshold for the sliding window.
#   The scan stops trimming once the window average reaches this value.
#
#   Phred scale reminder:
#       Q10 = 90% accuracy per base   (very lenient — keep almost everything)
#       Q13 = 95% accuracy per base   ← recommended for RSV Sanger data
#       Q15 = 97% accuracy per base   (moderate)
#       Q20 = 99% accuracy per base   (strict — may leave gaps between amplicons)
#
#   Typical read lengths after trimming (~1050–1254 bp raw input):
#
#       TRIM_WINDOW   TRIM_MIN_QUAL   kept FW     kept RV     notes
#       ───────────   ─────────────   ─────────   ─────────   ──────────────────
#            10             20          992 bp      903 bp    strict — too short,
#                                                             amplicons may not overlap
#            10             15         1070 bp     1027 bp    previous default
#            10             13         1217 bp     1045 bp    ← recommended
#            10             10         1239 bp     1049 bp    very lenient
#             5             15         1209 bp     1036 bp    narrow window
#
#   ↑ Higher  → more aggressive trimming, shorter reads, cleaner ends, more gaps
#   ↓ Lower   → longer reads, more overlap between amplicons, slightly noisier ends
#   Recommended: 13
#
TRIM_MIN_QUAL = 13


# ── Hard floor: never trim a read shorter than this ──────────────────────────
#
#   Even if quality trimming would produce a shorter result, the pipeline
#   will always keep at least this many bases. If the floor is hit, a WARNING
#   is printed so you know the read end may be lower quality than usual.
#
#   ↑ Higher  → more bases always kept, but poor-quality ends are retained
#   ↓ Lower   → allows very short reads through (useful if your amplicons are short)
#   Recommended: 500
#
MIN_KEEP_BP = 500


# ── Discard threshold: throw away reads shorter than this after trimming ─────
#
#   If a read is so degraded that fewer than this many bases survive trimming,
#   it is discarded entirely and a WARNING is logged.
#   Only relevant for severely degraded samples.
#
#   ↑ Higher  → discard more short reads (stricter)
#   ↓ Lower   → keep very short reads (may introduce noise)
#   Recommended: 100
#
MIN_USABLE_BP = 100


# ┌─────────────────────────────────────────────────────────────────────────┐
# │  SECTION 2 — PILEUP / ASSEMBLY QUALITY  (Step 2)                       │
# │                                                                         │
# │  After trimming, reads are stacked on top of each other (pileup) and   │
# │  a consensus base is chosen at each position. These two parameters     │
# │  control how the pipeline handles low-quality bases in that process.   │
# └─────────────────────────────────────────────────────────────────────────┘

# ── Below what quality should a base be treated as nearly unreadable? ────────
#
#   Bases with Phred below this value are deposited into the pileup as
#   lowercase 'n' — meaning "probably something, but too noisy to trust".
#   They can still contribute to the pileup vote, but they LOSE to any
#   confident uppercase base call from the other strand.
#
#   ↑ Higher  → more bases demoted to 'n', more positions rely on the other strand
#   ↓ Lower   → even very noisy bases keep their identity (more noise tolerated)
#   Recommended: 5  (only truly flat signal is demoted)
#
MIN_DEPOSIT_QUAL = 5


# ── Below what quality should a base be masked during alignment anchoring? ───
#
#   When anchoring a read onto the reference, bases below this Phred are
#   temporarily replaced with 'N' to give the aligner a cleaner sequence
#   to work with. This does NOT change the base that gets deposited in the
#   pileup — it only helps the aligner find the correct position.
#
#   Think of it as: "hide noisy bases from the aligner, but still use them
#   in the final consensus vote."
#
#   ↑ Higher  → more bases hidden from the aligner (cleaner anchor, less precise)
#   ↓ Lower   → aligner sees noisier bases (may anchor slightly off in very bad reads)
#   Recommended: 15
#
MASK_QUAL = 15


# ┌─────────────────────────────────────────────────────────────────────────┐
# │  SECTION 3 — ALIGNMENT GAP PENALTIES  (Steps 2–5)                      │
# │                                                                         │
# │  These control how the aligner handles gaps (insertions / deletions).   │
# │  A Sanger read is continuous DNA — gaps in the read vs the reference   │
# │  are almost always sequencing artefacts, not real biology.             │
# │  High gap penalties force the aligner to call a mismatch instead of    │
# │  opening a gap, which prevents false frameshifts.                      │
# └─────────────────────────────────────────────────────────────────────────┘

# ── Final contig vs reference: how hard to resist opening a gap? ─────────────
#
#   Used in Steps 3–5 when aligning the assembled contig against the reference
#   to call mutations.
#
#   Score is NEGATIVE — more negative = stronger penalty = fewer gaps.
#
#   ↑ Less negative (e.g. -10) → aligner opens gaps more freely
#                                 (may introduce phantom frameshifts)
#   ↓ More negative (e.g. -30) → aligner almost never opens gaps
#                                 (may call false mismatches near real indels)
#   Recommended: -20
#
ALIGN_GAP_OPEN   = -20   # gap open penalty for final contig alignment
ALIGN_GAP_EXTEND =  -1   # gap extend penalty (each extra base in a gap)
                          # Recommended: -1 — small extra cost per gap base


# ── Read anchoring: gap penalties used during Step 2 assembly ────────────────
#
#   When anchoring individual reads onto the reference (Step 2), slightly
#   different penalties are used because reads are shorter and their ends
#   can be ragged.
#
#   ANCHOR_REF_GAP_OPEN
#       Penalty for a gap in the reference side (= deletion in the read).
#       A real deletion in the sample is biologically possible, so this is
#       moderately penalised rather than forbidden.
#       Recommended: -5
#
#   ANCHOR_QUERY_GAP_OPEN
#       Penalty for a gap in the read side (= insertion in the read vs ref).
#       For a continuous Sanger read this is almost always a sequencing
#       artefact. Set extremely high to effectively forbid it.
#       Recommended: -1000  (do not change unless you are studying real insertions)
#
ANCHOR_REF_GAP_OPEN   =    -5   # deletion in read vs ref  (moderate penalty)
ANCHOR_QUERY_GAP_OPEN = -1000   # insertion in read vs ref (effectively forbidden)


# =============================================================================
#   ╔══════════════════════════════════════════════════════════════════════════╗
#   ║             END OF PARAMETERS — do not edit below this line            ║
#   ╚══════════════════════════════════════════════════════════════════════════╝
# =============================================================================


import os
import re
import glob
import csv
import logging
import shutil
import struct
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

SEP = "=" * 60   # report section separator

# =============================================================================
# IUPAC AMBIGUITY TABLE
# Maps each IUPAC ambiguity code to the set of bases it represents.
# Used to resolve which bases contribute to an uncertain call.
# =============================================================================
IUPAC_BASES = {
    'R': ('A', 'G'),  # puRine
    'Y': ('C', 'T'),  # pYrimidine
    'S': ('G', 'C'),  # Strong
    'W': ('A', 'T'),  # Weak
    'K': ('G', 'T'),  # Keto
    'M': ('A', 'C'),  # aMino
    'B': ('C', 'G', 'T'),  # not A
    'D': ('A', 'G', 'T'),  # not C
    'H': ('A', 'C', 'T'),  # not G
    'V': ('A', 'C', 'G'),  # not T/U
    'N': ('A', 'C', 'G', 'T'),  # aNy
}

# IUPAC codes that indicate genuine ambiguity (exclude N which means 'no call')
IUPAC_AMBIG = set(IUPAC_BASES.keys()) - {'N'}


# =============================================================================
# SHARED HELPERS
# =============================================================================

def parse_filename(filepath):
    """
    Parses .ab1 filenames. Handles several naming conventions seen in the wild:

        Old:  {sample}_{amplicon}-{direction}_{run_id}.ab1
              2025_Adult_083_1_A6-FW_6445495.ab1

        New:  {sample}_{amplicon}-{amplicon}-{direction}_{run_id}.ab1
              RSV_A2_A6-A6-RV_6449075.ab1
              RSV_B1_B6-B6-FW_6449073.ab1

        Short direction codes (F/R) are normalised to FW/RV.

    Strain is derived from the amplicon prefix:
        A<n>  ->  RSV_A
        B<n>  ->  RSV_B

    Returns a dict with keys: sample_key, amplicon, strain, direction, run_id.
    Returns None if the filename cannot be parsed.
    """
    base = os.path.splitext(os.path.basename(filepath))[0]
    match = re.match(
        r'^(?P<sample>.+)_(?P<amplicon>[AB]\d+)-(?:[AB]\d+-)?(?P<direction>FW|RV|F|R)_(?P<run_id>\d+)$',
        base, re.IGNORECASE
    )
    if match:
        amplicon  = match.group("amplicon").upper()
        raw_dir   = match.group("direction").upper()
        direction = {"F": "FW", "R": "RV"}.get(raw_dir, raw_dir)
        strain    = "RSV_A" if amplicon.startswith("A") else "RSV_B"
        return {
            "sample_key": match.group("sample"),
            "amplicon":   amplicon,
            "strain":     strain,
            "direction":  direction,
            "run_id":     match.group("run_id"),
        }
    log.warning(f"Could not parse filename: {base} — skipping.")
    return None


def is_reverse_amplicon(filepath):
    parsed = parse_filename(filepath)
    if parsed:
        return parsed["direction"] == "RV"
    return bool(re.search(r'-RV[_\-]', os.path.basename(filepath), re.IGNORECASE))


def resolve_iupac_best_base(iupac_char, ab1_record, position):
    """
    For an IUPAC ambiguous base (yellow/uncertain letter in the chromatogram
    viewer) resolve it to the single most-likely base using raw peak heights
    from the .ab1 file.

    HOW IT WORKS
    ------------
    The base caller already chose the primary base call (the letter before the
    IUPAC code was emitted).  The IUPAC code tells us *which two bases* the
    caller was uncertain between.  We read the raw chromatogram signal at the
    peak index for this base position from all four channels (DATA9–DATA12) and
    pick whichever of the two IUPAC candidates has the stronger signal.

    This is equivalent to asking: "given that the base caller couldn't decide
    between A and G (R), which one had the taller peak?"  That answer is our
    best single-base estimate.

    The resolved base is returned as **lowercase** so downstream pileup voting
    knows it is a low-confidence call.  If the opposite strand has a confident
    uppercase call at the same position, the confident call wins the vote
    automatically.

    .ab1 channel layout (BioPython abif_raw):
        DATA9  = G channel
        DATA10 = A channel
        DATA11 = T channel
        DATA12 = C channel
        PLOC2  = array of peak-scan indices, one per base call

    Parameters
    ----------
    iupac_char : str       The IUPAC ambiguity character (e.g. 'R', 'Y').
    ab1_record : SeqRecord The parsed .ab1 SeqRecord from BioPython.
    position   : int       0-based index in the base-called sequence.

    Returns
    -------
    str  Single lowercase base character (e.g. 'a', 'g').
         Falls back to the first base in the IUPAC pair if peak data is absent.
    """
    possible_bases = IUPAC_BASES.get(iupac_char.upper(), ('A',))

    # Channel map: base → abif_raw key
    channel_map = {
        'G': 'DATA9',
        'A': 'DATA10',
        'T': 'DATA11',
        'C': 'DATA12',
    }

    try:
        raw = ab1_record.annotations.get('abif_raw', {})

        # PLOC2: array giving the scan-index of each base-call peak
        peak_locs = raw.get('PLOC2', None)
        if peak_locs is None:
            raise ValueError("PLOC2 not found in .ab1 — cannot read peak positions")

        if position >= len(peak_locs):
            raise ValueError(f"position {position} out of range for PLOC2 (len={len(peak_locs)})")

        peak_scan_idx = peak_locs[position]

        best_base   = None
        best_signal = -1

        for base in possible_bases:
            ch_key = channel_map.get(base)
            if ch_key is None:
                continue
            channel_data = raw.get(ch_key)
            if channel_data is None:
                continue
            # Guard against peak index falling outside the trace array
            if peak_scan_idx >= len(channel_data):
                signal = 0
            else:
                signal = channel_data[peak_scan_idx]

            if signal > best_signal:
                best_signal = signal
                best_base   = base

        if best_base is not None:
            log.debug(
                f"    IUPAC '{iupac_char}' pos {position}: resolved to "
                f"'{best_base}' (signal={best_signal}) from candidates {possible_bases}"
            )
            return best_base.lower()

    except Exception as e:
        log.debug(
            f"    Peak resolution failed for IUPAC '{iupac_char}' at pos {position}: {e}"
        )

    # Fallback: use the first base in the IUPAC pair as lowercase
    fallback = possible_bases[0].lower()
    log.debug(f"    IUPAC '{iupac_char}' pos {position}: fallback to '{fallback}'")
    return fallback


def mask_ambiguous_bases(seq_str, ab1_record=None):
    """
    Process a raw sequence string for Step 1.

    For each position:
      - Uppercase ATGC  -> confident call, kept as-is.
      - IUPAC code      -> resolved to the highest-peak candidate base
                           (lowercase = low confidence) using chromatogram data.
                           The full candidate set is stored in iupac_candidates
                           so the pileup voter can use cross-strand intersection.
      - N / '-'         -> kept as-is (handled downstream as uncovered).

    Returns
    -------
    seq_out         : str  — sequence with IUPAC replaced by best-peak lowercase base
    iupac_positions : list — 0-based positions that contained IUPAC codes
    iupac_candidates: dict — {0-based pos: frozenset of candidate bases}
                             e.g. {42: frozenset({'A','G'})} for R at pos 42
    """
    result           = []
    iupac_positions  = []
    iupac_candidates = {}   # position -> frozenset of IUPAC-implied bases

    for i, base in enumerate(seq_str.upper()):
        if base in 'ATGC':
            result.append(base)
        elif base in IUPAC_AMBIG:
            candidates = frozenset(IUPAC_BASES[base])
            iupac_candidates[i] = candidates
            # Best-peak single base (lowercase) deposited into the sequence
            resolved = resolve_iupac_best_base(base, ab1_record, i) if ab1_record \
                       else list(candidates)[0].lower()
            result.append(resolved)
            iupac_positions.append(i)
        elif base in ('N', '-'):
            result.append(base)
        else:
            result.append(base)

    return "".join(result), iupac_positions, iupac_candidates


def make_aligner(mode, match=2, mismatch=-1, open_gap=-5, extend_gap=-1):
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    return aligner


def format_alignment_visual(ref_seq, samp_seq, line_len=60):
    """
    Formats a pairwise alignment for display with:
      - Position numbers at the start and end of each line (1-based, ref positions)
      - Match/mismatch/gap indicators between ref and sample rows
      - Stop codons (*) are excluded from mutation marking to avoid edge artifacts
    """
    lines = []
    ref_pos = 1

    for i in range(0, len(ref_seq), line_len):
        r = ref_seq[i:i + line_len]
        s = samp_seq[i:i + line_len]

        mid = "".join(
            "|" if rc == sc and rc not in ("-", "*")
            else " " if rc == "-" or sc == "-"
            else "|" if rc == "*" and sc == "*"
            else "*"
            for rc, sc in zip(r, s)
        )

        chunk_ref_start = ref_pos
        chunk_ref_end   = ref_pos + sum(1 for c in r if c != "-") - 1

        label_w = len(str(len(ref_seq))) + 1
        lines.append(
            f"  {chunk_ref_start:<{label_w}} Ref:  {r}  {chunk_ref_end}\n"
            f"  {'':<{label_w}}       {mid}\n"
            f"  {'':<{label_w}} Samp: {s}\n"
        )

        ref_pos += sum(1 for c in r if c != "-")

    return "\n".join(lines)


# =============================================================================
# STEP 1 — Read .ab1 files (replaces SnapGene)
# =============================================================================

def process_ab1_file(ab1_path):
    """
    Reads a .ab1 chromatogram, reverse-complements RV amplicons,
    marks ambiguous bases as lowercase, resolves IUPAC codes using peak
    heights, and returns a SeqRecord.

    Improvements over original:
      - IUPAC ambiguous bases (R/Y/S/W/K/M/B/D/H/V) are resolved to their
        best single-base estimate (lowercase) using raw chromatogram peak data.
      - Track IUPAC positions so they can be reported separately from N gaps.
      - The raw ab1_record is passed to mask_ambiguous_bases() for peak lookup.
    """
    filename = os.path.basename(ab1_path)
    try:
        record = SeqIO.read(ab1_path, "abi")
    except Exception as e:
        log.error(f"Could not read {filename}: {e}")
        return None, []

    seq = str(record.seq)
    if not seq:
        log.warning(f"{filename}: Empty sequence — skipping.")
        return None, []

    quals = list(record.letter_annotations.get("phred_quality", [20] * len(seq)))

    if is_reverse_amplicon(ab1_path):
        seq   = str(Seq(seq).reverse_complement())
        quals = quals[::-1]
        log.info(f"{filename}: RV amplicon — reverse complemented.")

    # Pass the raw record for IUPAC peak resolution.
    # iupac_candidates maps each 0-based read position to the frozenset of
    # bases that the IUPAC code implies, e.g. R -> frozenset({'A','G'}).
    # This is stored on the SeqRecord so the pileup voter can perform
    # cross-strand IUPAC intersection (see _deposit_read_into_pileup).
    final_seq, iupac_positions, iupac_candidates = mask_ambiguous_bases(
        seq, ab1_record=record
    )

    if iupac_positions:
        log.info(f"{filename}: {len(iupac_positions)} IUPAC ambiguous bases resolved "
                 f"(positions: {iupac_positions[:10]}{'...' if len(iupac_positions) > 10 else ''})")

    out = SeqRecord(
        Seq(final_seq),
        id=os.path.splitext(filename)[0],
        description=f"len={len(final_seq)}"
    )
    out.letter_annotations["phred_quality"] = quals
    out.annotations["iupac_positions"]  = iupac_positions
    out.annotations["iupac_candidates"] = iupac_candidates   # {read_pos: frozenset}
    return out, iupac_positions


def step1_process_ab1_folder(ab1_folder, output_folder_base):
    """
    Scans a single folder for all .ab1 files and automatically separates them
    by strain (RSV-A / RSV-B) based on the amplicon prefix in the filename.

    Returns: {"RSV_A": {sample_key: [("A6-FW", SeqRecord), ...]},
              "RSV_B": {sample_key: [("B6-FW", SeqRecord), ...]}}
    Also returns a summary of IUPAC positions per file for reporting.
    """
    ab1_files = glob.glob(os.path.join(ab1_folder, "**", "*.ab1"), recursive=True)

    if not ab1_files:
        log.warning(f"No .ab1 files found in {ab1_folder} (searched recursively)")
        return {"RSV_A": {}, "RSV_B": {}}

    by_strain = {"RSV_A": defaultdict(list), "RSV_B": defaultdict(list)}

    for ab1_path in sorted(ab1_files):
        parsed = parse_filename(ab1_path)
        if parsed is None:
            continue

        record, iupac_pos = process_ab1_file(ab1_path)
        if record is None:
            continue

        strain     = parsed["strain"]
        clean_name = f"{parsed['sample_key']}_{parsed['amplicon']}-{parsed['direction']}"
        record.id  = clean_name
        record.description = f"len={len(record.seq)}"

        out_folder = os.path.join(output_folder_base, strain)
        os.makedirs(out_folder, exist_ok=True)
        SeqIO.write(record, os.path.join(out_folder, f"{clean_name}.fa"), "fasta")

        if len(record.seq) < 100:
            log.warning(f"  SKIPPING {clean_name} ({len(record.seq)} bp) — too short.")
            continue

        iupac_note = f", {len(iupac_pos)} IUPAC bases resolved" if iupac_pos else ""
        log.info(f"  OK [{strain}] {clean_name}.fa ({len(record.seq)} bp{iupac_note})")
        amp_dir_key = f"{parsed['amplicon']}-{parsed['direction']}"
        by_strain[strain][parsed["sample_key"]].append((amp_dir_key, record))

    for strain, sample_amplicons in by_strain.items():
        for key in sample_amplicons:
            sample_amplicons[key].sort(key=lambda item: (
                int(re.search(r'\d+', item[0]).group()),
                0 if item[0].endswith('FW') else 1
            ))

    return {strain: dict(samples) for strain, samples in by_strain.items()}


# =============================================================================
# STEP 2 — Two-level Sanger Assembly (replaces BioEdit CAP)
#
# IUPAC-AWARE PILEUP VOTING:
#   Lowercase bases (resolved IUPAC or low-quality) lose the majority vote
#   to any uppercase (confident) base at the same position.
#   If FW says 'a' (IUPAC-resolved) and RV says 'A' (confident), 'A' wins.
#   If both strands are lowercase/ambiguous, the reference base is used as
#   placeholder and the position is flagged as 'AMBIGUOUS' in the report.
#
# TRIMMING STRATEGY:
#   All trimming parameters are defined in the PIPELINE PARAMETERS block
#   at the top of this file.  See that block to adjust strictness.
#   IUPAC bases are deposited as lowercase regardless of Phred so the
#   majority-vote preferentially uses the unambiguous strand.
# =============================================================================


def _make_anchor_aligner():
    """
    Smith-Waterman local aligner used to anchor reads onto the reference.

    GAP PENALTY ASYMMETRY
    ─────────────────────
    Sanger reads are physically continuous DNA.  A gap in the *query* (sample)
    side of the alignment (ref_step=0, qry_step>0) means the aligner thinks
    there is an insertion in the sample relative to the reference.  For a
    continuous Sanger read this is almost always a sequencing artefact caused
    by a noisy base cluster, not a real insertion.  Real insertions in RSV
    F-protein are extremely rare and would need to be confirmed on both strands.

    We therefore use strongly asymmetric gap penalties:
      • target (ref) gaps  → open -5  : moderate, allows the read to cover a
                                         genuine deletion in the sample
      • query (sample) gaps → open -1000 : effectively forbidden; the aligner
                                           will accept many mismatches rather
                                           than introduce a spurious insertion
    This prevents the cascade of downstream frameshifts caused by small
    phantom insertions at low-quality positions.
    """
    a = PairwiseAligner()
    a.mode                         = 'local'
    a.match_score                  =  2
    a.mismatch_score               = -1
    # ref-side gap (deletion in read): allowed but penalised
    a.target_internal_open_gap_score   = ANCHOR_REF_GAP_OPEN
    a.target_internal_extend_gap_score = -1
    a.target_end_gap_score             =  0
    # query-side gap (insertion in read): essentially forbidden
    a.query_internal_open_gap_score    = ANCHOR_QUERY_GAP_OPEN
    a.query_internal_extend_gap_score  = ANCHOR_QUERY_GAP_OPEN
    a.query_end_gap_score              =  0
    return a


def _sliding_window_trim(quals):
    """
    End-only quality trimming — equivalent to Trimmomatic SLIDINGWINDOW.

    ALGORITHM  (Mott / Trimmomatic SLIDINGWINDOW style)
    ─────────────────────────────────────────────────────────────────────────
    The scan is strictly end-only.  A quality dip in the middle of a read
    NEVER causes trimming; only the two ends are examined.

    5′ end (left trim):
        Walk a window of TRIM_WINDOW bases from position 0 inward,
        one base at a time.  As long as the window mean Phred < threshold,
        advance.  Stop (and set trim_start = window_left_edge) the moment
        the window mean first meets or exceeds the threshold.

    3′ end (right trim):
        Mirror image: walk from the last base inward.  Stop (set
        trim_end = window_right_edge) at the first window that meets
        the threshold.

    Once both anchors are found, everything outside them is discarded.
    The interior of the read is not examined again.

    This matches the behaviour of:
        Trimmomatic  SLIDINGWINDOW:10:15
        Phred/Phrap  -trim_alt  (Ewing & Green 1998)
        EMBOSS trimseq  -window / -percent

    PARAMETERS  (all set in the PIPELINE PARAMETERS block at top of file)
        TRIM_WINDOW   — window width in bases
        TRIM_MIN_QUAL — mean Phred threshold (lower = keep more)
        MIN_KEEP_BP   — hard floor in bases
        MIN_USABLE_BP — absolute discard threshold in bases

    FALLBACK
        If the TRIM_MIN_QUAL trim leaves fewer than MIN_USABLE_BP bases, the
        threshold is relaxed step-wise (Q10 → Q7 → Q5) and the trim is retried.

    HARD FLOOR  (MIN_KEEP_BP = 500 bp)
        Quality trimming will NEVER reduce a read below MIN_KEEP_BP bases.
        If the best quality-based trim would leave fewer than MIN_KEEP_BP,
        the boundaries are expanded outward from the anchor windows until
        MIN_KEEP_BP bases are retained, and a WARNING is logged.
        The read is only discarded when the raw read itself is shorter than
        MIN_USABLE_BP.
    """
    n = len(quals)
    if n < TRIM_WINDOW:
        return 0, 0

    def _trim_at(threshold):
        """
        Strict end-only scan.
        5′: walk left→right; return left edge of first good window.
        3′: walk right→left; return right edge of first good window.
        Returns (0, 0) if no qualifying window exists on either end.
        """
        # ── 5′ end ──────────────────────────────────────────────────────
        trim_start = None
        for i in range(n - TRIM_WINDOW + 1):
            if sum(quals[i : i + TRIM_WINDOW]) / TRIM_WINDOW >= threshold:
                trim_start = i      # left edge of first good window
                break               # stop here — never look further inward

        if trim_start is None:
            return 0, 0             # entire read is below threshold

        # ── 3′ end ──────────────────────────────────────────────────────
        trim_end = None
        for i in range(n - TRIM_WINDOW, -1, -1):
            if sum(quals[i : i + TRIM_WINDOW]) / TRIM_WINDOW >= threshold:
                trim_end = i + TRIM_WINDOW  # one-past-end of last good window
                break               # stop here — never look further inward

        if trim_end is None:
            return 0, 0

        if trim_end <= trim_start:  # degenerate: anchors crossed
            return 0, 0

        return trim_start, trim_end

    # ── Try thresholds most → least strict ───────────────────────────────
    best_start, best_end = None, None
    used_threshold = TRIM_MIN_QUAL
    for threshold in (TRIM_MIN_QUAL, 10, 7, 5):
        s, e = _trim_at(threshold)
        if e - s >= MIN_USABLE_BP:
            best_start, best_end = s, e
            used_threshold = threshold
            break

    if best_start is None:
        # No threshold produced any usable window — retain full read
        best_start, best_end = 0, n
        used_threshold = 0

    if 0 < used_threshold < TRIM_MIN_QUAL:
        log.warning(
            f"  Quality trim fallback to Q{used_threshold} — "
            f"read has very poor ends; {best_end - best_start} bp core retained."
        )

    # ── HARD FLOOR: never trim below MIN_KEEP_BP ─────────────────────────
    trimmed_len = best_end - best_start
    if trimmed_len < MIN_KEEP_BP and n >= MIN_KEEP_BP:
        # Expand outward from current anchor edges symmetrically
        need         = MIN_KEEP_BP - trimmed_len
        expand_left  = need // 2
        expand_right = need - expand_left
        new_start    = max(0, best_start - expand_left)
        new_end      = min(n, best_end   + expand_right)
        # If one side hits the read boundary, give the remainder to the other
        if new_start == 0:
            new_end = min(n, MIN_KEEP_BP)
        if new_end == n:
            new_start = max(0, n - MIN_KEEP_BP)
        log.warning(
            f"  Quality trim would leave only {trimmed_len} bp — "
            f"enforcing {MIN_KEEP_BP} bp floor: retaining bases "
            f"{new_start}–{new_end} ({new_end - new_start} bp). "
            f"End regions may contain degraded quality bases."
        )
        return new_start, new_end

    if best_end - best_start < MIN_USABLE_BP:
        return 0, 0

    return best_start, best_end


def _is_iupac_resolved(base):
    """Return True if base is a lowercase resolved IUPAC call (not raw low-qual)."""
    # We distinguish IUPAC-resolved (lowercase single base from IUPAC expansion)
    # from low-quality masked bases by checking the annotation on the record.
    # At this point in the pipeline both are lowercase; we treat them equally
    # in the vote (lowercase loses to uppercase).
    return base.islower() and base.upper() in 'ATGC'


def _deposit_read_into_pileup(label, record, ref_str, pileup, pileup_qual,
                               pileup_candidates, aligner, ambiguous_pileup_positions):
    """
    Anchor one read onto the reference and walk it base-by-base into the pileup.

    pileup[i]           = list of base chars deposited at reference position i
    pileup_qual[i]      = list of Phred scores for those bases
    pileup_candidates[i]= list of frozensets — for IUPAC bases the set of all
                          candidate bases implied by the code (e.g. frozenset({'A','G'})
                          for R), or None for confident/low-qual bases.
                          Used by _pileup_to_consensus for cross-strand intersection.

    Returns (ref_first_1based, ref_last_1based) or (None, None).
    """
    ref_len         = len(pileup)
    amp_seq         = str(record.seq)
    amp_qual        = list(record.letter_annotations.get("phred_quality", [20] * len(amp_seq)))
    # iupac_candidates: {read_pos_before_trim: frozenset}
    # Keys are indices into the *original* (pre-trim) sequence.
    read_candidates = record.annotations.get("iupac_candidates", {})

    trim_start, trim_end = _sliding_window_trim(amp_qual)
    usable_bp = trim_end - trim_start
    if usable_bp < MIN_USABLE_BP:
        mean_q = sum(amp_qual) / len(amp_qual) if amp_qual else 0
        log.warning(
            f"    SKIP {label} — only {usable_bp} bp survive quality trim "
            f"(read length {len(amp_seq)} bp, mean Q={mean_q:.1f}). "
            f"Read is too degraded to use."
        )
        return None, None

    trimmed_seq  = amp_seq[trim_start:trim_end]
    trimmed_qual = amp_qual[trim_start:trim_end]

    # Re-key iupac_candidates to trimmed coordinates
    trimmed_candidates = {
        (k - trim_start): v
        for k, v in read_candidates.items()
        if trim_start <= k < trim_end
    }

    masked_seq = "".join(
        b if (q >= MASK_QUAL and b == b.upper() and b.upper() in "ATGC") else "N"
        for b, q in zip(trimmed_seq, trimmed_qual)
    )

    best = next(iter(aligner.align(ref_str, masked_seq)), None)
    if best is None:
        log.warning(f"    SKIP {label} — could not anchor to reference.")
        return None, None

    coords     = best.coordinates
    ref_coords = coords[0]
    qry_coords = coords[1]

    ref_pos   = int(ref_coords[0])
    amp_pos   = int(qry_coords[0])
    ref_first = None
    ref_last  = None
    deposited = 0

    for seg in range(len(ref_coords) - 1):
        ref_step = int(ref_coords[seg + 1]) - int(ref_coords[seg])
        qry_step = int(qry_coords[seg + 1]) - int(qry_coords[seg])

        if ref_step > 0 and qry_step > 0:
            for _ in range(ref_step):
                if ref_pos < ref_len and amp_pos < len(trimmed_seq):
                    orig_base = trimmed_seq[amp_pos]
                    orig_qual = trimmed_qual[amp_pos]
                    cand_set  = trimmed_candidates.get(amp_pos, None)  # frozenset or None

                    if orig_qual < MIN_DEPOSIT_QUAL:
                        deposit_base = "n"
                        cand_set     = None   # too noisy to carry candidates
                    elif orig_base.islower():
                        # IUPAC-resolved: deposit lowercase, keep candidate set
                        deposit_base = orig_base.lower()
                        ambiguous_pileup_positions.add(ref_pos)
                    elif orig_qual < MASK_QUAL:
                        deposit_base = orig_base.lower()
                        cand_set     = None   # low-qual but unambiguous — no set needed
                    else:
                        deposit_base = orig_base   # uppercase = confident
                        cand_set     = None        # confident call, no set needed

                    pileup[ref_pos].append(deposit_base)
                    pileup_qual[ref_pos].append(orig_qual)
                    pileup_candidates[ref_pos].append(cand_set)
                    deposited += 1
                    if ref_first is None:
                        ref_first = ref_pos + 1
                    ref_last = ref_pos + 1
                ref_pos += 1
                amp_pos += 1
        elif ref_step > 0 and qry_step == 0:
            ref_pos += ref_step
        elif ref_step == 0 and qry_step > 0:
            # Insertion in sample relative to reference.
            # With query_internal_open_gap_score = -1000 this should essentially
            # never happen.  If it does, log it clearly so it's visible.
            ins_bases = trimmed_seq[amp_pos: amp_pos + qry_step]
            ins_quals = trimmed_qual[amp_pos: amp_pos + qry_step]
            mean_q    = sum(ins_quals) / len(ins_quals) if ins_quals else 0
            log.warning(
                f"    {label}: {qry_step}-bp insertion at ref pos {ref_pos} "
                f"(bases '{ins_bases}', mean Q={mean_q:.0f}) — dropped to preserve frame."
            )
            amp_pos += qry_step

    if ref_first is not None:
        log.info(f"    {label}: {deposited} bases -> ref {ref_first}-{ref_last}")
    else:
        log.warning(f"    {label}: deposited 0 bases.")
    return ref_first, ref_last


def _pileup_to_consensus(pileup, pileup_candidates, ambiguous_pileup_positions=None):
    """
    Collapse a pileup into a consensus string using IUPAC cross-strand intersection.

    VOTING HIERARCHY (per position)
    --------------------------------
    1. CONFIDENT (uppercase ATGC) calls exist
       → simple majority vote among uppercase calls.  Done.

    2. No confident calls, but IUPAC candidate sets are present
       → Cross-strand intersection:
           a. Compute the intersection of all candidate sets at this position.
              Example: FW says R = {A,G}, RV says M = {A,C}  → {A,G} ∩ {A,C} = {A}
              A single base in the intersection = CONFIRMED without ambiguity.
           b. If intersection has exactly 1 base → use it as UPPERCASE (confirmed).
           c. If intersection has 2+ bases → take the best-peak deposited base
              as lowercase (still ambiguous, but narrowed down).
           d. If intersection is empty (contradictory IUPAC codes) → treat as N.

    3. No confident calls, no candidate sets (all low-quality lowercase)
       → majority vote among lowercase calls; result is lowercase (ambiguous).
       Position flagged in ambiguous_consensus_pos.

    4. Empty position → 'N' (uncovered).

    Returns: (consensus_str, uncovered_regions, ambiguous_consensus_positions)
    """
    if ambiguous_pileup_positions is None:
        ambiguous_pileup_positions = set()

    consensus               = []
    uncovered               = []
    ambiguous_consensus_pos = set()
    in_gap                  = False
    gap_start               = 0

    for i, bases in enumerate(pileup):
        cand_sets = pileup_candidates[i]   # list of frozenset-or-None, parallel to bases

        if not bases:
            consensus.append("N")
            if not in_gap:
                in_gap, gap_start = True, i
            continue

        if in_gap:
            uncovered.append((gap_start + 1, i))
            in_gap = False

        # ── Tier 1: confident uppercase calls ────────────────────────────────
        confident = [b for b in bases if b in "ATGC"]   # uppercase only
        if confident:
            counts = {}
            for b in confident:
                counts[b] = counts.get(b, 0) + 1
            consensus.append(max(counts, key=counts.get))
            continue

        # ── Tier 2: IUPAC cross-strand intersection ───────────────────────────
        # Gather only the positions that actually have a candidate set
        # (i.e. came from an IUPAC base, not just low-qual noise)
        iupac_entries = [
            (b, cs)
            for b, cs in zip(bases, cand_sets)
            if cs is not None and len(cs) > 0
        ]
        if iupac_entries:
            # Start with the union of all candidates, then intersect
            intersection = iupac_entries[0][1]
            for _, cs in iupac_entries[1:]:
                intersection = intersection & cs

            if len(intersection) == 1:
                # Perfect cross-strand agreement: single base confirmed
                resolved_base = list(intersection)[0].upper()
                log.debug(
                    f"    Position {i+1}: IUPAC intersection resolved to "
                    f"'{resolved_base}' from {[str(cs) for _, cs in iupac_entries]}"
                )
                consensus.append(resolved_base)   # uppercase = confirmed
                continue
            elif len(intersection) > 1:
                # Narrowed but not fully resolved: use best-peak deposited base
                # (already lowercase from Step 1), mark as ambiguous
                counts = {}
                for b, _ in iupac_entries:
                    counts[b.lower()] = counts.get(b.lower(), 0) + 1
                best = max(counts, key=counts.get)
                consensus.append(best)
                ambiguous_consensus_pos.add(i + 1)
                log.debug(
                    f"    Position {i+1}: IUPAC intersection narrowed to "
                    f"{intersection}, using '{best}' (ambiguous)"
                )
                continue
            else:
                # Empty intersection = contradictory IUPAC codes across strands
                # This is very unusual; treat as uncovered
                consensus.append("N")
                log.warning(
                    f"    Position {i+1}: IUPAC codes contradict each other "
                    f"({[str(cs) for _, cs in iupac_entries]}) — treating as N"
                )
                continue

        # ── Tier 3: all low-quality lowercase (no candidate sets) ────────────
        ambiguous = [b for b in bases if b.islower() and b.upper() in "ATGC"]
        if ambiguous:
            counts = {}
            for b in ambiguous:
                counts[b] = counts.get(b, 0) + 1
            consensus.append(max(counts, key=counts.get))
            ambiguous_consensus_pos.add(i + 1)
            log.debug(f"    Position {i+1}: only low-qual bases {ambiguous} — ambiguous")
        else:
            consensus.append("N")
            if not in_gap:
                in_gap, gap_start = True, i

    if in_gap:
        uncovered.append((gap_start + 1, len(pileup)))

    return "".join(consensus), uncovered, ambiguous_consensus_pos


def build_pileup_consensus(amplicon_records, ref_record):
    """
    Two-level Sanger assembly with IUPAC-aware pileup voting.

    amplicon_records : list of ("A1-FW" | "A1-RV" | ..., SeqRecord)
                       RV reads have already been reverse-complemented in Step 1.

    Returns: (consensus_seq_str, coverage_stats_dict, warnings_list)
    """
    ref_len = len(ref_record.seq)
    ref_str = str(ref_record.seq)
    aligner = _make_anchor_aligner()

    # ── Group reads by amplicon number ────────────────────────────────────────
    amp_groups = {}
    for label, record in amplicon_records:
        amp_num = re.match(r'([AB]\d+)', label, re.IGNORECASE)
        if amp_num is None:
            log.warning(f"  Cannot parse amplicon number from '{label}' — skipping.")
            continue
        key = amp_num.group(1).upper()
        amp_groups.setdefault(key, []).append((label, record))

    sorted_amp_keys = sorted(
        amp_groups.keys(),
        key=lambda k: int(re.search(r'\d+', k).group())
    )

    log.info(f"  Amplicons detected: {sorted_amp_keys}")

    # ── LEVEL 1: per-amplicon FW+RV consensus ─────────────────────────────────
    amp_consensus_records = []
    all_ambiguous_positions = set()   # aggregate across all amplicons

    for amp_key in sorted_amp_keys:
        reads = amp_groups[amp_key]
        fw_reads = [(lbl, rec) for lbl, rec in reads if lbl.upper().endswith('FW')]
        rv_reads = [(lbl, rec) for lbl, rec in reads if lbl.upper().endswith('RV')]

        log.info(f"\n  [Amplicon {amp_key}]  "
                 f"FW reads: {[l for l,_ in fw_reads]}  "
                 f"RV reads: {[l for l,_ in rv_reads]}")

        if not fw_reads and not rv_reads:
            log.warning(f"  Amplicon {amp_key}: no usable reads — skipping.")
            continue
        if not fw_reads:
            log.warning(f"  Amplicon {amp_key}: no FW read — consensus from RV only.")
        if not rv_reads:
            log.warning(f"  Amplicon {amp_key}: no RV read — consensus from FW only.")

        amp_pileup            = [[] for _ in range(ref_len)]
        amp_pileup_qual       = [[] for _ in range(ref_len)]
        amp_pileup_candidates = [[] for _ in range(ref_len)]
        amp_ambig_positions   = set()

        for lbl, rec in fw_reads + rv_reads:
            _deposit_read_into_pileup(lbl, rec, ref_str,
                                      amp_pileup, amp_pileup_qual,
                                      amp_pileup_candidates, aligner,
                                      amp_ambig_positions)

        amp_seq, _, amp_only_ambig = _pileup_to_consensus(
            amp_pileup, amp_pileup_candidates, amp_ambig_positions
        )

        # Positions where BOTH strands were ambiguous → flag for report
        all_ambiguous_positions |= amp_only_ambig

        # Build per-position candidate set summary for the amplicon consensus.
        # At each position, take the intersection of all IUPAC candidate sets
        # that were deposited at that position (same logic as the voter above).
        # This is stored so level-2 can use it for cross-amplicon intersection.
        amp_cand_per_pos = []
        for pos_cands in amp_pileup_candidates:
            sets = [cs for cs in pos_cands if cs is not None]
            if sets:
                inter = sets[0]
                for cs in sets[1:]:
                    inter = inter & cs
                amp_cand_per_pos.append(inter if inter else None)
            else:
                amp_cand_per_pos.append(None)

        amp_rec = SeqRecord(
            Seq(amp_seq), id=amp_key,
            description=(f"amplicon_consensus "
                         f"fw={[l for l,_ in fw_reads]} "
                         f"rv={[l for l,_ in rv_reads]}")
        )
        amp_rec.letter_annotations["phred_quality"] = [
            int(sum(q) / len(q)) if q else 0
            for q in amp_pileup_qual
        ]
        # Store candidate sets as a plain annotation (not letter_annotations,
        # which only supports numeric types in BioPython).
        amp_rec.annotations["iupac_cand_per_pos"] = amp_cand_per_pos
        amp_consensus_records.append((amp_key, amp_rec))
        log.info(f"  Amplicon {amp_key} consensus built "
                 f"({amp_seq.count('N')} N positions, "
                 f"{len(amp_only_ambig)} ambiguous-only positions).")

    if not amp_consensus_records:
        log.error("  No amplicon consensuses could be built — returning all-N contig.")
        empty = "N" * ref_len
        stats = {
            "ref_length": ref_len, "covered_positions": 0,
            "coverage_pct": "0.0%", "multi_covered": 0,
            "uncovered_regions": [(1, ref_len)],
            "pileup_qual": [[] for _ in range(ref_len)],
            "ambiguous_positions": set(),
        }
        return empty, stats, ["No reads could be assembled."]

    # ── LEVEL 2: merge amplicon consensuses into final contig ─────────────────
    log.info(f"\n  [Level 2] Merging {len(amp_consensus_records)} amplicon consensuses...")

    final_pileup            = [[] for _ in range(ref_len)]
    final_pileup_qual       = [[] for _ in range(ref_len)]
    final_pileup_candidates = [[] for _ in range(ref_len)]
    final_ambig_positions   = set()

    for amp_key, amp_rec in amp_consensus_records:
        amp_seq  = str(amp_rec.seq)
        amp_qual = list(amp_rec.letter_annotations["phred_quality"])

        amp_cands = amp_rec.annotations.get("iupac_cand_per_pos", [None] * ref_len)
        for ref_pos, (base, qual) in enumerate(zip(amp_seq, amp_qual)):
            if base.upper() == 'N':
                continue   # uncovered by this amplicon
            final_pileup[ref_pos].append(base)
            final_pileup_qual[ref_pos].append(qual)
            final_pileup_candidates[ref_pos].append(
                amp_cands[ref_pos] if ref_pos < len(amp_cands) else None
            )
            if base.islower():
                final_ambig_positions.add(ref_pos)

        covered_by_amp = sum(1 for b in amp_seq if b.upper() != 'N')
        log.info(f"    {amp_key}: contributed {covered_by_amp} positions to final pileup.")

    final_consensus, uncovered_regions, ambig_only_final = _pileup_to_consensus(
        final_pileup, final_pileup_candidates, final_ambig_positions
    )

    # Union: positions ambiguous at level-1 OR still ambiguous after level-2 vote
    all_ambiguous_positions |= ambig_only_final

    covered      = sum(1 for b in final_pileup if b)
    multi        = sum(1 for b in final_pileup if len(b) > 1)
    coverage_pct = f"{covered / ref_len * 100:.1f}%"

    coverage_stats = {
        "ref_length":          ref_len,
        "covered_positions":   covered,
        "coverage_pct":        coverage_pct,
        "multi_covered":       multi,
        "uncovered_regions":   uncovered_regions,
        "ambiguous_positions": all_ambiguous_positions,  # NEW: IUPAC/ambiguous flagged positions
        "pileup_qual":         final_pileup_qual,
        "amplicons_used":      sorted_amp_keys,
        "amp_fw_rv": {
            k: {
                "fw": [l for l, _ in amp_groups[k] if l.upper().endswith('FW')],
                "rv": [l for l, _ in amp_groups[k] if l.upper().endswith('RV')],
            }
            for k in sorted_amp_keys
        },
    }

    log.info(f"\n  Final contig: {covered}/{ref_len} positions covered ({coverage_pct}), "
             f"{multi} with multi-amplicon confirmation, "
             f"{len(all_ambiguous_positions)} ambiguous-only positions.")

    warnings = [
        f"No coverage at reference positions {s}–{e} ({e - s + 1} bp)"
        for s, e in uncovered_regions
    ]
    if all_ambiguous_positions:
        warnings.append(
            f"{len(all_ambiguous_positions)} positions resolved from IUPAC ambiguous bases only "
            f"(reference used as placeholder; mutations cannot be confirmed at these positions)."
        )
    if len(final_consensus) % 3 != 0:
        warnings.append(f"Contig length {len(final_consensus)} bp not divisible by 3.")

    return final_consensus, coverage_stats, warnings


def step2_assemble_contigs(sample_amplicons, contig_folder, ref_record):
    """
    Builds two-level pileup consensus for all samples and saves each as a .fa file.
    Writes .qual and .meta.json sidecars to a hidden folder (.contig_sidecars/)
    so that the contigs/ directory contains only .fa files.
    Returns: {sample_name: contig_filepath}
    """
    import json
    os.makedirs(contig_folder, exist_ok=True)
    sidecar_folder = os.path.join("output", ".contig_sidecars",
                                  os.path.basename(contig_folder))
    os.makedirs(sidecar_folder, exist_ok=True)
    contig_paths = {}

    for sample_name, amp_tuples in sample_amplicons.items():
        amp_ids = [a for a, _ in amp_tuples]
        log.info(f"\nAssembling: {sample_name} | reads: {amp_ids}")

        contig_seq, cov_stats, warnings = build_pileup_consensus(amp_tuples, ref_record)

        for w in warnings:
            log.warning(f"  WARNING: {w}")

        contig_record = SeqRecord(
            Seq(contig_seq), id=sample_name,
            description=(f"2level_pileup_consensus "
                         f"amplicons={'+'.join(cov_stats.get('amplicons_used', amp_ids))} "
                         f"coverage={cov_stats['coverage_pct']} "
                         f"len={len(contig_seq)}")
        )
        contig_path = os.path.join(contig_folder, f"{sample_name}_contig.fa")
        SeqIO.write(contig_record, contig_path, "fasta")

        # Per-position mean Phred quality sidecar (written to hidden sidecar folder)
        qual_path = os.path.join(sidecar_folder, f"{sample_name}_contig.qual")
        mean_quals = [
            round(sum(q) / len(q), 1) if q else 0.0
            for q in cov_stats["pileup_qual"]
        ]
        with open(qual_path, "w") as qf:
            qf.write("\n".join(str(q) for q in mean_quals))

        # Assembly metadata sidecar (written to hidden sidecar folder)
        ambig_sorted = sorted(cov_stats.get("ambiguous_positions", set()))
        meta = {
            "amplicons_used":      cov_stats.get("amplicons_used", []),
            "amp_fw_rv":           cov_stats.get("amp_fw_rv", {}),
            "ambiguous_positions": ambig_sorted,   # NEW: 1-based ref positions
        }
        meta_path = os.path.join(sidecar_folder, f"{sample_name}_contig.meta.json")
        with open(meta_path, "w") as mf:
            json.dump(meta, mf, indent=2)

        log.info(f"  OK {sample_name}_contig.fa ({len(contig_seq)} bp, "
                 f"{cov_stats['coverage_pct']} covered, "
                 f"{len(cov_stats.get('amplicons_used', []))} amplicons, "
                 f"{len(ambig_sorted)} ambiguous positions)")
        contig_paths[sample_name] = contig_path

    return contig_paths


# =============================================================================
# STEPS 3–5 — Alignment, Mutation Calling, Reporting (replaces MEGA)
# =============================================================================

def analyze_single_sample(ref_record, contig_filepath, strain_label):
    filename = os.path.basename(contig_filepath)
    error_metrics = {
        "DNA_Identity(%)":        "ERROR",
        "AA_Identity(%)":         "ERROR",
        "Total_Mutations":        "ERROR",
        "Frameshifts_Fixed":      "ERROR",
        "Ambiguous_Positions":    "ERROR",
        "Mutation_List":          "ERROR",
    }

    try:
        records = list(SeqIO.parse(contig_filepath, "fasta"))
        if not records:
            raise ValueError("Empty FASTA file")
        contig_record = records[0]
    except Exception as e:
        return f"Error loading {filename}: {e}\n", error_metrics, None, None

    _contig_dir  = os.path.dirname(contig_filepath)
    _sidecar_dir = os.path.join("output", ".contig_sidecars",
                                 os.path.basename(_contig_dir))
    _base        = os.path.splitext(os.path.basename(contig_filepath))[0]
    qual_path = os.path.join(_sidecar_dir, f"{_base}.qual")
    if os.path.exists(qual_path):
        with open(qual_path) as qf:
            contig_qual = [float(x.strip()) for x in qf if x.strip()]
    else:
        contig_qual = [20.0] * len(contig_record.seq)
        log.warning(f"No .qual sidecar found for {filename} — using default Q20.")

    import json as _json
    meta_path = os.path.join(_sidecar_dir, f"{_base}.meta.json")
    if os.path.exists(meta_path):
        with open(meta_path) as mf:
            assembly_meta = _json.load(mf)
    else:
        assembly_meta = {"amplicons_used": [], "amp_fw_rv": {}, "ambiguous_positions": []}
        log.warning(f"No .meta.json sidecar found for {filename}.")

    # Load ambiguous positions (1-based) — these are where IUPAC bases were the
    # only evidence; reference was used as placeholder.
    ambiguous_ref_positions = set(assembly_meta.get("ambiguous_positions", []))

    global_aligner = make_aligner('global')
    local_aligner  = make_aligner('local', open_gap=-10)

    fwd_score = getattr(next(iter(local_aligner.align(ref_record.seq, contig_record.seq)), None), 'score', 0)
    rev_score = getattr(next(iter(local_aligner.align(ref_record.seq, contig_record.seq.reverse_complement())), None), 'score', 0)

    if fwd_score == 0 and rev_score == 0:
        return f"ERROR: Could not align {filename} to reference.\n", error_metrics, None, None

    use_revcomp = rev_score > fwd_score
    working_seq = contig_record.seq.reverse_complement() if use_revcomp else contig_record.seq

    if use_revcomp:
        log.info(f"{filename}: Better alignment on reverse strand.")

    working_str  = str(working_seq)
    seq_stripped = working_str.strip('Nn')
    leading_ns   = len(working_str) - len(working_str.lstrip('Nn'))
    trailing_ns  = len(working_str) - len(working_str.rstrip('Nn'))

    seq_to_align = seq_stripped   # no collapsing — full length preserved

    # ── Semi-global aligner: Sanger-appropriate gap penalties ────────────────
    # A Sanger read is a physically continuous piece of DNA. Any gap in the
    # *sample* side of the alignment is almost certainly an alignment artefact
    # caused by low-quality bases, not a real deletion. We therefore set a very
    # high internal gap open penalty (-20) so the aligner strongly prefers
    # accepting mismatches over introducing a gap.  End gaps are free (semi-
    # global behaviour) so N-padded ends don't distort the score.
    semi_global = PairwiseAligner()
    semi_global.mode             = 'global'
    semi_global.match_score      =  2
    semi_global.mismatch_score   = -1
    semi_global.open_gap_score   = ALIGN_GAP_OPEN    # see PIPELINE PARAMETERS
    semi_global.extend_gap_score = ALIGN_GAP_EXTEND
    semi_global.query_end_gap_score  = 0.0
    semi_global.target_end_gap_score = 0.0

    best_dna = next(iter(semi_global.align(ref_record.seq, Seq(seq_to_align))), None)
    if best_dna is None:
        return f"ERROR: Alignment failed for {filename}.\n", error_metrics, None, None

    aligned_ref_dna  = best_dna[0]
    aligned_samp_dna = best_dna[1]

    log.info(f"{filename}: aligned (stripped {leading_ns} leading / {trailing_ns} trailing Ns), "
             f"sequenced region = {len(seq_to_align)} bp vs ref {len(ref_record.seq)} bp")

    # ─── Reference-fill ──────────────────────────────────────────────────────
    # Replace N/gap positions in the aligned sample with the corresponding ref
    # base (lowercase = ref-filled, not sequenced).
    # NEW: Also replace positions that were flagged as 'ambiguous' in the pileup
    # (only IUPAC bases from both strands). These get a distinct tag so the
    # report can differentiate them from genuinely uncovered positions.
    filled_samp_chars   = []
    ref_filled_dna_regions   = []    # (start, end) 1-based — uncovered (N/gap)
    ambig_filled_dna_regions = []    # (start, end) 1-based — IUPAC-only, ref-filled

    in_n_run      = False
    in_ambig_run  = False
    n_run_start   = 0
    ambig_start   = 0
    _scan_pos     = 1

    for r, s in zip(aligned_ref_dna, aligned_samp_dna):
        if r == '-':
            filled_samp_chars.append(s)
            continue

        is_ambig_position = _scan_pos in ambiguous_ref_positions

        if s.upper() == 'N' or s == '-':
            # Genuinely uncovered position
            filled_samp_chars.append(r.lower())
            if in_ambig_run:
                ambig_filled_dna_regions.append((ambig_start, _scan_pos - 1))
                in_ambig_run = False
            if not in_n_run:
                in_n_run, n_run_start = True, _scan_pos
        elif is_ambig_position:
            # IUPAC-ambiguous position: use lowercase ref as placeholder
            # but track separately from N-gap positions
            filled_samp_chars.append(r.lower())
            if in_n_run:
                ref_filled_dna_regions.append((n_run_start, _scan_pos - 1))
                in_n_run = False
            if not in_ambig_run:
                in_ambig_run, ambig_start = True, _scan_pos
        else:
            if in_n_run:
                ref_filled_dna_regions.append((n_run_start, _scan_pos - 1))
                in_n_run = False
            if in_ambig_run:
                ambig_filled_dna_regions.append((ambig_start, _scan_pos - 1))
                in_ambig_run = False
            filled_samp_chars.append(s)
        _scan_pos += 1

    if in_n_run:
        ref_filled_dna_regions.append((n_run_start, _scan_pos - 1))
    if in_ambig_run:
        ambig_filled_dna_regions.append((ambig_start, _scan_pos - 1))

    aligned_samp_dna_filled = "".join(filled_samp_chars)

    non_gap_ref  = sum(1 for r in aligned_ref_dna if r != '-')
    dna_matches  = sum(
        1 for r, s in zip(aligned_ref_dna, aligned_samp_dna_filled)
        if r != '-' and s.upper() == r.upper() and s == s.upper() and s.upper() in 'ATGC'
    )
    dna_identity = (dna_matches / non_gap_ref * 100) if non_gap_ref > 0 else 0

    # ─── Build corrected CDS (frame-safe, gap-stripped) ──────────────────────
    #
    # Pass 1: walk the alignment and produce an initial corrected_chars list.
    #   - Matches / mismatches  → keep sample base as-is
    #   - Deletion in sample    → fill with lowercase ref base (length preserved)
    #   - Insertion in sample   → drop (would cause +1 frameshift)
    #
    corrected_chars = []
    deletions_log   = []
    insertions_log  = []
    dna_ref_pos     = 1

    for r, s in zip(aligned_ref_dna, aligned_samp_dna_filled):
        if r != '-' and s != '-':
            corrected_chars.append(s)
            dna_ref_pos += 1
        elif r != '-' and s == '-':
            corrected_chars.append(r.lower())
            deletions_log.append(
                f"  - Deletion at Reference DNA Position {dna_ref_pos} "
                f"(ref base '{r}' inserted as placeholder)"
            )
            dna_ref_pos += 1
        elif r == '-' and s != '-':
            pos_label = ("before Position 1" if dna_ref_pos == 1
                         else f"after Reference DNA Position {dna_ref_pos - 1}")
            insertions_log.append(
                f"  - Extra base '{s}' inserted {pos_label} — removed to preserve frame"
            )

    corrected_cds = Seq("".join(corrected_chars))

    # ─── Pass 2: quality-aware minimum-mutation frameshift repair ────────────
    #
    # SCORING  (lower tuple = better, compared lexicographically)
    # ─────────
    #   (stops, mutations, neg_quality)
    #
    #   stops       — internal stop codons  (primary: must reach 0)
    #   mutations   — AA mismatches vs ref  (secondary: minimise)
    #   neg_quality — negated mean Phred at operated positions so that
    #                 LOW-quality positions (read ends, likely artefacts)
    #                 are preferred when stops and mutations are tied.
    #
    # SPEED OPTIMISATION
    # ──────────────────
    # The naive pair search is O(N²·T) where T = translation time (~0.5 ms).
    # For a 1725 bp CDS that is ~3M trials × 0.5 ms = ~25 minutes — unusable.
    #
    # Key insight: a single-base frameshift shifts EVERY codon downstream of
    # the indel position.  Therefore only operations near the FIRST internal
    # stop codon can restore the frame.  We restrict the search window to
    # [0, first_stop_codon_dna_pos] — typically <700 bp — which cuts the
    # pair search to ~250k trials and completes in a few seconds.
    #
    # Additionally, for the pair search we only try positions that produced
    # a stop-count improvement in the single-op scan (candidate positions),
    # further pruning the search space.
    #
    frameshift_repair_log = []
    ref_str_full = str(ref_record.seq)
    ref_len_full = len(ref_str_full)
    ref_aa_full  = str(ref_record.seq.translate(to_stop=False))

    # Quality array: per-reference-position mean Phred (0-based).
    _cq = list(contig_qual) if len(contig_qual) >= ref_len_full \
          else list(contig_qual) + [0.0] * (ref_len_full - len(contig_qual))

    def _neg_qual(i, j=None):
        """Negated mean Phred at position i (and optionally j). Lower = worse quality."""
        if j is None:
            return -_cq[i] if i < len(_cq) else 0.0
        qi = _cq[i] if i < len(_cq) else 0.0
        qj = _cq[j] if j < len(_cq) else 0.0
        return -(qi + qj) / 2.0

    def _score_cds(chars):
        """Return (n_internal_stops, n_aa_mutations)."""
        s = "".join(chars).upper()
        rem = len(s) % 3
        if rem:
            s += ref_str_full[len(s): len(s) + (3 - rem)]
        aa = str(Seq(s[:ref_len_full]).translate(to_stop=False))
        stops = sum(1 for k, a in enumerate(aa[:-1]) if a == "*")
        muts  = sum(1 for ra, sa in zip(ref_aa_full, aa) if ra != sa and sa != "*")
        return stops, muts

    def _apply_del(chars, i):
        t = chars[:i] + chars[i+1:]
        pad = len(chars) - len(t)
        t += list(ref_str_full[len(t): len(t) + pad].lower())
        return (t + list(ref_str_full[len(t):].lower()))[:len(chars)]

    def _apply_ins(chars, i):
        if i >= ref_len_full:
            return None
        t = chars[:i] + [ref_str_full[i].lower()] + chars[i: len(chars) - 1]
        return (t + list(ref_str_full[len(t):].lower()))[:len(chars)]

    cds_len       = len(corrected_chars)
    stops0, muts0 = _score_cds(corrected_chars)

    if stops0 > 0:
        # ── Multi-strategy frameshift repair ─────────────────────────────────
        #
        # STRATEGY OVERVIEW
        # ─────────────────
        # A frameshift in a Sanger read is almost always caused by a single
        # indel artefact — but sometimes two nearby low-quality bases cause
        # a net +2 / -2 / +1-1 compound shift.  We therefore try three
        # strategies in order, stopping as soon as any one reaches 0 stops:
        #
        #   Strategy 1 — WINDOW SCAN  (fast, most targeted)
        #     Scan a window of ±60 bp around EACH internal stop codon.
        #     For each position in those windows, try DEL and INS.
        #     Rationale: the frameshift must be upstream of (or at) the stop.
        #     Searching the upstream window is far more targeted than the
        #     whole CDS and catches the vast majority of single-indel cases.
        #
        #   Strategy 2 — LOW-QUALITY PAIRS  (handles compound frameshifts)
        #     Take the 20 lowest-quality positions in the CDS.
        #     Try all pairs (DEL+DEL, INS+INS, DEL+INS) at those positions.
        #     This handles +2/-2 and ±1 compound shifts that Strategy 1 misses.
        #     Search space: 20×19/2 × 4 = ~760 trials — fast.
        #
        #   Strategy 3 — FULL SINGLE SCAN  (exhaustive fallback)
        #     Try DEL and INS at every position in the CDS.
        #     Only reached when Strategies 1 and 2 both fail.
        #     Completes in ~seconds for a 1725 bp CDS.
        #
        # SCORING  (lexicographic, lower = better)
        #   (n_stops, n_aa_mutations, neg_mean_quality_at_operated_positions)
        #   Stops are the primary objective — must reach 0.
        #   AA mutations are secondary — minimise.
        #   Quality is the tiebreaker — prefer operating on low-Q positions.
        #
        log.warning(
            f"{filename}: {stops0} internal stop(s), {muts0} AA mutation(s) — "
            f"starting multi-strategy frameshift repair."
        )

        # Find positions of all internal stop codons (0-based AA positions).
        def _stop_positions(chars):
            s = "".join(chars).upper()
            rem = len(s) % 3
            if rem:
                s += ref_str_full[len(s): len(s) + (3 - rem)]
            aa = str(Seq(s[:ref_len_full]).translate(to_stop=False))
            return [k for k, a in enumerate(aa[:-1]) if a == "*"]

        best = {"score": (stops0, muts0, 0.0), "chars": None, "ops": ""}

        def _try(trial, op_desc, q_score):
            if trial is None:
                return
            s, m = _score_cds(trial)
            score = (s, m, q_score)
            if score < best["score"]:
                best["score"] = score
                best["chars"] = trial
                best["ops"]   = op_desc

        # ── Strategy 1: window scan around each stop codon ──────────────────
        WINDOW = 60   # bp upstream of each stop to scan
        stop_aa_positions = _stop_positions(corrected_chars)
        window_positions  = set()
        for aa_pos in stop_aa_positions:
            dna_pos = aa_pos * 3
            lo = max(0, dna_pos - WINDOW)
            window_positions.update(range(lo, dna_pos + 1))
        window_positions = sorted(window_positions)

        log.info(f"  Strategy 1: scanning {len(window_positions)} positions "
                 f"around {len(stop_aa_positions)} stop codon(s).")
        for i in window_positions:
            q = _neg_qual(i)
            _try(_apply_del(corrected_chars, i), f"DEL@{i+1}(Q={_cq[i] if i<len(_cq) else 0:.0f})", q)
            _try(_apply_ins(corrected_chars, i), f"INS@{i+1}(Q={_cq[i] if i<len(_cq) else 0:.0f})", q)

        if best["score"][0] == 0:
            log.info(f"  Strategy 1 succeeded: {best['ops']}")
        else:
            # ── Strategy 2: low-quality pairs ────────────────────────────────
            N_LOW = 20
            ranked   = sorted(range(cds_len), key=lambda i: _cq[i] if i < len(_cq) else 0.0)
            low_q    = ranked[:N_LOW]
            log.info(f"  Strategy 2: trying pairs among {N_LOW} lowest-quality positions.")

            for idx_a, i in enumerate(low_q):
                for j in low_q[idx_a + 1:]:
                    q2 = _neg_qual(i, j)
                    # DEL + DEL  (-2 net)
                    t = _apply_del(corrected_chars, i)
                    if t: t = _apply_del(t, j if j < i else j - 1)
                    _try(t, f"DEL@{i+1}+DEL@{j+1}", q2)
                    # INS + INS  (+2 net)
                    t = _apply_ins(corrected_chars, i)
                    if t: t = _apply_ins(t, j if j > i else j)
                    _try(t, f"INS@{i+1}+INS@{j+1}", q2)
                    # DEL + INS  (net 0 — fixes a substitution-induced stop)
                    t = _apply_del(corrected_chars, i)
                    if t: t = _apply_ins(t, j if j < i else j - 1)
                    _try(t, f"DEL@{i+1}+INS@{j+1}", q2)
                    # INS + DEL  (net 0, opposite order)
                    t = _apply_ins(corrected_chars, i)
                    if t: t = _apply_del(t, j if j > i else j)
                    _try(t, f"INS@{i+1}+DEL@{j+1}", q2)

            if best["score"][0] == 0:
                log.info(f"  Strategy 2 succeeded: {best['ops']}")
            else:
                # ── Strategy 3: full CDS single scan ─────────────────────────
                log.info(f"  Strategy 3: exhaustive single-op scan of all {cds_len} positions.")
                for i in range(cds_len):
                    q = _neg_qual(i)
                    _try(_apply_del(corrected_chars, i),
                         f"DEL@{i+1}(Q={_cq[i] if i<len(_cq) else 0:.0f})", q)
                    _try(_apply_ins(corrected_chars, i),
                         f"INS@{i+1}(Q={_cq[i] if i<len(_cq) else 0:.0f})", q)

                if best["score"][0] == 0:
                    log.info(f"  Strategy 3 succeeded: {best['ops']}")
                else:
                    log.warning(f"  All three strategies failed — residual stops remain.")

        # ── Apply the best repair found (if it eliminates all stops) ─────────
        if best["chars"] is not None and best["score"][0] == 0:
            corrected_chars = best["chars"]
            corrected_cds   = Seq("".join(corrected_chars))
            desc = (
                f"  ✓ FRAMESHIFT REPAIRED  [{best['ops']}]  "
                f"→ {best['score'][1]} AA mutation(s) remain  — MANUAL CONFIRMATION RECOMMENDED"
            )
            frameshift_repair_log.append(desc)
            log.warning(
                f"{filename}: frameshift repaired ({best['ops']}) "
                f"→ {best['score'][1]} AA mutation(s). Manual confirmation recommended."
            )
        else:
            # Partial improvement: report best partial result if it reduced stops
            partial_desc = ""
            if best["chars"] is not None and best["score"][0] < stops0:
                corrected_chars = best["chars"]
                corrected_cds   = Seq("".join(corrected_chars))
                partial_desc = (
                    f"  ⚠ PARTIAL REPAIR  [{best['ops']}]  "
                    f"reduced stops {stops0} → {best['score'][0]}  "
                    f"but could not reach 0 — manual review required."
                )
                frameshift_repair_log.append(partial_desc)
                log.warning(
                    f"{filename}: partial repair only — stops reduced "
                    f"{stops0} → {best['score'][0]}. Manual review needed."
                )
            else:
                frameshift_repair_log.append(
                    f"  ! REPAIR FAILED: all three strategies exhausted "
                    f"({stops0} stop(s) remain). Manual review required.\n"
                    f"    Hint: check the chromatogram around DNA positions "
                    f"{stop_aa_positions[0]*3 - 60}–{stop_aa_positions[-1]*3 + 10} "
                    f"for a noisy insertion or deletion artefact."
                )
                log.warning(f"{filename}: all repair strategies failed — manual review needed.")
    else:
        log.info(f"{filename}: no internal stop codons — frame is clean.")


    # ─── Build filled-position tracking sets from (possibly repaired) chars ──
    ref_filled_dna_pos = {
        i + 1
        for i, ch in enumerate(corrected_chars)
        if ch.islower() and (i + 1) not in ambiguous_ref_positions
    }
    ambig_filled_dna_pos = {
        i + 1
        for i, ch in enumerate(corrected_chars)
        if ch.islower() and (i + 1) in ambiguous_ref_positions
    }

    all_filled_dna_pos   = ref_filled_dna_pos | ambig_filled_dna_pos
    ref_filled_aa_pos    = {(p - 1) // 3 + 1 for p in ref_filled_dna_pos}
    ambig_filled_aa_pos  = {(p - 1) // 3 + 1 for p in ambig_filled_dna_pos}
    all_filled_aa_pos    = ref_filled_aa_pos | ambig_filled_aa_pos

    try:
        ref_aa      = ref_record.seq.translate(to_stop=False)
        raw_samp_aa = corrected_cds.upper().translate(to_stop=False)
        samp_aa     = Seq(str(raw_samp_aa))
    except Exception as e:
        error_metrics["DNA_Identity(%)"] = f"{dna_identity:.2f}"
        return f"Warning during translation: {e}\n", error_metrics, Seq(str(corrected_cds)), None

    best_aa         = next(iter(global_aligner.align(ref_aa, samp_aa)))
    aligned_ref_aa  = best_aa[0]
    aligned_samp_aa = best_aa[1]

    non_gap_ref_aa = sum(1 for r in aligned_ref_aa if r != '-')
    aa_matches = 0
    _aa_pos = 1
    for r, s in zip(aligned_ref_aa, aligned_samp_aa):
        if r != '-':
            if r == s and _aa_pos not in all_filled_aa_pos:
                aa_matches += 1
            _aa_pos += 1
    aa_identity = (aa_matches / non_gap_ref_aa * 100) if non_gap_ref_aa > 0 else 0

    # Build mapping: AA position (1-based) -> mean Phred of its 3 codon bases
    # Also track which AA positions contain at least one IUPAC-resolved (lowercase)
    # base in the corrected CDS — these are single-strand-resolved calls and need
    # an extra manual-check warning even though they are not fully ref-filled.
    # aa_iupac_details maps AA pos -> list of (dna_pos_1based, original_base_in_contig)
    # for each IUPAC-resolved base in that codon.
    aa_qual_map      = {}
    aa_has_iupac_resolved = set()
    aa_iupac_details = {}   # AA pos -> [(dna_pos, contig_base), ...]
    _ref_dna_pos = 1
    for r, s in zip(aligned_ref_dna, aligned_samp_dna_filled):
        if r != '-':
            aa_pos_for_base = (_ref_dna_pos - 1) // 3 + 1
            if _ref_dna_pos not in all_filled_dna_pos and _ref_dna_pos - 1 < len(contig_qual):
                if aa_pos_for_base not in aa_qual_map:
                    aa_qual_map[aa_pos_for_base] = []
                aa_qual_map[aa_pos_for_base].append(contig_qual[_ref_dna_pos - 1])
            # Track IUPAC-resolved bases: these are lowercase in corrected_chars
            # but NOT in all_filled_dna_pos (which covers ref-filled / ambig-filled).
            # A lowercase base that is not ref-filled means it was resolved from a
            # single-strand IUPAC code (best-peak estimate, lower confidence).
            cds_idx = _ref_dna_pos - 1
            if (cds_idx < len(corrected_chars)
                    and corrected_chars[cds_idx].islower()
                    and _ref_dna_pos not in all_filled_dna_pos):
                aa_has_iupac_resolved.add(aa_pos_for_base)
                aa_iupac_details.setdefault(aa_pos_for_base, []).append(
                    (f"DNA pos {_ref_dna_pos}", corrected_chars[cds_idx])
                )
            _ref_dna_pos += 1

    def phred_to_confidence(phred_scores):
        import math
        if not phred_scores:
            return "UNKNOWN  (no quality data)"

        codon_probs = []
        for i in range(0, len(phred_scores), 3):
            codon_bases = phred_scores[i:i + 3]
            p_codon = 1.0
            for q in codon_bases:
                p_base_correct = 1.0 - 10 ** (-q / 10.0)
                p_codon *= p_base_correct
            codon_probs.append(p_codon)

        best_prob = max(codon_probs)
        pct = best_prob * 100.0
        n_reads = len(codon_probs)

        if pct >= 99.9:
            level = "VERY HIGH"
        elif pct >= 99.0:
            level = "HIGH"
        elif pct >= 95.0:
            level = "MODERATE"
        else:
            level = "LOW"

        reads_note = f", confirmed by {n_reads} reads" if n_reads > 1 else ""
        return f"{pct:.2f}%  [{level}{reads_note}]"

    aa_substitutions  = []
    concise_mutations = []
    ambiguous_flags   = []   # mutations at IUPAC-ambiguous positions
    residual_stops    = []   # any stop codons that survived repair (should be 0)
    aa_ref_pos        = 1

    for r, s in zip(aligned_ref_aa, aligned_samp_aa):
        if r != '-' and s != '-' and r != s:
            if aa_ref_pos not in ref_filled_aa_pos:
                is_last_pos = (aa_ref_pos == len(ref_aa))
                if s == '*' and not is_last_pos:
                    residual_stops.append(
                        f"  ! AA Position {aa_ref_pos}: internal stop codon NOT repaired — "
                        f"manual review required (check chromatogram around DNA pos "
                        f"{aa_ref_pos * 3 - 60}–{aa_ref_pos * 3})."
                    )
                elif aa_ref_pos in ambig_filled_aa_pos:
                    # Mutation at an IUPAC-ambiguous position — flag, don't report as confirmed.
                    # Compute which specific DNA positions in the codon are ambig-filled.
                    codon_dna_start = (aa_ref_pos - 1) * 3 + 1
                    ambig_dna = [
                        str(p) for p in range(codon_dna_start, codon_dna_start + 3)
                        if p in ambig_filled_dna_pos
                    ]
                    pos_str = f" (ambig DNA pos: {', '.join(ambig_dna)})" if ambig_dna else ""
                    ambiguous_flags.append(
                        f"  ? AA Position {aa_ref_pos}: '{r}' -> '{s}'{pos_str}  |  "
                        f"UNCERTAIN — both strands IUPAC-ambiguous, reference used as placeholder; "
                        f"cannot confirm this mutation — verify chromatogram"
                    )
                else:
                    quals_for_pos = aa_qual_map.get(aa_ref_pos, [])
                    confidence    = phred_to_confidence(quals_for_pos)
                    iupac_detail  = aa_iupac_details.get(aa_ref_pos, [])
                    if iupac_detail:
                        pos_str  = ", ".join(f"{dp} (base '{b}')" for dp, b in iupac_detail)
                        iupac_note = (
                            f"⚠ IUPAC-resolved base(s) in codon: {pos_str} — "
                            f"single-strand peak-height estimate; verify chromatogram"
                        )
                    else:
                        iupac_note = ""
                    entry = (
                        f"  - AA Position {aa_ref_pos}: '{r}' -> '{s}'  |  Confidence: {confidence}"
                    )
                    if iupac_note:
                        entry += f"\n    {iupac_note}"
                    aa_substitutions.append(entry)
                    concise_mutations.append(f"{r}{aa_ref_pos}{s}"
                                             + ("*" if aa_ref_pos in aa_has_iupac_resolved else ""))
        if r != '-':
            aa_ref_pos += 1

    # ─── Build report ─────────────────────────────────────────────────────────
    W  = 72
    HR = "═" * W   # thick rule — major section divider
    SH = "─" * W   # thin rule  — sub-section divider

    sample_id = os.path.splitext(filename)[0].replace("_contig", "")

    # ── Pre-compute display values ────────────────────────────────────────────
    n_uncov = len(ref_filled_dna_regions)
    n_ambig_gaps = len(ambig_filled_dna_regions)
    n_confirmed  = len(concise_mutations)
    n_uncertain  = len(ambiguous_flags)
    n_stops      = len(residual_stops)

    if n_uncov == 0 and n_ambig_gaps == 0:
        cov_status = "COMPLETE  ✓  all reference positions covered"
    else:
        parts_cov = []
        if n_uncov      > 0: parts_cov.append(f"{n_uncov} uncovered region(s)")
        if n_ambig_gaps > 0: parts_cov.append(f"{n_ambig_gaps} IUPAC-only region(s)")
        cov_status = "PARTIAL   —  " + ",  ".join(parts_cov)

    mut_list_str = ", ".join(
        m.replace("*", "⚠") for m in concise_mutations
    ) if concise_mutations else "—  none detected"

    del_rows = [l for l in deletions_log  if l.strip()]
    ins_rows = [l for l in insertions_log if l.strip()]

    amp_fw_rv      = assembly_meta.get("amp_fw_rv", {})
    amplicons_used = assembly_meta.get("amplicons_used", [])

    csv_metrics = {
        "DNA_Identity(%)":      f"{dna_identity:.2f}",
        "AA_Identity(%)":       f"{aa_identity:.2f}",
        "Total_Mutations":      n_confirmed,
        "Frameshifts_Fixed":    len(frameshift_repair_log),
        "Residual_Stops":       n_stops,
        "Ambiguous_Positions":  len(ambiguous_ref_positions),
        "Mutation_List":        ", ".join(concise_mutations) if concise_mutations else "None",
    }

    # ═════════════════════════════════════════════════════════════════════════
    # SECTION 1 — SAMPLE OVERVIEW
    # ═════════════════════════════════════════════════════════════════════════
    overview = [
        HR,
        f"  SAMPLE   :  {sample_id}",
        f"  STRAIN   :  {strain_label}",
        f"  REFERENCE:  {ref_record.id}  ({len(ref_record.seq)} bp)",
        SH,
        f"  DNA identity  :  {dna_identity:.2f}%",
        f"  AA  identity  :  {aa_identity:.2f}%",
        f"  Coverage      :  {cov_status}",
        SH,
    ]

    # Quick-read status flags
    if n_stops > 0:
        overview.append(f"  ⛔  {n_stops} unresolved internal stop codon(s) — manual review required")
    if frameshift_repair_log and any("✓" in l for l in frameshift_repair_log):
        overview.append(f"  ⚠   Frameshift repair was applied — manual confirmation recommended")
    elif frameshift_repair_log and any("PARTIAL" in l for l in frameshift_repair_log):
        overview.append(f"  ⚠   Partial frameshift repair only — manual review required")
    elif frameshift_repair_log:
        overview.append(f"  ⚠   Frameshift repair attempted but failed — manual review required")
    if n_uncertain > 0:
        overview.append(f"  ❓  {n_uncertain} mutation(s) could NOT be confirmed (IUPAC-ambiguous bases)")
    if n_confirmed == 0 and n_stops == 0:
        overview.append(f"  ✓   No amino acid mutations detected vs. reference")
    else:
        overview.append(
            f"  ►   AA mutations (confirmed): {n_confirmed}   —   {mut_list_str}"
        )

    overview.append(HR)
    overview_block = "\n".join(overview)

    # ═════════════════════════════════════════════════════════════════════════
    # SECTION 2 — ASSEMBLY DETAILS
    # ═════════════════════════════════════════════════════════════════════════
    asm_lines = [f"\n  [ 1 ]  ASSEMBLY DETAILS\n  {SH}"]

    if amp_fw_rv:
        asm_lines.append("  Amplicons used (FW ▶ forward read  /  RV ◀ reverse read):\n")
        for amp_key in amplicons_used:
            fw     = amp_fw_rv[amp_key]["fw"]
            rv     = amp_fw_rv[amp_key]["rv"]
            fw_str = ", ".join(fw) if fw else "—  (no FW read)"
            rv_str = ", ".join(rv) if rv else "—  (no RV read)"
            asm_lines.append(f"    {amp_key:<6}  ▶  {fw_str:<35}  ◀  {rv_str}")
    else:
        asm_lines.append("  FW + RV pileup consensus  (no per-amplicon metadata available)")

    asm_lines.append(f"\n  Alignment score vs. reference : {best_dna.score}")
    asm_lines.append(
        f"  IUPAC-ambiguous positions     : {len(ambiguous_ref_positions)}  "
        f"(both strands unresolvable — ref used as placeholder)"
    )
    asm_block = "\n".join(asm_lines)

    # ═════════════════════════════════════════════════════════════════════════
    # SECTION 3 — COVERAGE GAPS  (only if gaps exist)
    # ═════════════════════════════════════════════════════════════════════════
    gap_block = ""
    if ref_filled_dna_regions or ambig_filled_dna_regions:
        gl = [f"\n  [ 2 ]  COVERAGE GAPS\n  {SH}",
              "  The reference sequence was used as a placeholder in these regions.",
              "  Mutations inside these regions CANNOT be confirmed.\n",
              "  Symbol key:  ✗ = no reads at all"
              "                ? = reads present but both strands IUPAC-ambiguous\n"]
        if ref_filled_dna_regions:
            gl.append("  Uncovered regions (no reads):")
            for s, e in ref_filled_dna_regions:
                gl.append(f"    ✗  DNA pos {s:>5}–{e:<5}  ({e-s+1:>3} bp)")
        if ambig_filled_dna_regions:
            if ref_filled_dna_regions:
                gl.append("")
            gl.append("  IUPAC-ambiguous regions (reads present but unresolvable):")
            for s, e in ambig_filled_dna_regions:
                gl.append(f"    ?  DNA pos {s:>5}–{e:<5}  ({e-s+1:>3} bp)")
        gap_block = "\n".join(gl)

    # ═════════════════════════════════════════════════════════════════════════
    # SECTION 4 — FRAMESHIFT REPAIR  (only if attempted)
    # ═════════════════════════════════════════════════════════════════════════
    repair_block = ""
    if frameshift_repair_log:
        rl = [f"\n  [ 3 ]  FRAMESHIFT REPAIR\n  {SH}",
              "  A frameshift is a shift in the reading frame caused by an insertion",
              "  or deletion artefact in the sequencing read.  The pipeline attempts",
              "  to detect and correct these automatically.\n"]
        rl += frameshift_repair_log
        if any("✓" in l or "REPAIRED" in l for l in frameshift_repair_log):
            rl.append("\n  ⚠  Any auto-repair is a best-guess correction based on quality scores.")
            rl.append("     Always verify the corrected sequence in the chromatogram.")
        repair_block = "\n".join(rl)

    # ═════════════════════════════════════════════════════════════════════════
    # SECTION 5 — UNRESOLVED STOP CODONS  (only if present)
    # ═════════════════════════════════════════════════════════════════════════
    stop_block = ""
    if residual_stops:
        sl = [f"\n  [ 4 ]  UNRESOLVED INTERNAL STOP CODONS  ⛔\n  {SH}",
              "  These stop codons could NOT be removed by the repair algorithm.",
              "  They almost certainly represent sequencing artefacts, not real biology.",
              "  The AA sequence downstream of each stop is unreliable.\n"]
        sl += residual_stops
        sl.append("\n  Action required: inspect the chromatogram at the DNA positions")
        sl.append("  listed above and check for a noisy insertion or deletion.")
        stop_block = "\n".join(sl)

    # ═════════════════════════════════════════════════════════════════════════
    # SECTION 6 — DNA INDELS  (only if present)
    # ═════════════════════════════════════════════════════════════════════════
    indel_block = ""
    if del_rows or ins_rows:
        il = [f"\n  [ 5 ]  DNA-LEVEL INDELS  (detected during alignment)\n  {SH}",
              "  These insertions / deletions were found when aligning your contig",
              "  to the reference.  Each one has been corrected to preserve the",
              "  reading frame.  If any of these are real (not artefacts), the",
              "  downstream AA sequence may be incorrect — verify in chromatogram.\n"]
        if del_rows:
            il.append("  Deletions  (reference base inserted as placeholder to fill the gap):")
            il += del_rows
        if ins_rows:
            if del_rows: il.append("")
            il.append("  Insertions  (extra base removed to preserve reading frame):")
            il += ins_rows
        indel_block = "\n".join(il)

    # ═════════════════════════════════════════════════════════════════════════
    # SECTION 7 — CONFIRMED AA MUTATIONS
    # ═════════════════════════════════════════════════════════════════════════
    mut_sec_num = sum(1 for b in [gap_block, repair_block, stop_block, indel_block] if b) + 2
    mut_lines = [
        f"\n  [ {mut_sec_num} ]  CONFIRMED AMINO ACID MUTATIONS\n  {SH}",
        "  These positions differ from the reference and are covered by confident reads.",
        "  Confidence = probability that all 3 bases of the codon are correctly read",
        "  (calculated from Phred quality scores of the sequencing reads).",
        "  ⚠ after a mutation code = codon contains an IUPAC-resolved base; verify.\n",
    ]
    if aa_substitutions:
        mut_lines += aa_substitutions
    else:
        mut_lines.append("  — No confirmed amino acid mutations vs. reference.")
    mut_block = "\n".join(mut_lines)

    # ═════════════════════════════════════════════════════════════════════════
    # SECTION 8 — UNCERTAIN AA MUTATIONS  (only if present)
    # ═════════════════════════════════════════════════════════════════════════
    unc_block = ""
    if ambiguous_flags:
        unc_sec_num = mut_sec_num + 1
        ul = [
            f"\n  [ {unc_sec_num} ]  UNCERTAIN AMINO ACID MUTATIONS  (require manual verification)\n  {SH}",
            "  These positions appear to differ from the reference, BUT at least one",
            "  base in the codon was IUPAC-ambiguous on both strands.  The reference",
            "  base was used as a placeholder, so the exact amino acid is unknown.",
            "  DO NOT report these mutations without first checking the chromatogram.\n",
        ]
        ul += ambiguous_flags
        unc_block = "\n".join(ul)

    # ═════════════════════════════════════════════════════════════════════════
    # SECTION 9 — ALIGNMENTS  (always shown, at the end)
    # ═════════════════════════════════════════════════════════════════════════
    aln_sec_num = mut_sec_num + (2 if ambiguous_flags else 1)
    aa_aln_block = (
        f"\n  [ {aln_sec_num} ]  AMINO ACID ALIGNMENT\n  {SH}\n"
        f"  Legend:  |  = identical   *  = mismatch   (space) = gap\n"
        f"  Note: lowercase sample residues = translated from ref-filled DNA positions.\n\n"
        + format_alignment_visual(aligned_ref_aa, aligned_samp_aa, 60)
    )
    dna_aln_block = (
        f"\n  [ {aln_sec_num + 1} ]  DNA ALIGNMENT\n  {SH}\n"
        f"  Legend:  |  = identical   *  = mismatch   (space) = gap\n"
        f"  Note: lowercase sample bases = ref-filled (uncovered or IUPAC-ambiguous).\n\n"
        + format_alignment_visual(aligned_ref_dna, aligned_samp_dna_filled, 60)
    )

    # ── Assemble all sections in logical order ────────────────────────────────
    all_parts = [overview_block, asm_block]
    if gap_block:    all_parts.append(gap_block)
    if repair_block: all_parts.append(repair_block)
    if stop_block:   all_parts.append(stop_block)
    if indel_block:  all_parts.append(indel_block)
    all_parts.append(mut_block)
    if unc_block:    all_parts.append(unc_block)
    all_parts.append(aa_aln_block)
    all_parts.append(dna_aln_block)
    all_parts.append(f"\n{HR}\n")

    report = "\n".join(all_parts)
    return report, csv_metrics, corrected_cds, samp_aa


def process_batch(ref_record, contig_paths, results_folder, fasta_folder, strain_label, master_summary):
    os.makedirs(results_folder, exist_ok=True)
    os.makedirs(fasta_folder, exist_ok=True)

    for sample_name, contig_path in contig_paths.items():
        report_text, csv_metrics, extr_dna, extr_aa = analyze_single_sample(
            ref_record, contig_path, strain_label
        )

        report_path = os.path.join(results_folder, f"{sample_name}_{strain_label}_report.txt")
        with open(report_path, "w") as f:
            f.write(report_text)
        log.info(f"  OK Report -> {sample_name}_{strain_label}_report.txt")

        if extr_dna is not None and extr_aa is not None:
            fasta_path = os.path.join(fasta_folder, f"{sample_name}_{strain_label}_extracted.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{sample_name}_F_protein_DNA\n{extr_dna}\n")
                f.write(f">{sample_name}_F_protein_AA\n{extr_aa}\n")

        master_summary.append({"Sample_Name": sample_name, "Strain": strain_label, **csv_metrics})


def write_pipeline_summary(summary_path, master_summary, ref_a, ref_b):
    """
    Writes a short, schematic plain-text report to output/pipeline_summary.txt
    describing what the pipeline did and which tools / algorithms it used.
    """
    import datetime

    n_samples  = len(master_summary)
    n_rsv_a    = sum(1 for r in master_summary if r.get("Strain") == "RSV_A")
    n_rsv_b    = sum(1 for r in master_summary if r.get("Strain") == "RSV_B")
    total_muts = sum(int(r.get("Total_Mutations", 0)) for r in master_summary)
    total_fix  = sum(int(r.get("Frameshifts_Fixed", 0)) for r in master_summary)
    total_stop = sum(int(r.get("Residual_Stops", 0)) for r in master_summary)

    HR  = "=" * 68
    SH  = "-" * 68
    now = datetime.datetime.now().strftime("%Y-%m-%d  %H:%M")

    lines = [
        HR,
        "  SANGER SEQUENCING PIPELINE — RUN SUMMARY",
        f"  Generated : {now}",
        HR,
        "",
        "  PIPELINE OVERVIEW",
        "  " + SH,
        "  Input  : Sanger .ab1 chromatogram files  (RSV-A and/or RSV-B)",
        "  Output : Amplicon FASTAs · Assembled contigs · Mutation reports · F-protein",
        "",
        "  ┌─────────────────────────────────────────────────────────────────┐",
        "  │  STEP 1 — Chromatogram reading  (replaces SnapGene)            │",
        "  │  Tool   : BioPython  SeqIO.read(ab1, 'abi')                    │",
        "  │  Action : Reads raw .ab1 files; reverse-complements RV reads.  │",
        "  │           IUPAC ambiguous bases resolved via peak heights       │",
        "  │           (DATA9–DATA12 channels, PLOC2 index).                │",
        "  │  Trim   : Q15 sliding window (10 bp); hard floor 500 bp.       │",
        "  │  Output : output/01_amplicon_reads/  — one .fa per amplicon    │",
        "  └─────────────────────────────────────────────────────────────────┘",
        "           │",
        "  ┌─────────────────────────────────────────────────────────────────┐",
        "  │  STEP 2 — Contig assembly  (replaces BioEdit CAP)              │",
        "  │  Tool   : BioPython  PairwiseAligner  (Smith-Waterman local)   │",
        "  │  Action : Anchors each trimmed read onto the reference.        │",
        "  │           Pileup voting:  uppercase > lowercase;               │",
        "  │           cross-strand IUPAC intersection to resolve ambiguity.│",
        "  │           Unresolved positions → N (flagged in report).        │",
        "  │  Output : output/02_assembled_contigs/  — one .fa per sample   │",
        "  └─────────────────────────────────────────────────────────────────┘",
        "           │",
        "  ┌─────────────────────────────────────────────────────────────────┐",
        "  │  STEPS 3–5 — Alignment & mutation calling  (replaces MEGA)     │",
        "  │  Tool   : BioPython  PairwiseAligner  (semi-global, gap=-20)   │",
        "  │  Action : Aligns contig to RSV reference (A or B).             │",
        "  │           Extracts F-protein CDS; translates AA sequence.      │",
        "  │           Conservative frameshift repair (alignment-gap based).│",
        "  │           Reports AA/DNA identity, mutations, indels.          │",
        "  │  Output : output/03_analysis_reports/  — per-sample .txt + CSV │",
        "  │           output/04_extracted_F_protein/ — DNA + AA FASTA      │",
        "  └─────────────────────────────────────────────────────────────────┘",
        "",
        "  REFERENCES USED",
        "  " + SH,
        f"  RSV-A : {ref_a.id}  ({len(ref_a.seq)} bp)",
        f"  RSV-B : {ref_b.id}  ({len(ref_b.seq)} bp)",
        "",
        "  RUN STATISTICS",
        "  " + SH,
        f"  Samples processed : {n_samples}  "
        f"(RSV-A: {n_rsv_a}  |  RSV-B: {n_rsv_b})",
        f"  Total AA mutations : {total_muts}",
        f"  Frameshifts repaired: {total_fix}",
        f"  Residual stop codons: {total_stop}"
        + ("  ← manual review recommended" if total_stop else ""),
        "",
    ]

    if master_summary:
        lines += [
            "  PER-SAMPLE RESULTS",
            "  " + SH,
            f"  {'Sample':<30} {'Strain':<7} {'DNA%':>6} {'AA%':>6} {'Muts':>5} {'Stops':>6}",
            "  " + "-" * 60,
        ]
        for r in master_summary:
            lines.append(
                f"  {r['Sample_Name']:<30} {r['Strain']:<7} "
                f"{float(r.get('DNA_Identity(%)', 0)):>6.2f} "
                f"{float(r.get('AA_Identity(%)', 0)):>6.2f} "
                f"{int(r.get('Total_Mutations', 0)):>5} "
                f"{int(r.get('Residual_Stops', 0)):>6}"
            )

    lines += ["", HR, "  END OF SUMMARY", HR, ""]

    with open(summary_path, "w") as f:
        f.write("\n".join(lines))
    print(f"  Saved -> {summary_path}")


def write_master_csv(master_summary, path):
    fieldnames = [
        "Sample_Name", "Strain", "DNA_Identity(%)",
        "AA_Identity(%)", "Total_Mutations", "Frameshifts_Fixed",
        "Residual_Stops", "Ambiguous_Positions", "Mutation_List"
    ]
    with open(path, mode="w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(master_summary)
    print(f"  Saved -> {path}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    # =========================================================
    # CONFIGURE PATHS HERE
    # =========================================================
    RAW_AB1    = os.path.join("input", "Raw_Data")
    REF_A_FILE = os.path.join("input", "references", "ref_RSV_A.fasta")
    REF_B_FILE = os.path.join("input", "references", "ref_RSV_B.fasta")

    FASTAS_DIR  = os.path.join("output", "01_amplicon_reads")
    CONTIGS_DIR = os.path.join("output", "02_assembled_contigs")
    RESULTS_DIR = os.path.join("output", "03_analysis_reports")
    FASTA_DIR   = os.path.join("output", "04_extracted_F_protein")
    MASTER_CSV  = os.path.join("output", "03_analysis_reports", "master_mutation_summary.csv")
    # =========================================================

    print("=" * 60)
    print("  Sanger Sequencing Pipeline — Full Automation")
    print("=" * 60)

    try:
        ref_A = SeqIO.read(REF_A_FILE, "fasta")
        ref_B = SeqIO.read(REF_B_FILE, "fasta")
        log.info(f"Loaded reference A: {ref_A.id} ({len(ref_A.seq)} bp)")
        log.info(f"Loaded reference B: {ref_B.id} ({len(ref_B.seq)} bp)")
    except Exception as e:
        log.error(f"Could not load reference files: {e}")
        exit(1)

    refs = {"RSV_A": ref_A, "RSV_B": ref_B}

    os.makedirs(os.path.join("input", "references"), exist_ok=True)
    os.makedirs("output", exist_ok=True)
    for folder in [RESULTS_DIR, FASTA_DIR, FASTAS_DIR, CONTIGS_DIR]:
        if os.path.exists(folder):
            shutil.rmtree(folder)
            log.info(f"Cleared old outputs: {folder}/")
        os.makedirs(folder, exist_ok=True)

    master_summary = []

    print("\n[STEP 1] Reading chromatograms (auto-detecting RSV-A / RSV-B)...")
    all_amplicons = step1_process_ab1_folder(RAW_AB1, FASTAS_DIR)

    for strain, sample_amplicons in all_amplicons.items():
        if not sample_amplicons:
            log.info(f"No {strain} samples found — skipping.")
            continue

        ref       = refs[strain]
        contigs_d = os.path.join(CONTIGS_DIR, strain)
        results_d = os.path.join(RESULTS_DIR, strain)
        fasta_d   = os.path.join(FASTA_DIR,   strain)

        print(f"\n[STEP 2] Assembling {strain} contigs ({len(sample_amplicons)} sample(s))...")
        contig_paths = step2_assemble_contigs(sample_amplicons, contigs_d, ref)

        print(f"\n[STEPS 3–5] Aligning and calling mutations — {strain}...")
        process_batch(ref, contig_paths, results_d, fasta_d, strain, master_summary)

    if master_summary:
        print("\n[OUTPUT] Writing master CSV...")
        write_master_csv(master_summary, MASTER_CSV)

    print("\n[OUTPUT] Writing pipeline summary report...")
    write_pipeline_summary(
        summary_path=os.path.join("output", "pipeline_summary.txt"),
        master_summary=master_summary,
        ref_a=ref_A, ref_b=ref_B,
    )

    print("\nPipeline complete.")