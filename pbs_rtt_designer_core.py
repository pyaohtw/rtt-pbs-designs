from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Optional

import pandas as pd

import pbs_tm

PBS_MIN_LEN = 5
PBS_MAX_LEN = 25
TARGET_PBS_TM = 37.0

RTT_SELECTOR_SPACER_CONTEXT = 10
RTT_SELECTOR_POST_SPACER_ROW1 = 10
RTT_SELECTOR_POST_SPACER_ROW2 = 20

# Insertion orientation modes (addition mode only).
INSERTION_ORIENTATION_AUTO = "auto"
INSERTION_ORIENTATION_FORWARD = "forward"
INSERTION_ORIENTATION_REVERSE = "reverse"


@dataclass
class SpacerMatch:
    strand: str
    spacer_start_plus: int
    spacer_start_target: int
    spacer_end_plus: int
    spacer_end_target: int
    nick_plus: int
    nick_target: int


@dataclass
class RTTSelectorBase:
    target_index: int
    base: str
    is_nick: bool
    is_in_spacer: bool
    clickable: bool


@dataclass
class RTTStartSelector:
    match: SpacerMatch
    target_orientation: str
    default_rtt_start_target: int
    row1: list[RTTSelectorBase]
    row2: list[RTTSelectorBase]


@dataclass
class DesignResult:
    match: SpacerMatch
    target_orientation: str
    default_pbs_len: Optional[int]
    default_pbs_tm: Optional[float]
    rtt_df: pd.DataFrame
    pbs_df: pd.DataFrame
    rtt_text: str
    pbs_text: str
    rtt_pbs_text: str
    rtt_insert_pbs_text: str
    summary_line: str
    warnings: list[str]
    # Addition-mode extras (None / defaults in legacy mode).
    addition_mode: bool = False
    addition_insert_rna: Optional[str] = None
    addition_retained_dna: Optional[str] = None
    insertion_orientation_applied: Optional[str] = None
    detected_strand: Optional[str] = None


def clean_dna(seq: str) -> str:
    seq = re.sub(r"\s+", "", str(seq or "")).upper().replace("U", "T")
    if not seq:
        raise ValueError("Sequence is empty.")
    bad = sorted({base for base in seq if base not in "ACGT"})
    if bad:
        raise ValueError(f"Sequence contains invalid character(s): {', '.join(bad)}")
    return seq


def clean_spacer(spacer: str) -> str:
    spacer = re.sub(r"\s+", "", str(spacer or "")).upper().replace("U", "T")
    if not spacer:
        raise ValueError("Spacer is empty.")
    bad = sorted({base for base in spacer if base not in "ACGT"})
    if bad:
        raise ValueError(f"Spacer contains invalid character(s): {', '.join(bad)}")
    return spacer


def validate_insertion_sequence(seq: str, allow_empty: bool = True) -> str:
    seq = re.sub(r"\s+", "", str(seq or ""))
    if not seq:
        if allow_empty:
            return ""
        raise ValueError("Insertion sequence is empty.")
    bad = sorted({base for base in seq if base not in "ATCGUatcgu"})
    if bad:
        raise ValueError(f"Insertion contains invalid character(s): {', '.join(bad)}")
    return seq


def revcomp_dna(seq: str) -> str:
    return clean_dna(seq).translate(str.maketrans("ACGT", "TGCA"))[::-1]


def dna_to_rna(seq_dna: str) -> str:
    return str(seq_dna or "").upper().replace("T", "U")


def insertion_to_combo_output(seq: str) -> str:
    return str(seq or "").replace("T", "U").replace("t", "u")


def score_rnadna_tm_from_pbs_dna(pbs_dna: str) -> float:
    return float(pbs_tm.score_one_pbs(dna_to_rna(pbs_dna))["tm"])


def find_spacer_matches(genomic_seq: str, spacer: str, nick_offset: int = 3) -> list[SpacerMatch]:
    seq_plus = clean_dna(genomic_seq)
    seq_minus = revcomp_dna(seq_plus)
    spacer_dna = clean_spacer(spacer)
    n = len(seq_plus)
    spacer_len = len(spacer_dna)
    nick_offset = int(nick_offset)
    if nick_offset < 0:
        raise ValueError("Nick offset must be 0 or greater.")
    pattern = re.compile(rf"(?=({re.escape(spacer_dna)}))", flags=re.IGNORECASE)
    matches: list[SpacerMatch] = []
    for match in pattern.finditer(seq_plus):
        spacer_start_plus = match.start()
        spacer_end_plus = spacer_start_plus + spacer_len
        nick_plus = spacer_end_plus - nick_offset
        if nick_plus < 0 or nick_plus > n:
            continue
        matches.append(
            SpacerMatch(
                strand="+",
                spacer_start_plus=spacer_start_plus,
                spacer_start_target=spacer_start_plus,
                spacer_end_plus=spacer_end_plus,
                spacer_end_target=spacer_end_plus,
                nick_plus=nick_plus,
                nick_target=nick_plus,
            )
        )
    for match in pattern.finditer(seq_minus):
        spacer_start_minus = match.start()
        spacer_end_minus = spacer_start_minus + spacer_len
        nick_minus = spacer_end_minus - nick_offset
        if nick_minus < 0 or nick_minus > len(seq_minus):
            continue
        spacer_start_plus = n - spacer_end_minus
        spacer_end_plus = n - spacer_start_minus
        nick_plus = n - nick_minus
        matches.append(
            SpacerMatch(
                strand="-",
                spacer_start_plus=spacer_start_plus,
                spacer_start_target=spacer_start_minus,
                spacer_end_plus=spacer_end_plus,
                spacer_end_target=spacer_end_minus,
                nick_plus=nick_plus,
                nick_target=nick_minus,
            )
        )
    return matches


def resolve_unique_match(genomic_seq: str, spacer: str, nick_offset: int = 3) -> SpacerMatch:
    matches = find_spacer_matches(genomic_seq, spacer, nick_offset=nick_offset)
    if len(matches) == 0:
        raise ValueError("Spacer was not found in the genomic sequence on either strand.")
    if len(matches) > 1:
        summary = "; ".join(
            f"{m.strand} strand at position {m.spacer_start_plus + 1}" for m in matches[:10]
        )
        raise ValueError(
            f"Spacer matched multiple sites ({len(matches)} total). "
            f"Please disambiguate. Matches: {summary}"
        )
    return matches[0]


def get_target_sequence_and_match(genomic_seq: str, spacer: str, nick_offset: int = 3) -> tuple[str, SpacerMatch]:
    seq_plus = clean_dna(genomic_seq)
    seq_minus = revcomp_dna(seq_plus)
    match = resolve_unique_match(seq_plus, spacer, nick_offset=nick_offset)
    target_seq = seq_plus if match.strand == "+" else seq_minus
    return target_seq, match


def _make_selector_bases(target_seq: str, match: SpacerMatch, start: int, end: int) -> list[RTTSelectorBase]:
    out: list[RTTSelectorBase] = []
    start = max(0, start)
    end = min(len(target_seq), end)
    for idx in range(start, end):
        out.append(
            RTTSelectorBase(
                target_index=idx,
                base=target_seq[idx],
                is_nick=(idx == match.nick_target),
                is_in_spacer=(match.spacer_start_target <= idx < match.spacer_end_target),
                clickable=(idx >= match.nick_target),
            )
        )
    return out


def build_rtt_start_selector(genomic_seq: str, spacer: str, nick_offset: int = 3) -> RTTStartSelector:
    target_seq, match = get_target_sequence_and_match(
        genomic_seq=genomic_seq,
        spacer=spacer,
        nick_offset=nick_offset,
    )
    row1_start = match.spacer_end_target - RTT_SELECTOR_SPACER_CONTEXT
    row1_end = match.spacer_end_target + RTT_SELECTOR_POST_SPACER_ROW1
    row2_start = row1_end
    row2_end = row2_start + RTT_SELECTOR_POST_SPACER_ROW2
    row1 = _make_selector_bases(target_seq, match, row1_start, row1_end)
    row2 = _make_selector_bases(target_seq, match, row2_start, row2_end)
    return RTTStartSelector(
        match=match,
        target_orientation=("sense/+" if match.strand == "+" else "reverse-complement/-"),
        default_rtt_start_target=match.nick_target,
        row1=row1,
        row2=row2,
    )


def parse_length_list(text: str, minimum: int, maximum: int, label: str) -> list[int]:
    values: list[int] = []
    for part in re.split(r"[\s,;]+", str(text or "").strip()):
        if not part:
            continue
        try:
            value = int(part)
        except ValueError as exc:
            raise ValueError(f"{label} contains a non-integer value: {part}") from exc
        if value < minimum or value > maximum:
            raise ValueError(f"{label} length must be between {minimum} and {maximum}: {value}")
        values.append(value)
    values = sorted(set(values))
    if not values:
        raise ValueError(f"Please provide at least one {label} length.")
    return values


def evenly_spaced_lengths(min_len: int, max_len: int, count: int) -> list[int]:
    min_len = int(min_len)
    max_len = int(max_len)
    count = int(count)
    if min_len > max_len:
        raise ValueError("RTT min length cannot be greater than RTT max length.")
    if count < 2:
        raise ValueError("RTT count must be at least 2 in range mode.")
    if count > (max_len - min_len + 1):
        raise ValueError("RTT count is too large for the selected min/max range to produce unique integer lengths.")
    if count == 2:
        return [min_len, max_len]
    step = (max_len - min_len) / (count - 1)
    values = [math.ceil(min_len + i * step) for i in range(count)]
    values[0] = min_len
    values[-1] = max_len
    values = sorted(set(values))
    if len(values) != count:
        raise ValueError(
            "RTT count and range produced duplicate lengths after rounding. "
            "Please reduce count or widen the range."
        )
    return values


def build_pbs_candidate(target_seq_dna: str, nick_target: int, pbs_len: int) -> Optional[dict]:
    if nick_target - pbs_len < 0:
        return None
    template = target_seq_dna[nick_target - pbs_len:nick_target]
    if len(template) != pbs_len:
        return None
    pbs_dna = revcomp_dna(template)
    pbs_rna = dna_to_rna(pbs_dna)
    tm = score_rnadna_tm_from_pbs_dna(pbs_dna)
    return {
        "Sequence": pbs_rna,
        "Length": pbs_len,
        "Tm": tm,
    }


def build_rtt_candidate(target_seq_dna: str, rtt_start_target: int, rtt_len: int) -> Optional[dict]:
    if rtt_start_target < 0:
        return None
    if rtt_start_target + rtt_len > len(target_seq_dna):
        return None
    template = target_seq_dna[rtt_start_target:rtt_start_target + rtt_len]
    if len(template) != rtt_len:
        return None
    rtt_dna = revcomp_dna(template)
    rtt_rna = dna_to_rna(rtt_dna)
    return {
        "Sequence": rtt_rna,
        "Length": rtt_len,
    }


def orient_insertion_for_flap(typed_ins_dna: str, strand: str, orientation: str) -> str:
    """Return the insertion in the flap (matched/nicked-strand) orientation.

    The user types the insertion in their input-DNA (sense/plus) orientation. The flap is
    synthesized in the matched-strand orientation (= input orientation on the + strand,
    = reverse complement of the input on the - strand).

    orientation:
      - "auto"/"forward": the insertion should read forward (as typed) in the final input DNA.
      - "reverse": the insertion should read reverse-complemented in the final input DNA.

    Net effect (default forward): the edited input strand always carries the insertion exactly
    as typed, regardless of which strand the spacer matched.
    """
    ins = clean_dna(typed_ins_dna)
    if strand == "-":
        ins = revcomp_dna(ins)  # convert input-orientation -> minus-strand (flap) orientation
    if orientation == INSERTION_ORIENTATION_REVERSE:
        ins = revcomp_dna(ins)
    return ins


def build_addition_homology_rtt(target_seq_dna: str, insert_after_target: int, homology_len: int) -> Optional[dict]:
    """Addition-mode RTT column = reverse-complement RNA of the 3' homology arm only.

    The homology arm is the genomic stretch immediately downstream of the insertion site.
    The retained genomic (nick -> insertion site) and the insertion itself live in the
    'insert' column, so that RTT + insert + PBS concatenates into the pegRNA 3' extension.
    """
    if homology_len < 1:
        return None
    downstream_start = insert_after_target + 1
    if downstream_start < 0:
        return None
    homology = target_seq_dna[downstream_start:downstream_start + homology_len]
    if len(homology) != homology_len:
        return None
    rtt_dna = revcomp_dna(homology)
    return {
        "Sequence": dna_to_rna(rtt_dna),
        "Length": homology_len,
    }


def pick_tm_optimal_pbs(target_seq_dna: str, nick_target: int) -> tuple[Optional[int], Optional[float]]:
    best_len = None
    best_tm = None
    best_score = None
    for length in range(PBS_MIN_LEN, PBS_MAX_LEN + 1):
        candidate = build_pbs_candidate(target_seq_dna, nick_target, length)
        if candidate is None:
            continue
        tm = float(candidate["Tm"])
        score = (abs(tm - TARGET_PBS_TM), length)
        if best_score is None or score < best_score:
            best_score = score
            best_len = length
            best_tm = tm
    return best_len, best_tm


def pbs_lengths_from_mode(target_seq_dna: str, nick_target: int, shorter: int = 0, longer: int = 0) -> tuple[list[int], Optional[int], Optional[float]]:
    default_len, default_tm = pick_tm_optimal_pbs(target_seq_dna, nick_target)
    if default_len is None:
        raise ValueError("No valid PBS could be generated at this nick site.")
    shorter = max(0, int(shorter))
    longer = max(0, int(longer))
    lengths = list(range(max(PBS_MIN_LEN, default_len - shorter), min(PBS_MAX_LEN, default_len + longer) + 1))
    return lengths, default_len, default_tm


def rtt_lengths_from_mode(rtt_mode: str, manual_lengths_text: str = "", min_len: int = 10, max_len: int = 20, count: int = 3) -> list[int]:
    if rtt_mode == "manual":
        return parse_length_list(manual_lengths_text, 1, 500, "RTT")
    return evenly_spaced_lengths(min_len=min_len, max_len=max_len, count=count)


def format_plain_sequences(sequences: list[str]) -> str:
    return "\n".join(sequences)


def format_delimited_rows(rows: list[tuple[str, ...]]) -> str:
    return "\n".join("\t".join(row) for row in rows)


def design_pbs_rtt(
    genomic_seq: str,
    spacer: str,
    nick_offset: int = 3,
    pbs_shorter: int = 0,
    pbs_longer: int = 0,
    rtt_mode: str = "range",
    rtt_manual_lengths: str = "",
    rtt_min: int = 10,
    rtt_max: int = 20,
    rtt_count: int = 3,
    rtt_start_mode: str = "selected",
    rtt_start_target: Optional[int] = None,
    include_insertion: bool = True,
    insertion_sequence: str = "",
    addition_mode: bool = False,
    insertion_orientation: str = INSERTION_ORIENTATION_AUTO,
    insert_at_nick: bool = False,
) -> DesignResult:
    target_seq, match = get_target_sequence_and_match(
        genomic_seq=genomic_seq,
        spacer=spacer,
        nick_offset=nick_offset,
    )

    insertion_sequence = validate_insertion_sequence(insertion_sequence, allow_empty=not addition_mode)
    if addition_mode and not insertion_sequence:
        raise ValueError("Addition (pure insertion) mode requires a non-empty insertion sequence.")

    pbs_lengths, default_pbs_len, default_pbs_tm = pbs_lengths_from_mode(
        target_seq_dna=target_seq,
        nick_target=match.nick_target,
        shorter=pbs_shorter,
        longer=pbs_longer,
    )
    rtt_lengths = rtt_lengths_from_mode(
        rtt_mode=rtt_mode,
        manual_lengths_text=rtt_manual_lengths,
        min_len=rtt_min,
        max_len=rtt_max,
        count=rtt_count,
    )

    # Resolve the site.
    #  - legacy mode: first base included in RTT.
    #  - addition mode: the base AFTER which the insertion is placed.
    if rtt_start_mode == "nick":
        effective_site_target = match.nick_target
    else:
        effective_site_target = match.nick_target if rtt_start_target is None else int(rtt_start_target)

    if addition_mode and insert_at_nick:
        # Insert right at the nick junction: retain zero genomic bases.
        effective_site_target = match.nick_target - 1

    if not (addition_mode and insert_at_nick):
        if effective_site_target < match.nick_target:
            raise ValueError("Selected site must be at or downstream of the nick site in the matched orientation.")
    if effective_site_target >= len(target_seq):
        raise ValueError("Selected site is outside the target sequence.")

    pbs_rows: list[dict] = []
    rtt_rows: list[dict] = []
    warnings: list[str] = []
    skipped_pbs: list[int] = []
    skipped_rtt: list[int] = []

    for pbs_len in pbs_lengths:
        candidate = build_pbs_candidate(target_seq, match.nick_target, pbs_len)
        if candidate is None:
            skipped_pbs.append(pbs_len)
        else:
            pbs_rows.append(candidate)

    addition_insert_rna: Optional[str] = None
    addition_retained_dna: Optional[str] = None

    if addition_mode:
        typed_ins_dna = clean_dna(insertion_sequence)
        ins_flap = orient_insertion_for_flap(typed_ins_dna, match.strand, insertion_orientation)
        # Retained genomic (matched-strand orientation) from the nick up to and including the site.
        retained = target_seq[match.nick_target:effective_site_target + 1]
        addition_retained_dna = retained
        # 'insert' column = reverse-complement RNA of (retained + insertion), so that
        # RTT (revcomp homology) + insert + PBS = the pegRNA 3' extension (5'->3').
        insert_region_dna = retained + ins_flap
        addition_insert_rna = dna_to_rna(revcomp_dna(insert_region_dna))
        for homology_len in rtt_lengths:
            candidate = build_addition_homology_rtt(target_seq, effective_site_target, homology_len)
            if candidate is None:
                skipped_rtt.append(homology_len)
            else:
                rtt_rows.append(candidate)
    else:
        for rtt_len in rtt_lengths:
            candidate = build_rtt_candidate(target_seq, effective_site_target, rtt_len)
            if candidate is None:
                skipped_rtt.append(rtt_len)
            else:
                rtt_rows.append(candidate)

    if not pbs_rows:
        raise ValueError("None of the requested PBS lengths could be generated at this nick site.")
    if not rtt_rows:
        if addition_mode:
            raise ValueError("None of the requested homology-arm lengths fit downstream of the insertion site.")
        raise ValueError("None of the requested RTT lengths could be generated at this RTT start site.")

    if skipped_pbs:
        warnings.append("Skipped PBS lengths outside sequence bounds at this nick site: " + ", ".join(map(str, skipped_pbs)))
    if skipped_rtt:
        if addition_mode:
            warnings.append("Skipped homology-arm lengths that run past the end of the sequence downstream of the insertion site: " + ", ".join(map(str, skipped_rtt)))
        else:
            warnings.append("Skipped RTT lengths outside sequence bounds at this RTT start site: " + ", ".join(map(str, skipped_rtt)))

    pbs_df = pd.DataFrame(pbs_rows).sort_values(["Length", "Sequence"]).reset_index(drop=True)
    pbs_df["Tm"] = pbs_df["Tm"].map(lambda x: round(float(x), 2))

    rtt_df = pd.DataFrame(rtt_rows).sort_values(["Length", "Sequence"]).reset_index(drop=True)

    rtt_sequences = rtt_df["Sequence"].tolist()
    pbs_sequences = pbs_df["Sequence"].tolist()

    rtt_pbs_rows = [(rtt_seq, pbs_seq) for pbs_seq in pbs_sequences for rtt_seq in rtt_sequences]
    rtt_pbs_text = format_delimited_rows(rtt_pbs_rows)

    rtt_insert_pbs_text = ""
    if addition_mode:
        # 3-column pegRNA-extension pieces: RTT (revcomp homology) | insert (revcomp retained+insertion) | PBS.
        rtt_insert_pbs_rows = [
            (rtt_seq, addition_insert_rna, pbs_seq)
            for pbs_seq in pbs_sequences
            for rtt_seq in rtt_sequences
        ]
        rtt_insert_pbs_text = format_delimited_rows(rtt_insert_pbs_rows)
    elif include_insertion and insertion_sequence:
        # Legacy behavior: insertion block placed between RTT and PBS (forward, not templated into RTT).
        insert_for_output = insertion_to_combo_output(insertion_sequence)
        rtt_insert_pbs_rows = [
            (rtt_seq, insert_for_output, pbs_seq)
            for pbs_seq in pbs_sequences
            for rtt_seq in rtt_sequences
        ]
        rtt_insert_pbs_text = format_delimited_rows(rtt_insert_pbs_rows)

    summary_line = f"Spacer matching strand: {match.strand}; Default PBS: {default_pbs_len} nt, Tm = {round(float(default_pbs_tm), 2)} \u00b0C"
    if addition_mode:
        if insert_at_nick or effective_site_target < match.nick_target:
            site_desc = "at the nick junction (0 retained bases)"
        else:
            offset = effective_site_target - match.nick_target
            site_desc = f"after site +{offset} (base '{target_seq[effective_site_target]}')"
        orient_desc = {
            INSERTION_ORIENTATION_AUTO: "auto / strand-aware (forward in input DNA)",
            INSERTION_ORIENTATION_FORWARD: "forced forward in input DNA",
            INSERTION_ORIENTATION_REVERSE: "forced reverse-complement in input DNA",
        }.get(insertion_orientation, insertion_orientation)
        summary_line += (
            f"; Mode: Addition (pure insertion) on the {match.strand} strand \u2014 "
            f"{len(clean_dna(insertion_sequence))} nt inserted {site_desc}; "
            f"orientation: {orient_desc}. RTT column = reverse-complement RNA of the 3\u2032 homology arm; "
            f"insert (revcomp of retained+insertion) = {addition_insert_rna}. "
            f"RTT + insert + PBS = pegRNA 3\u2032 extension."
        )

    return DesignResult(
        match=match,
        target_orientation=("sense/+" if match.strand == "+" else "reverse-complement/-"),
        default_pbs_len=default_pbs_len,
        default_pbs_tm=(round(float(default_pbs_tm), 2) if default_pbs_tm is not None else None),
        rtt_df=rtt_df,
        pbs_df=pbs_df,
        rtt_text=format_plain_sequences(rtt_sequences),
        pbs_text=format_plain_sequences(pbs_sequences),
        rtt_pbs_text=rtt_pbs_text,
        rtt_insert_pbs_text=rtt_insert_pbs_text,
        summary_line=summary_line,
        warnings=warnings,
        addition_mode=addition_mode,
        addition_insert_rna=addition_insert_rna,
        addition_retained_dna=addition_retained_dna,
        insertion_orientation_applied=(insertion_orientation if addition_mode else None),
        detected_strand=match.strand,
    )


# ---------------------------------------------------------------------------
# Batch pure-insertion mode
# ---------------------------------------------------------------------------

BATCH_PLACEMENT_NICK = "nick"          # insert exactly at the nick (retain 0 bases)
BATCH_PLACEMENT_POSITION = "position"  # insert right after a 1-based input-DNA position

QWC_HALF_WINDOW = 22

BATCH_COLUMNS = [
    "spacer ID",
    "spacer seq",
    "insertion ID",
    "insertion seq",
    "RTT",
    "insert",
    "PBS",
    "input DNA",
    "edited DNA",
    "Quantification_Window_Coordinates",
]


def parse_two_column_table(text: str, label: str, validate_dna: bool = True) -> list[tuple[str, str]]:
    """Parse a 2-column (ID, sequence) block.

    Accepts tab / comma / multi-space separated columns, one entry per line.
    Blank lines and lines starting with '#' are ignored. Returns [(id, seq), ...]
    preserving input order. Duplicate IDs raise an error.
    """
    rows: list[tuple[str, str]] = []
    seen_ids: set[str] = set()
    for lineno, raw in enumerate(str(text or "").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"[\t,]+|\s{2,}|\s+", line)
        parts = [p for p in parts if p != ""]
        if len(parts) < 2:
            raise ValueError(f"{label} line {lineno} needs two columns (ID and sequence): {raw!r}")
        entry_id = parts[0].strip()
        seq = "".join(parts[1:]).strip()
        if not entry_id:
            raise ValueError(f"{label} line {lineno} is missing an ID.")
        if entry_id in seen_ids:
            raise ValueError(f"{label} has a duplicate ID: {entry_id}")
        if validate_dna:
            seq = clean_spacer(seq)  # uppercases, U->T, validates ACGT
        else:
            seq = clean_dna(seq)
        seen_ids.add(entry_id)
        rows.append((entry_id, seq))
    if not rows:
        raise ValueError(f"No valid {label.lower()} entries were provided.")
    return rows


def quantification_window_coordinates(nick_plus: int, seq_len: int, half: int = QWC_HALF_WINDOW) -> str:
    """CRISPResso -qwc string, 0-based, clipped to [0, seq_len-1].

    Window center is the nick site expressed in input-DNA (plus-strand) coordinates,
    which equals the 0-based index of the base immediately 3' of the cut on the plus
    strand, regardless of which strand the spacer matched.
    """
    start = max(0, nick_plus - half)
    stop = min(seq_len - 1, nick_plus + half)
    return f"{start}-{stop}"


def _insertion_in_input_orientation(typed_ins_dna: str, orientation: str) -> str:
    """Insertion as it appears in the edited input (plus) DNA.

    auto/forward -> as typed; reverse -> reverse complement of typed.
    """
    ins = clean_dna(typed_ins_dna)
    if orientation == INSERTION_ORIENTATION_REVERSE:
        ins = revcomp_dna(ins)
    return ins


def _resolve_batch_site(match: SpacerMatch, seq_len: int, placement_mode: str, position_1based: Optional[int]):
    """Map a placement request to (rtt_start_target, insert_at_nick, insert_pos_plus).

    insert_pos_plus is the 0-based plus-strand slice index where the insertion is
    spliced into the input DNA (edited = input[:insert_pos_plus] + ins + input[insert_pos_plus:]).

    Returns None if the requested position is not reachable from this spacer's nick
    on the RTT-templated (downstream) side.
    """
    if placement_mode == BATCH_PLACEMENT_NICK:
        return {"rtt_start_target": None, "insert_at_nick": True, "insert_pos_plus": match.nick_plus}

    if position_1based is None:
        raise ValueError("Position mode requires a 1-based input-DNA position.")
    P = int(position_1based)
    if P < 1 or P > seq_len:
        raise ValueError(f"Position must be between 1 and {seq_len} (got {P}).")

    # Insertion junction sits at plus coordinate P (insert right AFTER plus index P-1).
    insert_pos_plus = P

    if match.strand == "+":
        # Reachable when junction is at or downstream (>=) of the nick.
        if P < match.nick_plus:
            return None
        if P == match.nick_plus:
            return {"rtt_start_target": None, "insert_at_nick": True, "insert_pos_plus": insert_pos_plus}
        # retain (P - nick_plus) genomic bases; last retained plus index = P-1
        return {"rtt_start_target": P - 1, "insert_at_nick": False, "insert_pos_plus": insert_pos_plus}
    else:
        # Minus strand: reachable when junction is at or downstream on the minus strand,
        # i.e. P <= nick_plus (lower plus coordinate).
        if P > match.nick_plus:
            return None
        if P == match.nick_plus:
            return {"rtt_start_target": None, "insert_at_nick": True, "insert_pos_plus": insert_pos_plus}
        # last retained minus index = n - P - 1
        return {"rtt_start_target": seq_len - P - 1, "insert_at_nick": False, "insert_pos_plus": insert_pos_plus}


def design_batch_insertion(
    genomic_seq: str,
    spacers: list[tuple[str, str]],
    insertions: list[tuple[str, str]],
    placement_mode: str = BATCH_PLACEMENT_NICK,
    position_1based: Optional[int] = None,
    nick_offset: int = 3,
    pbs_shorter: int = 0,
    pbs_longer: int = 0,
    rtt_mode: str = "range",
    rtt_manual_lengths: str = "",
    rtt_min: int = 10,
    rtt_max: int = 20,
    rtt_count: int = 3,
    insertion_orientation: str = INSERTION_ORIENTATION_AUTO,
    qwc_half: int = QWC_HALF_WINDOW,
) -> tuple[pd.DataFrame, list[str]]:
    """Enumerate pure-insertion pegRNA designs for every spacer x insertion combination.

    Rows are ordered spacer-major then insertion-major; within a combination all
    RTT x PBS length combinations are expanded (each becomes one row). Returns the
    result DataFrame (BATCH_COLUMNS) and a list of human-readable warnings/skips.
    """
    input_dna = clean_dna(genomic_seq)
    seq_len = len(input_dna)
    warnings: list[str] = []
    out_rows: list[dict] = []

    if placement_mode not in (BATCH_PLACEMENT_NICK, BATCH_PLACEMENT_POSITION):
        raise ValueError(f"Unknown placement mode: {placement_mode}")

    for spacer_id, spacer_seq in spacers:
        # Resolve the (unique) match for this spacer.
        try:
            matches = find_spacer_matches(input_dna, spacer_seq, nick_offset=nick_offset)
        except ValueError as exc:
            warnings.append(f"[{spacer_id}] invalid spacer skipped: {exc}")
            continue
        if len(matches) == 0:
            warnings.append(f"[{spacer_id}] not found in the input DNA on either strand - all combinations skipped.")
            continue
        if len(matches) > 1:
            sites = "; ".join(f"{m.strand} strand @ {m.spacer_start_plus + 1}" for m in matches[:6])
            warnings.append(f"[{spacer_id}] matched {len(matches)} sites ({sites}) - ambiguous, all combinations skipped.")
            continue
        match = matches[0]
        target_seq = input_dna if match.strand == "+" else revcomp_dna(input_dna)

        qwc = quantification_window_coordinates(match.nick_plus, seq_len, half=qwc_half)

        # Resolve placement once per spacer (position mode) or trivially (nick mode).
        try:
            placement = _resolve_batch_site(match, seq_len, placement_mode, position_1based)
        except ValueError as exc:
            warnings.append(f"[{spacer_id}] {exc}")
            continue
        if placement is None:
            side = "downstream (higher position)" if match.strand == "+" else "downstream (lower position)"
            warnings.append(
                f"[{spacer_id}] position {position_1based} is not reachable from this spacer's nick "
                f"(nick at input-DNA index {match.nick_plus}, {match.strand} strand extends {side}) - "
                f"all its insertion combinations skipped."
            )
            continue

        for ins_id, ins_seq in insertions:
            try:
                result = design_pbs_rtt(
                    genomic_seq=input_dna,
                    spacer=spacer_seq,
                    nick_offset=nick_offset,
                    pbs_shorter=pbs_shorter,
                    pbs_longer=pbs_longer,
                    rtt_mode=rtt_mode,
                    rtt_manual_lengths=rtt_manual_lengths,
                    rtt_min=rtt_min,
                    rtt_max=rtt_max,
                    rtt_count=rtt_count,
                    rtt_start_mode="selected",
                    rtt_start_target=placement["rtt_start_target"],
                    include_insertion=False,
                    insertion_sequence=ins_seq,
                    addition_mode=True,
                    insertion_orientation=insertion_orientation,
                    insert_at_nick=placement["insert_at_nick"],
                )
            except ValueError as exc:
                warnings.append(f"[{spacer_id} x {ins_id}] skipped: {exc}")
                continue

            # Edited DNA (input/plus orientation), insertion forward-as-typed unless forced reverse.
            ins_in_input = _insertion_in_input_orientation(ins_seq, insertion_orientation)
            pos = placement["insert_pos_plus"]
            edited_dna = input_dna[:pos] + ins_in_input + input_dna[pos:]

            insert_rna = result.addition_insert_rna
            rtt_sequences = result.rtt_df.sort_values(["Length", "Sequence"])["Sequence"].tolist()
            pbs_sequences = result.pbs_df.sort_values(["Length", "Sequence"])["Sequence"].tolist()

            for rtt_seq in rtt_sequences:
                for pbs_seq in pbs_sequences:
                    out_rows.append({
                        "spacer ID": spacer_id,
                        "spacer seq": spacer_seq,
                        "insertion ID": ins_id,
                        "insertion seq": ins_seq,
                        "RTT": rtt_seq,
                        "insert": insert_rna,
                        "PBS": pbs_seq,
                        "input DNA": input_dna,
                        "edited DNA": edited_dna,
                        "Quantification_Window_Coordinates": qwc,
                    })

    df = pd.DataFrame(out_rows, columns=BATCH_COLUMNS)
    return df, warnings
