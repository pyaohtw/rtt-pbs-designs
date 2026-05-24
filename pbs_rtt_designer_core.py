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
    summary_line: str
    warnings: list[str]



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



def revcomp_dna(seq: str) -> str:
    return clean_dna(seq).translate(str.maketrans("ACGT", "TGCA"))[::-1]



def dna_to_rna(seq_dna: str) -> str:
    return str(seq_dna or "").upper().replace("T", "U")



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
            f"Spacer matched multiple sites ({len(matches)} total). Please disambiguate. Matches: {summary}"
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
        raise ValueError("RTT count and range produced duplicate lengths after rounding. Please reduce count or widen the range.")
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
) -> DesignResult:
    target_seq, match = get_target_sequence_and_match(
        genomic_seq=genomic_seq,
        spacer=spacer,
        nick_offset=nick_offset,
    )

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

    if rtt_start_mode == "nick":
        effective_rtt_start_target = match.nick_target
    else:
        effective_rtt_start_target = match.nick_target if rtt_start_target is None else int(rtt_start_target)
        if effective_rtt_start_target < match.nick_target:
            raise ValueError("Selected RTT start must be at or downstream of the nick site in the matched orientation.")

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

    for rtt_len in rtt_lengths:
        candidate = build_rtt_candidate(target_seq, effective_rtt_start_target, rtt_len)
        if candidate is None:
            skipped_rtt.append(rtt_len)
        else:
            rtt_rows.append(candidate)

    if not pbs_rows:
        raise ValueError("None of the requested PBS lengths could be generated at this nick site.")
    if not rtt_rows:
        raise ValueError("None of the requested RTT lengths could be generated at this RTT start site.")

    if skipped_pbs:
        warnings.append("Skipped PBS lengths outside sequence bounds at this nick site: " + ", ".join(map(str, skipped_pbs)))
    if skipped_rtt:
        warnings.append("Skipped RTT lengths outside sequence bounds at this RTT start site: " + ", ".join(map(str, skipped_rtt)))

    pbs_df = pd.DataFrame(pbs_rows).sort_values(["Length", "Sequence"]).reset_index(drop=True)
    pbs_df["Tm"] = pbs_df["Tm"].map(lambda x: round(float(x), 2))
    rtt_df = pd.DataFrame(rtt_rows).sort_values(["Length", "Sequence"]).reset_index(drop=True)

    summary_line = f"Spacer matching strand: {match.strand}; Default PBS: {default_pbs_len} nt, Tm = {round(float(default_pbs_tm), 2)} °C"

    return DesignResult(
        match=match,
        target_orientation=("sense/+" if match.strand == "+" else "reverse-complement/-"),
        default_pbs_len=default_pbs_len,
        default_pbs_tm=(round(float(default_pbs_tm), 2) if default_pbs_tm is not None else None),
        rtt_df=rtt_df,
        pbs_df=pbs_df,
        rtt_text=format_plain_sequences(rtt_df["Sequence"].tolist()),
        pbs_text=format_plain_sequences(pbs_df["Sequence"].tolist()),
        summary_line=summary_line,
        warnings=warnings,
    )
