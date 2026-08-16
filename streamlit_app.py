from __future__ import annotations

import hashlib
from html import escape

import streamlit as st

from pbs_rtt_designer_core import (
    BATCH_COLUMNS,
    BATCH_PLACEMENT_NICK,
    BATCH_PLACEMENT_POSITION,
    INSERTION_ORIENTATION_AUTO,
    INSERTION_ORIENTATION_FORWARD,
    INSERTION_ORIENTATION_REVERSE,
    build_rtt_start_selector,
    design_batch_insertion,
    design_pbs_rtt,
    parse_two_column_table,
    validate_insertion_sequence,
)

st.set_page_config(page_title="PBS/RTT Sub-tool", layout="wide")
st.title("PBS / RTT designer")
st.caption("Minimal Streamlit sub-tool for spacer QC, PBS-by-Tm selection, and RTT length design.")

st.markdown(
    """
    <style>
    .rtt-seq-wrap {
        --nt-size: 1.55rem;
        --nt-gap: 0.18rem;
        margin-top: 0.5rem;
        margin-bottom: 0.5rem;
        display: inline-block;
        background: #ffffff;
        border: 1px solid #d9d9d9;
        border-radius: 0.6rem;
        padding: 0.65rem 0.7rem;
        box-shadow: 0 1px 2px rgba(0, 0, 0, 0.04);
    }
    .rtt-seq-row {
        display: grid;
        grid-template-columns: repeat(40, var(--nt-size));
        gap: var(--nt-gap);
        margin-bottom: 0.3rem;
        align-items: center;
    }
    .rtt-seq-row:last-child {
        margin-bottom: 0;
    }
    .rtt-nt {
        width: var(--nt-size);
        height: var(--nt-size);
        display: flex;
        align-items: center;
        justify-content: center;
        border-radius: 0.25rem;
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
        font-size: 0.92rem;
        line-height: 1;
        text-decoration: none;
        color: #111111;
        background: #ffffff;
        border: 1px solid transparent;
        box-sizing: border-box;
    }
    .rtt-nt.spacer {
        color: #d62828;
        font-weight: 700;
    }
    .rtt-nt.selected {
        background: #ffeb3b;
        color: #111111;
    }
    .rtt-nt.nick {
        box-shadow: inset 0 -3px 0 #4c78ff;
    }
    .rtt-legend {
        font-size: 0.9rem;
        color: #555555;
        margin-top: 0.25rem;
    }
    .rtt-offset-row {
        display: grid;
        grid-template-columns: repeat(40, var(--nt-size));
        gap: var(--nt-gap);
        margin-top: 0.35rem;
        margin-bottom: 0.12rem;
    }
    .rtt-offset-cell {
        text-align: center;
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
        font-size: 0.65rem;
        line-height: 1;
        color: #aeb6c2;
        min-height: 0.95rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


def _selector_signature(genomic_seq: str, spacer: str, nick_offset: int) -> str:
    raw = f"{genomic_seq}|{spacer}|{nick_offset}"
    return hashlib.md5(raw.encode("utf-8")).hexdigest()


def _render_selector_row(row, selected_target: int) -> str:
    cells: list[str] = []
    for base in row:
        classes = ["rtt-nt"]
        if base.is_in_spacer:
            classes.append("spacer")
        if base.is_nick:
            classes.append("nick")
        if base.target_index == selected_target:
            classes.append("selected")
        cells.append(f'<span class="{" ".join(classes)}">{escape(base.base)}</span>')
    return '<div class="rtt-seq-row">' + "".join(cells) + "</div>"


def _offset_label(offset: int) -> str:
    if offset == 0:
        return "0"
    return f"{offset:+d}"


def render_rtt_start_selector(selector, signature: str, show_buttons: bool = True, addition_mode: bool = False) -> int:
    all_bases = selector.row1 + selector.row2
    selectable_bases = [base for base in all_bases if base.clickable]
    if not selectable_bases:
        st.warning("No selectable positions are available in the preview window.")
        return selector.default_rtt_start_target
    valid_targets = {base.target_index for base in selectable_bases}
    if st.session_state.get("rtt_selector_signature") != signature:
        st.session_state["rtt_selector_signature"] = signature
        st.session_state["selected_rtt_start_target"] = selector.default_rtt_start_target
    selected_target = int(st.session_state.get("selected_rtt_start_target", selector.default_rtt_start_target))
    if selected_target not in valid_targets:
        selected_target = selector.default_rtt_start_target
        st.session_state["selected_rtt_start_target"] = selected_target
    if not show_buttons:
        selected_target = selector.default_rtt_start_target
        st.session_state["selected_rtt_start_target"] = selected_target

    if addition_mode:
        st.markdown("### Insertion site selector")
        st.caption("Matched orientation is shown left-to-right. Choose the site AFTER which the insertion is placed.")
    else:
        st.markdown("### RTT start selector")
        st.caption("Matched orientation is shown left-to-right. Choose the first base included in RTT.")

    html = '<div class="rtt-seq-wrap">'
    offset_labels = []
    for base in all_bases:
        offset = base.target_index - selector.match.nick_target
        offset_labels.append(f'<div class="rtt-offset-cell">{escape(_offset_label(offset))}</div>')
    html += '<div class="rtt-offset-row">' + ''.join(offset_labels) + '</div>'
    html += _render_selector_row(all_bases, selected_target)
    html += '</div>'
    st.markdown(html, unsafe_allow_html=True)

    if addition_mode:
        legend = "Red = spacer &middot; Yellow = insertion site (insertion goes right after it) &middot; Blue underline = nick site"
    else:
        legend = "Red = spacer &middot; Yellow = selected RTT start &middot; Blue underline = nick site"
    st.markdown(f"<div class='rtt-legend'>{legend}</div>", unsafe_allow_html=True)

    if show_buttons:
        label_cols = st.columns(len(all_bases), gap="small")
        for col, base in zip(label_cols, all_bases):
            with col:
                st.caption(_offset_label(base.target_index - selector.match.nick_target))
        button_cols = st.columns(len(all_bases), gap="small")
        button_help = (
            "Insert the insertion sequence right after this base."
            if addition_mode
            else "Choose this base as the first nucleotide included in RTT."
        )
        for col, base in zip(button_cols, all_bases):
            with col:
                if base.clickable:
                    clicked = st.button(
                        base.base,
                        key=f"rtt_start_btn_{signature}_{base.target_index}",
                        type=("primary" if base.target_index == selected_target else "secondary"),
                        width="stretch",
                        help=button_help,
                    )
                    if clicked:
                        st.session_state["selected_rtt_start_target"] = base.target_index
                        st.rerun()
                else:
                    st.button(
                        base.base,
                        key=f"rtt_start_btn_disabled_{signature}_{base.target_index}",
                        disabled=True,
                        width="stretch",
                    )
    else:
        if addition_mode:
            st.caption("Insertion site is locked to the nick site (offset 0) in this mode.")
        else:
            st.caption("RTT start is locked to the nick site in this mode.")
    return int(st.session_state.get("selected_rtt_start_target", selected_target))



def render_batch_mode() -> None:
    st.markdown("## Batch pure-insertion designer")
    st.caption(
        "Enter a list of spacers and a list of insertions. Every spacer \u00d7 insertion "
        "combination is designed in pure-insertion (Addition) mode. For each design, "
        "RTT + insert + PBS concatenate into that pegRNA's 3\u2032 extension. Output is a table "
        "(sorted spacer-major, then insertion-major) with a CSV download."
    )

    genomic_seq = st.text_area(
        "Input DNA (amplicon / reference)",
        height=150,
        placeholder="Paste the reference DNA (A/C/G/T; U accepted and converted to T).",
        key="batch_input_dna",
    )

    col_sp, col_ins = st.columns(2)
    with col_sp:
        spacers_text = st.text_area(
            "Spacers  \u2014  ID <tab/space/comma> sequence, one per line",
            height=170,
            value=(
                "c9sg066\tGCAGCGGCCGGGGCTGGCCACC\n"
                "c9sg068\tCTGGCAGCAGCGGCCGGGGCTG\n"
                "c9sg069\tTGCTGGGTAGAGGTGGCCAGCC"
            ),
            key="batch_spacers",
        )
    with col_ins:
        insertions_text = st.text_area(
            "Insertions  \u2014  ID <tab/space/comma> sequence, one per line",
            height=170,
            value=("ETS-ETS\tACCGGAAGTAGCACCGGAAGTA\nETS\tACCGGAAGTA"),
            key="batch_insertions",
        )

    st.markdown("### Insertion placement")
    place_col, pos_col = st.columns([6, 4])
    with place_col:
        placement_label = st.radio(
            "Where to place each insertion",
            options=[
                "At the nick (retain 0 bases)",
                "After a specific input-DNA position",
            ],
            index=0,
            key="batch_placement",
        )
    placement_mode = (
        BATCH_PLACEMENT_NICK
        if placement_label.startswith("At the nick")
        else BATCH_PLACEMENT_POSITION
    )
    position_1based = None
    with pos_col:
        if placement_mode == BATCH_PLACEMENT_POSITION:
            position_1based = st.number_input(
                "Insert AFTER this position (1-based)",
                min_value=1,
                max_value=1000000,
                value=1,
                step=1,
                key="batch_position",
                help="1-based position in the input DNA; the insertion goes right after this base. "
                "Example: AAAAATTTTT, position 5, insertion CCCC \u2192 AAAAACCCCTTTTT. Spacers whose nick "
                "cannot reach this position (wrong side / too far) are skipped and listed in a warning.",
            )
        else:
            st.caption("Insertion goes right at each spacer's nick junction (0 retained bases).")

    orient_col, off_col = st.columns([6, 4])
    with orient_col:
        orientation_label = st.selectbox(
            "Insertion orientation (in your input DNA)",
            options=[
                "Auto \u2014 strand-aware (forward)",
                "Force forward",
                "Force reverse-complement",
            ],
            index=0,
            key="batch_orient",
            help="Auto/forward: the insertion reads forward (exactly as typed) in the edited input DNA "
            "regardless of which strand the spacer matched. Force reverse-complement flips it in the input DNA.",
        )
    insertion_orientation = {
        "Auto \u2014 strand-aware (forward)": INSERTION_ORIENTATION_AUTO,
        "Force forward": INSERTION_ORIENTATION_FORWARD,
        "Force reverse-complement": INSERTION_ORIENTATION_REVERSE,
    }[orientation_label]
    with off_col:
        nick_offset = st.number_input(
            "Nick offset",
            min_value=0,
            max_value=50,
            value=3,
            key="batch_nick_offset",
            help="Nick placed this many nt upstream of the spacer 3\u2032 end in the matched orientation.",
        )

    st.markdown("### PBS & RTT (3\u2032 homology arm) settings")
    pcol1, pcol2 = st.columns(2)
    with pcol1:
        pbs_shorter = st.slider("Include shorter PBS lengths", 0, 20, 0, key="batch_pbs_short")
    with pcol2:
        pbs_longer = st.slider("Include longer PBS lengths", 0, 20, 0, key="batch_pbs_long")

    rmode_col, r1, r2, r3 = st.columns([3, 2, 2, 2])
    with rmode_col:
        rtt_mode_label = st.radio(
            "RTT / homology-arm lengths",
            options=["Min / max / count", "Manual exact lengths"],
            index=0,
            key="batch_rtt_mode",
        )
    rtt_mode = "range" if rtt_mode_label == "Min / max / count" else "manual"
    rtt_manual_lengths = ""
    rtt_min, rtt_max, rtt_count = 10, 20, 3
    if rtt_mode == "range":
        with r1:
            rtt_min = st.number_input("Arm min", 1, 500, 10, key="batch_rtt_min")
        with r2:
            rtt_max = st.number_input("Arm max", 1, 500, 20, key="batch_rtt_max")
        with r3:
            rtt_count = st.number_input("How many", 2, 100, 3, key="batch_rtt_count")
    else:
        with r1:
            rtt_manual_lengths = st.text_input(
                "Manual arm lengths", value="10,15,20", key="batch_rtt_manual",
                help="Comma- or space-separated homology-arm lengths.",
            )

    run = st.button("Run batch design", type="primary", width="stretch")
    if not run:
        return

    if not genomic_seq.strip():
        st.error("Please paste the input DNA.")
        return
    try:
        spacers = parse_two_column_table(spacers_text, "Spacers", validate_dna=True)
        insertions = parse_two_column_table(insertions_text, "Insertions", validate_dna=False)
    except Exception as exc:
        st.error(str(exc))
        return

    try:
        df, batch_warnings = design_batch_insertion(
            genomic_seq=genomic_seq,
            spacers=spacers,
            insertions=insertions,
            placement_mode=placement_mode,
            position_1based=(int(position_1based) if position_1based is not None else None),
            nick_offset=int(nick_offset),
            pbs_shorter=int(pbs_shorter),
            pbs_longer=int(pbs_longer),
            rtt_mode=rtt_mode,
            rtt_manual_lengths=rtt_manual_lengths,
            rtt_min=int(rtt_min),
            rtt_max=int(rtt_max),
            rtt_count=int(rtt_count),
            insertion_orientation=insertion_orientation,
        )
    except Exception as exc:
        st.error(str(exc))
        return

    if batch_warnings:
        for warning in batch_warnings:
            st.warning(warning)

    if df.empty:
        st.info("No valid designs were produced \u2014 see the warnings above.")
        return

    st.success(
        f"Generated {len(df)} row(s) across {df['spacer ID'].nunique()} spacer(s) "
        f"and {df['insertion ID'].nunique()} insertion(s)."
    )
    st.dataframe(df, width="stretch", hide_index=True)

    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download CSV",
        data=csv_bytes,
        file_name="batch_pegRNA_insertion_designs.csv",
        mime="text/csv",
        type="primary",
        width="stretch",
    )


app_mode = st.radio(
    "Tool mode",
    options=["Single design", "Batch pure-insertion"],
    index=0,
    horizontal=True,
    key="app_mode",
)

if app_mode == "Batch pure-insertion":
    render_batch_mode()
    st.stop()


left_col, right_col = st.columns(2)

with left_col:
    genomic_seq = st.text_area(
        "Input DNA",
        height=220,
        placeholder="Paste genomic DNA sequence here (A/C/G/T; U is also accepted and converted to T).",
    )
    spacer = st.text_input(
        "pegRNA spacer",
        placeholder="Paste spacer here (DNA or RNA alphabet accepted).",
    )

    st.markdown("**Insertion / addition**")
    add_toggle_col, add_help_col = st.columns([2.4, 7])
    with add_toggle_col:
        addition_mode = st.toggle("Addition (pure insertion) mode", value=False)
    with add_help_col:
        with st.popover("What does Addition mode do?"):
            st.write(
                "OFF (default): 'site included / rewrite' behavior. The RTT templates the "
                "genomic sequence starting at the selected base, and you edit it yourself to "
                "encode deletions, substitutions, or insertions."
            )
            st.write(
                "ON: pure-insertion helper. The genomic bases from the nick up to and INCLUDING the "
                "selected site plus your insertion form the 'insert' column, and a 3\u2032 homology arm "
                "forms the 'RTT' column. RTT + insert + PBS concatenated = the pegRNA 3\u2032 extension. "
                "The homology-arm length(s) come from the RTT-length control below."
            )
            st.write(
                "Example (+ strand): downstream CTAGCTAG (offsets 0..7), insertion AA, site +3 (G). "
                "New DNA flap = CTAG + AA + homology. insert column = revcomp(CTAGAA) = UUCUAG; "
                "RTT column = revcomp(homology arm)."
            )

    insert_input_col, insert_help_col = st.columns([8, 0.8])
    with insert_input_col:
        insertion_sequence = st.text_input(
            "Insertion sequence",
            value="aaaaggggttttcccc",
            placeholder="Insertion sequence (A/T/C/G/U accepted).",
        )
    with insert_help_col:
        with st.popover("?"):
            st.write(
                "Type the insertion in your input-DNA (sense) orientation. In Addition mode it is "
                "templated into the pegRNA (U/u treated as T/t); in legacy mode it is only placed "
                "between RTT and PBS in the combination output."
            )

    # Addition-mode specific controls.
    insertion_orientation = INSERTION_ORIENTATION_AUTO
    insert_at_nick = False
    if addition_mode:
        orient_col, atnick_col = st.columns([6, 5])
        with orient_col:
            orientation_label = st.selectbox(
                "Insertion orientation (in your input DNA)",
                options=[
                    "Auto \u2014 strand-aware (forward)",
                    "Force forward",
                    "Force reverse-complement",
                ],
                index=0,
                help="Auto: detect the spacer strand and place the insertion so it reads forward (as typed) "
                "in your input DNA \u2014 sense spacer keeps it forward-in-input (reverse-complemented inside the RTT), "
                "antisense spacer keeps the insertion forward inside the RTT. Force forward / Force "
                "reverse-complement override that so the insertion reads forward or reverse-complemented "
                "in your input DNA regardless of strand.",
            )
        insertion_orientation = {
            "Auto \u2014 strand-aware (forward)": INSERTION_ORIENTATION_AUTO,
            "Force forward": INSERTION_ORIENTATION_FORWARD,
            "Force reverse-complement": INSERTION_ORIENTATION_REVERSE,
        }[orientation_label]
        with atnick_col:
            insert_at_nick = st.checkbox(
                "Insert exactly at the nick (retain 0 bases)",
                value=False,
                help="If ON, the insertion goes right at the nick junction and the selected site is ignored; "
                "e.g. downstream CTAG gives new DNA AA + CTAG.",
            )

    if addition_mode:
        include_insertion = False
        if not insertion_sequence.strip():
            st.warning("Addition mode requires an insertion sequence.")
        else:
            try:
                validate_insertion_sequence(insertion_sequence, allow_empty=False)
            except ValueError as exc:
                st.warning(str(exc))
            else:
                st.caption(
                    "Addition mode is ON: RTT column = 3\u2032 homology arm (reverse-complement RNA); "
                    "insert column = reverse-complement of (retained genomic + insertion); "
                    "RTT + insert + PBS = pegRNA 3\u2032 extension. RTT-length values set the homology-arm length."
                )
    else:
        include_insertion = st.checkbox(
            "Include insertion in RTT \u00d7 insert \u00d7 PBS combo output",
            value=True,
            help="Legacy option: places the insertion between RTT and PBS in the combination output only "
            "(it is NOT templated into the RTT). Ignored while Addition mode is ON.",
        )
        if insertion_sequence.strip():
            try:
                validate_insertion_sequence(insertion_sequence, allow_empty=True)
            except ValueError as exc:
                st.warning(str(exc))

with right_col:
    st.markdown("### PBS settings")
    pbs_col1, pbs_col2 = st.columns(2)
    with pbs_col1:
        pbs_shorter = st.slider(
            "Include shorter PBS lengths",
            min_value=0,
            max_value=20,
            value=0,
            help="The default PBS is the design with RNA:DNA Tm closest to 37 \u00b0C. "
            "This slider expands the output set to include shorter PBS lengths around that default. "
            "Minimum PBS length=5.",
        )
    with pbs_col2:
        pbs_longer = st.slider(
            "Include longer PBS lengths",
            min_value=0,
            max_value=20,
            value=0,
            help="The default PBS is the design with RNA:DNA Tm closest to 37 \u00b0C. "
            "This slider expands the output set to include longer PBS lengths around that default. "
            "Maximum PBS length=25.",
        )

    if addition_mode:
        st.markdown("### RTT settings \u2014 homology arm (Addition mode)")
        st.caption(
            "In Addition mode these length values set the 3\u2032 homology-arm length (the RTT column). "
            "The insertion and retained genomic live in the 'insert' column."
        )
    else:
        st.markdown("### RTT settings")

    rtt_mode_label_col, rtt_mode_help_col, rtt_mode_widget_col = st.columns([2, 0.8, 6])
    with rtt_mode_label_col:
        st.markdown("**RTT mode**")
    with rtt_mode_help_col:
        with st.popover("?"):
            st.write(
                "Min / max / count returns evenly spaced lengths including the minimum and maximum. "
                "Manual exact lengths returns only the lengths you specify. In Addition mode these are "
                "3\u2032 homology-arm lengths."
            )
    with rtt_mode_widget_col:
        rtt_mode = st.radio(
            "RTT mode",
            options=["Min / max / count", "Manual exact lengths"],
            index=0,
            horizontal=True,
            label_visibility="collapsed",
        )

    rtt_start_label_col, rtt_start_help_col, rtt_start_widget_col = st.columns([2, 0.8, 6])
    with rtt_start_label_col:
        st.markdown("**Insertion site**" if addition_mode else "**RTT start**")
    with rtt_start_help_col:
        with st.popover("?"):
            if addition_mode:
                st.write(
                    "Selected site: the insertion is placed right after the highlighted base you choose "
                    "in the viewer below. Start from nick site: insertion is placed after the first base "
                    "downstream of the nick (offset 0). Use the 'Insert exactly at the nick' checkbox to "
                    "retain zero bases."
                )
            else:
                st.write(
                    "Selected RTT start uses the highlighted base you choose in the viewer below as the "
                    "first base included in RTT. Start from nick site locks RTT to begin at the "
                    "nick-defined start."
                )
    with rtt_start_widget_col:
        rtt_start_mode = st.radio(
            "RTT start",
            options=["Selected site", "Start from nick site"],
            index=0,
            horizontal=True,
            label_visibility="collapsed",
        )

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        nick_offset = st.number_input(
            "Nick offset",
            min_value=0,
            max_value=50,
            value=3,
            help="Nick is placed this many nucleotides upstream of the spacer 3\u2032 end in the matched orientation.",
        )
    if rtt_mode == "Min / max / count":
        with c2:
            rtt_min = st.number_input("RTT min", min_value=1, max_value=500, value=10)
        with c3:
            rtt_max = st.number_input("RTT max", min_value=1, max_value=500, value=20)
        with c4:
            rtt_count = st.number_input("How many to design", min_value=2, max_value=100, value=3)
        rtt_manual_lengths = ""
    else:
        with c2:
            rtt_manual_lengths = st.text_input(
                "Manual RTT lengths",
                value="10,15,20",
                placeholder="e.g. 10,15,20",
                help="Comma- or space-separated lengths.",
            )
        with c3:
            st.empty()
        with c4:
            st.empty()
        rtt_min = 10
        rtt_max = 20
        rtt_count = 3

selected_rtt_start_target = None
if genomic_seq.strip() and spacer.strip():
    try:
        selector = build_rtt_start_selector(
            genomic_seq=genomic_seq,
            spacer=spacer,
            nick_offset=int(nick_offset),
        )
    except Exception as exc:
        st.info(f"Selector preview unavailable: {exc}")
    else:
        signature = _selector_signature(genomic_seq, spacer, int(nick_offset))
        selected_rtt_start_target = render_rtt_start_selector(
            selector,
            signature,
            show_buttons=(rtt_start_mode == "Selected site"),
            addition_mode=addition_mode,
        )

run = st.button("Run design", type="primary")

if run:
    try:
        result = design_pbs_rtt(
            genomic_seq=genomic_seq,
            spacer=spacer,
            nick_offset=int(nick_offset),
            pbs_shorter=int(pbs_shorter),
            pbs_longer=int(pbs_longer),
            rtt_mode=("range" if rtt_mode == "Min / max / count" else "manual"),
            rtt_manual_lengths=rtt_manual_lengths,
            rtt_min=int(rtt_min),
            rtt_max=int(rtt_max),
            rtt_count=int(rtt_count),
            rtt_start_mode=("selected" if rtt_start_mode == "Selected site" else "nick"),
            rtt_start_target=selected_rtt_start_target,
            include_insertion=bool(include_insertion),
            insertion_sequence=insertion_sequence,
            addition_mode=bool(addition_mode),
            insertion_orientation=insertion_orientation,
            insert_at_nick=bool(insert_at_nick),
        )
    except Exception as exc:
        st.error(str(exc))
    else:
        st.success("Design completed.")
        st.write(result.summary_line)
        if result.warnings:
            for warning in result.warnings:
                st.warning(warning)
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("### RTT" + (" (3\u2032 homology arm)" if result.addition_mode else ""))
            st.code(result.rtt_text, language=None)
            st.dataframe(result.rtt_df, width="stretch")
        with c2:
            st.markdown("### PBS")
            st.code(result.pbs_text, language=None)
            st.dataframe(result.pbs_df, width="stretch")
        if result.addition_mode and result.addition_insert_rna:
            st.markdown("### insert (retained genomic + insertion, reverse-complemented)")
            st.code(result.addition_insert_rna, language=None)
        st.markdown("### RTT x PBS")
        st.code(result.rtt_pbs_text, language=None)
        if result.rtt_insert_pbs_text:
            st.markdown("### RTT x insert x PBS")
            if result.addition_mode:
                st.caption("Columns are RTT (homology) \u00b7 insert (retained+insertion, revcomp) \u00b7 PBS. "
                           "Concatenate a row left-to-right to get that pegRNA's 3\u2032 extension.")
            st.code(result.rtt_insert_pbs_text, language=None)
