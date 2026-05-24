
# PBS/RTT Designer Streamlit App

A small Streamlit app for designing PBS and RTT sequences from a target DNA sequence and a pegRNA spacer.

## What this app does

This app takes:
- a genomic DNA sequence
- a pegRNA spacer

It then:
- finds the spacer in the input DNA on either strand
- requires a unique match
- places the nick site using a configurable offset from the spacer 3′ end
- designs PBS candidates around the default PBS closest to 37 °C
- designs RTT candidates from either:
  - a selected RTT start site
  - the nick site
- returns output sequences in RNA format (T converted to U)

## Main design logic

### Spacer QC
- The spacer is matched against the input DNA in both forward and reverse-complement orientation
- No PAM is required
- If there is no match, the app returns an error
- If there is more than one match, the app returns an error and asks the user to disambiguate

### Nick site
- The nick site is defined relative to the matched spacer
- Default behavior: nick is 3 nt upstream of the spacer 3′ end
- This offset is user-adjustable

### PBS design
- PBS lengths are evaluated from 5 to 25 nt
- The default PBS is the sequence whose RNA:DNA Tm is closest to 37 °C
- Users can expand the PBS output to include shorter and longer PBS designs around that default

### RTT design
Two RTT modes are supported:

1. Min / max / count
- Input a minimum RTT length
- Input a maximum RTT length
- Input how many RTT designs to return
- The app generates evenly spaced RTT lengths including the minimum and maximum
- If spacing is not an integer, values are rounded up

2. Manual exact lengths
- Input specific RTT lengths directly

Two RTT start modes are supported:

1. Selected RTT start
- User chooses the first base included in the RTT from the displayed sequence ribbon

2. Start from nick site
- RTT starts directly from the nick-defined position

### Output
The app returns:
- PBS sequences in RNA format
- RTT sequences in RNA format
- table output for visualization
- copy-friendly sequence output

## UI notes

### Layout
The app uses a 2-column input layout:
- left column: input DNA and pegRNA spacer
- right column: PBS settings and RTT settings

### RTT selector
The RTT selector shows:
- a one-row sequence ribbon
- spacer bases highlighted in red
- selected RTT start highlighted in yellow
- nick site marked by blue underline

When RTT start mode is set to Start from nick site, the selector buttons are hidden and the ribbon remains visible for reference.

## File structure

- `streamlit_app.py` or versioned app file such as `streamlit_app.py`: Streamlit UI
- `pbs_rtt_designer_core.py`: core PBS/RTT design logic
- `pbs_tm.py`: PBS Tm calculation
- `requirements.txt`: Python dependencies

## Run locally

```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

If using a versioned app file:

```bash
streamlit run streamlit_app.py
```

## Input expectations

### DNA sequence
- DNA alphabet expected
- Case-insensitive
- Non-sequence formatting should be removed before input

### pegRNA spacer
- DNA alphabet expected
- The app handles forward and reverse-complement matching internally

## Important assumptions
- Spacer matching does not require PAM
- Spacer match must be unique
- PBS and RTT outputs are reported as RNA
- RTT start selection defines the first base included in RTT

## Typical workflow
1. Paste the genomic DNA sequence
2. Paste the pegRNA spacer
3. Adjust PBS shorter/longer range if needed
4. Choose RTT mode
5. Choose RTT start mode
6. If using Selected RTT start, pick the first RTT base from the selector
7. Run the design
8. Copy PBS and RTT outputs for downstream use

## Notes
This app is intended as a focused sub-tool for PBS/RTT enumeration and visualization. It does not model edit type directly and does not infer edit boundaries from edited sequence input.
