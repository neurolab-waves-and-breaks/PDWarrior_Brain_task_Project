# BRAIN Task – Python Implementation

This folder contains the **Python / Pygame implementation** of the
Bradykinesia Akinesia Incoordination (BRAIN) alternating finger-tapping task,
designed to assess upper limb motor function.

The task records high-resolution timing of key presses during a fixed time window
and exports raw keystroke data for downstream preprocessing and analysis.

---

## Task description

Participants are instructed to:

- Press the **“S” and “;” keys in alternation** using a **single finger**
- Start tapping immediately after the **auditory “GO” cue**
- Stop tapping when the **auditory stop cue (trumpet sound)** is played

Participants alternately strike the “S” and “;” keys on a laptop keyboard as
quickly and accurately as possible for **30 seconds**, using a single finger of
their choice. The centres of the two keys are separated by **15.5 cm**, and the
keys are highlighted with green tape to ensure clear visual distinction from
surrounding keys.

### Task structure
- Optional training block (15 s)
- Recording block (30 s)
- Left and right hands tested separately

---

## Purpose

- Present a timed finger-tapping task to participants
- Record:
  - Key press onset times
  - Key dwell times
  - Pressed key identity
- Save raw keystroke data in `.csv` format for subsequent preprocessing in MATLAB

---

## Folder contents

01_brain-task/
├── brain_task.py # Main task script
├── Go.wav # Auditory go cue
├── Trumpet.wav # Auditory stop cue
└── README.md


---

## Requirements

- Python ≥ 3.8
- Python packages:
  - pygame
  - pandas
  - numpy

---

## Running the task

Before running the task, set participant-specific parameters at the top of
`brain_task.py`:

```python
ID = 'P1'              # Participant ID
TRAINING_ON = 0        # 1 = training (15 s), 0 = recording (30 s)
hand = 'right'         # 'right' or 'left'
exerc = 'preExerc'     # 'preExerc' or 'postExerc'

Running the script will start the task and automatically save a .csv file
containing raw keystroke timing data for that session.