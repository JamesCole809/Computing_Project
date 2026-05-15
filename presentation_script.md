# Presentation plan (~8-9 min, two speakers)

His /pop/ advice: it's a performance, slides are just visual aids, less
is more, rehearse it, have the tech tested. The brief says "15 minute"
but he also says shorter is better — we aim for a tight 8-9 min that
actually shows the thing working, rather than padding.

These are talking points, NOT a word-for-word script. Say them in your
own words so it sounds like us, not a read-out.

## Before recording

- [ ] Slides done (7 slides max, see below).
- [ ] `pytest test_project.py -q` green, `doctest README.md` silent.
- [ ] `demo_project.ipynb` open on a fresh kernel, font zoomed up so
      text and cursor are clearly visible on the recording.
- [ ] Mic + screen record + webcam tested on a 20-sec sample.
- [ ] One full practice run with a stopwatch.

Record in sections and stitch in a video editor. Each person records
their own parts whenever — no need to be in the same room. Because the
demo is pre-recorded, nothing can crash live, so that doubles as the
"backup plan" /pop/ asks for. If a take is fluffed, redo just that
section.

## Visuals / editing rule

His rule: visual aids must SUPPORT the performance, not replace it.
A muted screen recording with no voice = a "manuscript" = the failure
case he names. So: your voice runs on top of every visual, the whole
time, explaining what and why.

Three shot types:

- **Talking head, full screen** — face only. Used for intro and close
  (human connection, "spark interest").
- **Slide, full screen** — only while you're actively pointing at it.
  Cut away the second you're done; never leave a dead slide up.
- **Screen capture + webcam picture-in-picture** — demo/terminal fills
  the screen, your face is a small box in a corner so you're still
  "present". This is the main demo shot.

Match the visual to the sentence: never cut to a visual before you
start talking about it, and don't read a slide aloud (that's the
"adds nothing" trap — the slide shows it, you explain it).

## Slides (keep them near-empty)

1. Title — names, module, date.
2. The problem — one image (revising / textbook).
3. The seven ODE kinds — bullet list.
4. `generate -> solve -> Solution` sketch.
5. "Demo" — placeholder, nothing else.
6. Testing — three numbers: 118 pytest, 19 doctest, 7 kinds.
7. Comparison table + GitHub URL + thanks.

## Running order

### 0:00-1:00 — Intro (James)
**On screen:** talking head, full screen (no slide).
- Who we are, module, one-line: a Python library that generates and
  solves ODE practice problems for first year Cardiff Maths students.
- The gap: notes/textbook have a fixed number of examples; once you've
  done them you're stuck; Wolfram Alpha hides working behind a paywall.

### 1:00-2:30 — What it does (Gruff)
**On screen:** slide 3 then slide 4, full screen. Cut back to talking
head between points so it isn't a static slide for 90 sec.
- Seven ODE kinds (point at slide 3).
- Each kind has a generator (random valid question, 3 difficulties) and
  a solver (walks the lecture method, returns each step).
- Top level: `generate` -> Problem, `solve` -> Solution with answer +
  steps + a verified flag.

### 2:30-5:30 — Demo (James, in the notebook)
**On screen:** full-screen screen capture of the notebook, webcam as a
small picture-in-picture in a corner. Voice runs the whole time.
- Import as `ml`.
- `ml.generate("ode1_bernoulli", difficulty="medium")` — show LaTeX question.
- `ml.solve(problem)` — show answer.
- `result.show_steps()` — this is the point: full working, same shape
  as on paper.
- `result.verified` -> True (plugged back into the ODE and checked).
- `ml.problem_sheet(n=10, kind="mixed", export=True)` — open the PDF for
  two seconds.
- Mention `ml.launch_web()` exists (only click through it if time).

### 5:30-6:30 — How a solver works (Gruff)
**On screen:** the one Bernoulli-formula slide. Cut away from it the
moment the point is made.
- Take the Bernoulli one just shown. Substitution
  `z = y^(1-n) e^((1-n) int a dx)`, y terms cancel, integrate in x.
- Key point: we don't just call `sympy.dsolve` and dump the answer —
  each step is a separate object so the student sees the method.

### 6:30-7:30 — Testing + docs (James)
**On screen:** screen capture of the terminal (`118 passed`, then the
silent doctest run), webcam PiP in a corner.
- 118 pytest covering every kind/difficulty + edge cases (resonance,
  repeated roots, Bernoulli domain).
- 19 doctests: every README example is runnable, docs can't go stale.
- README follows Diataxis (Tutorial / How-to / Explanation / Reference).

### 7:30-8:30 — Comparison + close (Gruff)
**On screen:** slide 7 (comparison table) while comparing, then cut to
talking head full screen for the close.
- Not replacing SymPy (we use it). SymPy gives the answer, no working.
  Wolfram Alpha paywalls steps. Textbooks have a fixed pool. None
  generate random practice questions with full working.
- Future: more ODE classes, fuller step derivations.
- GitHub link, leaving it public for future first years. Thanks.

## On the day
- Record each section separately; assemble in the editor afterwards.
- Stay roughly to timings; if you overrun, cut the web UI bit first.
- Talk to the points in your own words — don't read this aloud.
- In the edit: keep voice over every visual, match the cut to the
  sentence, never leave a dead slide up.
