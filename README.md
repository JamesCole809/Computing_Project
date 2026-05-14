# Differential Equation Practice Generator

A small Python library that makes random differential equation questions and
solves them step by step. Built for first year Cardiff Mathematics students
who want extra practice past what their textbook gives them.

The library covers seven kinds of ODE:

- First order linear
- Separable
- Bernoulli
- Homogeneous (using the substitution $z = y/x$)
- Second order with constant coefficients (homogeneous)
- Second order with constant coefficients (inhomogeneous)
- Coupled first order linear ODEs

## Installation

The library is a plain Python package, no install step. Just clone the repo
and make sure the dependencies are there:

```
pip install sympy numpy matplotlib flask pytest
```

Then run Python from the repo root so that `project.py` is on the path.

## Tutorial

This section walks through the main things the library can do. The examples
below are real Python sessions and can be checked automatically with
`python -m doctest README.md`.

### Importing

```pycon
>>> import project as ml

```

### Listing what's available

```pycon
>>> ml.list_topics()
['ode']
>>> len(ml.list_kinds())
7
>>> 'ode1_linear' in ml.list_kinds()
True

```

### Generating a question

`ml.generate` builds a random problem of a given kind. We pass `seed` here
so the result is reproducible for the doctest.

```pycon
>>> p = ml.generate("ode1_linear", difficulty="easy", with_ics=True, seed=0)
>>> p.kind
'ode1_linear'
>>> bool(p.prompt)
True

```

In a Jupyter notebook the problem will render as nice LaTeX automatically.
Outside a notebook you can read `p.prompt` (plain text) or `p.prompt_latex`
(the LaTeX version).

### Solving

`ml.solve` does the working and returns a `Solution` object. It has the
final answer, the list of steps and a flag saying whether the answer was
verified by plugging it back into the original ODE.

```pycon
>>> r = ml.solve(p)
>>> r.kind
'ode1_linear'
>>> r.verified
True
>>> r.verify_msg
'passed'
>>> len(r.steps) > 0
True

```

If you are in Jupyter, `r.show_steps()` displays each step with rendered
LaTeX. Outside Jupyter it falls back to plain text.

### Quick verification

`ml.check` is a shortcut that just confirms the answer is correct.

```pycon
>>> rep = ml.check(p, r)
>>> rep['verified']
True
>>> rep['message']
'passed'

```

### Making a problem sheet

`ml.problem_sheet` produces a list of questions plus their solutions. With
`export=True` it also writes a PDF of the sheet.

```pycon
>>> qs = ml.problem_sheet(n=3, kind="ode1_linear", difficulty="easy", seed=0)
>>> len(qs)
3
>>> qs[0]['problem'].kind
'ode1_linear'
>>> qs[0]['solution'].verified
True

```

## How to guides

### How to ask for a specific kind of ODE

Pass any of the names returned by `list_kinds()` as the `kind`. For a random
choice across all kinds pass `"mixed"`:

```python
p = ml.generate("ode1_bernoulli", difficulty="medium")
p = ml.generate("mixed")
```

### How to control the difficulty

Each generator supports `difficulty="easy"`, `"medium"` or `"hard"`. Easy
sticks to integer coefficients, medium adds things like $1/x$ and trig,
hard mixes them all together.

```python
p = ml.generate("ode2_cc_hom", difficulty="hard")
```

### How to generate a problem without an initial condition

```python
p = ml.generate("ode1_separable", with_ics=False)
```

### How to make the same problem every time

Pass a `seed`:

```python
p = ml.generate("ode1_linear", seed=42)
```

### How to export a problem sheet as a PDF

```python
ml.problem_sheet(n=10, kind="mixed", difficulty="mixed", export=True)
```

By default the PDF is written next to the library. To send it somewhere
else pass `saveto`:

```python
ml.problem_sheet(n=10, kind="ode2_cc_hom", export=True, saveto="my_sheets")
```

### How to launch the browser UI

```python
ml.launch_web()
```

This starts a small Flask server (default `http://127.0.0.1:8080`) which
gives you a button-driven page for generating, solving and plotting.

### How to get a LaTeX block of a worked solution

```python
out = ml.export_latex(result, problem)
```

The returned string can be pasted straight into a `.tex` file.

## Explanation

### Why this library and not something else

There are already plenty of tools that solve ODEs, so it is worth saying
upfront what makes this one different and what it is *not* trying to
replace.

**SymPy** can solve a huge range of ODEs through `sympy.dsolve`, and we
use it under the hood for the actual symbolic work. The problem is that
`dsolve` only gives you the final answer. It does not show the method,
which is the bit that students actually need to learn. It also still
expects you to come up with the question yourself.

**Wolfram Alpha** will solve an ODE if you type it in, and the website
will offer step-by-step working too, but the steps sit behind a Pro
subscription. It also is not scriptable and is not built around
generating fresh questions in bulk.

**Maxima**, **SageMath**, **Maple** and **Mathematica** are full
computer algebra systems. They will solve an ODE and most of them can
show some workings, but they are heavyweight tools aimed at researchers
and they take a lot of setup. None of them are designed around a first
year student wanting to grind through twenty separable equations the
night before an exam.

**Symbolab** and similar online tutors come closest in spirit but they
are proprietary, you cannot ask them for a *random* question, and you
cannot export a problem sheet PDF in one line.

**Textbooks** such as Boyce and DiPrima are still the standard reference
but each chapter only has a fixed number of practice problems. Once you
have done them once you have to re-read the same questions, and there
is no way to get a brand new ODE in the same difficulty band.

This library deliberately covers far less ground than any of the systems
above, only the seven ODE types taught in our first year module. In
return it does two things none of them do well:

1. It **generates** a valid random question of a chosen kind and
   difficulty, so practice problems do not run out.
2. It **shows the method**, not just the answer, in the same form that
   the lecture notes use.

It is meant to sit next to those bigger tools rather than compete with
them. If you want to solve one specific arbitrary ODE, use SymPy or
Wolfram Alpha. If you want twenty fresh practice problems with full
working at exam revision time, use this.

### How the solver works internally

Each ODE kind has its own module inside `diff_eq/`. Every module has a
**generator** that picks random coefficients of a known solvable form,
and a **solver** that walks through the standard method we were taught
in lectures. The solver uses SymPy under the hood for the actual symbolic
work (integration, root finding, substitution) but the high level method
is written out by hand so that the steps shown to the user match the
method a student would follow on paper.

After the answer is found the solver plugs it back into the original ODE
and checks that the residual simplifies to zero. If an initial condition
was given it also checks that. The result of the check is stored in
`Solution.verified`.

### Why these seven kinds

These are the seven types covered in our first year Differential Equations
module so the library lines up with what a student would actually be asked
to revise. Picking too many types would have made the codebase huge for
this project, picking fewer would have left obvious gaps.

### Limitations

- The step explanations are short. They name the method and show the
  intermediate expression but do not derive each line in full.
- Only the seven ODE kinds listed above are supported. Things like exact
  equations, Laplace transform methods and higher order systems are not.
- The system solver assumes the matrix is non-degenerate enough that you
  can eliminate one variable. It does not fall back to eigenvector methods.

## Reference

### Public functions

| Function | What it does |
|----------|--------------|
| `generate(kind, difficulty, with_ics, seed)` | Build a random Problem of the given kind. |
| `solve(problem, want_steps, want_verify)` | Solve a Problem and return a Solution. |
| `check(problem, result)` | Quick correctness check; returns `{"verified": bool, "message": str}`. |
| `export_latex(result, problem=None)` | Return a LaTeX block of the worked solution. |
| `problem_sheet(n, kind, difficulty, export, saveto, seed)` | Build a list of questions, optionally export a PDF. |
| `list_kinds()` | All seven ODE kind names. |
| `list_topics()` | Topic groupings (just `'ode'` for now). |
| `launch_web(host, port, open_browser)` | Start the Flask front end. |
| `help()` | Print a short summary of the API. |

### Problem attributes

| Attribute | Meaning |
|-----------|---------|
| `.kind` | The ODE kind, one of `list_kinds()`. |
| `.prompt` | Plain text version of the question. |
| `.prompt_latex` | LaTeX version of the question. |
| `.data` | Internal dictionary used by the solvers. |

### Solution attributes

| Attribute | Meaning |
|-----------|---------|
| `.kind` | Same kind as the problem. |
| `.final_answer` | String form of the answer. |
| `.answer_expr` | SymPy expression of the answer. |
| `.steps` | List of `Step` objects. |
| `.verified` | True if the answer satisfies the ODE and ICs. |
| `.verify_msg` | `'passed'` or the reason it failed. |
| `.warnings` | Any domain or branch warnings. |
| `.show_steps()` | Display each step rendered in LaTeX (Jupyter). |

### ODE kinds

| Name | Standard form |
|------|---------------|
| `ode1_linear` | $y' + p(x) y = q(x)$ |
| `ode1_separable` | $g(y) y' = f(x)$ |
| `ode1_bernoulli` | $y' + a(x) y = b(x) y^n$ |
| `ode1_homogeneous_sub` | $y' = f(y/x)$ |
| `ode2_cc_hom` | $a y'' + b y' + c y = 0$ |
| `ode2_cc_inhom` | $a y'' + b y' + c y = f(x)$ |
| `ode_sys2_linear` | $u' = a u + b v,\; v' = c u + d v$ |

### Tests

Two test commands:

```
pytest test_project.py
python -m doctest README.md
```

`pytest` runs the full unit test suite (118 tests covering generators,
solvers, verification, edge cases). `doctest` re-runs every `>>>` example
in this file and checks the output still matches.

### Bibliography

1. Meurer, A. et al. *SymPy: symbolic computing in Python*, PeerJ Computer
   Science, 2017. <https://www.sympy.org>
2. Maxima, a Computer Algebra System. <https://maxima.sourceforge.io>
3. The Sage Developers. *SageMath, the Sage Mathematics Software System*.
   <https://www.sagemath.org>
4. Wolfram Research. *Mathematica*. <https://www.wolfram.com/mathematica>
5. Wolfram Alpha LLC. *Wolfram Alpha computational engine*.
   <https://www.wolframalpha.com>
6. Boyce, W. E. and DiPrima, R. C. *Elementary Differential Equations and
   Boundary Value Problems*, 11th ed., Wiley, 2017.
7. Procida, D. *The Diataxis documentation framework*.
   <https://diataxis.fr>
8. Knight, V. *Nashpy: A Python library for the computation of Nash
   equilibria*, Journal of Open Source Software, 3(30), 904, 2018.
9. Wilde, H. and Knight, V. *Matching: A Python library for solving
   matching games*, Journal of Open Source Software, 5(48), 2169, 2020.
