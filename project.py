import sympy as sp
import re

class TutorError(Exception):
    pass

class Problem:
    def __init__(self, kind, prompt, data=None, metadata=None):
        self.kind = kind
        self.prompt = prompt
        self.data = data or {}
        self.metadata = metadata or {}

    def _repr_latex_(self):
        return f"${getattr(self, 'prompt_latex', self.prompt)}$"


class Step:
    def __init__(self, title, math, explanation):
        self.title = title
        self.math = math
        self.explanation = explanation

    def _repr_latex_(self):
        return f"${self.math}$"


class Solution:
    def __init__(self, kind, final_answer, answer_expr=None,
                 steps=None, verified=None, verify_msg="",
                 warnings=None):
        self.kind = kind
        self.final_answer = final_answer
        self.answer_expr = answer_expr
        self.steps = steps or []
        self.verified = verified
        self.verify_msg = verify_msg
        self.warnings = warnings or []

    def _repr_latex_(self):
        if self.answer_expr is not None:
            return f"${clean_latex(sp.latex(self.answer_expr))}$"
        return f"${self.final_answer}$"

    def show_steps(self):
        try:
            from IPython.display import display, Math, Markdown
            for s in self.steps:
                display(Markdown(f"**{s.title}**"))
                display(Math(s.math))
                if s.explanation:
                    display(Markdown(f"*{s.explanation}*"))
        except ImportError:
            for s in self.steps:
                print(f"{s.title}: {s.math}")
                if s.explanation:
                    print(f"  {s.explanation}")


def make_ic(func, order, x0, value):
    return {"func": func, "order": order, "x0": x0, "value": value}

def make_problem_data(x_symbol, functions, equations, ics=None):
    return {
        "x_symbol": x_symbol,
        "functions": functions,
        "equations": equations,
        "ics": ics or [],
    }

def extract_symbols(problem):
    d = problem.data
    x = sp.Symbol(d["x_symbol"])
    if len(d["functions"]) == 1:
        y = sp.Function(d["functions"][0])(x)
        return x, y
    else:
        funcs = [sp.Function(f)(x) for f in d["functions"]]
        return x, funcs

def apply_ics(general_sol, ics, x, constants):
    if not ics:
        return general_sol, None
    if len(ics) == 1:
        ic = ics[0]
        eq = general_sol.subs(x, ic["x0"]) - ic["value"]
        vals = sp.solve(eq, constants[0])
        if vals:
            return sp.simplify(general_sol.subs(constants[0], vals[0])), {constants[0]: vals[0]}
    elif len(ics) >= 2:
        eq1 = general_sol.subs(x, ics[0]["x0"]) - ics[0]["value"]
        eq2 = sp.diff(general_sol, x).subs(x, ics[1]["x0"]) - ics[1]["value"]
        vals = sp.solve([eq1, eq2], list(constants))
        if vals:
            return sp.simplify(general_sol.subs(vals)), vals
    return general_sol, None

def apply_ics_system(u_sol, v_sol, ics, x, func_names, constants):
    if not ics or len(ics) < 2:
        return u_sol, v_sol, None
    ic_u = None
    ic_v = None
    for ic in ics:
        if ic["func"] == func_names[0]:
            ic_u = ic
        elif ic["func"] == func_names[1]:
            ic_v = ic
    if ic_u and ic_v:
        eq_a = u_sol.subs(x, ic_u["x0"]) - ic_u["value"]
        eq_b = v_sol.subs(x, ic_v["x0"]) - ic_v["value"]
        vals = sp.solve([eq_a, eq_b], list(constants))
        if vals:
            return sp.simplify(u_sol.subs(vals)), sp.simplify(v_sol.subs(vals)), vals
    return u_sol, v_sol, None

def verify_scalar(problem, sol):
    d = problem.data
    x = sp.Symbol(d["x_symbol"])
    eq = d["equations"]
    ics = d.get("ics", [])
    func_name = d["functions"][0]

    if isinstance(sol, sp.Eq):
        y_sym = sp.Symbol("y")
        F = sol.lhs - sol.rhs
        Fy = sp.diff(F, y_sym)
        if Fy == 0:
            return False, "Can't verify implicit solution"
        dy_dx = -sp.diff(F, x) / Fy
        y_func = sp.Function(func_name)
        ode_expr = eq.lhs - eq.rhs if isinstance(eq, sp.Eq) else eq
        residual = sp.simplify(ode_expr.subs(y_func(x).diff(x), dy_dx).subs(y_func(x), y_sym))
        ode_ok = (residual == 0)
        for ic in ics:
            if ic["order"] == 0:
                if sp.simplify(F.subs({x: ic["x0"], y_sym: ic["value"]})) != 0:
                    return False, "IC check failed"
    else:
        check = sp.checkodesol(eq, sp.Eq(sp.Function(func_name)(x), sol))
        ode_ok = check[0]
        for ic in ics:
            if ic["order"] == 0:
                if sp.simplify(sol.subs(x, ic["x0"]) - ic["value"]) != 0:
                    return False, "IC check failed"
            elif ic["order"] == 1:
                if sp.simplify(sp.diff(sol, x).subs(x, ic["x0"]) - ic["value"]) != 0:
                    return False, "IC check failed"

    return ode_ok, "passed" if ode_ok else "ODE check failed"

def verify_system(problem, u_sol, v_sol):
    d = problem.data
    x = sp.Symbol(d["x_symbol"])
    funcs = d["functions"]
    eqs = d["equations"]
    ics = d.get("ics", [])

    u_func = sp.Function(funcs[0])(x)
    v_func = sp.Function(funcs[1])(x)

    res1 = sp.simplify(sp.diff(u_sol, x) - eqs[0].rhs.subs(u_func, u_sol).subs(v_func, v_sol))
    res2 = sp.simplify(sp.diff(v_sol, x) - eqs[1].rhs.subs(u_func, u_sol).subs(v_func, v_sol))
    if res1 != 0 or res2 != 0:
        res1 = sp.simplify(res1.rewrite(sp.exp))
        res2 = sp.simplify(res2.rewrite(sp.exp))
        if res1 != 0 or res2 != 0:
            return False, "ODE check failed"

    for ic in ics:
        if ic["func"] == funcs[0]:
            if sp.simplify(u_sol.subs(x, ic["x0"]) - ic["value"]) != 0:
                return False, "IC check failed"
        elif ic["func"] == funcs[1]:
            if sp.simplify(v_sol.subs(x, ic["x0"]) - ic["value"]) != 0:
                return False, "IC check failed"

    return True, "passed"

def maximal_domain(problem, sol):
    d = problem.data
    x = sp.Symbol(d["x_symbol"])
    ics = d.get("ics", [])
    excludes = d.get("domain_excludes", [])

    if isinstance(sol, sp.Eq):
        return None, ""

    try:
        full_domain = sp.calculus.util.continuous_domain(sol, x, sp.S.Reals)
    except Exception:
        return None, ""

    for pt in excludes:
        full_domain = full_domain - sp.FiniteSet(pt)

    if ics:
        x0 = ics[0]["x0"]
        if isinstance(full_domain, sp.Union):
            for part in full_domain.args:
                if isinstance(part, sp.Interval) and part.contains(x0):
                    return part, "Interval containing IC"
        if full_domain.contains(x0):
            return full_domain, "Where the solution is defined"

    return full_domain, "Where the solution is defined"

def add_domain_and_verify(steps, final_sol, want_steps, want_verify, problem):
    verified = None
    verify_msg = ""

    if final_sol is not None and not isinstance(final_sol, sp.Eq):
        domain, reason = maximal_domain(problem, final_sol)
        if domain is not None and want_steps:
            steps.append(Step("Maximal domain",
                rf"\text{{Domain: }} {sp.latex(domain)}", reason))

    if want_verify and final_sol is not None:
        verified, verify_msg = verify_scalar(problem, final_sol)

    return verified, verify_msg

def clean_latex(s):
    for f in "yuv":
        s = s.replace(rf'\frac{{d^{{2}}}}{{d x^{{2}}}} {f}{{\left(x \right)}}', rf'\frac{{d^{{2}}{f}}}{{dx^{{2}}}}')
        s = s.replace(rf'\frac{{d}}{{d x}} {f}{{\left(x \right)}}', rf'\frac{{d{f}}}{{dx}}')
    s = re.sub(r'([yuvz])(\^\{[^}]*\})?\{\\left\(x \\right\)\}', r'\1\2', s)
    return s

def to_latex(expr):
    return clean_latex(sp.latex(expr))

def latex_prompt(problem):
    d = problem.data
    eqs = d["equations"]
    if isinstance(eqs, list):
        latex = ", \\quad ".join(to_latex(eq) for eq in eqs)
    else:
        latex = to_latex(eqs)
    ics = d.get("ics", [])
    if ics:
        ic_parts = []
        for ic in ics:
            f, x0, val = ic["func"], sp.latex(ic["x0"]), sp.latex(ic["value"])
            if ic["order"] == 0:
                ic_parts.append(f"{f}({x0}) = {val}")
            else:
                ic_parts.append(f"{f}'({x0}) = {val}")
        latex += ", \\quad " + ", \\; ".join(ic_parts)
    return latex

from diff_eq import TOPICS, SOLVERS, GENERATORS

def generate(kind, difficulty="easy", with_ics=True, seed=None):
    import random
    if seed is not None:
        random.seed(seed)
    if kind == "mixed":
        kind = random.choice(list_kinds())
    if kind not in GENERATORS:
        raise TutorError(f"Unknown: {kind}")
    problem = GENERATORS[kind](difficulty, with_ics)
    problem.prompt_latex = latex_prompt(problem)
    return problem

def solve(problem, want_steps=True, want_verify=True):
    if problem.kind not in SOLVERS:
        raise TutorError(f"Unknown: {problem.kind}")
    return SOLVERS[problem.kind](problem, want_steps, want_verify)

def check(problem, result):
    d = problem.data
    if len(d["functions"]) == 1:
        ok, msg = verify_scalar(problem, result.answer_expr)
    else:
        parts = result.answer_expr
        ok, msg = verify_system(problem, parts[0], parts[1])
    return {"verified": ok, "message": msg}

def export_latex(result, problem=None):
    lines = []
    if problem is not None:
        lines.append(r"\textbf{Problem}")
        lines.append(getattr(problem, "prompt_latex", problem.prompt))
        lines.append("")
    lines.append(r"\textbf{Solution}")
    for s in result.steps:
        lines.append(rf"\textbf{{{s.title}}}")
        lines.append(s.math)
        if s.explanation:
            lines.append(rf"\text{{{s.explanation}}}")
        lines.append("")
    if result.answer_expr is not None:
        lines.append(r"\textbf{Answer}")
        lines.append(clean_latex(sp.latex(result.answer_expr)))
        lines.append("")
    if result.verified is not None:
        status = "Verified" if result.verified else "Not verified"
        lines.append(rf"\text{{{status}: {result.verify_msg}}}")
    return "\n".join(lines)

def list_topics():
    return list(TOPICS.keys())

def list_kinds():
    kinds = []
    for topic in TOPICS:
        kinds.extend(TOPICS[topic])
    return kinds

def help(): #Had to use AI to make this as it kept printing horribly when I tried to do it
    print("""
ml.generate(kind, difficulty, with_ics, seed) -> Problem
    .kind           - ODE type e.g. "ode1_linear"
    .prompt         - plain text question
    .prompt_latex   - LaTeX formatted question

ml.solve(problem) -> Solution
    .final_answer   - answer as a string
    .answer_expr    - answer as a SymPy expression
    .steps          - list of Steps (each has .title, .math, .explanation)
    .verified       - True/False
    .verify_msg     - "passed" or reason it failed
    .warnings       - list of warnings (e.g. domain issues)
    .kind           - ODE type that was solved
    .show_steps()   - display steps with rendered LaTeX (Jupyter)

ml.check(problem, result) -> {"verified": True/False, "message": "..."}
ml.export_latex(result, problem=None) -> LaTeX string of full solution

ml.problem_sheet(n, kind, difficulty, export, saveto, seed) -> list of questions
ml.list_kinds()   -> all 7 ODE type names (+ "mixed" in problem_sheet)
ml.list_topics()  -> topic groupings
ml.launch_web()   -> start browser UI
ml.help()         -> this message

Jupyter: Problems and Solutions auto-render LaTeX when displayed in cells.
         Use result.show_steps() to see all solution steps rendered.
""")

def problem_sheet(n=5, kind="ode1_linear", difficulty="mixed", export=False, saveto=None, seed=None):
    import random
    if seed is not None:
        random.seed(seed)
    questions = []
    for i in range(n):
        diff = random.choice(["easy", "medium", "hard"]) if difficulty == "mixed" else difficulty
        k = random.choice([x for x in list_kinds() if x != "ode_sys2_linear"]) if kind == "mixed" else kind
        prob = generate(k, difficulty=diff)
        questions.append({"number": i + 1, "difficulty": diff,
                          "problem": prob, "solution": solve(prob)})

    if export:
        import pathlib
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        from matplotlib.mathtext import MathTextParser
        parser = MathTextParser("path")

        def safe_math(s):
            try:
                parser.parse(f"${s}$")
                return f"${s}$", 15
            except ValueError:
                return str(s), 10

        def write_page(pdf, title, rows):
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(0.5, 0.95, title, ha="center", fontsize=16, weight="bold")
            y = 0.90
            for label, math in rows:
                if y < 0.08:
                    pdf.savefig(fig); plt.close(fig)
                    fig = plt.figure(figsize=(8.27, 11.69)); y = 0.95
                fig.text(0.06, y, label, fontsize=11, weight="bold")
                text, size = safe_math(math)
                fig.text(0.08, y - 0.03, text, fontsize=size)
                y -= 0.10
            pdf.savefig(fig); plt.close(fig)

        if saveto:
            folder = pathlib.Path(saveto)
        else:
            folder = pathlib.Path(__file__).parent
        folder.mkdir(parents=True, exist_ok=True)
        path = str(folder / f"problem_sheet_{kind}.pdf")

        with PdfPages(path) as pdf:
            write_page(pdf, f"Problem Sheet — {kind.replace('_', ' ')}",
                [(f"Q{q['number']}.", q["problem"].prompt_latex) for q in questions])
            write_page(pdf, "Answers",
                [(f"Q{q['number']}:", to_latex(q["solution"].answer_expr)) for q in questions])

        print(f"Saved to {path}")

    return questions

def launch_web(host="127.0.0.1", port=8080, open_browser=True):
    import threading, webbrowser, logging
    from frontend.app import app
    log = logging.getLogger("werkzeug")
    log.setLevel(logging.ERROR)
    if open_browser:
        threading.Timer(1.0, lambda: webbrowser.open(f"http://{host}:{port}")).start()
    try:
        get_ipython()
        t = threading.Thread(target=lambda: app.run(host=host, port=port, debug=False, use_reloader=False), daemon=True)
        t.start()
        print(f"Server running http://{host}:{port}")
    except NameError:
        app.run(host=host, port=port, debug=False, use_reloader=False)
