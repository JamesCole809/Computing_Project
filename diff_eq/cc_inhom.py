import sympy as sp
import random
from project import (
    Problem, Step, Solution, TutorError,
    extract_symbols, to_latex, apply_ics, add_domain_and_verify,
    make_ic, make_problem_data,
)

#This one too way too long

def format_dict(d):
    parts = sorted(d.items(), key=lambda kv: str(kv[0]))
    return ", \\; ".join(f"{sp.latex(k)} = {sp.latex(v)}" for k, v in parts)


def get_complementary(a_coeff, b_coeff, c_coeff, x):
    mu = sp.Symbol("mu")
    aux_eq = a_coeff * mu**2 + b_coeff * mu + c_coeff
    roots = sp.solve(aux_eq, mu)
    disc = sp.simplify(b_coeff**2 - 4 * a_coeff * c_coeff)

    A = sp.Symbol("A")
    B = sp.Symbol("B")

    disc_pos = disc.is_positive
    disc_zero = disc.is_zero
    disc_neg = disc.is_negative

    #idk if this is the right way tbh, sympy kept returning None on the sign so just did it numerically. seems to work
    if disc_pos is None or disc_zero is None or disc_neg is None:
        try:
            disc_num = float(sp.N(disc))
            disc_pos = disc_num > 1e-12
            disc_zero = abs(disc_num) <= 1e-12
            disc_neg = disc_num < -1e-12
        except (TypeError, ValueError):
            raise TutorError("Can't evaluate discriminant")

    if disc_pos:
        if len(roots) != 2:
            raise TutorError("Expected two distinct roots but got unexpected root")
        # sort just so output is consistent, not really needed
        try:
            mu1, mu2 = sorted(roots, key=lambda r: float(sp.N(r)))
        except (TypeError, ValueError):
            mu1, mu2 = roots[0], roots[1]
        yc = A * sp.exp(mu1 * x) + B * sp.exp(mu2 * x)
        case = "two_real"
    elif disc_zero:
        if len(roots) < 1:
            raise TutorError("Expected repeated root but got no roots")
        mu1 = roots[0]
        yc = A * sp.exp(mu1 * x) + B * x * sp.exp(mu1 * x)
        case = "repeated"
    elif disc_neg:
        p = -b_coeff / (2 * a_coeff)
        q = sp.Abs(sp.sqrt(4 * a_coeff * c_coeff - b_coeff**2) / (2 * a_coeff))
        yc = sp.exp(p * x) * (A * sp.cos(q * x) + B * sp.sin(q * x))
        case = "complex"
    else:
        raise TutorError("Can't classify discriminant")

    return yc, roots, case


def term_trial(rest, x, idx):
    if rest == 1:
        c = sp.Symbol(f"K{idx}")
        return c, [c], idx + 1

    if rest.is_Pow and rest.base == x:
        terms, coeffs = 0, []
        for i in range(int(rest.exp) + 1):
            c = sp.Symbol(f"K{idx}")
            terms += c * x**i
            coeffs.append(c)
            idx += 1
        return terms, coeffs, idx

    if x in rest.free_symbols and rest.func != sp.exp and not rest.has(sp.sin, sp.cos, sp.exp):
        deg = sp.Poly(rest, x).degree()
        terms, coeffs = 0, []
        for i in range(deg + 1):
            c = sp.Symbol(f"K{idx}")
            terms += c * x**i
            coeffs.append(c)
            idx += 1
        return terms, coeffs, idx

    if rest.has(sp.exp) and not rest.has(sp.sin, sp.cos):
        exp_term = list(rest.atoms(sp.exp))[0]
        c = sp.Symbol(f"K{idx}")
        return c * exp_term, [c], idx + 1

    if rest.has(sp.sin) or rest.has(sp.cos):
        terms, coeffs = 0, []
        for arg in {a.args[0] for a in rest.atoms(sp.sin, sp.cos)}:
            c1 = sp.Symbol(f"K{idx}")
            c2 = sp.Symbol(f"K{idx+1}")
            terms += c1 * sp.cos(arg) + c2 * sp.sin(arg)
            coeffs.extend([c1, c2])
            idx += 2
        return terms, coeffs, idx

    c = sp.Symbol(f"K{idx}")
    return c * rest, [c], idx + 1


def check_overlap(sub_trial, yc, x):
    #multiply by x if test clashes with yc, 3 loops should be enough hopefully

    for _ in range(3):
        overlap = False
        for t in sp.Add.make_args(yc):
            _, yc_core = t.as_independent(x, as_Add=False)
            yc_core = sp.expand(yc_core)
            for tr in sp.Add.make_args(sub_trial):
                _, tr_core = tr.as_independent(x, as_Add=False)
                tr_core = sp.expand(tr_core)
                if sp.simplify(yc_core - tr_core) == 0:
                    overlap = True
                    break
            if overlap:
                break
        if overlap:
            sub_trial = sub_trial * x
        else:
            break
    return sub_trial


def build_trial(f_x, x, yc):
    undetermined = []
    idx = 0
    total = 0
    #avoid duplicate test terms as I was getting double unknowns without these
    seen_trig_args = set()
    seen_exp_bases = set()

    for term in sp.Add.make_args(f_x):
        _, rest = term.as_coeff_Mul()

        trig_atoms = rest.atoms(sp.sin, sp.cos)
        if trig_atoms:
            arg = next(iter({a.args[0] for a in trig_atoms}))
            if arg in seen_trig_args:
                continue
            seen_trig_args.add(arg)

        exp_atoms = rest.atoms(sp.exp)
        if exp_atoms and not trig_atoms and len(exp_atoms) == 1:
            exp_term = next(iter(exp_atoms))
            coeff_part, dep_part = rest.as_independent(x, as_Add=False)
            if dep_part == exp_term:
                if exp_term in seen_exp_bases:
                    continue
                seen_exp_bases.add(exp_term)

        sub_trial, coeffs, idx = term_trial(rest, x, idx)
        sub_trial = check_overlap(sub_trial, yc, x)
        total += sub_trial
        undetermined.extend(coeffs)

    if total == 0:
        c = sp.Symbol("K0")
        return c * f_x, [c]

    return total, undetermined


def solve_ode2_cc_inhom(problem, want_steps, want_verify):
    d = problem.data
    x, y = extract_symbols(problem)
    eq = d["equations"]
    steps = []

    dy = y.diff(x)
    d2y = y.diff(x, 2)
    expr = eq.lhs - eq.rhs if isinstance(eq, sp.Eq) else eq
    expr = sp.expand(expr)

    a_coeff = expr.coeff(d2y)
    b_coeff = expr.coeff(dy)
    c_coeff = expr.coeff(y.func(x))
    f_x = -sp.simplify(expr - a_coeff * d2y - b_coeff * dy - c_coeff * y)

    if a_coeff == 0:
        raise TutorError("No y'' term found")

    remainder = sp.simplify(expr - (a_coeff * d2y + b_coeff * dy + c_coeff * y - f_x))
    if remainder != 0 and remainder.equals(0) is False:
        raise TutorError("Couldn't put nto standard inhomogeneous form")

    for name, coeff in [("a", a_coeff), ("b", b_coeff), ("c", c_coeff)]:
        if x in coeff.free_symbols:
            raise TutorError(f"Coefficient {name} depends on x and needs constant coefficients")

    if want_steps:
        steps.append(Step(
            "Identify standard form",
            rf"{to_latex(a_coeff)} \frac{{d^2 y}}{{dx^2}} + \left({to_latex(b_coeff)}\right) \frac{{dy}}{{dx}} + \left({to_latex(c_coeff)}\right) y = {to_latex(f_x)}",
            rf"$a = {to_latex(a_coeff)}, \; b = {to_latex(b_coeff)}, \; c = {to_latex(c_coeff)}, \; f(x) = {to_latex(f_x)}$"
        ))

    A = sp.Symbol("A")
    B = sp.Symbol("B")
    yc, roots, case = get_complementary(a_coeff, b_coeff, c_coeff, x)

    if want_steps:
        steps.append(Step(
            "Find complementary function",
            rf"y_c = {to_latex(yc)}",
            f"Auxiliary roots: ${', \\; '.join(to_latex(r) for r in roots)}$ ({case})"
        ))

    trial, undetermined = build_trial(f_x, x, yc)

    if want_steps:
        steps.append(Step(
            "Choose test particular solution",
            rf"y_p = {to_latex(trial)}",
            "Form chosen from $f(x)$ type with check"
        ))

    lhs_sub = a_coeff * sp.diff(trial, x, 2) + b_coeff * sp.diff(trial, x) + c_coeff * trial
    lhs_sub = sp.expand(lhs_sub)
    rhs_expanded = sp.expand(f_x)

    diff_expr = sp.expand(lhs_sub - rhs_expanded)

    #In lectures we just equate coefficients of like terms but couldn't get
    # that working with sympy, so plugging in numbers for x to make a system
    # of linear equations instead. probably not the proper way but it works
    n_unknowns = len(undetermined)
    coeff_vals = None
    for offset in range(3):
        eval_pts = [sp.Rational(k + offset * 10, 4) for k in range(n_unknowns)]
        eqs = [diff_expr.subs(x, pt) for pt in eval_pts]
        coeff_vals = sp.solve(eqs, undetermined)
        if coeff_vals:
            break
    if not coeff_vals:
        coeff_vals = sp.solve(diff_expr, undetermined)

    def _normalise(cv):
        if isinstance(cv, dict):
            return cv
        if isinstance(cv, list):
            if len(cv) == 1 and isinstance(cv[0], dict):
                return cv[0]
            if len(cv) == 1 and isinstance(cv[0], tuple):
                return dict(zip(undetermined, cv[0]))
            if cv and not isinstance(cv[0], (dict, tuple)):
                return {undetermined[0]: cv[0]}
        return None

    coeff_vals = _normalise(coeff_vals)

    if coeff_vals is None or not all(k in coeff_vals for k in undetermined):
        fallback = sp.solve(diff_expr, undetermined)
        coeff_vals = _normalise(fallback)
        if coeff_vals is None or not all(k in coeff_vals for k in undetermined):
            raise TutorError("Could not determine coefficients")

    yp = trial.subs(coeff_vals)
    yp = sp.simplify(yp)

    if want_steps:
        steps.append(Step(
            "Determine coefficients",
            rf"{format_dict(coeff_vals)}",
            rf"$y_p = {to_latex(yp)}$"
        ))

    general_sol = sp.expand(yc + yp)

    if want_steps:
        steps.append(Step(
            "General solution",
            rf"y = {to_latex(general_sol)}",
            "Complementary + particular"
        ))

    ics = d.get("ics", [])
    sorted_ics = sorted(ics, key=lambda ic: ic["order"]) if len(ics) >= 2 else ics
    final_sol, ic_vals = apply_ics(general_sol, sorted_ics, x, (A, B))
    if ic_vals and len(sorted_ics) >= 2 and want_steps:
        ic0 = sorted_ics[0]
        ic1 = sorted_ics[1]
        steps.append(Step(
            "Apply initial conditions",
            rf"y({ic0['x0']}) = {ic0['value']}, \quad y'({ic1['x0']}) = {ic1['value']}",
            rf"$A = {to_latex(ic_vals[A])}, \; B = {to_latex(ic_vals[B])}$"))
        steps.append(Step("Particular solution",
            rf"y = {to_latex(final_sol)}", ""))

    verified, verify_msg = add_domain_and_verify(
        steps, final_sol, want_steps, want_verify, problem)

    return Solution(
        kind="ode2_cc_inhom", final_answer=str(final_sol),
        answer_expr=final_sol, steps=steps, verified=verified,
        verify_msg=verify_msg, warnings=[])


def gen_ode2_cc_inhom(difficulty, with_ics=True):
    x = sp.Symbol("x")
    y = sp.Function("y")(x)

    if difficulty == "easy":
        mu1 = random.choice([-3, -2, -1, 1, 2, 3])
        mu2 = random.choice([-3, -2, -1, 1, 2, 3])
        while mu2 == mu1:
            mu2 = random.choice([-3, -2, -1, 1, 2, 3])
        b = -(mu1 + mu2)
        c = mu1 * mu2
        f_x = random.choice([random.choice([-3, -2, -1, 1, 2, 3]) * sp.S.One,
                              random.choice([-3, -2, -1, 1, 2, 3]) * x,
                              sp.exp(random.choice([-4, -3, 3, 4]) * x)])

    elif difficulty == "medium":
        b = random.choice([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
        c = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        w = random.choice([1, 2, 3])
        f_x = random.choice([sp.sin(w * x), sp.cos(w * x),
                              random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]) * sp.sin(w * x),
                              sp.exp(random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]) * x)])

    else:
        mu1 = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        mu2 = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        b = -(mu1 + mu2)
        c = mu1 * mu2
        clash_root = random.choice([mu1, mu2])
        f_x = sp.exp(clash_root * x)

    eq = sp.Eq(y.diff(x, 2) + b * y.diff(x) + c * y, f_x)
    ics = []
    if with_ics:
        x0 = 0
        if difficulty == "easy":
            y0 = random.choice([-3, -2, -1, 1, 2, 3])
            y1 = random.choice([-3, -2, -1, 1, 2, 3])
        elif difficulty == "medium":
            y0 = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
            y1 = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        else:
            y0 = random.choice(list(range(-10, 0)) + list(range(1, 11)))
            y1 = random.choice(list(range(-10, 0)) + list(range(1, 11)))
        ics = [make_ic("y", 0, x0, y0), make_ic("y", 1, x0, y1)]

    data = make_problem_data("x", ["y"], eq, ics)
    prompt = f"Solve y'' + ({b})*y' + ({c})*y = {f_x}"
    if ics:
        prompt += f", y({x0}) = {y0}, y'({x0}) = {y1}"

    return Problem("ode2_cc_inhom", prompt, data,
                   metadata={"difficulty": difficulty})
