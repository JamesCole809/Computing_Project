import sympy as sp
import random
from project import (
    Problem, Step, Solution, TutorError,
    extract_symbols, to_latex, apply_ics, add_domain_and_verify,
    make_ic, make_problem_data,
)

#File for constant coefficient homogeneous ODEs

def solve_ode2_cc_hom(problem, want_steps, want_verify):
    """Second order constant coefficient homogeneous ODE via the auxiliary equation"""
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

    if a_coeff == 0:
        raise TutorError("No y'' term found")

    remainder = sp.simplify(expr - a_coeff * d2y - b_coeff * dy - c_coeff * y.func(x))
    if remainder != 0 and remainder.equals(0) is False:
        raise TutorError("Not a CC hom ODE in standard form")

    for name, coeff in [("a", a_coeff), ("b", b_coeff), ("c", c_coeff)]:
        if x in coeff.free_symbols:
            raise TutorError(f"Coefficient {name} depends on x and needs constant coefficients")

    if want_steps:
        steps.append(Step("Identify standard form",
            rf"{to_latex(a_coeff)} \frac{{d^2 y}}{{dx^2}} + \left({to_latex(b_coeff)}\right) \frac{{dy}}{{dx}} + \left({to_latex(c_coeff)}\right) y = 0",
            rf"$a = {to_latex(a_coeff)}, \; b = {to_latex(b_coeff)}, \; c = {to_latex(c_coeff)}$"))

    mu = sp.Symbol("mu")
    aux_eq = a_coeff * mu**2 + b_coeff * mu + c_coeff
    roots = sp.solve(aux_eq, mu)

    if want_steps:
        steps.append(Step("make auxiliary equation",
            rf"{to_latex(aux_eq)} = 0",
            f"Roots: ${', \\; '.join(to_latex(r) for r in roots)}$"))

    disc = sp.simplify(b_coeff**2 - 4 * a_coeff * c_coeff)
    A = sp.Symbol("A")
    B = sp.Symbol("B")

    disc_pos = disc.is_positive
    disc_zero = disc.is_zero
    disc_neg = disc.is_negative
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
        try:
            mu1, mu2 = sorted(roots, key=lambda r: float(sp.N(r)))
        except (TypeError, ValueError):
            mu1, mu2 = roots[0], roots[1]
        general_sol = A * sp.exp(mu1 * x) + B * sp.exp(mu2 * x)
        if want_steps:
            steps.append(Step("Two distinct real roots",
                rf"\mu_1 = {to_latex(mu1)}, \quad \mu_2 = {to_latex(mu2)}",
                rf"$y = A e^{{{to_latex(mu1)} x}} + B e^{{{to_latex(mu2)} x}}$"))
    elif disc_zero:
        if len(roots) < 1:
            raise TutorError("Expected repeated root but got no roots")
        mu1 = roots[0]
        general_sol = A * sp.exp(mu1 * x) + B * x * sp.exp(mu1 * x)
        if want_steps:
            steps.append(Step("Repeated root",
                rf"\mu = {to_latex(mu1)}",
                rf"$y = A e^{{{to_latex(mu1)} x}} + B x e^{{{to_latex(mu1)} x}}$"))
    elif disc_neg:
        p = -b_coeff / (2 * a_coeff)
        q_raw = sp.sqrt(4 * a_coeff * c_coeff - b_coeff**2) / (2 * a_coeff)
        q_disp = sp.Abs(q_raw)
        general_sol = sp.exp(p * x) * (A * sp.cos(q_raw * x) + B * sp.sin(q_raw * x))
        if want_steps:
            steps.append(Step("Complex roots",
                rf"p \pm iq = {to_latex(p)} \pm {to_latex(q_disp)} i",
                rf"$y = e^{{{to_latex(p)} x}} \left(A \cos({to_latex(q_disp)} x) + B \sin({to_latex(q_disp)} x)\right)$"))
    else:
        raise TutorError("Can't classify discriminant")

    if want_steps:
        steps.append(Step("General solution",
            rf"y = {to_latex(general_sol)}",
            "General solution with two constants"))

    ics = d.get("ics", [])
    sorted_ics = sorted(ics, key=lambda ic: ic["order"]) if len(ics) >= 2 else ics
    final_sol, ic_vals = apply_ics(general_sol, sorted_ics, x, (A, B))
    if ic_vals and len(sorted_ics) >= 2 and want_steps:
        ic0 = sorted_ics[0]
        ic1 = sorted_ics[1]
        steps.append(Step("Apply initial conditions",
            rf"y({ic0['x0']}) = {ic0['value']}, \quad y'({ic1['x0']}) = {ic1['value']}",
            rf"$A = {to_latex(ic_vals[A])}, \; B = {to_latex(ic_vals[B])}$"))
        steps.append(Step("Particular solution",
            rf"y = {to_latex(final_sol)}", ""))

    verified, verify_msg = add_domain_and_verify(
        steps, final_sol, want_steps, want_verify, problem)

    return Solution(
        kind="ode2_cc_hom", final_answer=str(final_sol),
        answer_expr=final_sol, steps=steps, verified=verified,
        verify_msg=verify_msg, warnings=[])


def gen_ode2_cc_hom(difficulty, with_ics=True):
    """Make a random second order CC homogeneous ODE"""
    x = sp.Symbol("x")
    y = sp.Function("y")(x)

    if difficulty == "easy":
        mu1 = random.choice([-4, -3, -2, -1, 1, 2, 3, 4])
        mu2 = random.choice([-4, -3, -2, -1, 1, 2, 3, 4])
        while mu2 == mu1:
            mu2 = random.choice([-4, -3, -2, -1, 1, 2, 3, 4])
        b = -(mu1 + mu2)
        c = mu1 * mu2

    elif difficulty == "medium":
        case = random.choice(["complex", "repeated"])
        if case == "complex":
            p = random.choice([-3, -2, -1, 0, 1, 2, 3])
            q = random.choice([1, 2, 3])
            b = -2 * p
            c = p**2 + q**2
        else:
            mu1 = random.choice([-4, -3, -2, -1, 0, 1, 2, 3, 4])
            b = -2 * mu1
            c = mu1**2

    else:
        case = random.choice(["complex", "repeated"])
        if case == "complex":
            p = random.choice(list(range(-5, 0)) + list(range(1, 6)))
            q = random.choice([1, 2, 3, 4, 5])
            b = -2 * p
            c = p**2 + q**2
        else:
            mu1 = random.choice(list(range(-5, 0)) + list(range(1, 6)))
            b = -2 * mu1
            c = mu1**2

    eq = sp.Eq(y.diff(x, 2) + b * y.diff(x) + c * y, 0)
    ics = []
    if with_ics:
        x0 = 0
        if difficulty == "easy":
            pool = [-3, -2, -1, 1, 2, 3]
        elif difficulty == "medium":
            pool = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
        else:
            pool = list(range(-10, 0)) + list(range(1, 11))
        y0 = random.choice(pool)
        y1 = random.choice(pool)
        ics = [make_ic("y", 0, x0, y0), make_ic("y", 1, x0, y1)]

    data = make_problem_data("x", ["y"], eq, ics)
    prompt = f"Solve y'' + ({b})*y' + ({c})*y = 0"
    if ics:
        prompt += f", y({x0}) = {y0}, y'({x0}) = {y1}"

    return Problem("ode2_cc_hom", prompt, data,
                   metadata={"difficulty": difficulty})
