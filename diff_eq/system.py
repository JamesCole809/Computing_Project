import sympy as sp
import random
from project import (
    Problem, Step, Solution, TutorError,
    extract_symbols, to_latex, apply_ics_system, verify_system,
    make_ic, make_problem_data,
)


def solve_ode_sys2_linear(problem, want_steps, want_verify):
    """Solve a coupled first order linear system by eliminating one variable into a second order ODE"""
    d = problem.data
    x, (u, v) = extract_symbols(problem)
    funcs = d["functions"]
    eqs = d["equations"]

    steps = []

    eq1 = eqs[0]
    eq2 = eqs[1]

    if want_steps:
        steps.append(Step("System of equations",
            str(rf"\frac{{d{funcs[0]}}}{{dx}} = {to_latex(eq1.rhs)}, \quad \frac{{d{funcs[1]}}}{{dx}} = {to_latex(eq2.rhs)}"),
            str("Two coupled FO linear ODEs")))

    rhs1 = eq1.rhs
    u_double_prime = sp.diff(rhs1, x)
    u_double_prime = u_double_prime.subs(v.diff(x), eq2.rhs)
    u_double_prime = u_double_prime.subs(u.diff(x), eq1.rhs)

    if want_steps:
        steps.append(Step(f"Differentiate first equation",
            str(rf"\frac{{d^2 {funcs[0]}}}{{dx^2}} = {to_latex(u_double_prime)}"),
            str(rf"Substituted ${funcs[1]}'$ from the second equation")))

    v_expr_list = sp.solve(eq1.rhs - u.diff(x), v)
    if not v_expr_list:
        raise TutorError("Cannot isolate v from the first equation")
    v_expr = v_expr_list[0]

    u_double_prime = sp.simplify(u_double_prime.subs(v, v_expr))

    if want_steps:
        steps.append(Step(f"Eliminate {funcs[1]}",
            str(rf"{funcs[1]} = {to_latex(v_expr)}"),
            str(rf"Substituted into ${funcs[0]}''$ to get a single second-order ODE")))

    ode2_expr = u.diff(x, 2) - u_double_prime
    ode2_eq = sp.Eq(ode2_expr, 0)

    if want_steps:
        steps.append(Step(f"Second-order ODE in {funcs[0]}",
            str(to_latex(ode2_eq)),
            str("Now solve this using the SO methods")))

    ode2_data = make_problem_data(d["x_symbol"], [funcs[0]], ode2_eq)

    dy = u.diff(x)
    d2y = u.diff(x, 2)
    a_c = ode2_expr.coeff(d2y)
    b_c = ode2_expr.coeff(dy)
    c_c = ode2_expr.coeff(u)
    f_x = -sp.simplify(ode2_expr - a_c * d2y - b_c * dy - c_c * u)

    kind_2 = "ode2_cc_hom" if f_x == 0 else "ode2_cc_inhom"
    import project
    ode2_report = project.solve(Problem(kind_2, "", ode2_data),
                        want_steps=False, want_verify=False)
    u_sol = ode2_report.answer_expr

    if want_steps:
        steps.append(Step(f"Solve for {funcs[0]}",
            str(rf"{funcs[0]} = {to_latex(u_sol)}"),
            str("General solution of the second-order ODE")))

    v_sol = sp.simplify(v_expr.subs(u, u_sol).subs(u.diff(x), sp.diff(u_sol, x)))

    if want_steps:
        steps.append(Step(f"Recover {funcs[1]}",
            str(rf"{funcs[1]} = {to_latex(v_sol)}"),
            str(rf"From ${funcs[1]} = {to_latex(v_expr)}$")))

    A = sp.Symbol("A")
    B = sp.Symbol("B")
    final_u, final_v, ic_vals = apply_ics_system(
        u_sol, v_sol, d.get("ics", []), x, funcs, (A, B))
    if ic_vals and want_steps:
        steps.append(Step("Apply initial conditions",
            str(rf"{funcs[0]}(0) = \ldots, \; {funcs[1]}(0) = \ldots"),
            str(rf"$A = {to_latex(ic_vals[A])}, \; B = {to_latex(ic_vals[B])}$")))
        steps.append(Step("Particular solution",
            str(rf"{funcs[0]} = {to_latex(final_u)}, \quad {funcs[1]} = {to_latex(final_v)}"),
            str("")))

    verified = None
    verify_msg = ""
    if want_verify:
        verify_problem = Problem("ode_sys2_linear", "", {
            "x_symbol": d["x_symbol"], "functions": funcs,
            "equations": eqs, "ics": d.get("ics", [])
        })
        verified, verify_msg = verify_system(verify_problem, final_u, final_v)

    answer = f"{funcs[0]} = {final_u}, {funcs[1]} = {final_v}"

    return Solution(
        kind="ode_sys2_linear", final_answer=answer,
        answer_expr={"u": final_u, "v": final_v}, steps=steps,
        verified=verified, verify_msg=verify_msg,
        warnings=[])


def gen_ode_sys2_linear(difficulty, with_ics=True):
    """Make a random coupled first order linear system"""
    x = sp.Symbol("x")
    u = sp.Function("u")(x)
    v = sp.Function("v")(x)

    if difficulty == "easy":
        a = random.choice([-3, -2, -1, 1, 2, 3])
        b = random.choice([-3, -2, -1, 1, 2, 3])
        c = random.choice([-3, -2, -1, 1, 2, 3])
        d = random.choice([-3, -2, -1, 1, 2, 3])
        while b == 0:
            b = random.choice([-3, -2, -1, 1, 2, 3])
        f_x = 0
        g_x = 0

    elif difficulty == "medium":
        a = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        b = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        c = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        d = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        while b == 0:
            b = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        f_x = 0
        g_x = 0

    else:
        a = random.choice(list(range(-10, 0)) + list(range(1, 11)))
        b = random.choice(list(range(-10, 0)) + list(range(1, 11)))
        c = random.choice(list(range(-10, 0)) + list(range(1, 11)))
        d = random.choice(list(range(-10, 0)) + list(range(1, 11)))
        while b == 0:
            b = random.choice(list(range(-10, 0)) + list(range(1, 11)))
        forcing = [sp.exp(x), random.choice(list(range(-10, 0)) + list(range(1, 11))) * sp.exp(x),
                    x, random.choice(list(range(-10, 0)) + list(range(1, 11))) * x]
        f_x = random.choice(forcing)
        g_x = random.choice([0] + forcing)

    rhs1 = a * u + b * v
    rhs2 = c * u + d * v
    if f_x != 0:
        rhs1 = rhs1 + f_x
    if g_x != 0:
        rhs2 = rhs2 + g_x

    eq1 = sp.Eq(u.diff(x), rhs1)
    eq2 = sp.Eq(v.diff(x), rhs2)

    ics = []
    if with_ics:
        if difficulty == "easy":
            u0 = random.choice([-3, -2, -1, 1, 2, 3])
            v0 = random.choice([-3, -2, -1, 1, 2, 3])
        elif difficulty == "medium":
            u0 = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
            v0 = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        else:
            u0 = random.choice(list(range(-10, 0)) + list(range(1, 11)))
            v0 = random.choice(list(range(-10, 0)) + list(range(1, 11)))
        ics = [make_ic("u", 0, 0, u0), make_ic("v", 0, 0, v0)]

    data = make_problem_data("x", ["u", "v"], [eq1, eq2], ics)
    prompt = f"Solve u' = {rhs1}, v' = {rhs2}"
    if ics:
        prompt += f", u(0) = {u0}, v(0) = {v0}"

    return Problem("ode_sys2_linear", prompt, data,
                   metadata={"difficulty": difficulty})
