#imports
import sympy as sp
import random
from project import (
    Problem, Step, Solution, TutorError,
    extract_symbols, to_latex, add_domain_and_verify,
    make_ic, make_problem_data,
)

def solve_ode1_separable(problem, want_steps, want_verify):
    """Separable ODE: split variables and integrate both sides"""
    d = problem.data
    x, y = extract_symbols(problem)
    eq = d["equations"]

    steps = []

    dy = y.diff(x)
    expr = eq.lhs - eq.rhs if isinstance(eq, sp.Eq) else eq
    g_y = expr.coeff(dy)
    f_x = -sp.simplify(expr - g_y * dy)

    if g_y == 0:
        raise TutorError("No y' term for separable form")

    if want_steps:
        steps.append(Step("Rearrange for separable form",
            str(rf"\left({to_latex(g_y)}\right) \frac{{dy}}{{dx}} = {to_latex(f_x)}"),
            str(rf"$g(y) = {to_latex(g_y)}, \; f(x) = {to_latex(f_x)}$")))

    y_sym = sp.Symbol("y")
    g_y_sym = g_y.subs(y.func(x), y_sym)
    lhs_integral = sp.integrate(g_y_sym, y_sym)
    rhs_integral = sp.integrate(f_x, x)
    C1 = sp.Symbol("C1")

    if want_steps:
        steps.append(Step("Integrate both sides",
            str(rf"\int {to_latex(g_y_sym)} \, dy = \int {to_latex(f_x)} \, dx"),
            str(rf"${to_latex(lhs_integral)} = {to_latex(rhs_integral)} + C_1$")))

    implicit_sol = sp.Eq(lhs_integral, rhs_integral + C1)

    if want_steps:
        steps.append(Step("Implicit solution", str(to_latex(implicit_sol)),
            str("Relationship between x and y")))

    explicit_sols = sp.solve(implicit_sol, y_sym)
    general_sol = explicit_sols[0] if explicit_sols else None

    if general_sol is not None and want_steps:
        steps.append(Step("Solve for y",
            str(rf"y = {to_latex(general_sol)}"), str("General solution")))

    final_sol = general_sol
    final_implicit = implicit_sol
    if d.get("ics"):
        ic = d["ics"][0]
        c_val = sp.solve(implicit_sol.subs({y_sym: ic["value"], x: ic["x0"]}), C1)
        if c_val:
            if general_sol is not None:
                found = False
                for sol in explicit_sols:
                    test = sol.subs(C1, c_val[0]).subs(x, ic["x0"])
                    if sp.simplify(test - ic["value"]) == 0:
                        general_sol = sol
                        found = True
                        break
                if found:
                    final_sol = sp.simplify(general_sol.subs(C1, c_val[0]))
                else:
                    final_sol = None
            final_implicit = implicit_sol.subs(C1, c_val[0])
            if want_steps:
                steps.append(Step("Apply initial condition",
                    str(rf"y({ic['x0']}) = {ic['value']} \Rightarrow C_1 = {to_latex(c_val[0])}"),
                    str("Substituting into general solution")))
                steps.append(Step("Particular solution",
                    str(to_latex(final_implicit)),
                    str("")))

    result = final_sol if final_sol is not None else final_implicit
    verified, verify_msg = add_domain_and_verify(
        steps, result, want_steps, want_verify, problem)

    answer = str(result)

    return Solution(
        kind="ode1_separable", final_answer=answer,
        answer_expr=result,
        steps=steps, verified=verified, verify_msg=verify_msg,
        warnings=[])


def gen_ode1_separable(difficulty, with_ics=True):
    """Make a random separable ODE"""
    x = sp.Symbol("x")
    y = sp.Function("y")(x)

    if difficulty == "easy":
        g = random.choice([1, y, y**2, sp.exp(y)])
        f = random.choice([x, 1, sp.exp(x), sp.cos(x), x**2])
        x0 = random.choice([0, 1])
        y0 = random.choice([-3, -2, -1, 1, 2, 3])

    elif difficulty == "medium":
        g = random.choice([1/(1 + y**2), y/(1 + y**2), sp.cos(y), y])
        f = random.choice([x, sp.exp(x), sp.cos(x), x**2,
                         random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]) * x])
        x0 = random.choice([0, 1])
        y0 = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])

    else:
        g = random.choice([sp.exp(-y), sp.exp(y), 1/y**2])
        f = random.choice([x, x**2, sp.sin(x)])
        x0 = 0
        if g == 1/y**2:
            y0 = random.choice([1, 2, -1, -2])
        else:
            y0 = random.choice([0, 1, -1])

    eq = sp.Eq(g * y.diff(x), f)
    ics = []
    if with_ics:
        ics = [make_ic("y", 0, x0, y0)]

    data = make_problem_data("x", ["y"], eq, ics)
    prompt = f"Solve ({g})*y' = {f}"
    if ics:
        prompt += f", y({x0}) = {y0}"

    return Problem("ode1_separable", prompt, data,
                   metadata={"difficulty": difficulty})
