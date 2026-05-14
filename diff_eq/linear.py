#imports
import sympy as sp
import random
from project import (
    Problem, Step, Solution, TutorError,
    extract_symbols, to_latex, apply_ics, add_domain_and_verify,
    make_ic, make_problem_data,
)

#function that solves ode
def solve_ode1_linear(problem, want_steps, want_verify):
    """First order linear ODE using the integrating factor method"""
    #shortcuts
    d = problem.data
    x, y = extract_symbols(problem)
    eq = d["equations"]

    #empty lists
    steps = []
    warnings = []

    #shorthand
    dy = y.diff(x)
    #moves everything for = 0
    expr = eq.lhs - eq.rhs if isinstance(eq, sp.Eq) else eq
    #reads coeff infront of dy/dx
    coeff_dy = expr.coeff(dy)

    #guard
    if coeff_dy == 0:
        raise TutorError("No y' term found")
    #if not 1 divide through for standard form
    if coeff_dy != 1:
        warnings.append(f"Divided by {coeff_dy} exclude zeros from domain")
        expr = sp.simplify(expr / coeff_dy)

    #coeff of y is p(x)
    p = expr.coeff(y.func(x))
    #strips away to get q(x)
    q = -sp.simplify(expr - dy - p * y)

    #steps of course
    if want_steps:
        steps.append(Step("Rearrange for standard form",
            str(rf"\frac{{dy}}{{dx}} + \left({to_latex(p)}\right) y = {to_latex(q)}"),
            str(rf"$p(x) = {to_latex(p)}, \; q(x) = {to_latex(q)}$")))

    #computes integral of p(x)
    int_p = sp.integrate(p, x)
    #integrating factor
    mu = sp.simplify(sp.exp(int_p))

    #steps :)
    if want_steps:
        steps.append(Step("Compute integrating factor",
            str(rf"I(x) = e^{{\int {to_latex(p)} \, dx}} = {to_latex(mu)}"),
            str("Multiply both sides by this to make LHS a product rule derivative")))
        steps.append(Step("Multiply entire equation by I",
            str(rf"{to_latex(mu)} \cdot \frac{{dy}}{{dx}} + {to_latex(mu)} \left({to_latex(p)}\right) y = {to_latex(mu)} \left({to_latex(q)}\right)"),
            str("LHS is now d/dx[I(x)*y]")))
        steps.append(Step("Rewrite LHS as a single derivative",
            str(rf"\frac{{d}}{{dx}}\!\left[{to_latex(mu)} \cdot y\right] = {to_latex(mu)} \left({to_latex(q)}\right)"),
            str("This is product rule in reverse")))

    #integrates right hand side
    rhs_integral = sp.integrate(mu * q, x)
    #symbolic integration constant
    C1 = sp.Symbol("C1")

    #steps again
    if want_steps:
        steps.append(Step("Integrate both sides",
            str(rf"{to_latex(mu)} \cdot y = {to_latex(rhs_integral)} + C_1"),
            str("Add the integration constant")))

    #gives general solution
    general_sol = sp.simplify((rhs_integral + C1) / mu)

    #yup steps
    if want_steps:
        steps.append(Step("Solve for y",
            str(rf"y = {to_latex(general_sol)}"), str("This is the general solution")))

    #final solution
    final_sol, ic_vals = apply_ics(general_sol, d.get("ics", []), x, (C1,))
    if ic_vals and want_steps:
        steps.append(Step("Apply initial condition",
            str(rf"y({d['ics'][0]['x0']}) = {d['ics'][0]['value']} \Rightarrow C_1 = {to_latex(ic_vals[C1])}"),
            str(rf"Substituting into the general solution")))
        steps.append(Step("Particular solution",
            str(rf"y = {to_latex(final_sol)}"),
            str("")))

    #if statement for warnings
    if warnings and d.get("ics") and want_steps:
        steps.append(Step("Domain check",
            str(rf"\text{{Excluded zeros from division; IC at }} x = {d['ics'][0]['x0']}"),
            str("Pick the interval containing the IC point")))

    #verify and plug in answer to original ODE
    verified, verify_msg = add_domain_and_verify(
        steps, final_sol, want_steps, want_verify, problem)

    #returns the solution
    return Solution(
        kind="ode1_linear", final_answer=str(final_sol),
        answer_expr=final_sol, steps=steps, verified=verified,
        verify_msg=verify_msg, warnings=warnings)

#function to generate random problems
def gen_ode1_linear(difficulty, with_ics=True):
    """Make a random first order linear ODE"""
    x = sp.Symbol("x")
    y = sp.Function("y")(x)

    #difficulty
    if difficulty == "easy":
        p = random.choice([1, 2, -1, 3, -2])
        q = random.choice([x, 1, sp.exp(x), sp.sin(x), x**2])
        x0 = random.choice([0, 1])
        y0 = random.choice([-3, -2, -1, 1, 2, 3])

    elif difficulty == "medium":
        k = random.choice([1, 2, 3, -1, -2])
        p = k / x
        q = random.choice([x, x**2, sp.exp(x), sp.sin(x),
                         random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]) * x])
        x0 = random.choice([1, 2])
        y0 = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])

    else:
        p = random.choice([1/x, 2/x, 3/x, x, sp.cos(x)])
        q = random.choice([x**2, sp.exp(x), x * sp.exp(x),
                         sp.sin(x), sp.cos(x)])
        has_x_denom = p in [1/x, 2/x, 3/x]
        x0 = random.choice([1, 2]) if has_x_denom else random.choice([0, 1])
        y0 = random.choice(list(range(-10, 0)) + list(range(1, 11)))

    eq = sp.Eq(y.diff(x) + p * y, q)
    ics = []
    if with_ics:
        ics = [make_ic("y", 0, x0, y0)]

    data = make_problem_data("x", ["y"], eq, ics)
    prompt = f"Solve y' + ({p})*y = {q}"
    if ics:
        prompt += f", y({x0}) = {y0}"

    #returns the problem
    return Problem("ode1_linear", prompt, data,
                   metadata={"difficulty": difficulty})
