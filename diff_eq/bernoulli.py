# Imports everything
import sympy as sp
import random
from project import (
    Problem, Step, Solution, TutorError,
    extract_symbols, to_latex, apply_ics, add_domain_and_verify,
    make_ic, make_problem_data,
)


#Function that solves the ODE (first order bernoulli)
def solve_ode1_bernoulli(problem, want_steps, want_verify):
    #shortcut to problem data
    d = problem.data
    x, y = extract_symbols(problem)

    #the actual ODE
    eq = d["equations"]
    #bernoulli exponent from RHS of ODE
    n = d.get("n")

    #empty lists
    steps = []
    warnings = []
    domain_excludes = []

    #shorthand
    dy = y.diff(x)
    #moves everything to one side so it = 0
    expr = eq.lhs - eq.rhs if isinstance(eq, sp.Eq) else eq

    #reads coeff
    coeff_dy = expr.coeff(dy)
    #if statement for coeff != 1 so we can divide through for standard form
    if coeff_dy != 1:
        roots = sp.solve(coeff_dy, x)
        if isinstance(roots, list) and roots:
            warnings.append(f"Divided by {coeff_dy}, exclude x = {roots} from domain")
            domain_excludes.extend(roots)
        expr = sp.simplify(expr / coeff_dy)

    #checks if bernoulli
    if n is None:
        raise TutorError("No 'n' in problem data")

    #symbolics
    yn = y ** n
    b = -expr.coeff(yn)
    remainder = expr - dy + b * yn
    a = sp.simplify(remainder.coeff(y.func(x)))
    b = sp.simplify(b)

    #if statement for if user wants steps
    if want_steps:
        steps.append(Step("Put into correct bernoulli form",
            rf"\frac{{dy}}{{dx}} + \left({to_latex(a)}\right) y = \left({to_latex(b)}\right) y^{{{n}}}",
            rf"$a(x) = {to_latex(a)}, \; b(x) = {to_latex(b)}, \; n = {n}$"))

    #exponent for LHS
    m = 1 - n
    #computes integral of a(x)
    int_a = sp.integrate(a, x)
    #integrating factor
    exp_part = sp.simplify(sp.exp(m * int_a))

    #steps once again
    if want_steps:
        steps.append(Step(r"Form $z = y^{1-n} \cdot e^{(1-n)\int a\,dx}$",
            rf"z = y^{{{m}}} \cdot {to_latex(exp_part)}", rf"with $1-n = {m}$"))

    rhs_part = sp.simplify(m * exp_part * b)

    if want_steps:
        steps.append(Step("Differentiate z and substitute the ODE",
            rf"\frac{{dz}}{{dx}} = {to_latex(rhs_part)}",
            "The y terms cancel after substituting"))

    #integration constant symbol
    C1 = sp.Symbol("C1")
    #integrate for z
    z_expr = sp.integrate(rhs_part, x)

    #steps...
    if want_steps:
        steps.append(Step("Integrate to find z",
            rf"z = {to_latex(z_expr)} + C_1",
            "Add constant"))

    #gives y^(1-n)
    y_m_expr = sp.simplify((z_expr + C1) / exp_part)
    general_sol = sp.simplify(y_m_expr ** sp.Rational(1, m))

    #guess what this does
    if want_steps:
        steps.append(Step("Substitute back for y",
            rf"y = {to_latex(general_sol)}", "General solution"))


    #plugs IC into GC and solves for c1
    final_sol, ic_vals = apply_ics(general_sol, d.get("ics", []), x, (C1,))
    if ic_vals and want_steps:
        steps.append(Step("Apply initial condition",
            rf"y({d['ics'][0]['x0']}) = {d['ics'][0]['value']} \Rightarrow C_1 = {to_latex(ic_vals[C1])}",
            "Substituting into the general solution"))
        steps.append(Step("Particular solution",
            rf"y = {to_latex(final_sol)}", ""))

    #warnings if statement
    if sp.simplify(a * x).is_number:
        if 0 not in domain_excludes:
            domain_excludes.append(0)
        ics = d.get("ics", [])
        if ics and ics[0]["x0"] > 0:
            warnings.append("a(x) has 1/x and x0 > 0, so domain restricted to x > 0")
        else:
            warnings.append("a(x) has 1/x, so x = 0 removed from domain")

    d["domain_excludes"] = domain_excludes

    verified, verify_msg = add_domain_and_verify(
        steps, final_sol, want_steps, want_verify, problem)

    return Solution(
        kind="ode1_bernoulli", final_answer=str(final_sol),
        answer_expr=final_sol, steps=steps, verified=verified,
        verify_msg=verify_msg, warnings=warnings)


#deals with sympy functions we dont want in problems we generate
SPECIAL_FUNCS = (sp.Ei, sp.erf, sp.erfi, sp.Si, sp.Ci, sp.Shi, sp.Chi,
                  sp.gamma, sp.uppergamma, sp.meijerg, sp.hyper, sp.expint)

#returns false if the expression contains unevaluated integral or any of the above functions otherwise true
def is_elementary(expr):
    if expr.has(sp.Integral):
        return False
    for f in SPECIAL_FUNCS:
        if expr.has(f):
            return False
    return True

#random problem gen, pretty self explanitory code
def gen_ode1_bernoulli(difficulty, with_ics=True):
    x = sp.Symbol("x")
    y = sp.Function("y")(x)

    if difficulty == "easy":
        n = random.choice([2, 3])
        a = random.choice([1, -1, 2, -2])
        b = random.choice([1, -1, x])

    elif difficulty == "medium":
        family = random.choice(["const", "const", "invx"])
        n = random.choice([2, 3, -1])
        if family == "invx":
            a = random.choice([sp.Rational(1, 1)/x, sp.Rational(2, 1)/x])
            b = random.choice([1, -1, 2, x])
        else:
            a = random.choice([1, -1, 2, 3])
            b = random.choice([1, -1, 2, x])

    else:
        family = random.choice(["const", "const", "invx"])
        n = random.choice([2, 3, -1, 4])
        if family == "invx":
            a = random.choice([sp.Rational(1, 1)/x, sp.Rational(2, 1)/x,
                               sp.Rational(3, 1)/x])
            b = random.choice([1, -1, x, 2*x, x**2])
        else:
            a = random.choice([1, -1, 2, -2, 3, -3])
            b = random.choice([1, -1, x, 2*x, sp.exp(x)])

    if difficulty == "easy":
        retry_a = [1, -1, 2, -2]
        retry_b = [1, -1, x]
        retry_n = [2, 3]
    elif difficulty == "medium":
        retry_a = [1, -1, 2, 3]
        retry_b = [1, -1, 2, x]
        retry_n = [2, 3, -1]
    else:
        retry_a = [1, -1, 2, -2, 3, -3]
        retry_b = [1, -1, x, 2*x]
        retry_n = [2, 3, -1, 4]

    for _attempt in range(5):
        m = 1 - n
        int_a = sp.integrate(a, x)
        exp_part = sp.exp(m * int_a)
        rhs = sp.simplify(m * exp_part * b)
        z_integral = sp.integrate(rhs, x)
        if is_elementary(z_integral):
            break
        a = random.choice(retry_a)
        b = random.choice(retry_b)
        n = random.choice(retry_n)
    else:
        a, b, n = 1, 1, 2
        m = 1 - n

    eq = sp.Eq(y.diff(x) + a * y, b * y**n)
    ics = []
    if with_ics:
        x0 = random.choice([1, 2, 3])
        y0 = random.choice([1, 2, sp.Rational(1, 2)])
        ics = [make_ic("y", 0, x0, y0)]

    data = make_problem_data("x", ["y"], eq, ics)
    data["n"] = n
    prompt = f"Solve y' + ({a})*y = ({b})*y^{n}"
    if ics:
        prompt += f", y({x0}) = {y0}"

    #returns the problem
    return Problem("ode1_bernoulli", prompt, data,
                   metadata={"difficulty": difficulty})
