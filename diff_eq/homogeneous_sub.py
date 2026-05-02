import sympy as sp
import random
from project import (
    Problem, Step, Solution, TutorError,
    extract_symbols, to_latex, add_domain_and_verify,
    make_ic, make_problem_data,
)


def solve_ode1_homogeneous_sub(problem, want_steps, want_verify):
    d = problem.data
    x, y = extract_symbols(problem)
    eq = d["equations"]

    steps = []
    warnings = []
    domain_excludes = []

    dy = y.diff(x)
    expr = eq.lhs - eq.rhs if isinstance(eq, sp.Eq) else eq

    coeff_dy = expr.coeff(dy)
    if coeff_dy == 0:
        raise TutorError("No y' term found")

    if coeff_dy != 1:
        raw_coeff_roots = sp.solve(coeff_dy, x)
        try:
            coeff_roots = list(raw_coeff_roots)
        except TypeError:
            coeff_roots = []
            warnings.append(f"Exclude where {coeff_dy} = 0")
        if coeff_roots:
            domain_excludes.extend(coeff_roots)
            warnings.append(f"Divided by {coeff_dy} so exclude x = {coeff_roots} from domain")

    rhs_orig = sp.simplify(-(expr - coeff_dy * dy) / coeff_dy)

    z = sp.Symbol("z")
    f_z = rhs_orig.subs(y.func(x), z * x)
    f_z = sp.simplify(f_z)

    if x in f_z.free_symbols:
        raise TutorError("Not homogeneous because RHS is not a function of y/x only")

    if 0 not in domain_excludes:
        domain_excludes.append(0)
        warnings.append("Substitution z = y/x requires x != 0") #my keyboard doesn't have the not equal to symbol

    if want_steps:
        steps.append(Step(
            "Identify homogeneous form",
            rf"\frac{{dy}}{{dx}} = {to_latex(rhs_orig)}",
            rf"Substituting $y = zx$ gives $f(z) = {to_latex(f_z)}$"))
        steps.append(Step(
            "Let z = y/x",
            r"z = \frac{y}{x}, \quad y = zx",
            "We will rewrite the ODE in terms of z"))
        steps.append(Step(
            "Differentiate z",
            r"\frac{dy}{dx} = x \frac{dz}{dx} + z",
            "By the product rule on y = zx"))

    sep_rhs = sp.simplify(f_z - z)
    sep_rhs_is_zero = (sep_rhs == 0) or (sep_rhs.equals(0) is True)

    singular_z_vals = []
    if not sep_rhs_is_zero:
        raw_roots = sp.solve(sep_rhs, z)
        try:
            singular_z_vals = [v for v in list(raw_roots)
                               if x not in v.free_symbols]
        except (TypeError, AttributeError):
            warnings.append("f(z) - z has roots so constant solutions may exist")
        for z0 in singular_z_vals:
            warnings.append(f"Constant solution y = {z0}*x (from f(z) - z = 0 at z = {z0})")
        if singular_z_vals and want_steps:
            sol_strs = ", \\; ".join(rf"y = {to_latex(z0)} x" for z0 in singular_z_vals)
            steps.append(Step(
                "Additional constant solutions",
                sol_strs,
                "From roots of f(z) - z = 0 (lost when separating variables)"))

    if want_steps:
        steps.append(Step(
            "Rewrite ODE in terms of z",
            rf"x \frac{{dz}}{{dx}} + z = {to_latex(f_z)} \Rightarrow x \frac{{dz}}{{dx}} = {to_latex(sep_rhs)}",
            "This is now separable in z and x"
        ))

    C1 = sp.Symbol("C1")
    if sep_rhs_is_zero:
        pass
    else:
        lhs_integrand = 1 / sep_rhs
        lhs_integral = sp.integrate(lhs_integrand, z)

        if want_steps:
            steps.append(Step(
                "Solve ODE",
                rf"\int \frac{{1}}{{{to_latex(sep_rhs)}}} \, dz = \int \frac{{1}}{{x}} \, dx",
                rf"${to_latex(lhs_integral)} = \ln|x| + C_1$"
            ))

    y_sym = sp.Symbol("y")
    if sep_rhs_is_zero:
        implicit_in_y = sp.Eq(y_sym, C1 * x)
    else:
        lhs_in_y = lhs_integral.subs(z, y_sym / x)
        implicit_in_y = sp.Eq(lhs_in_y, sp.log(sp.Abs(x)) + C1)

    if want_steps:
        steps.append(Step(
            "Substitute back z = y/x",
            to_latex(implicit_in_y),
            "Solution in terms of x and y"
        ))

    explicit_sols = sp.solve(implicit_in_y, y_sym)
    general_sol = explicit_sols[0] if explicit_sols else None

    if general_sol is not None and want_steps:
        steps.append(Step(
            "Solve for y explicitly",
            rf"y = {to_latex(general_sol)}",
            "General solution"
        ))

    final_sol = general_sol
    final_implicit = implicit_in_y
    if d.get("ics"):
        ic = d["ics"][0]
        x0 = ic["x0"]
        y0 = ic["value"]

        ic_on_singular = None
        for z0 in singular_z_vals:
            if sp.simplify(z0 * x0 - y0) == 0:
                ic_on_singular = z0
                break

        c_val = sp.solve(implicit_in_y.subs({y_sym: y0, x: x0}), C1)
        if c_val:
            if general_sol is not None:
                best = general_sol
                for sol in explicit_sols:
                    test = sol.subs(C1, c_val[0]).subs(x, x0)
                    if sp.simplify(test - y0) == 0:
                        best = sol
                        break
                general_sol = best
                final_sol = sp.simplify(general_sol.subs(C1, c_val[0]))
            final_implicit = implicit_in_y.subs(C1, c_val[0])
            if want_steps:
                steps.append(Step(
                    "Apply initial condition",
                    rf"y({x0}) = {y0} \Rightarrow C_1 = {to_latex(c_val[0])}",
                    "Substituting general solution"))
                steps.append(Step("Particular solution",
                    to_latex(final_implicit), ""))
        elif ic_on_singular is not None:
            final_sol = ic_on_singular * x
            final_implicit = sp.Eq(y_sym, ic_on_singular * x)
            if want_steps:
                steps.append(Step(
                    "Apply IC",
                    rf"y({x0}) = {y0} \Rightarrow y = {to_latex(ic_on_singular)} x",
                    "IC matches constant solution"))

        if ic_on_singular is not None and c_val:
            warnings.append(f"IC also satisfied by constant solution y = {ic_on_singular}*x")

    d["domain_excludes"] = domain_excludes

    result = final_sol if final_sol is not None else final_implicit
    verified, verify_msg = add_domain_and_verify(
        steps, result, want_steps, want_verify, problem)

    answer = str(final_sol) if final_sol is not None else str(final_implicit)

    return Solution(
        kind="ode1_homogeneous_sub", final_answer=answer,
        answer_expr=result, steps=steps,
        verified=verified, verify_msg=verify_msg,
        warnings=warnings,
    )


def gen_ode1_homogeneous_sub(difficulty, with_ics=True):
    x = sp.Symbol("x")
    y = sp.Function("y")(x)
    z = sp.Symbol("z")

    if difficulty == "easy":
        f_z = random.choice([z + 1, z - 1, z + 2, z - 2])

    elif difficulty == "medium":
        f_z = random.choice([z + 1/z, (z + 1)/(z - 1)])

    else:
        f_z = random.choice([z**2 + z, z + z**2/2, z + 1/(2*z)])

    rhs = f_z.subs(z, y / x)

    eq = sp.Eq(y.diff(x), rhs)
    ics = []
    if with_ics:
        x0 = random.choice([1, 2])
        if difficulty == "easy":
            y0 = random.choice([-3, -2, -1, 1, 2, 3])
        elif difficulty == "medium":
            y0 = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        else:
            y0 = random.choice(list(range(-10, 0)) + list(range(1, 11)))
        ics = [make_ic("y", 0, x0, y0)]

    data = make_problem_data("x", ["y"], eq, ics)
    prompt = f"Solve y' = {rhs}"
    if ics:
        prompt += f", y({x0}) = {y0}"

    return Problem("ode1_homogeneous_sub", prompt, data,
                   metadata={"difficulty": difficulty})
