import pytest
import sympy as sp
import project as ml
from project import Problem, make_problem_data, make_ic


KINDS = ml.list_kinds()


#parametrize runs the test once per item
@pytest.mark.parametrize("kind", KINDS)
def test_generate(kind):
    p = ml.generate(kind, difficulty="easy", with_ics=True)
    assert p.kind == kind
    assert p.prompt
    assert p.prompt_latex


@pytest.mark.parametrize("kind", KINDS)
def test_solve(kind):
    p = ml.generate(kind, difficulty="easy", with_ics=True)
    res = ml.solve(p)
    assert res.final_answer
    assert res.answer_expr is not None
    assert len(res.steps) > 0


#split into 3 cause nested was being weird
@pytest.mark.parametrize("kind", KINDS)
def test_verify_easy(kind):
    p = ml.generate(kind, difficulty="easy", with_ics=True)
    res = ml.solve(p)
    assert res.verified is True, res.verify_msg


@pytest.mark.parametrize("kind", KINDS)
def test_verify_medium(kind):
    p = ml.generate(kind, difficulty="medium", with_ics=True)
    res = ml.solve(p)
    assert res.verified is True, res.verify_msg


@pytest.mark.parametrize("kind", KINDS)
def test_verify_hard(kind):
    p = ml.generate(kind, difficulty="hard", with_ics=True)
    res = ml.solve(p)
    assert res.verified is True, res.verify_msg


@pytest.mark.parametrize("kind", KINDS)
def test_no_ics(kind):
    p = ml.generate(kind, difficulty="easy", with_ics=False)
    res = ml.solve(p)
    assert res.final_answer


def test_mixed():
    p = ml.generate("mixed")
    assert p.kind in KINDS


def test_check():
    p = ml.generate("ode1_linear", with_ics=True)
    res = ml.solve(p)
    rep = ml.check(p, res)
    assert rep["verified"] is True
    assert rep["message"] == "passed"


def test_export_latex_with_problem():
    p = ml.generate("ode1_linear")
    res = ml.solve(p)
    out = ml.export_latex(res, p)
    assert r"\textbf{Problem}" in out
    assert r"\textbf{Solution}" in out
    assert r"\textbf{Answer}" in out


def test_export_latex_solo():
    p = ml.generate("ode1_linear")
    res = ml.solve(p)
    out = ml.export_latex(res)
    assert r"\textbf{Problem}" not in out
    assert r"\textbf{Solution}" in out


def test_kinds_count():
    kinds = ml.list_kinds()
    assert len(kinds) == 7
    assert "ode1_linear" in kinds


def test_topics():
    assert "ode" in ml.list_topics()


def test_unknown_kind_generate():
    with pytest.raises(ml.TutorError):
        ml.generate("not_a_real_kind")


def test_unknown_kind_solve():
    p = ml.generate("ode1_linear")
    p.kind = "fake"
    with pytest.raises(ml.TutorError):
        ml.solve(p)


def test_problem_sheet():
    qs = ml.problem_sheet(n=3, kind="ode1_linear", difficulty="easy")
    assert len(qs) == 3
    for q in qs:
        assert q["problem"].kind == "ode1_linear"
        assert q["solution"].verified is True


#edge cases
def test_resonance():
    x = sp.Symbol("x")
    y = sp.Function("y")
    eq = sp.Eq(y(x).diff(x, 2) + y(x), sp.sin(x))
    data = make_problem_data("x", ["y"], eq, [
        make_ic("y", 0, 0, 0),
        make_ic("y", 1, 0, 0),
    ])
    p = Problem("ode2_cc_inhom", "y'' + y = sin(x), y(0)=0, y'(0)=0", data)
    res = ml.solve(p)
    assert res.verified is True, res.verify_msg
    assert sp.simplify(res.answer_expr.subs(x, 0)) == 0


def test_repeated_root():
    x = sp.Symbol("x")
    y = sp.Function("y")
    eq = sp.Eq(y(x).diff(x, 2) - 2*y(x).diff(x) + y(x), 0)
    data = make_problem_data("x", ["y"], eq, [
        make_ic("y", 0, 0, 1),
        make_ic("y", 1, 0, 3),
    ])
    p = Problem("ode2_cc_hom", "y'' - 2y' + y = 0, y(0)=1, y'(0)=3", data)
    res = ml.solve(p)
    assert res.verified is True, res.verify_msg
    assert sp.simplify(res.answer_expr.subs(x, 0) - 1) == 0


def test_bernoulli_domain():
    x = sp.Symbol("x")
    y = sp.Function("y")
    eq = sp.Eq(y(x).diff(x) + y(x), y(x)**2)
    data = make_problem_data("x", ["y"], eq, [
        make_ic("y", 0, 0, 2),
    ])
    data["n"] = 2
    p = Problem("ode1_bernoulli", "y' + y = y^2, y(0) = 2", data)
    res = ml.solve(p)
    assert res.verified is True, res.verify_msg
    dom, _ = ml.maximal_domain(p, res.answer_expr)
    assert dom is not None
    assert dom.contains(0)
    assert not dom.contains(sp.log(2))


SPECIAL = (sp.Ei, sp.erf, sp.erfi, sp.Si, sp.Ci, sp.Shi, sp.Chi,
           sp.gamma, sp.uppergamma, sp.meijerg, sp.hyper, sp.expint)


#these slip in randomly so spam test
@pytest.mark.parametrize("i", range(60))
def test_bernoulli_no_special(i):
    import random
    random.seed(1000 + i)
    diff = ["easy", "medium", "hard"][i % 3]
    p = ml.generate("ode1_bernoulli", difficulty=diff, with_ics=True)
    res = ml.solve(p)
    ans = res.answer_expr
    assert not ans.has(sp.Integral)
    for f in SPECIAL:
        assert not ans.has(f)
    assert res.verified is True


@pytest.mark.parametrize("difficulty", ["easy", "medium", "hard"])
def test_bernoulli_n_valid(difficulty):
    import random
    random.seed(999)
    for _ in range(20):
        p = ml.generate("ode1_bernoulli", difficulty=difficulty, with_ics=True)
        n = p.data["n"]
        assert n not in (0, 1)


def test_bernoulli_ic_positive():
    import random
    random.seed(777)
    for _ in range(50):
        p = ml.generate("ode1_bernoulli", difficulty="hard", with_ics=True)
        ics = p.data.get("ics", [])
        if ics:
            assert ics[0]["x0"] > 0
