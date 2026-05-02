from .linear import solve_ode1_linear, gen_ode1_linear
from .bernoulli import solve_ode1_bernoulli, gen_ode1_bernoulli
from .separable import solve_ode1_separable, gen_ode1_separable
from .homogeneous_sub import solve_ode1_homogeneous_sub, gen_ode1_homogeneous_sub
from .cc_hom import solve_ode2_cc_hom, gen_ode2_cc_hom
from .cc_inhom import solve_ode2_cc_inhom, gen_ode2_cc_inhom
from .system import solve_ode_sys2_linear, gen_ode_sys2_linear

TOPICS = {
    "ode": [
        "ode1_linear",
        "ode1_bernoulli",
        "ode1_separable",
        "ode1_homogeneous_sub",
        "ode2_cc_hom",
        "ode2_cc_inhom",
        "ode_sys2_linear",
    ],
}

SOLVERS = {
    "ode1_linear": solve_ode1_linear,
    "ode1_bernoulli": solve_ode1_bernoulli,
    "ode1_separable": solve_ode1_separable,
    "ode1_homogeneous_sub": solve_ode1_homogeneous_sub,
    "ode2_cc_hom": solve_ode2_cc_hom,
    "ode2_cc_inhom": solve_ode2_cc_inhom,
    "ode_sys2_linear": solve_ode_sys2_linear,
}

GENERATORS = {
    "ode1_linear": gen_ode1_linear,
    "ode1_bernoulli": gen_ode1_bernoulli,
    "ode1_separable": gen_ode1_separable,
    "ode1_homogeneous_sub": gen_ode1_homogeneous_sub,
    "ode2_cc_hom": gen_ode2_cc_hom,
    "ode2_cc_inhom": gen_ode2_cc_inhom,
    "ode_sys2_linear": gen_ode_sys2_linear,
}
