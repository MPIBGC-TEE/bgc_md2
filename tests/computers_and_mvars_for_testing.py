# This module is just providing fake classes (=mvars) and functions(=computers) to test
# the algoritms that do only use the information stored in the annotation. E.g. the computatioof the computability graph
class A_minus_1:
    pass


class A_minus_2:
    pass


class A:
    """A variable we assume to be given"""

    pass


class B:
    pass


class A_minus_2:
    pass


class A_minus_1:
    pass


class A0:
    pass


class A1:
    pass


class A2:
    pass


class A3:
    pass


class B_minus_2:
    pass


class B_minus_1:
    pass


class B0:
    pass


class B1:
    pass


class B2:
    pass


class B3:
    pass


class C:
    pass


class D:
    pass


class E:
    pass


class F:
    pass


class G:
    pass


class H:
    pass


class I:
    pass


class X:
    pass


class Y:
    pass


class Z:
    pass


def a_from_b_c(b: B, c: C) -> A:
    return A()


def a_from_b_d(b: B, d: D) -> A:
    return A()


def a_from_x(x: X) -> A:
    return A()


def b_from_x(x: X) -> B:
    return B()


def a_from_y(y: Y) -> A:
    return A()


def b_from_y(y: Y) -> B:
    return B()


def a_from_z(z: Z) -> A:
    return A()


def b_from_z(z: Z) -> B:
    return B()


def c_from_z(z: Z) -> C:
    return C()


def a_from_i(i: I) -> A:
    """Computes a from i"""
    return A()


def b_from_c_d(c: C, d: D) -> B:
    return B()


def b_from_d_e(d: D, e: E) -> B:
    return B()


def b_from_e_f(e: E, f: F) -> B:
    return B()


def c_from_b(b: B) -> C:
    """Computes c from b"""
    return C()


def d_from_a(a: A) -> D:
    return D()


def d_from_b(b: B) -> D:
    """Computes d from b"""
    return D()


def d_from_g_h(g: G, h: H) -> D:
    """Computes d from g and h"""
    return D()


def e_from_b(b: B) -> E:
    """Computes e from b"""
    return E()


def f_from_b(b: B) -> F:
    """Computes f from b"""
    return F()


def a_minus_1_from_a_minus_2(x: A_minus_2) -> A_minus_1:
    return A_minus_1()


def a0_from_a_minus_1(x: A_minus_1) -> A0:
    return A0()


def a1_from_a0(a0: A0) -> A1:
    return A1()


def a2_from_a1(a1: A1) -> A2:
    return A2()


def a3_from_a2(a2: A2) -> A3:
    return A3()


def a0_from_b0(x: B0) -> A0:
    return A0()


def b_minus_1_from_b_minus_2(x: B_minus_2) -> B_minus_1:
    return B_minus_1()


def b0_from_b_minus_1(x: B_minus_1) -> B0:
    return B0()


def b1_from_b0(b0: B0) -> B1:
    return B1()


def b2_from_b1(b1: B1) -> B2:
    return B2()


def b3_from_b2(b2: B2) -> B3:
    return B3()
