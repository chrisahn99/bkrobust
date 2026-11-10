"""E3S: mechanical translations of the submitted E1 encoding (CP-SAT) to MILP and CNF.

The E1 model (``bkrobust.sat.e1._build``, criterion "gac") is exported to its
protocol buffer and translated constraint by constraint, with no
reformulation, so that every solver sees the same problem:

    minimise the number of retracted closure orientations, subject to the
    state being Meek-closed and Z failing (GAC) in a witness DAG of it.

Constraint types handled (the only ones the E1 model emits; anything else
raises): bool_or, bool_and, exactly_one, at_most_one, lin_max over 0/1
expressions, linear over Booleans with +-1 coefficients, linear over two
integer variables with coefficients (1, -1) (acyclicity by positions), each
possibly enforced by literals.

MILP: literal l -> x or 1 - x; enforced linear constraints by big-M from the
bounds of the expression; lin_max over 0/1 as t >= e_i, t <= sum e_i.

CNF: Booleans as variables; integer positions in the order encoding
(o[v][i] <-> pos_v >= i); Boolean linear constraints by sequential counters
(PySAT CardEnc), the enforcement literals added to every clause of the
encoding (sound: the auxiliary variables are fresh); position differences by
the order encoding.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field

from google.protobuf import text_format
from ortools.sat import cp_model_pb2


def export_proto(model) -> cp_model_pb2.CpModelProto:
    """The CpModelProto of an OR-Tools 9.15 CpModel (via its text export)."""
    with tempfile.NamedTemporaryFile(suffix=".pb.txt", delete=True) as f:
        ok = model.ExportToFile(f.name)
        assert ok, "ExportToFile failed"
        pb = cp_model_pb2.CpModelProto()
        text_format.Parse(open(f.name).read(), pb)
    return pb


def _dom(v):
    d = list(v.domain)
    return d[0], d[-1]


def _is_bool(pb, i):
    lo, hi = _dom(pb.variables[i])
    return lo >= 0 and hi <= 1


def _ref(lit):
    """(var index, positive?) for a literal reference."""
    return (lit, True) if lit >= 0 else (-lit - 1, False)


def _intervals(dom):
    return [(dom[i], dom[i + 1]) for i in range(0, len(dom), 2)]


# ------------------------------------------------------------------------ MILP


@dataclass
class MilpModel:
    solver: object
    x: list
    n_rows: int = 0


def to_milp(pb, objective_vars, backend="SCIP"):
    """pywraplp model minimising sum(objective_vars) (proto indices of Booleans)."""
    from ortools.linear_solver import pywraplp

    s = pywraplp.Solver.CreateSolver(backend)
    assert s is not None, backend
    inf = s.infinity()
    x = []
    for i, v in enumerate(pb.variables):
        lo, hi = _dom(v)
        if lo >= 0 and hi <= 1:
            x.append(s.BoolVar(f"b{i}") if lo != hi else s.NumVar(lo, hi, f"c{i}"))
        else:
            x.append(s.IntVar(lo, hi, f"i{i}"))
    bounds = [_dom(v) for v in pb.variables]

    def lit(l):
        i, pos = _ref(l)
        return x[i] if pos else 1 - x[i]

    def relax(enf):
        """sum over enforcement literals of (1 - e): 0 iff all enforced."""
        return sum((1 - lit(e)) for e in enf) if enf else 0

    rows = 0
    for c in pb.constraints:
        kind = c.WhichOneof("constraint")
        enf = list(c.enforcement_literal)
        r = relax(enf)
        if kind == "bool_or":
            s.Add(sum(lit(l) for l in c.bool_or.literals) + r >= 1)
            rows += 1
        elif kind == "bool_and":
            for l in c.bool_and.literals:
                s.Add(lit(l) + r >= 1)
                rows += 1
        elif kind in ("exactly_one", "at_most_one"):
            lits = c.exactly_one.literals if kind == "exactly_one" else c.at_most_one.literals
            e = sum(lit(l) for l in lits)
            n = len(lits)
            s.Add(e <= 1 + n * r)
            rows += 1
            if kind == "exactly_one":
                s.Add(e >= 1 - r)
                rows += 1
        elif kind == "linear":
            iv = _intervals(list(c.linear.domain))
            assert len(iv) == 1, "disjoint linear domain"
            lo, hi = iv[0]
            terms = list(zip(c.linear.vars, c.linear.coeffs))
            e = sum(k * lit(v) for v, k in terms)
            emin = sum(k * (bounds[_ref(v)[0]][0] if k > 0 else bounds[_ref(v)[0]][1]) for v, k in terms)
            emax = sum(k * (bounds[_ref(v)[0]][1] if k > 0 else bounds[_ref(v)[0]][0]) for v, k in terms)
            if hi < emax:
                s.Add(e <= hi + (emax - hi) * r)
                rows += 1
            if lo > emin:
                s.Add(e >= lo - (lo - emin) * r)
                rows += 1
        elif kind == "lin_max":
            assert not enf, "enforced lin_max"
            def expr(le):
                assert len(le.vars) <= 1 and all(_is_bool(pb, _ref(v)[0]) for v in le.vars)
                return le.offset + sum(k * lit(v) for v, k in zip(le.vars, le.coeffs))
            t = expr(c.lin_max.target)
            es = [expr(e) for e in c.lin_max.exprs]
            if not es:
                continue
            for e in es:
                s.Add(t >= e)
            s.Add(t <= sum(es))
            rows += len(es) + 1
        else:
            raise NotImplementedError(kind)
    s.Minimize(sum(x[i] for i in objective_vars))
    return MilpModel(s, x, rows)


# ------------------------------------------------------------------------- CNF


@dataclass
class CnfModel:
    hard: list = field(default_factory=list)
    top: int = 0
    var_of: dict = field(default_factory=dict)   # proto Boolean index -> CNF var
    soft: list = field(default_factory=list)     # CNF literals that should be FALSE (drops)

    def new(self):
        self.top += 1
        return self.top


def to_cnf(pb, objective_vars):
    from pysat.card import CardEnc, EncType

    m = CnfModel()
    order = {}  # int var index -> (lo, hi, {i: cnf var for pos >= i, lo < i <= hi})
    for i, v in enumerate(pb.variables):
        lo, hi = _dom(v)
        if lo >= 0 and hi <= 1:
            m.var_of[i] = m.new()
            if lo == hi:
                m.hard.append([m.var_of[i] if lo == 1 else -m.var_of[i]])
        else:
            o = {j: m.new() for j in range(lo + 1, hi + 1)}
            for j in range(lo + 2, hi + 1):
                m.hard.append([-o[j], o[j - 1]])
            order[i] = (lo, hi, o)

    def lit(l):
        i, pos = _ref(l)
        return m.var_of[i] if pos else -m.var_of[i]

    def ge(i, j):
        """Literal for pos_i >= j, or True/False as Python bools."""
        lo, hi, o = order[i]
        if j <= lo:
            return True
        if j > hi:
            return False
        return o[j]

    def add(clause, neg_enf):
        m.hard.append(neg_enf + clause)

    def card(lits, lo, hi, neg_enf):
        n = len(lits)
        if hi < n:
            if hi < 0:
                add([], neg_enf)
            else:
                enc = CardEnc.atmost(lits=lits, bound=hi, top_id=m.top, encoding=EncType.seqcounter)
                m.top = max(m.top, enc.nv)
                for cl in enc.clauses:
                    add(list(cl), neg_enf)
        if lo > 0:
            if lo > n:
                add([], neg_enf)
            elif lo == 1:
                add(list(lits), neg_enf)
            else:
                enc = CardEnc.atleast(lits=lits, bound=lo, top_id=m.top, encoding=EncType.seqcounter)
                m.top = max(m.top, enc.nv)
                for cl in enc.clauses:
                    add(list(cl), neg_enf)

    for c in pb.constraints:
        kind = c.WhichOneof("constraint")
        ne = [-lit(e) for e in c.enforcement_literal]
        if kind == "bool_or":
            add([lit(l) for l in c.bool_or.literals], ne)
        elif kind == "bool_and":
            for l in c.bool_and.literals:
                add([lit(l)], ne)
        elif kind == "exactly_one":
            card([lit(l) for l in c.exactly_one.literals], 1, 1, ne)
        elif kind == "at_most_one":
            card([lit(l) for l in c.at_most_one.literals], 0, 1, ne)
        elif kind == "linear":
            iv = _intervals(list(c.linear.domain))
            assert len(iv) == 1, "disjoint linear domain"
            lo, hi = iv[0]
            vs, ks = list(c.linear.vars), list(c.linear.coeffs)
            if all(_ref(v)[0] in m.var_of for v in vs):
                # sum k*l over literals, k in {+1,-1}: -l = (not l) - 1
                lits, off = [], 0
                for v, k in zip(vs, ks):
                    assert k in (1, -1), k
                    if k == 1:
                        lits.append(lit(v))
                    else:
                        lits.append(-lit(v))
                        off -= 1
                n = len(lits)
                card(lits, max(lo - off, 0) if lo > -(2 ** 62) else 0,
                     min(hi - off, n) if hi < 2 ** 62 else n, ne)
            else:
                assert len(vs) == 2 and sorted(ks) == [-1, 1] and all(_ref(v)[0] in order for v in vs) \
                    and all(v >= 0 for v in vs), (vs, ks)
                a, b = (vs[0], vs[1]) if ks[0] == 1 else (vs[1], vs[0])
                # a - b <= hi :  for every j, (a >= j) -> (b >= j - hi)
                def diff_le(a, b, c):
                    la, ha, _ = order[a]
                    for j in range(la + 1, ha + 2):
                        pa = ge(a, j) if j <= ha else False
                        if pa is False:
                            continue
                        pb_ = ge(b, j - c)
                        if pb_ is True:
                            continue
                        cl = ([] if pa is True else [-pa]) + ([] if pb_ is False else [pb_])
                        add(cl, ne)
                    # j = la: a >= la always true -> b >= la - c
                    pb_ = ge(b, la - c)
                    if pb_ is not True:
                        add([] if pb_ is False else [pb_], ne)
                if hi < 2 ** 62:
                    diff_le(a, b, hi)
                if lo > -(2 ** 62):
                    diff_le(b, a, -lo)
        elif kind == "lin_max":
            assert not c.enforcement_literal
            def one(le):
                assert len(le.vars) == 1 and _ref(le.vars[0])[0] in m.var_of
                k, off = le.coeffs[0], le.offset
                if k == 1 and off == 0:
                    return lit(le.vars[0])
                assert k == -1 and off == 1, (k, off)
                return -lit(le.vars[0])
            t = one(c.lin_max.target)
            es = [one(e) for e in c.lin_max.exprs]
            for e in es:
                m.hard.append([-e, t])
            m.hard.append([-t] + es)
        else:
            raise NotImplementedError(kind)
    m.soft = [m.var_of[i] for i in objective_vars]
    return m
