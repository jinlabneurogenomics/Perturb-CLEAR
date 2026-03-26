#!/usr/bin/env python3
"""
Morphometric_feature_extraction.py

Morphometric analysis from SWC reconstructions.

"""

import argparse
import io
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, Delaunay, cKDTree


# SWC input

def read_swc(path: Path) -> pd.DataFrame:
    cols = ["id", "type", "x", "y", "z", "radius", "parent"]
    try:
        df = pd.read_csv(
            path,
            comment="#",
            sep=r"\s+",
            names=cols,
            usecols=range(7),
            engine="python",
        )
    except Exception:
        text = Path(path).read_text(errors="ignore")
        lines = [ln for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
        if not lines:
            raise ValueError(f"No SWC data rows found in {path}")
        buf = io.StringIO("\n".join(lines))
        tmp = pd.read_csv(buf, sep=r"\s+", header=None, engine="python")
        if tmp.shape[1] < 7:
            raise ValueError(
                f"SWC file appears malformed (found {tmp.shape[1]} columns, need >=7): {path}"
            )
        df = tmp.iloc[:, :7].copy()
        df.columns = cols

    df["id"] = df["id"].astype(int)
    df["type"] = df["type"].astype(int)
    df["parent"] = df["parent"].astype(int)
    for c in ["x", "y", "z", "radius"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["x", "y", "z", "radius"]).copy()
    df["radius"] = df["radius"].astype(float).fillna(0.0)
    return df


# Core helpers

def soma_center(df: pd.DataFrame) -> np.ndarray:
    roots = df[df["parent"] == -1]
    if len(roots):
        r = roots.iloc[0]
        return np.array([r.x, r.y, r.z], float)
    soma = df[df["type"] == 1]
    if len(soma):
        return soma[["x", "y", "z"]].mean().to_numpy(float)
    return df[["x", "y", "z"]].mean().to_numpy(float)


def build_edges(df: pd.DataFrame):
    by = df.set_index("id")
    edges = []
    for _, ch in df.iterrows():
        pid = int(ch.parent)
        if pid == -1 or pid not in by.index:
            continue
        pa = by.loc[pid]
        p0 = np.array([pa.x, pa.y, pa.z], float)
        p1 = np.array([ch.x, ch.y, ch.z], float)
        L = float(np.linalg.norm(p1 - p0))
        edges.append(
            dict(
                pid=int(pid),
                cid=int(ch.id),
                type=int(ch.type),
                p0=p0,
                p1=p1,
                L=L,
            )
        )
    return edges


def selected_edges(edges, type_filter=None):
    if type_filter is None:
        return list(edges)
    t = int(type_filter)
    return [e for e in edges if e["type"] == t]


def typed_node_ids(df: pd.DataFrame, type_filter=None):
    if type_filter is None:
        return set(df["id"].astype(int).tolist())
    t = int(type_filter)
    return set(df[df["type"] == t]["id"].astype(int).tolist())


def analysis_node_ids(df: pd.DataFrame, edges, type_filter=None):
    """Return the node set used by geometric outputs.

    For ALL, use every node.
    For a filtered type, use all nodes touched by selected edges.
    This includes support parents at the starts of selected segments.
    Isolated typed nodes are added back as singletons.
    """
    if type_filter is None:
        return set(df["id"].astype(int).tolist())

    nodes = set()
    for e in selected_edges(edges, type_filter):
        nodes.update([e["pid"], e["cid"]])

    nodes.update(typed_node_ids(df, type_filter) - nodes)
    return nodes


def build_subgraph(df, edges, type_filter=None):
    """Build an undirected graph from the selected edges.

    Edge type follows the child node type.
    With a type filter, parent support nodes are kept if they anchor selected edges.
    """
    E = selected_edges(edges, type_filter)
    adj = {}
    elen = {}
    nodes = set()
    for e in E:
        u, v = e["pid"], e["cid"]
        L = e["L"]
        nodes.update([u, v])
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
        elen[frozenset((u, v))] = L
    return nodes, adj, elen


def components_for_type(df, edges, type_filter=None):
    """Return connected components for the selected graph.

    For type-specific runs, isolated typed nodes are added as singleton components.
    Components are returned in a stable order by minimum node id.
    """
    nodes, adj, _ = build_subgraph(df, edges, type_filter)
    comps = []
    seen = set()

    for n in sorted(nodes):
        if n in seen:
            continue
        stack = [n]
        seen.add(n)
        comp = {n}
        while stack:
            u = stack.pop()
            for v in sorted(adj.get(u, [])):
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
                    comp.add(v)
        comps.append(comp)

    if type_filter is not None:
        present = set().union(*comps) if comps else set()
        isolated = sorted(typed_node_ids(df, type_filter) - present)
        comps.extend([{nid} for nid in isolated])

    comps.sort(key=lambda c: (min(c), len(c)))
    return comps


def node_coords(df, nid):
    row = df.loc[df["id"] == nid]
    if row.empty:
        return None
    r = row.iloc[0]
    return np.array([r.x, r.y, r.z], float)


# Terminals and stems

def count_terminals(df, edges, type_filter=None):
    """Count terminal nodes in the selected graph.

    A node is terminal if it appears as a child but never as a parent.
    Isolated typed nodes count as one terminal in type-specific outputs.
    """
    if type_filter is None:
        parents = {e["pid"] for e in edges}
        children = {e["cid"] for e in edges}
        soma_ids = set(df[df["type"] == 1]["id"].astype(int).tolist())
        isolated = set(df["id"].astype(int).tolist()) - (parents | children | soma_ids)
    else:
        E = selected_edges(edges, type_filter)
        parents = {e["pid"] for e in E}
        children = {e["cid"] for e in E}
        isolated = typed_node_ids(df, type_filter) - (parents | children)
    terminals = children - parents
    return len(terminals) + len(isolated)


def count_primary_stems(df, type_filter=None):
    """Count soma-to-dendrite stems.

    For ALL, this counts soma-to-{3,4,9} edges.
    For a filtered type, it counts soma-to-that-type edges.
    """
    by = df.set_index("id")
    stems = 0
    target_types = {3, 4, 9} if type_filter is None else {int(type_filter)}
    for _, ch in df.iterrows():
        pid = int(ch.parent)
        ctype = int(ch.type)
        if pid in by.index and int(by.loc[pid]["type"]) == 1 and ctype in target_types:
            stems += 1
    return stems


def trunk_total_length(edges):
    return float(sum(e["L"] for e in edges if e["type"] == 5))


def dendrite_total_length(edges, type_filter=None):
    dend_types = {3, 4, 9}
    if type_filter is None:
        return float(sum(e["L"] for e in edges if e["type"] in dend_types))
    t = int(type_filter)
    if t in dend_types:
        return float(sum(e["L"] for e in edges if e["type"] == t))
    return float("nan")


# Branch order

def _max_branch_order_from_component_edges(comp_edges, comp_nodes):
    if not comp_nodes:
        return float("nan")
    if not comp_edges:
        return 0.0

    by_parent = {}
    parents = set()
    children = set()
    for e in comp_edges:
        by_parent.setdefault(e["pid"], []).append(e["cid"])
        parents.add(e["pid"])
        children.add(e["cid"])

    roots = sorted(parents - children)
    if not roots:
        roots = sorted(children - parents)
    if not roots:
        roots = sorted(comp_nodes)

    bif = {p: (len(cs) >= 2) for p, cs in by_parent.items()}
    local_max = 0
    visited = set()
    stack = [(r, 0) for r in roots]

    while stack:
        u, order = stack.pop()
        if u in visited:
            continue
        visited.add(u)
        local_max = max(local_max, order)
        inc = 1 if bif.get(u, False) else 0
        for v in by_parent.get(u, []):
            stack.append((v, order + inc))

    return float(local_max)


def max_branch_order(df, edges, type_filter=None, comps=None):
    """Return overall max order and per-component orders.

    Per-component orders follow the same component order as components_for_type().
    """
    if comps is None:
        comps = components_for_type(df, edges, type_filter)

    E = selected_edges(edges, type_filter)
    per_component = []
    for comp in comps:
        comp_edges = [e for e in E if e["pid"] in comp and e["cid"] in comp]
        per_component.append(_max_branch_order_from_component_edges(comp_edges, comp))

    if not per_component:
        return float("nan"), []
    return float(np.nanmax(per_component)), [float(x) for x in per_component]


# Sholl analysis

def segment_sphere_hits(p0, p1, c, r, eps=1e-9):
    v = p1 - p0
    a = float(np.dot(v, v))
    if a <= eps:
        return 0
    u = p0 - c
    b = 2.0 * float(np.dot(u, v))
    c0 = float(np.dot(u, u) - r * r)
    D = b * b - 4 * a * c0
    if D < -1e-12:
        return 0
    if abs(D) <= 1e-12:
        t = -b / (2 * a)
        return 1 if (eps < t < 1 - eps) else 0
    sD = math.sqrt(max(0.0, D))
    t1 = (-b - sD) / (2 * a)
    t2 = (-b + sD) / (2 * a)
    return int(eps < t1 < 1 - eps) + int(eps < t2 < 1 - eps)


def sholl_counts(edges, center, radii, type_filter=None):
    E = selected_edges(edges, type_filter)
    out = np.zeros(len(radii), dtype=int)
    for i, r in enumerate(radii):
        s = 0
        for e in E:
            s += segment_sphere_hits(e["p0"], e["p1"], center, r)
        out[i] = s
    return out


def sholl_scalar_metrics(r, c):
    """Return scalar summaries of a Sholl profile.

    The semilog fit uses log(count + 1) versus radius on nonzero bins.
    """
    r = np.asarray(r, float)
    c = np.asarray(c, float)
    out = {}
    auc = float(np.trapezoid(c, r)) if len(r) >= 2 else float(c.sum())
    out["sholl_auc"] = auc

    if c.size:
        i_peak = int(np.argmax(c))
        out["sholl_peak"] = int(c[i_peak])
        out["sholl_radius_at_peak"] = float(r[i_peak])
    else:
        out["sholl_peak"] = 0
        out["sholl_radius_at_peak"] = float("nan")

    mask = c > 0
    if mask.sum() >= 2:
        y = np.log(c[mask] + 1.0)
        x = r[mask]
        X = np.vstack([np.ones_like(x), x]).T
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        a, k = beta
        yhat = X @ beta
        ss_res = float(((y - yhat) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        out["sholl_semilog_slope"] = float(k)
        out["sholl_semilog_intercept"] = float(a)
        out["sholl_semilog_r2"] = float(r2)
    else:
        out["sholl_semilog_slope"] = float("nan")
        out["sholl_semilog_intercept"] = float("nan")
        out["sholl_semilog_r2"] = float("nan")

    S = c.sum()
    if S > 0:
        out["sholl_radius_com"] = float((r * c).sum() / S)
        p = c.cumsum() / S

        def interp_at(q):
            return float(np.interp(q, p, r))

        out["sholl_r10"] = interp_at(0.10)
        out["sholl_r50"] = interp_at(0.50)
        out["sholl_r90"] = interp_at(0.90)
        out["sholl_width_r90_r10"] = out["sholl_r90"] - out["sholl_r10"]
        p_nonzero = c[c > 0] / S
        out["sholl_entropy"] = float(-(p_nonzero * np.log(p_nonzero)).sum())
    else:
        for k in [
            "sholl_radius_com",
            "sholl_r10",
            "sholl_r50",
            "sholl_r90",
            "sholl_width_r90_r10",
            "sholl_entropy",
        ]:
            out[k] = float("nan")
    return out


# Territory volumes

def hull3d_volume_and_faces(XYZ: np.ndarray):
    if XYZ.shape[0] < 4:
        return float("nan"), None
    try:
        hull = ConvexHull(XYZ)
        faces = hull.simplices
        return float(hull.volume), faces
    except Exception:
        return float("nan"), None


def alphashape3d_from_delaunay(tri: Delaunay, XYZ: np.ndarray, alpha: float):
    if not np.isfinite(alpha) or alpha <= 0 or tri is None:
        return float("nan"), np.zeros((0, 3), dtype=int)
    tets = tri.simplices
    if tets.size == 0:
        return float("nan"), np.zeros((0, 3), dtype=int)

    def _circ_r(P):
        A = 2.0 * (P[1:] - P[0])
        b = np.sum(P[1:] ** 2, axis=1) - np.sum(P[0] ** 2)
        try:
            c = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return np.inf
        return float(np.linalg.norm(c - P[0]))

    kept = []
    vol = 0.0
    for tet in tets:
        P = XYZ[tet]
        r = _circ_r(P)
        if r <= alpha:
            kept.append(tet)
            vol += abs(np.dot(P[3] - P[0], np.cross(P[1] - P[0], P[2] - P[0]))) / 6.0

    if not kept:
        return float(vol), np.zeros((0, 3), dtype=int)

    from collections import Counter

    faces = []
    for tet in kept:
        a, b, c, d = tet
        faces += [
            tuple(sorted((a, b, c))),
            tuple(sorted((a, b, d))),
            tuple(sorted((a, c, d))),
            tuple(sorted((b, c, d))),
        ]
    cnt = Counter(faces)
    boundary = [f for f, c in cnt.items() if c == 1]
    return float(vol), np.array(boundary, dtype=int)


def auto_alpha_for_type(df: pd.DataFrame, type_filter=None, method="seglen", factor=2.0):
    pts = []
    seg_lens = []
    by = df.set_index("id")
    for _, ch in df.iterrows():
        pid = int(ch["parent"])
        if pid == -1:
            continue
        if (type_filter is not None) and (int(ch["type"]) != type_filter):
            continue
        p0 = np.array([by.loc[pid]["x"], by.loc[pid]["y"], by.loc[pid]["z"]], float)
        p1 = np.array([ch["x"], ch["y"], ch["z"]], float)
        seg_lens.append(float(np.linalg.norm(p1 - p0)))
        pts.append(p0)
        pts.append(p1)
    if not pts:
        return float("nan")
    XYZ = np.unique(np.array(pts), axis=0)
    if method == "nn":
        if XYZ.shape[0] < 3:
            L = np.median(seg_lens) if len(seg_lens) else 0.0
            return factor * float(L)
        tree = cKDTree(XYZ)
        d, _ = tree.query(XYZ, k=2)
        nn = d[:, 1]
        return factor * float(np.median(nn))
    L = np.median(seg_lens) if len(seg_lens) else 0.0
    return factor * float(L)


# Branch paths and straightness

def branch_paths(df, edges, type_filter=None):
    nodes, adj, elen = build_subgraph(df, edges, type_filter)
    if not nodes:
        return []
    deg = {n: len(adj.get(n, [])) for n in nodes}
    endpoints = {n for n in nodes if deg.get(n, 0) != 2}
    visited = set()
    branches = []
    for a in sorted(endpoints):
        for b in sorted(adj.get(a, [])):
            e = frozenset((a, b))
            if e in visited:
                continue
            path = [a, b]
            Lsum = elen[e]
            visited.add(e)
            prev, cur = a, b
            while True:
                if deg.get(cur, 0) != 2:
                    break
                nbrs = sorted(adj[cur])
                nxt = nbrs[0] if nbrs[1] == prev else nbrs[1]
                e2 = frozenset((cur, nxt))
                if e2 in visited:
                    break
                visited.add(e2)
                Lsum += elen[e2]
                path.append(nxt)
                prev, cur = cur, nxt
            branches.append(dict(path=path, length=Lsum))
    return branches


def tortuosity_contraction_for_branch(df, branch):
    path = branch["path"]
    pts = np.array([node_coords(df, nid) for nid in path])
    if len(pts) < 2:
        return float("nan"), float("nan")
    chord = float(np.linalg.norm(pts[-1] - pts[0]))
    path_len = float(branch["length"])
    if chord <= 0 or path_len <= 0:
        return float("nan"), float("nan")
    tort = path_len / chord
    contr = chord / path_len
    return tort, contr


# Bifurcation geometry

def unit_vec(a, b):
    v = b - a
    n = np.linalg.norm(v)
    if n == 0:
        return None
    return v / n


def angle_between(u, v):
    if u is None or v is None:
        return float("nan")
    d = np.clip(float(np.dot(u, v)), -1.0, 1.0)
    return math.degrees(math.acos(d))


def bifurcation_nodes_for_type(df, edges, type_filter=None):
    by_parent = {}
    for e in selected_edges(edges, type_filter):
        by_parent.setdefault(e["pid"], []).append(e["cid"])
    bif_nodes = [nid for nid, chs in by_parent.items() if len(chs) >= 2]
    return bif_nodes, by_parent


def parent_of(df, nid):
    row = df.loc[df["id"] == nid]
    if row.empty:
        return None
    pid = int(row.iloc[0]["parent"])
    return pid if pid != -1 else None


def upstream_bifurcation(df, edges, node_id, type_filter=None, max_hops=50):
    bif_nodes, _ = bifurcation_nodes_for_type(df, edges, type_filter)
    is_bif = set(bif_nodes)
    cur = parent_of(df, node_id)
    hops = 0
    while cur is not None and hops < max_hops:
        if cur in is_bif:
            return cur
        cur = parent_of(df, cur)
        hops += 1
    return None


def bifurcation_metrics(df, edges, type_filter=None):
    bif_nodes, by_parent = bifurcation_nodes_for_type(df, edges, type_filter)
    rows = []
    for nid in bif_nodes:
        P = node_coords(df, nid)
        child_ids = by_parent.get(nid, [])
        child_vecs = [unit_vec(P, node_coords(df, cid)) for cid in child_ids]

        angle_min = float("nan")
        angle_max = float("nan")
        if len(child_vecs) >= 2:
            angles = []
            for i in range(len(child_vecs)):
                for j in range(i + 1, len(child_vecs)):
                    ang = angle_between(child_vecs[i], child_vecs[j])
                    if np.isfinite(ang):
                        angles.append(ang)
            if angles:
                angle_min = float(np.min(angles))
                angle_max = float(np.max(angles))

        pid = parent_of(df, nid)
        tilt_min = float("nan")
        if pid is not None and len(child_vecs) >= 1:
            U = unit_vec(node_coords(df, pid), P)
            t_angles = [angle_between(U, cv) for cv in child_vecs if cv is not None]
            t_angles = [a for a in t_angles if np.isfinite(a)]
            if t_angles:
                tilt_min = float(np.min(t_angles))

        n_curr = None
        if len(child_vecs) >= 2:
            max_pair = None
            max_ang = -1
            for i in range(len(child_vecs)):
                for j in range(i + 1, len(child_vecs)):
                    ang = angle_between(child_vecs[i], child_vecs[j])
                    if np.isfinite(ang) and ang > max_ang:
                        max_ang = ang
                        max_pair = (child_vecs[i], child_vecs[j])
            if max_pair is not None:
                n = np.cross(max_pair[0], max_pair[1])
                if np.linalg.norm(n) > 0:
                    n_curr = n / np.linalg.norm(n)

        torque = float("nan")
        prev_bif = upstream_bifurcation(df, edges, nid, type_filter=type_filter)
        if prev_bif is not None:
            Pprev = node_coords(df, prev_bif)
            child_ids_prev = by_parent.get(prev_bif, [])
            child_vecs_prev = [unit_vec(Pprev, node_coords(df, cid)) for cid in child_ids_prev]
            n_prev = None
            if len(child_vecs_prev) >= 2:
                max_pair = None
                max_ang = -1
                for i in range(len(child_vecs_prev)):
                    for j in range(i + 1, len(child_vecs_prev)):
                        ang = angle_between(child_vecs_prev[i], child_vecs_prev[j])
                        if np.isfinite(ang) and ang > max_ang:
                            max_ang = ang
                            max_pair = (child_vecs_prev[i], child_vecs_prev[j])
                if max_pair is not None:
                    n = np.cross(max_pair[0], max_pair[1])
                    if np.linalg.norm(n) > 0:
                        n_prev = n / np.linalg.norm(n)
            if n_prev is not None and n_curr is not None:
                torque = angle_between(n_prev, n_curr)

        rows.append(
            dict(
                node_id=int(nid),
                n_children=len(child_vecs),
                bif_angle_min_deg=angle_min,
                bif_angle_max_deg=angle_max,
                bif_tilt_min_deg=tilt_min,
                bif_torque_deg=torque,
            )
        )
    return rows


# Output names

def type_prefix(t):
    if t == 4:
        return "Tuft_"
    if t == 3:
        return "Basal_"
    if t == 9:
        return "Oblique_"
    return ""


# Summary stats

def series_stats(arr):
    a = np.asarray(arr, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return dict(mean=np.nan, median=np.nan, sd=np.nan, min=np.nan, max=np.nan, n=0)
    sd = float(np.std(a, ddof=1)) if a.size >= 2 else np.nan
    return dict(
        mean=float(np.mean(a)),
        median=float(np.median(a)),
        sd=sd,
        min=float(np.min(a)),
        max=float(np.max(a)),
        n=int(a.size),
    )


# Global axes for width, height, and depth

def pca3d_axes(XYZ: np.ndarray):
    if XYZ.shape[0] < 2:
        return np.eye(3)
    C = np.cov((XYZ - XYZ.mean(axis=0)).T)
    w, V = np.linalg.eigh(C)
    idx = np.argsort(w)[::-1]
    V = V[:, idx]
    for k in range(3):
        n = np.linalg.norm(V[:, k])
        if n > 0:
            V[:, k] /= n
    return V


def trunk_direction_3d_length_weighted(df: pd.DataFrame):
    by = df.set_index("id")
    acc = np.array([0.0, 0.0, 0.0], float)
    nseg = 0
    for _, ch in df.iterrows():
        if int(ch["type"]) != 5:
            continue
        pid = int(ch["parent"])
        if pid == -1 or pid not in by.index:
            continue
        p0 = np.array([by.loc[pid]["x"], by.loc[pid]["y"], by.loc[pid]["z"]], float)
        p1 = np.array([ch["x"], ch["y"], ch["z"]], float)
        acc += p1 - p0
        nseg += 1
    nrm = np.linalg.norm(acc)
    if nseg == 0 or nrm <= 1e-9:
        return None, nseg
    return acc / nrm, nseg


def compute_global_axes(df: pd.DataFrame):
    """Set the frame for width, height, and depth.

    If trunk segments are present, use trunk direction as y.
    Otherwise, fall back to PCA on all node coordinates.
    """
    XYZ_all = df[["x", "y", "z"]].to_numpy(float)
    V_all = pca3d_axes(XYZ_all)
    df_trunk = df[df["type"] == 5]
    if len(df_trunk) >= 2:
        ydir, _ = trunk_direction_3d_length_weighted(df)
        if ydir is None:
            V_tr = pca3d_axes(df_trunk[["x", "y", "z"]].to_numpy(float))
            ydir = V_tr[:, 0]
        PCs = [V_all[:, 0], V_all[:, 1], V_all[:, 2]]
        dots = [abs(float(np.dot(ydir, pk))) for pk in PCs]
        xcand = PCs[int(np.argmin(dots))]
        x = xcand - float(np.dot(xcand, ydir)) * ydir
        nx = np.linalg.norm(x)
        if nx <= 1e-9:
            for pk in PCs:
                x = pk - float(np.dot(pk, ydir)) * ydir
                nx = np.linalg.norm(x)
                if nx > 1e-9:
                    break
            if nx <= 1e-9:
                arb = np.array([1.0, 0.0, 0.0])
                if abs(float(np.dot(arb, ydir))) > 0.9:
                    arb = np.array([0.0, 1.0, 0.0])
                x = arb - float(np.dot(arb, ydir)) * ydir
                nx = np.linalg.norm(x)
        x /= nx + 1e-12
        z = np.cross(x, ydir)
        z /= np.linalg.norm(z) + 1e-12
        x = np.cross(ydir, z)
        x /= np.linalg.norm(x) + 1e-12
        return x, ydir, z

    e1, e2, e3 = V_all[:, 0], V_all[:, 1], V_all[:, 2]
    e1 = e1 / (np.linalg.norm(e1) + 1e-12)
    e2 = e2 - float(np.dot(e2, e1)) * e1
    e2 /= np.linalg.norm(e2) + 1e-12
    e3 = np.cross(e1, e2)
    e3 /= np.linalg.norm(e3) + 1e-12
    return e1, e2, e3


def extents_along_axes(XYZ: np.ndarray, ex, ey, ez):
    if XYZ.size == 0:
        return float("nan"), float("nan"), float("nan")
    px = XYZ @ ex
    py = XYZ @ ey
    pz = XYZ @ ez
    return float(px.max() - px.min()), float(py.max() - py.min()), float(pz.max() - pz.min())


# Main workflow

def main():
    ap = argparse.ArgumentParser(description="Morphometrics for SWC files.")
    ap.add_argument("input", type=str, help="SWC file or directory")
    ap.add_argument("-o", "--outdir", type=str, default="morph_out_v4_8", help="Output directory")
    ap.add_argument("--sholl-r0", type=float, default=0.0)
    ap.add_argument("--sholl-step", type=float, default=10.0)
    ap.add_argument("--sholl-rmax", type=float, default=-1.0, help="<=0 means auto per file")
    ap.add_argument("--alpha", type=float, default=-1.0, help="alpha (µm) for alpha-shape; if <=0 use auto")
    ap.add_argument(
        "--alpha-from",
        choices=["seglen", "nn"],
        default="seglen",
        help="auto alpha base (median segment length or nearest-neighbor)",
    )
    ap.add_argument("--alpha-factor", type=float, default=2.0, help="multiplier for auto alpha")
    ap.add_argument("--types", type=str, default=None, help="comma-separated SWC types to analyze; default: all present")
    args = ap.parse_args()

    inpath = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    files = [inpath] if (inpath.is_file() and inpath.suffix.lower() == ".swc") else sorted(inpath.rglob("*.swc"))

    types_filter = None
    if args.types:
        types_filter = [int(x.strip()) for x in args.types.split(",") if x.strip()]

    rows_summary = []
    rows_sholl = []
    rows_branches = []
    rows_branch_geom = []
    rows_bif = []
    rows_components = []

    for f in files:
        df = read_swc(f)
        edges = build_edges(df)
        center = soma_center(df)
        ex, ey, ez = compute_global_axes(df)

        if args.sholl_rmax is None or args.sholl_rmax <= 0:
            dists = np.linalg.norm(df[["x", "y", "z"]].to_numpy(float) - center[None, :], axis=1)
            rmax = float(dists.max()) if dists.size else 0.0
        else:
            rmax = float(args.sholl_rmax)
        r0 = float(args.sholl_r0)
        dr = float(args.sholl_step) if args.sholl_step > 0 else 10.0
        nsteps = max(1, int(math.floor((rmax - r0) / dr)) + 1)
        radii = r0 + dr * np.arange(nsteps, dtype=float)

        types_present = sorted(set(int(e["type"]) for e in edges))
        if types_filter:
            types_present = sorted(set(types_present).intersection(types_filter))

        trunk_len = trunk_total_length(edges)

        for t in [None] + types_present:
            tlabel = "ALL" if t is None else str(t)

            counts = sholl_counts(edges, center, radii, type_filter=t)
            scalar = sholl_scalar_metrics(radii, counts)
            for r, c in zip(radii, counts):
                rows_sholl.append(dict(file=f.name, type=tlabel, radius=r, count=int(c)))

            branches = branch_paths(df, edges, t)
            blens = [b["length"] for b in branches]
            branch_torts = []
            branch_contrs = []
            for i, b in enumerate(branches):
                rows_branches.append(dict(file=f.name, type=tlabel, branch_id=i, length=float(b["length"])))
                tort, contr = tortuosity_contraction_for_branch(df, b)
                rows_branch_geom.append(
                    dict(file=f.name, type=tlabel, branch_id=i, tortuosity=tort, contraction=contr)
                )
                branch_torts.append(tort)
                branch_contrs.append(contr)
            stats_len = series_stats(blens)
            stats_tort = series_stats(branch_torts)
            stats_contr = series_stats(branch_contrs)

            bif_rows = bifurcation_metrics(df, edges, t)
            for rowb in bif_rows:
                rowb["file"] = f.name
                rowb["type"] = tlabel
                rows_bif.append(rowb)
            ang_arr = [row["bif_angle_min_deg"] for row in bif_rows]
            tilt_arr = [row["bif_tilt_min_deg"] for row in bif_rows]
            torq_arr = [row["bif_torque_deg"] for row in bif_rows]
            stats_ang = series_stats(ang_arr)
            stats_tilt = series_stats(tilt_arr)
            stats_torq = series_stats(torq_arr)

            metric_node_ids = analysis_node_ids(df, edges, t)
            nodes_df = df[df["id"].isin(metric_node_ids)] if t is not None else df
            XYZ = nodes_df[["x", "y", "z"]].to_numpy(float)
            vol_hull3d, _ = hull3d_volume_and_faces(XYZ)
            alpha_val = (
                float(args.alpha)
                if (args.alpha and args.alpha > 0)
                else float(auto_alpha_for_type(df, type_filter=t, method=args.alpha_from, factor=args.alpha_factor))
            )
            tri = None
            if XYZ.shape[0] >= 4:
                try:
                    tri = Delaunay(XYZ)
                except Exception:
                    tri = None
            if tri is not None:
                vol_alpha3d, _ = alphashape3d_from_delaunay(tri, XYZ, alpha_val)
            else:
                vol_alpha3d = float("nan")

            n_terms = count_terminals(df, edges, t)
            n_stems = count_primary_stems(df, t)

            comps = components_for_type(df, edges, t)
            max_order, per_comp_orders = max_branch_order(df, edges, t, comps=comps)
            n_comps = len(comps)

            W, H, D = extents_along_axes(XYZ, ex, ey, ez)

            vol_hull_sum = float("nan")
            vol_alpha_sum = float("nan")
            if t == 9:
                hsum = 0.0
                asum = 0.0
                for cid, comp in enumerate(comps):
                    sub_df = df[df["id"].isin(comp)]
                    XYZc = sub_df[["x", "y", "z"]].to_numpy(float)
                    vh, _ = hull3d_volume_and_faces(XYZc)
                    tri_c = None
                    if XYZc.shape[0] >= 4:
                        try:
                            tri_c = Delaunay(XYZc)
                        except Exception:
                            tri_c = None
                    if tri_c is not None:
                        va, _ = alphashape3d_from_delaunay(tri_c, XYZc, alpha_val)
                    else:
                        va = float("nan")
                    rows_components.append(
                        dict(
                            file=f.name,
                            type=tlabel,
                            component_id=cid,
                            n_nodes=len(sub_df),
                            convex_3D_volume=vh,
                            alpha_3D_volume=va,
                            max_branch_order=per_comp_orders[cid],
                        )
                    )
                    if np.isfinite(vh):
                        hsum += vh
                    if np.isfinite(va):
                        asum += va
                vol_hull_sum = hsum
                vol_alpha_sum = asum

            dend_len = dendrite_total_length(edges, t)

            row = dict(
                file=f.name,
                type=tlabel,
                **scalar,
                convex_3D_volume=vol_hull3d,
                alpha_3D_volume=vol_alpha3d,
                n_components=n_comps,
                convex_3D_volume_component_sum=vol_hull_sum,
                alpha_3D_volume_component_sum=vol_alpha_sum,
                n_branches=int(stats_len["n"]),
                branch_len_mean=stats_len["mean"],
                branch_len_median=stats_len["median"],
                branch_len_sd=stats_len["sd"],
                branch_len_min=stats_len["min"],
                branch_len_max=stats_len["max"],
                branch_tortuosity_mean=stats_tort["mean"],
                branch_tortuosity_median=stats_tort["median"],
                branch_tortuosity_sd=stats_tort["sd"],
                branch_tortuosity_max=stats_tort["max"],
                branch_contraction_mean=stats_contr["mean"],
                branch_contraction_median=stats_contr["median"],
                branch_contraction_sd=stats_contr["sd"],
                branch_contraction_min=stats_contr["min"],
                n_bifurcations=int(stats_ang["n"]),
                bif_childpair_angle_min_mean=stats_ang["mean"],
                bif_childpair_angle_min_median=stats_ang["median"],
                bif_childpair_angle_min_sd=stats_ang["sd"],
                bif_childpair_angle_min_min=stats_ang["min"],
                bif_childpair_angle_min_max=stats_ang["max"],
                bif_parent_child_tilt_min_mean=stats_tilt["mean"],
                bif_parent_child_tilt_min_median=stats_tilt["median"],
                bif_parent_child_tilt_min_sd=stats_tilt["sd"],
                bif_parent_child_tilt_min_min=stats_tilt["min"],
                bif_parent_child_tilt_min_max=stats_tilt["max"],
                bif_plane_torque_mean=stats_torq["mean"],
                bif_plane_torque_median=stats_torq["median"],
                bif_plane_torque_sd=stats_torq["sd"],
                bif_plane_torque_min=stats_torq["min"],
                bif_plane_torque_max=stats_torq["max"],
                n_terminals=int(n_terms),
                N_stems=int(n_stems),
                trunk_total_length=(trunk_len if t is None else float("nan")),
                max_branch_order=float(max_order),
                extent_width=float(W),
                extent_height=float(H),
                extent_depth=float(D),
                dendrite_total_length=float(dend_len),
            )
            rows_summary.append(row)

    pd.DataFrame(rows_summary).to_csv(outdir / "summary_metrics_tidy.csv", index=False)
    pd.DataFrame(rows_sholl).to_csv(outdir / "sholl_counts.csv", index=False)
    pd.DataFrame(rows_branches).to_csv(outdir / "branch_lengths.csv", index=False)
    pd.DataFrame(rows_branch_geom).to_csv(outdir / "branches_geom.csv", index=False)
    pd.DataFrame(rows_bif).to_csv(outdir / "bifurcations.csv", index=False)
    pd.DataFrame(rows_components).to_csv(outdir / "components.csv", index=False)

    tidy = pd.DataFrame(rows_summary)
    if not tidy.empty:
        all_rows = tidy[tidy["type"] == "ALL"].copy()
        typed = tidy[tidy["type"] != "ALL"].copy()
        metric_cols = [c for c in typed.columns if c not in ["file", "type"]]

        pieces = []
        for _, r in typed.iterrows():
            prefix = type_prefix(int(r["type"]))
            out = {"file": r["file"]}
            for c in metric_cols:
                if c in ("convex_3D_volume_component_sum", "alpha_3D_volume_component_sum") and int(r["type"]) != 9:
                    continue
                out[prefix + c] = r[c]
            pieces.append(out)
        wide_typed = pd.DataFrame(pieces)
        if not wide_typed.empty:
            wide_typed = wide_typed.groupby("file", as_index=False).first()

        g = all_rows.drop(columns=["type"]).copy()
        rename_map = {c: ("Global_" + c) for c in g.columns if c != "file"}
        g = g.rename(columns=rename_map)

        wide = g.merge(wide_typed, on="file", how="left") if not wide_typed.empty else g

        if wide.shape[0] > 0:
            mask_any = wide.drop(columns=["file"]).notna().any(axis=0)
            keep_cols = ["file"] + [c for c in wide.columns if c != "file" and bool(mask_any.get(c, False))]
            wide = wide[keep_cols]

        wide.to_csv(outdir / "summary_metrics_prefixed_wide.csv", index=False)

    print(f"[OK] Wrote outputs to {outdir}")


if __name__ == "__main__":
    main()
