#!/usr/bin/env python3
"""
swc_branch_cleanup.py

The script fix basic parent-graph issues:
missing parents, self-parent nodes, parent cycles, and zero-length edges.
As well as tracing artifacts: over-branch repair and tip pruning.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set, Optional

EPS = 1e-9


@dataclass
class Node:
    nid: int
    ntype: int
    x: float
    y: float
    z: float
    r: float
    parent: int


# ------------------------------ I/O ------------------------------


def read_swc(path: str) -> Tuple[Dict[int, Node], List[str]]:
    nodes: Dict[int, Node] = {}
    comments: List[str] = []

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for ln, line in enumerate(f, 1):
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                comments.append(s)
                continue

            parts = s.split()
            if len(parts) < 7:
                raise ValueError(
                    f"Malformed SWC row at line {ln} in {path}: expected at least 7 fields, got {len(parts)}"
                )

            try:
                nid = int(float(parts[0]))
                ntype = int(float(parts[1]))
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                r = float(parts[5])
                parent = int(float(parts[6]))
            except Exception as e:
                raise ValueError(f"Parse error at line {ln} in {path}: {e}") from e

            if nid in nodes:
                raise ValueError(f"Duplicate node ID {nid} at line {ln} in {path}")
            if not all(math.isfinite(v) for v in (x, y, z, r)):
                raise ValueError(f"Non-finite coordinate or radius at line {ln} in {path}")

            nodes[nid] = Node(nid, ntype, x, y, z, r, parent)

    if not nodes:
        raise ValueError(f"No SWC nodes found in {path}")

    return nodes, comments


def write_swc(path: str, nodes: Dict[int, Node], comments: List[str], edits_log: List[str]) -> None:
    header = [
        "# swc_branch_cleanup.py",
        "# " + (" | ".join(edits_log) if edits_log else "(no changes)"),
    ]
    lines: List[str] = []
    lines.extend(comments)
    lines.extend(header)
    for nid in sorted(nodes):
        n = nodes[nid]
        lines.append(f"{n.nid} {n.ntype} {n.x:.6f} {n.y:.6f} {n.z:.6f} {n.r:.6f} {n.parent}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")




def clone_nodes(nodes: Dict[int, Node]) -> Dict[int, Node]:
    return {nid: Node(n.nid, n.ntype, n.x, n.y, n.z, n.r, n.parent) for nid, n in nodes.items()}

# ----------------------------- Helpers -----------------------------


def build_children(nodes: Dict[int, Node]) -> Dict[int, List[int]]:
    children: Dict[int, List[int]] = defaultdict(list)
    for n in nodes.values():
        if n.parent != -1 and n.parent in nodes:
            children[n.parent].append(n.nid)
    return children


def degree(children: Dict[int, List[int]], nid: int) -> int:
    return len(children.get(nid, []))


def neighbors(nodes: Dict[int, Node], children: Dict[int, List[int]], nid: int) -> List[int]:
    out = list(children.get(nid, []))
    p = nodes[nid].parent
    if p in nodes:
        out.append(p)
    return out


def leaves(nodes: Dict[int, Node], children: Dict[int, List[int]]) -> List[int]:
    return [nid for nid in nodes if len(children.get(nid, [])) == 0]


def aeq(a: float, b: float, eps: float = EPS) -> bool:
    return abs(a - b) <= eps


def euclid(a: Node, b: Node) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def zero_length_edges(nodes: Dict[int, Node]) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []
    for n in nodes.values():
        if n.parent != -1 and n.parent in nodes:
            p = nodes[n.parent]
            if aeq(n.x, p.x) and aeq(n.y, p.y) and aeq(n.z, p.z):
                out.append((n.parent, n.nid))
    return out


def overall_components(nodes: Dict[int, Node]) -> List[Set[int]]:
    adj: Dict[int, List[int]] = defaultdict(list)
    for nid, n in nodes.items():
        if n.parent in nodes:
            adj[nid].append(n.parent)
            adj[n.parent].append(nid)

    comps: List[Set[int]] = []
    seen: Set[int] = set()
    for nid in nodes:
        if nid in seen:
            continue
        comp: Set[int] = set()
        dq = deque([nid])
        seen.add(nid)
        while dq:
            u = dq.popleft()
            comp.add(u)
            for v in adj.get(u, []):
                if v not in seen:
                    seen.add(v)
                    dq.append(v)
        comps.append(comp)
    return comps


def subtree_nodes(nodes: Dict[int, Node], children: Dict[int, List[int]], nid: int) -> Set[int]:
    acc: Set[int] = set()
    stack = [nid]
    while stack:
        u = stack.pop()
        if u in acc:
            continue
        acc.add(u)
        stack.extend(children.get(u, []))
    return acc


def subtree_length(nodes: Dict[int, Node], children: Dict[int, List[int]], nid: int) -> float:
    length = 0.0
    sub = subtree_nodes(nodes, children, nid)
    for u in sub:
        p = nodes[u].parent
        if p != -1 and p in sub:
            length += euclid(nodes[u], nodes[p])
    return length


# ------------------------ Parent-graph checks ------------------------


def find_parent_cycles(nodes: Dict[int, Node]) -> List[List[int]]:
    state: Dict[int, int] = {}
    seen_keys: Set[Tuple[int, ...]] = set()
    cycles: List[List[int]] = []

    for start in nodes:
        if state.get(start, 0) == 2:
            continue

        trail: List[int] = []
        pos: Dict[int, int] = {}
        cur = start

        while True:
            if cur not in nodes:
                break

            st = state.get(cur, 0)
            if st == 2:
                break
            if st == 1:
                if cur in pos:
                    cyc = trail[pos[cur]:]
                    key = tuple(sorted(cyc))
                    if key not in seen_keys:
                        seen_keys.add(key)
                        cycles.append(cyc)
                break

            state[cur] = 1
            pos[cur] = len(trail)
            trail.append(cur)

            parent = nodes[cur].parent
            if parent == -1 or parent not in nodes:
                break
            cur = parent

        for nid in trail:
            state[nid] = 2

    return cycles


def fix_missing_parents(nodes: Dict[int, Node], edits_log: Optional[List[str]] = None) -> int:
    changed = 0
    for n in nodes.values():
        if n.parent != -1 and n.parent not in nodes:
            n.parent = -1
            changed += 1
    if edits_log is not None and changed:
        edits_log.append(f"missing parents: promoted {changed} nodes to roots")
    return changed


def fix_self_parent(nodes: Dict[int, Node], edits_log: Optional[List[str]] = None) -> int:
    changed = 0
    for n in nodes.values():
        if n.parent == n.nid:
            n.parent = -1
            changed += 1
    if edits_log is not None and changed:
        edits_log.append(f"self-parent: fixed {changed} nodes")
    return changed


def fix_parent_cycles(nodes: Dict[int, Node], edits_log: Optional[List[str]] = None) -> int:
    broken = 0
    examples: List[str] = []
    while True:
        cycles = find_parent_cycles(nodes)
        if not cycles:
            break
        for cyc in cycles:
            breaker = min(cyc)
            nodes[breaker].parent = -1
            broken += 1
            if len(examples) < 10:
                examples.append(f"{cyc} -> root {breaker}")
    if edits_log is not None and broken:
        msg = f"parent cycles: broke {broken} cycle(s)"
        if examples:
            msg += f" ({'; '.join(examples)})"
        edits_log.append(msg)
    return broken


def collapse_zero_length_edges(nodes: Dict[int, Node], edits_log: Optional[List[str]] = None) -> int:
    removed = 0
    while True:
        changed = False
        for p_id, c_id in zero_length_edges(nodes):
            if p_id not in nodes or c_id not in nodes:
                continue
            for ch in [nid for nid, n in nodes.items() if n.parent == c_id]:
                nodes[ch].parent = p_id
            nodes.pop(c_id, None)
            removed += 1
            changed = True
        if not changed:
            break
    if edits_log is not None and removed:
        edits_log.append(f"zero-length edges: collapsed {removed} nodes")
    return removed


# -------------------------- Tip pruning --------------------------


def twig_length_to_leaf(nodes: Dict[int, Node], children: Dict[int, List[int]], leaf_id: int) -> Tuple[float, int]:
    length = 0.0
    curr = leaf_id
    while True:
        parent = nodes[curr].parent
        if parent == -1 or parent not in nodes:
            return length, parent
        length += euclid(nodes[curr], nodes[parent])
        if len(children.get(parent, [])) != 1 or nodes[parent].parent == -1:
            return length, parent
        curr = parent


def plan_prune_short_terminal_branches(nodes: Dict[int, Node], threshold: float, tip_mode: str) -> dict:
    children = build_children(nodes)
    leaf_ids = leaves(nodes, children)
    twig_infos = []
    nodes_to_delete: Set[int] = set()
    leaves_flagged: Set[int] = set()

    for leaf in leaf_ids:
        if nodes[leaf].parent == -1:
            continue
        parent = nodes[leaf].parent
        if tip_mode == "segment":
            seglen = euclid(nodes[leaf], nodes[parent])
            if seglen < threshold - EPS:
                leaves_flagged.add(leaf)
                nodes_to_delete.add(leaf)
                twig_infos.append(
                    {"leaf": leaf, "measure": "segment", "value": seglen, "cut_parent": parent, "nodes_on_path": [leaf]}
                )
        else:
            twiglen, cut_parent = twig_length_to_leaf(nodes, children, leaf)
            if twiglen < threshold - EPS and cut_parent in nodes:
                path = []
                curr = leaf
                while nodes[curr].parent != cut_parent:
                    path.append(curr)
                    curr = nodes[curr].parent
                    if curr not in nodes or curr == -1:
                        break
                if curr in nodes:
                    path.append(curr)
                    for nid in path:
                        nodes_to_delete.add(nid)
                    leaves_flagged.add(leaf)
                    twig_infos.append(
                        {
                            "leaf": leaf,
                            "measure": "twig",
                            "value": twiglen,
                            "cut_parent": cut_parent,
                            "nodes_on_path": list(reversed(path)),
                        }
                    )

    return {
        "twig_count": len(leaves_flagged),
        "node_count": len(nodes_to_delete),
        "leaves": sorted(leaves_flagged),
        "nodes": sorted(nodes_to_delete),
        "examples": twig_infos[:20],
    }


def prune_short_terminal_branches(nodes: Dict[int, Node], threshold: float, tip_mode: str, edits_log: Optional[List[str]] = None) -> int:
    removed = 0
    while True:
        children = build_children(nodes)
        leaf_ids = leaves(nodes, children)
        to_delete: Set[int] = set()

        for leaf in leaf_ids:
            if nodes[leaf].parent == -1:
                continue
            parent = nodes[leaf].parent
            if tip_mode == "segment":
                if euclid(nodes[leaf], nodes[parent]) < threshold - EPS:
                    to_delete.add(leaf)
            else:
                twiglen, cut_parent = twig_length_to_leaf(nodes, children, leaf)
                if twiglen < threshold - EPS and cut_parent in nodes:
                    curr = leaf
                    while nodes[curr].parent != cut_parent:
                        curr = nodes[curr].parent
                        if curr not in nodes or curr == -1:
                            break
                    if curr in nodes:
                        path: Set[int] = set()
                        u = leaf
                        while u != curr:
                            path.add(u)
                            u = nodes[u].parent
                        path.add(curr)
                        to_delete |= path

        if not to_delete:
            break

        current_children = build_children(nodes)
        for nid in to_delete:
            for ch in current_children.get(nid, []):
                if ch in nodes:
                    nodes[ch].parent = -1
            nodes.pop(nid, None)

        removed += len(to_delete)
        if edits_log is not None:
            edits_log.append(f"tip pruning: removed {len(to_delete)} nodes (<{threshold} µm, mode={tip_mode})")

    return removed


# ----------------------- Over-branch repair -----------------------


def plan_overbranching(nodes: Dict[int, Node], keep_strategy: str, allow_cross_type: bool) -> dict:
    children = build_children(nodes)
    over_nodes = []
    extras_total = 0
    extras_movable = 0
    unresolved_nodes: Set[int] = set()
    details = []

    for bn, ch in children.items():
        if nodes[bn].ntype == 1:
            continue
        if len(ch) <= 2:
            continue

        over_nodes.append(bn)

        if keep_strategy == "longest_subtree":
            scores = [(subtree_length(nodes, children, d), d) for d in ch]
        elif keep_strategy == "longest_segment":
            scores = [(euclid(nodes[d], nodes[bn]), d) for d in ch]
        elif keep_strategy == "largest_radius":
            scores = [(nodes[d].r, d) for d in ch]
        else:
            scores = [(subtree_length(nodes, children, d), d) for d in ch]

        scores.sort(reverse=True)
        keep = {d for _, d in scores[:2]}
        extras = [d for d in ch if d not in keep]
        extras_total += len(extras)

        movable_here = 0
        for d in extras:
            forb = subtree_nodes(nodes, children, d)
            forb.add(bn)
            req_type = None if allow_cross_type else nodes[d].ntype
            cands = [
                nid
                for nid, node in nodes.items()
                if nid not in forb
                and (degree(children, nid) < 2 or node.ntype == 1)
                and (req_type is None or node.ntype == req_type or node.ntype == 1)
            ]
            if cands:
                movable_here += 1

        extras_movable += movable_here
        if movable_here < len(extras):
            unresolved_nodes.add(bn)

        details.append({"branch_node": bn, "extras": extras, "movable": movable_here})

    return {
        "over_nodes": sorted(over_nodes),
        "over_count": len(over_nodes),
        "extras_total": extras_total,
        "extras_movable": extras_movable,
        "unresolved_nodes": sorted(unresolved_nodes),
        "details": details[:50],
    }


def fix_overbranching(nodes: Dict[int, Node], keep_strategy: str, allow_cross_type: bool, edits_log: Optional[List[str]] = None) -> Tuple[int, int, List[int]]:
    total_considered = 0
    total_moves = 0
    unresolved: Set[int] = set()

    while True:
        children = build_children(nodes)
        over_list = [nid for nid, ch in children.items() if nodes[nid].ntype != 1 and len(ch) > 2]
        if not over_list:
            break

        moved = 0
        for bn in over_list:
            children = build_children(nodes)
            dtrs = list(children[bn])
            if len(dtrs) <= 2:
                continue

            total_considered += 1

            if keep_strategy == "longest_subtree":
                scores = [(subtree_length(nodes, children, d), d) for d in dtrs]
            elif keep_strategy == "longest_segment":
                scores = [(euclid(nodes[d], nodes[bn]), d) for d in dtrs]
            elif keep_strategy == "largest_radius":
                scores = [(nodes[d].r, d) for d in dtrs]
            else:
                scores = [(subtree_length(nodes, children, d), d) for d in dtrs]

            scores.sort(reverse=True)
            keep = {d for _, d in scores[:2]}
            extras = [d for d in dtrs if d not in keep]

            if edits_log is not None:
                edits_log.append(f"over-branch at {bn}: kept {sorted(keep)}, reattach {sorted(extras)}")

            for d in extras:
                children = build_children(nodes)
                forb = subtree_nodes(nodes, children, d)
                forb.add(bn)
                req_type = None if allow_cross_type else nodes[d].ntype
                cands = [
                    nid
                    for nid, node in nodes.items()
                    if nid not in forb
                    and (degree(children, nid) < 2 or node.ntype == 1)
                    and (req_type is None or node.ntype == req_type or node.ntype == 1)
                ]
                if not cands:
                    unresolved.add(bn)
                    continue
                cands.sort(key=lambda nid: euclid(nodes[d], nodes[nid]))
                nodes[d].parent = cands[0]
                total_moves += 1
                moved += 1

        if moved == 0:
            break

    if edits_log is not None:
        edits_log.append(
            f"over-branch repair: considered={total_considered}, reattachments={total_moves}, unresolved={sorted(unresolved)}"
        )
    return total_considered, total_moves, sorted(unresolved)


# ----------------------------- Checks -----------------------------


def check_sanity(nodes: Dict[int, Node]) -> dict:
    missing_parents = [nid for nid, n in nodes.items() if n.parent != -1 and n.parent not in nodes]
    self_parent = [nid for nid, n in nodes.items() if n.parent == n.nid]
    zedges = zero_length_edges(nodes)
    cycles = find_parent_cycles(nodes)
    comps = [sorted(c) for c in overall_components(nodes)]
    return {
        "missing_parent_nodes": missing_parents,
        "self_parent_nodes": self_parent,
        "parent_cycles": cycles,
        "zero_length_edges": zedges,
        "overall_components": comps,
        "num_components": len(comps),
    }


# -------------------------- Orchestration --------------------------


def analyze_only(path: str, args) -> dict:
    raw_nodes, _ = read_swc(path)
    sanity = check_sanity(raw_nodes)

    nodes = clone_nodes(raw_nodes)
    structure_preview = {
        "missing_parent_fixed": fix_missing_parents(nodes),
        "self_parent_fixed": fix_self_parent(nodes),
        "parent_cycles_fixed": fix_parent_cycles(nodes),
        "zero_length_collapsed": collapse_zero_length_edges(nodes),
    }

    ob = plan_overbranching(nodes, args.keep_strategy, args.allow_cross_type)
    tp = plan_prune_short_terminal_branches(nodes, args.tip_threshold, args.tip_mode)

    summary = {
        "overbranch_nodes": ob["over_count"],
        "extras_movable": ob["extras_movable"],
        "extras_total": ob["extras_total"],
        "short_twigs": tp["twig_count"],
        "short_twig_nodes": tp["node_count"],
        "missing_parents": len(sanity["missing_parent_nodes"]),
        "parent_cycles": len(sanity["parent_cycles"]),
        "zero_len_edges": len(sanity["zero_length_edges"]),
        "overall_comps": sanity["num_components"],
    }

    return {
        "summary": summary,
        "details": {
            "structure_preview_fixes": structure_preview,
            "overbranch": ob,
            "short_twigs": tp,
            "sanity": sanity,
        },
        "config": {
            "keep_strategy": args.keep_strategy,
            "allow_cross_type_reattach": args.allow_cross_type,
            "tip_threshold": args.tip_threshold,
            "tip_mode": args.tip_mode,
        },
    }


def apply_fixes(path: str, out_path: Optional[str], args) -> dict:
    nodes, comments = read_swc(path)
    edits_log: List[str] = []

    missing_fixed = fix_missing_parents(nodes, edits_log)
    self_fixed = fix_self_parent(nodes, edits_log)
    cycle_fixed = fix_parent_cycles(nodes, edits_log)
    zero_fixed = collapse_zero_length_edges(nodes, edits_log)

    over_considered = over_moves = 0
    over_unresolved: List[int] = []
    if not args.skip_overbranch:
        over_considered, over_moves, over_unresolved = fix_overbranching(
            nodes, args.keep_strategy, args.allow_cross_type, edits_log
        )

    pruned = 0
    if not args.skip_tip_prune and args.tip_threshold > 0:
        pruned = prune_short_terminal_branches(nodes, args.tip_threshold, args.tip_mode, edits_log)

    if args.second_overbranch and not args.skip_overbranch:
        _c2, m2, u2 = fix_overbranching(nodes, args.keep_strategy, args.allow_cross_type, edits_log)
        over_moves += m2
        over_unresolved = sorted(set(over_unresolved) | set(u2))

    if args.reindex:
        nodes = reindex_bfs(nodes, edits_log=edits_log)

    sanity = check_sanity(nodes)
    summary = {
        "status": "changes_made"
        if any([missing_fixed, self_fixed, cycle_fixed, zero_fixed, over_moves, pruned])
        else "no_changes",
        "counters": {
            "missing_parent_fixed": missing_fixed,
            "self_parent_fixed": self_fixed,
            "parent_cycles_fixed": cycle_fixed,
            "zero_length_collapsed": zero_fixed,
            "overbranch_nodes_considered": over_considered,
            "overbranch_reattachments": over_moves,
            "overbranch_unresolved_nodes": over_unresolved,
            "pruned_tip_nodes": pruned,
            "final_num_nodes": len(nodes),
        },
    }
    report = {
        "summary": summary,
        "sanity": sanity,
        "edits": edits_log,
    }

    if out_path:
        write_swc(out_path, nodes, comments, edits_log)
        base = os.path.splitext(out_path)[0]
        rep_json = base + ".branch_cleanup_report.json"
        with open(rep_json, "w", encoding="utf-8") as f:
            json.dump({path: report}, f, indent=2)
        rep_csv = base + ".branch_cleanup_report.csv"
        write_apply_csv({path: report}, rep_csv)
        print(f"Wrote fixed SWC: {out_path}\nWrote cleanup report: {rep_json}\nWrote CSV summary: {rep_csv}")

    return report


# ---------------------------- Utilities ----------------------------


def reindex_bfs(nodes: Dict[int, Node], edits_log: Optional[List[str]] = None) -> Dict[int, Node]:
    roots = [nid for nid, n in nodes.items() if n.parent == -1 or n.parent not in nodes]
    start_id = min(roots) if roots else min(nodes)
    children = build_children(nodes)

    order: List[int] = []
    seen: Set[int] = set()
    dq = deque([start_id])
    seen.add(start_id)

    while dq:
        u = dq.popleft()
        order.append(u)
        for v in children.get(u, []):
            if v not in seen:
                seen.add(v)
                dq.append(v)

    for nid in sorted(nodes):
        if nid not in seen:
            order.append(nid)
            seen.add(nid)

    idmap = {old: i + 1 for i, old in enumerate(order)}
    new_nodes: Dict[int, Node] = {}
    for old in order:
        n = nodes[old]
        new_id = idmap[old]
        new_parent = -1 if n.parent == -1 else idmap.get(n.parent, -1)
        new_nodes[new_id] = Node(new_id, n.ntype, n.x, n.y, n.z, n.r, new_parent)

    if edits_log is not None:
        edits_log.append(f"reindexed to BFS IDs 1..{len(order)}")
    return new_nodes


# ---------------------------- Renderers ----------------------------


def render_table(all_reports: Dict[str, dict]) -> str:
    cols = [
        "File",
        "OverBranches",
        "Extras(movable/total)",
        "ShortTwigs(twigs/nodes)",
        "MissingParents",
        "ParentCycles",
        "ZeroLenEdges",
        "OverallComps",
    ]
    rows = []
    for fp, rep in all_reports.items():
        s = rep["summary"]
        rows.append(
            [
                os.path.basename(fp),
                str(s["overbranch_nodes"]),
                f"{s['extras_movable']}/{s['extras_total']}",
                f"{s['short_twigs']}/{s['short_twig_nodes']}",
                str(s["missing_parents"]),
                str(s["parent_cycles"]),
                str(s["zero_len_edges"]),
                str(s["overall_comps"]),
            ]
        )
    widths = [max(len(cols[i]), max((len(r[i]) for r in rows), default=0)) for i in range(len(cols))]
    def fmt_row(r: List[str]) -> str:
        return " | ".join(cell.ljust(widths[i]) for i, cell in enumerate(r))
    sep = "-+-".join("-" * w for w in widths)
    return "\n".join([fmt_row(cols), sep] + [fmt_row(r) for r in rows])


def write_preview_csv(reports: Dict[str, dict], path: str) -> None:
    cols = [
        "File",
        "OverBranches",
        "ExtrasMovable",
        "ExtrasTotal",
        "ShortTwigs",
        "ShortTwigNodes",
        "MissingParents",
        "ParentCycles",
        "ZeroLenEdges",
        "OverallComps",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for fp, rep in reports.items():
            s = rep["summary"]
            w.writerow(
                [
                    os.path.basename(fp),
                    s["overbranch_nodes"],
                    s["extras_movable"],
                    s["extras_total"],
                    s["short_twigs"],
                    s["short_twig_nodes"],
                    s["missing_parents"],
                    s["parent_cycles"],
                    s["zero_len_edges"],
                    s["overall_comps"],
                ]
            )


def write_apply_csv(reports: Dict[str, dict], path: str) -> None:
    cols = [
        "File",
        "Status",
        "MissingParentFixed",
        "SelfParentFixed",
        "ParentCyclesFixed",
        "ZeroLenCollapsed",
        "OverBranchConsidered",
        "OverBranchReattached",
        "OverBranchUnresolvedCount",
        "PrunedTipNodes",
        "FinalNumNodes",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for fp, rep in reports.items():
            s = rep["summary"]
            c = s["counters"]
            unresolved_count = len(c.get("overbranch_unresolved_nodes", []))
            w.writerow(
                [
                    os.path.basename(fp),
                    s["status"],
                    c.get("missing_parent_fixed", 0),
                    c.get("self_parent_fixed", 0),
                    c.get("parent_cycles_fixed", 0),
                    c.get("zero_length_collapsed", 0),
                    c.get("overbranch_nodes_considered", 0),
                    c.get("overbranch_reattachments", 0),
                    unresolved_count,
                    c.get("pruned_tip_nodes", 0),
                    c.get("final_num_nodes", 0),
                ]
            )


# ------------------------------- CLI -------------------------------


def infer_output_path(in_path: str, out_dir: Optional[str] = None) -> str:
    base = os.path.basename(in_path)
    root, ext = os.path.splitext(base)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        return os.path.join(out_dir, f"{root}{ext or '.swc'}")
    return os.path.join(os.path.dirname(in_path), f"{root}.branch_clean{ext or '.swc'}")


def run_cli() -> None:
    p = argparse.ArgumentParser(description="Standalone over-branch repair and tip pruning for SWC files.")
    p.add_argument("input", help="Path to one .swc file or a directory of .swc files.")
    p.add_argument("-o", "--output", help="Output .swc path in single-file apply mode.")
    p.add_argument("--out-dir", help="Output folder. Keeps the original filename.")

    p.add_argument("--dry-run", action="store_true", help="Preview only.")
    p.add_argument("--report-format", choices=["json", "table", "csv"], default="json")

    p.add_argument("--skip-overbranch", action="store_true", help="Skip over-branch repair.")
    p.add_argument("--skip-tip-prune", action="store_true", help="Skip tip pruning.")
    p.add_argument("--keep-strategy", choices=["longest_subtree", "longest_segment", "largest_radius"], default="longest_subtree")
    p.add_argument(
        "--allow-cross-type-reattach",
        dest="allow_cross_type",
        action="store_true",
        help="Allow reattachment across types. Default is off.",
    )
    p.add_argument("--second-overbranch", action="store_true", help="Run a second over-branch pass after pruning.")
    p.add_argument("--tip-threshold", type=float, default=2.0, help="Threshold in µm for removing terminal branches.")
    p.add_argument("--tip-mode", choices=["twig", "segment"], default="twig")
    p.add_argument("--reindex", action="store_true", help="Reindex nodes to BFS order after fixes.")

    args = p.parse_args()
    in_path = args.input

    if not os.path.exists(in_path):
        print(json.dumps({"error": "input path does not exist", "input": in_path}, indent=2))
        sys.exit(1)

    if args.dry_run:
        reports: Dict[str, dict] = {}
        file_list = (
            [os.path.join(in_path, f) for f in os.listdir(in_path) if f.lower().endswith(".swc")]
            if os.path.isdir(in_path)
            else [in_path]
        )
        if not file_list:
            print(json.dumps({"error": "no .swc files found in directory", "directory": in_path}, indent=2))
            sys.exit(1)

        for fp in file_list:
            reports[fp] = analyze_only(fp, args)

        if args.report_format == "json":
            print(json.dumps(reports, indent=2))
        elif args.report_format == "table":
            print(render_table(reports))
        else:
            csv_path = os.path.join(in_path, "branch_cleanup_report.csv") if os.path.isdir(in_path) else os.path.splitext(in_path)[0] + ".branch_cleanup_report.csv"
            write_preview_csv(reports, csv_path)
            print(f"Wrote {csv_path}")
        return

    if os.path.isdir(in_path):
        swc_files = [os.path.join(in_path, f) for f in os.listdir(in_path) if f.lower().endswith(".swc")]
        if not swc_files:
            print(json.dumps({"error": "no .swc files found in directory", "directory": in_path}, indent=2))
            sys.exit(1)
        reports_apply: Dict[str, dict] = {}
        for fp in swc_files:
            outp = infer_output_path(fp, out_dir=args.out_dir)
            reports_apply[fp] = apply_fixes(fp, outp, args)
        agg_dir = args.out_dir if args.out_dir else in_path
        agg_csv = os.path.join(agg_dir, "branch_cleanup.aggregate.csv")
        write_apply_csv(reports_apply, agg_csv)
        print(f"Wrote aggregate CSV: {agg_csv}")
    else:
        outp = args.output if args.output else infer_output_path(in_path, out_dir=args.out_dir)
        apply_fixes(in_path, outp, args)


if __name__ == "__main__":
    run_cli()
