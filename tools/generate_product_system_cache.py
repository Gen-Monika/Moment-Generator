#!/usr/bin/env python3
"""Generate exact product-system caches for Moment Generator.

The cache stores the matrix A_k(n) in the exact term form

    c * (n)_b * n^p

where c is an arbitrary-size integer, b is the falling-factor length,
(n)_b = n(n-1)...(n-b+1), (n)_0 = 1, and p may be negative. This avoids expensive symbolic simplification while keeping every entry
exact. Qt reads these JSON files for fast UI rendering.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from functools import lru_cache
from math import comb
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Tuple

Counts = Tuple[int, ...]
Partition = Tuple[int, ...]
TermKey = Tuple[int, int]  # (n_power, falling_blocks)


def partitions_without_ones(total: int, max_part: int | None = None) -> Iterable[Partition]:
    if max_part is None or max_part > total:
        max_part = total
    if total == 0:
        yield ()
        return
    for part in range(max_part, 1, -1):
        if part <= total:
            next_max = min(part, total - part) if total - part else 0
            for rest in partitions_without_ones(total - part, next_max):
                yield (part,) + rest


def add_part(partition: Partition, moment_order: int) -> Partition:
    parts = list(partition)
    parts.append(moment_order)
    parts.sort(reverse=True)
    return tuple(parts)


def block_compositions(counts: Counts, anchor_order: int) -> Iterable[Counts]:
    max_order = len(counts) - 1
    block = [0] * (max_order + 1)

    def rec(order: int) -> Iterable[Counts]:
        if order > max_order:
            yield tuple(block)
            return
        lower = 1 if order == anchor_order else 0
        for value in range(lower, counts[order] + 1):
            block[order] = value
            yield from rec(order + 1)
        block[order] = 0

    yield from rec(1)


@lru_cache(maxsize=None)
def expectation_of_power_counts(counts: Counts) -> Dict[Partition, int]:
    """Return coefficients in E[prod_s p_s^counts[s]].

    The coefficient for a partition lambda already includes the number of labeled
    set partitions producing lambda. The falling factor (n)_{len(lambda)} is added
    later, because the UI/cache keeps it symbolic.
    """
    if sum(counts) == 0:
        return {(): 1}

    max_order = len(counts) - 1
    anchor_order = next(order for order in range(1, max_order + 1) if counts[order] > 0)
    out: DefaultDict[Partition, int] = defaultdict(int)

    for block in block_compositions(counts, anchor_order):
        moment_order = sum(order * block[order] for order in range(1, max_order + 1))
        if moment_order == 1:
            # v_1 = 0, so this collision type contributes nothing.
            continue

        ways = comb(counts[anchor_order] - 1, block[anchor_order] - 1)
        for order in range(1, max_order + 1):
            if order != anchor_order:
                ways *= comb(counts[order], block[order])

        remainder = list(counts)
        for order in range(1, max_order + 1):
            remainder[order] -= block[order]

        for partition, coeff in expectation_of_power_counts(tuple(remainder)).items():
            out[add_part(partition, moment_order)] += ways * coeff

    return dict(out)


def single_moment_terms(moment_order: int, total_order: int) -> List[Tuple[Counts, int, int]]:
    terms: List[Tuple[Counts, int, int]] = []
    for mean_power in range(moment_order + 1):
        counts = [0] * (total_order + 1)
        counts[1] += mean_power
        n_power = -(mean_power + 1)
        raw_power = moment_order - mean_power
        if raw_power == 0:
            n_power += 1  # p_0 = n
        else:
            counts[raw_power] += 1
        coeff = (-1) ** mean_power * comb(moment_order, mean_power)
        terms.append((tuple(counts), n_power, coeff))
    return terms


def sample_product_expansion(row: Partition, total_order: int) -> Dict[Tuple[Counts, int], int]:
    zero_counts: Counts = tuple([0] * (total_order + 1))
    poly: Dict[Tuple[Counts, int], int] = {(zero_counts, 0): 1}

    for moment_order in row:
        new_poly: DefaultDict[Tuple[Counts, int], int] = defaultdict(int)
        for (counts, n_power), coeff in poly.items():
            for term_counts, term_n_power, term_coeff in single_moment_terms(moment_order, total_order):
                combined = list(counts)
                for order in range(1, total_order + 1):
                    combined[order] += term_counts[order]
                new_poly[(tuple(combined), n_power + term_n_power)] += coeff * term_coeff
        poly = {key: value for key, value in new_poly.items() if value}

    return poly


def row_coefficients(row: Partition, total_order: int) -> Dict[Partition, Dict[TermKey, int]]:
    out: DefaultDict[Partition, DefaultDict[TermKey, int]] = defaultdict(lambda: defaultdict(int))

    for (counts, n_power), coeff in sample_product_expansion(row, total_order).items():
        for column, ways in expectation_of_power_counts(counts).items():
            out[column][(n_power, len(column))] += coeff * ways

    return {
        column: {key: value for key, value in terms.items() if value}
        for column, terms in out.items()
    }


def generate_order(order: int) -> dict:
    expectation_of_power_counts.cache_clear()
    partitions = list(partitions_without_ones(order))
    index = {partition: idx for idx, partition in enumerate(partitions)}
    entries = []
    term_count = 0

    started = time.time()
    for row_index, row in enumerate(partitions):
        coeffs = row_coefficients(row, order)
        for column, terms in coeffs.items():
            normalized_terms = [
                {
                    "coeff": str(coeff),
                    "n_power": n_power,
                    "falling": falling_blocks,
                }
                for (n_power, falling_blocks), coeff in sorted(terms.items())
            ]
            term_count += len(normalized_terms)
            entries.append({
                "row": row_index,
                "col": index[column],
                "terms": normalized_terms,
            })

    elapsed_ms = int(round((time.time() - started) * 1000))
    return {
        "version": 1,
        "order": order,
        "partitions": [list(partition) for partition in partitions],
        "dimension": len(partitions),
        "max_factors": max(len(partition) for partition in partitions),
        "entry_count": len(entries),
        "term_count": term_count,
        "term_form": "coeff * falling(n, falling) * n^n_power",
        "elapsed_ms": elapsed_ms,
        "entries": entries,
    }


def write_order(order: int, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    data = generate_order(order)
    output_path = output_dir / f"product_system_k{order}.json"
    output_path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(
        f"k={order}: dimension={data['dimension']} entries={data['entry_count']} "
        f"terms={data['term_count']} elapsed={data['elapsed_ms']}ms -> {output_path}"
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Moment Generator product-system caches.")
    parser.add_argument("--order", type=int, help="Generate one order.")
    parser.add_argument("--min-order", type=int, default=2, help="First order for a range generation.")
    parser.add_argument("--max-order", type=int, default=20, help="Last order for a range generation.")
    parser.add_argument("--output-dir", type=Path, default=Path("cache/product_systems"))
    args = parser.parse_args()

    if args.order is not None:
        orders = [args.order]
    else:
        orders = list(range(args.min_order, args.max_order + 1))

    for order in orders:
        if order < 2 or order > 25:
            raise SystemExit(f"order {order} is outside the supported range 2..25")
        write_order(order, args.output_dir)


if __name__ == "__main__":
    main()

