#!/usr/bin/env python3
# genvec_and_incidence.py
# Find a generating vector (binary) for strength-3 CAOA and build k x n circulant incidence matrix.
# - Uses reduced search over triplet frequencies (floor/ceil n/8)
# - Constructs Eulerian circuit on De Bruijn frequency graph
# - Outputs generating vector and k x n incidence matrix (rows = cyclic shifts)
# - Optional verification of strength-3 on all 3-row choices (costly for large k)

from itertools import combinations_with_replacement, combinations
from collections import deque
import math
import sys

TRIPLETS = ['000','001','010','011','100','101','110','111']
VERTS = ['00','01','10','11']

def all_compositions(rem, bins):
    if rem == 0:
        yield [0]*bins
        return
    # stars-and-bars via combinations_with_replacement
    from itertools import combinations_with_replacement
    for comb in combinations_with_replacement(range(bins), rem):
        arr = [0]*bins
        for idx in comb: arr[idx]+=1
        yield arr

def check_eulerian_balance(counts):
    out_deg = {v:0 for v in VERTS}
    in_deg  = {v:0 for v in VERTS}
    for trip,c in counts.items():
        u = trip[:2]; v = trip[1:]
        out_deg[u]+=c; in_deg[v]+=c
    return all(out_deg[v]==in_deg[v] for v in VERTS)

def build_multigraph(counts):
    adj = {v: deque() for v in VERTS}
    # deterministic ordering: TRIPLETS order
    for trip in TRIPLETS:
        c = counts[trip]
        if c<=0: continue
        u = trip[:2]; v = trip[1:]
        for _ in range(c):
            adj[u].append((v, trip))
    return adj

def hierholzer(adj):
    # copy adjacency to avoid mutating original
    adjc = {v: deque(adj[v]) for v in adj}
    start = next((v for v in VERTS if adjc[v]), None)
    if start is None: return None
    stack=[start]; edge_stack=[]; path=[]
    while stack:
        v = stack[-1]
        if adjc[v]:
            to,label = adjc[v].popleft()
            stack.append(to); edge_stack.append(label)
        else:
            stack.pop()
            if edge_stack:
                lab = edge_stack.pop()
                path.append(lab)
    path.reverse()
    return path

def make_generating_vector_from_triplet_path(path):
    if not path: return ''
    seq = [path[0][0], path[0][1]]
    for t in path: seq.append(t[2])
    return ''.join(seq)  # length edges + 2

def circular_triplet_counts_from_g(g):
    n = len(g)
    counts = {t:0 for t in TRIPLETS}
    for i in range(n):
        a = g[i % n]; b = g[(i+1) % n]; c = g[(i+2) % n]
        counts[a+b+c] += 1
    return counts

def is_uniform_counts(n, counts):
    vals = list(counts.values())
    return sum(vals)==n and (max(vals)-min(vals) <= 1)

def find_generating_vector(n, require_all_positive=False, verbose=False):
    if n < 8:
        if verbose: print("n must be >= 8 for t=3 designs.")
        return None
    base = n // 8
    rem = n - base*8
    if require_all_positive and base == 0:
        if verbose: print("Impossible: cannot force every triplet >=1 when n < 8.")
        return None
    if verbose:
        print(f"[search] n={n}, base={base}, remainder={rem}")

    # enumerate ways to distribute rem across 8 triplets
    for extra in all_compositions(rem, 8):
        counts = {TRIPLETS[i]: base + extra[i] for i in range(8)}
        if require_all_positive and any(v<1 for v in counts.values()):
            continue
        if sum(counts.values()) != n:
            continue
        if not check_eulerian_balance(counts):
            continue
        adj = build_multigraph(counts)
        path = hierholzer(adj)
        if not path: continue
        seq = make_generating_vector_from_triplet_path(path)
        gv = seq[:n]
        obs_counts = circular_triplet_counts_from_g(gv)
        if obs_counts != counts:
            if verbose:
                print("constructed g had mismatch with assigned counts; skipping")
            continue
        if verbose:
            print("Found generating vector with counts:", counts)
        return gv, counts
    if verbose: print("No feasible generating vector found under these uniform constraints.")
    return None, None

def incidence_matrix_from_genvec(g, k):
    n = len(g)
    if k > n:
        # you can have k > n but rows repeat after n shifts; warn the user
        print(f"Warning: k={k} > n={n}. Rows will repeat every n shifts (periodic).")
    A = []
    for i in range(k):
        row = [g[(i+j) % n] for j in range(n)]
        A.append(row)
    return A

def verify_strength3(A, n, verbose=False):
    # Check every combination of 3 rows for triplet counts equal floor/ceil pattern
    from itertools import combinations
    target_floor = n//8
    target_ceil = math.ceil(n/8)
    bad = []
    k = len(A)
    for (i,j,l) in combinations(range(k), 3):
        counts = {t:0 for t in TRIPLETS}
        for col in range(n):
            trip = str(A[i][col]) + str(A[j][col]) + str(A[l][col])
            counts[trip] += 1
        vals = list(counts.values())
        if not (sum(vals)==n and max(vals)-min(vals) <= 1):
            bad.append(((i,j,l), counts))
            if verbose:
                print("Violation for rows", (i,j,l), "counts:", counts)
                break
    return bad

def pretty_print_matrix(A):
    for row in A:
        print(''.join(str(x) for x in row))

if __name__ == "__main__":
    # Option A: set n,k here
    #n = 8; k = 3
    #n = 20; k = 6
    # Option B: simple CLI
    if len(sys.argv) >= 3:
        n = int(sys.argv[1]); k = int(sys.argv[2])
    else:
        # default interactive prompt
        try:
            n = int(input("Enter n (length of generating vector / #columns) >=8: ").strip())
        except:
            print("Invalid n"); sys.exit(1)
        try:
            k = int(input("Enter k (number of rows / factors to produce, k<=n recommended): ").strip())
        except:
            print("Invalid k"); sys.exit(1)

    gv, counts = find_generating_vector(n, require_all_positive=True, verbose=True)
    if not gv:
        print("No generating vector found for n=", n)
        sys.exit(1)
    print("\nGenerating vector g (length n):")
    print(gv)
    print("\nTriplet counts (observed) for g:")
    print(counts)

    # build incidence matrix for given k
    A = incidence_matrix_from_genvec(gv, k)
    print(f"\nIncidence matrix A ({k} x {n}):")
    pretty_print_matrix(A)

    # optional verification (can be slow if k large)
    do_verify = input("\nVerify strength-3 for all 3-row combinations? (y/N): ").strip().lower()
    if do_verify == 'y':
        print("Verifying... (may take time if k large)")
        bad = verify_strength3(A, n, verbose=True)
        if not bad:
            print("All 3-row submatrices pass the uniformity check (counts differ by at most 1).")
        else:
            print("Found violations (first shown above).")
