// Strength-4 CAOA search with BANDWIDTH = 1 checks (2-symbols).
// Scans n = 8..20, brute-force partitions with 1 fixed in V0.
// Writes results to "results_strength4_bw1.txt".
#include <bits/stdc++.h>
using namespace std;
using Mat = vector<vector<int>>;
using namespace std::chrono;

////////////////////////////////////////////////////////////////////////////////
// Difference matrices and order-statistics (first/second/third) -- residues
// considered are 1..n-1 (we keep 0 in matrix but never count r==0).
////////////////////////////////////////////////////////////////////////////////

vector<vector<vector<vector<int>>>> generate_difference_matrix(const vector<int> &V0, const vector<int> &V1, int n)
{
    vector<vector<vector<vector<int>>>> diff(2, vector<vector<vector<int>>>(2));
    vector<vector<int>> P = {V0, V1};
    for (int a = 0; a < 2; a++)
    {
        for (int b = 0; b < 2; b++)
        {
            vector<vector<int>> sub(P[a].size(), vector<int>(P[b].size()));
            for (size_t i = 0; i < P[a].size(); ++i)
            {
                for (size_t j = 0; j < P[b].size(); ++j)
                {
                    int d = (P[a][i] - P[b][j]) % n;
                    if (d < 0)
                        d += n;
                    // Keep d in 0..n-1, ignore r=0 when computing statistics
                    sub[i][j] = d;
                }
            }
            diff[a][b] = move(sub);
        }
    }
    return diff;
}

// first-order: for r = 1..n-1, count occurrences in each block (a,b)
// Note: We map 0 to n in the diff matrix but only count residues 1..n-1 for strength
vector<Mat> first_order(const vector<vector<vector<vector<int>>>> &diff, int n)
{
    int s = 2;
    int maxR = n - 1; // Count only non-zero residues 1..n-1
    vector<Mat> res(maxR + 1, Mat(s, vector<int>(s, 0)));
    for (int r = 1; r <= maxR; ++r)
    {
        for (int a = 0; a < s; ++a)
        {
            for (int b = 0; b < s; ++b)
            {
                int cnt = 0;
                const auto &sub = diff[a][b];
                for (size_t i = 0; i < sub.size(); ++i)
                    for (size_t j = 0; j < sub[i].size(); ++j)
                        if (sub[i][j] == r)
                            ++cnt;
                res[r][a][b] = cnt;
            }
        }
    }
    return res;
}

// compute second-order matrix for a single pair (r1,r2) (rows that have both residues)
Mat second_onepair(const vector<vector<vector<vector<int>>>> &diff, int r1, int r2)
{
    int s = 2;
    Mat M(s, vector<int>(s, 0));
    for (int a = 0; a < s; ++a)
    {
        for (int b = 0; b < s; ++b)
        {
            int rowsWithBoth = 0;
            const auto &sub = diff[a][b];
            for (size_t i = 0; i < sub.size(); ++i)
            {
                bool h1 = false, h2 = false;
                for (size_t j = 0; j < sub[i].size(); ++j)
                {
                    if (sub[i][j] == r1)
                        h1 = true;
                    if (sub[i][j] == r2)
                        h2 = true;
                }
                if (h1 && h2)
                    ++rowsWithBoth;
            }
            M[a][b] = rowsWithBoth;
        }
    }
    return M;
}

// compute third-order matrix for a single triple (r1,r2,r3)
Mat third_onetriple(const vector<vector<vector<vector<int>>>> &diff, int r1, int r2, int r3)
{
    int s = 2;
    Mat M(s, vector<int>(s, 0));
    for (int a = 0; a < s; ++a)
    {
        for (int b = 0; b < s; ++b)
        {
            int rowsWithAll = 0;
            const auto &sub = diff[a][b];
            for (size_t i = 0; i < sub.size(); ++i)
            {
                bool h1 = false, h2 = false, h3 = false;
                for (size_t j = 0; j < sub[i].size(); ++j)
                {
                    if (sub[i][j] == r1)
                        h1 = true;
                    if (sub[i][j] == r2)
                        h2 = true;
                    if (sub[i][j] == r3)
                        h3 = true;
                }
                if (h1 && h2 && h3)
                    ++rowsWithAll;
            }
            M[a][b] = rowsWithAll;
        }
    }
    return M;
}

bool mat_within(const Mat &A, const Mat &B, int tol)
{
    if (A.size() != B.size())
        return false;
    for (size_t i = 0; i < A.size(); ++i)
    {
        if (A[i].size() != B[i].size())
            return false;
        for (size_t j = 0; j < A[i].size(); ++j)
        {
            if (abs(A[i][j] - B[i][j]) > tol)
                return false;
        }
    }
    return true;
}

// Check if matrix is approximately uniform (all elements within tol of each other)
// This is crucial for CAOA - we need uniform distribution, not arbitrary patterns
// For true CAOA: tol should be 0 (exact uniformity within each matrix)
// Bandwidth tolerance applies to COMPARING different residues' uniform values
bool mat_is_uniform(const Mat &M, int tol)
{
    if (M.empty() || M[0].empty())
        return true;
    int first_val = M[0][0];
    for (size_t i = 0; i < M.size(); ++i)
    {
        for (size_t j = 0; j < M[i].size(); ++j)
        {
            if (abs(M[i][j] - first_val) > tol)
                return false;
        }
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////
// prefix p: largest r such that first[r] is within tol of first[1]
// Following paper's definition of Λ and ID (Identity Distance)
// Per paper: compare matrices to reference, NO internal uniformity requirement
////////////////////////////////////////////////////////////////////////////////
int prefix_equal_first(const vector<Mat> &first, int n, int tol)
{
    int maxR = n - 1; // Only consider non-zero residues (1..n-1)
    int p = 1;
    for (int r = 2; r <= maxR; ++r)
    {
        // Compare first[r] to first[1] with bandwidth tolerance
        // Paper's ID definition: matrices must match element-wise within tol
        if (!mat_within(first[r], first[1], tol))
            break;
        p = r;
    }
    return p;
} ////////////////////////////////////////////////////////////////////////////////
// binary search for max t such that all second/triple among 1..t match
// the reference pair {1,2} and triple {1,2,3} within tol.
////////////////////////////////////////////////////////////////////////////////
bool second_all_equal_up_to_ref(const vector<vector<vector<vector<int>>>> &diff, int upto, int tol)
{
    if (upto <= 1)
        return true;
    Mat ref = second_onepair(diff, 1, 2);
    for (int i = 1; i <= upto; ++i)
    {
        for (int j = i + 1; j <= upto; ++j)
        {
            Mat cur = second_onepair(diff, i, j);
            if (!mat_within(cur, ref, tol))
                return false;
        }
    }
    return true;
}
bool third_all_equal_up_to_ref(const vector<vector<vector<vector<int>>>> &diff, int upto, int tol)
{
    if (upto <= 2)
        return true;
    Mat ref = third_onetriple(diff, 1, 2, 3);
    for (int i = 1; i <= upto; ++i)
    {
        for (int j = i + 1; j <= upto; ++j)
        {
            for (int k = j + 1; k <= upto; ++k)
            {
                Mat cur = third_onetriple(diff, i, j, k);
                if (!mat_within(cur, ref, tol))
                    return false;
            }
        }
    }
    return true;
}

int max_t_binary(const vector<vector<vector<vector<int>>>> &diff, int p, int tol)
{
    if (p <= 1)
        return p;
    int lo = 1, hi = p, ans = 1;
    while (lo <= hi)
    {
        int mid = (lo + hi) / 2;
        if (second_all_equal_up_to_ref(diff, mid, tol) && third_all_equal_up_to_ref(diff, mid, tol))
        {
            ans = mid;
            lo = mid + 1;
        }
        else
            hi = mid - 1;
    }
    return ans;
}

////////////////////////////////////////////////////////////////////////////////
// find k for given partition using BANDWIDTH tol (here tol=1)
////////////////////////////////////////////////////////////////////////////////
int find_k_for_partition_strength4(const vector<int> &V0, const vector<int> &V1, int n, int tol)
{
    auto diff = generate_difference_matrix(V0, V1, n);
    auto first = first_order(diff, n);
    int p = prefix_equal_first(first, n, tol);
    int t = max_t_binary(diff, p, tol);
    int k = t + 1;
    // k cannot exceed n in a valid CAOA
    if (k > n)
        k = n;
    return k;
}

////////////////////////////////////////////////////////////////////////////////
// Theorem-2 expected values (as in paper); for BANDWIDTH=1 we allow +/-1
// Returns lambda1, lambda2, lambda3_low, lambda3_high and an invalid flag
////////////////////////////////////////////////////////////////////////////////
struct Expected
{
    bool invalid_by_theorem;
    long lambda1;
    long lambda2;
    long lambda3_low;
    long lambda3_high;
};
Expected expected_theorem2(int n)
{
    Expected e;
    e.invalid_by_theorem = false;
    int r = n % 16;
    auto fl = [&](int d) -> long
    { return n / d; };
    auto ce = [&](int d) -> long
    { return (n + d - 1) / d; };
    if (r == 0)
    {
        e.lambda1 = fl(4);
        e.lambda2 = fl(8);
        e.lambda3_low = e.lambda3_high = fl(16);
    }
    else if (r == 1 || r == 6 || r == 7 || r == 9 || r == 14 || r == 15)
    {
        e.lambda1 = fl(4);
        e.lambda2 = fl(8);
        e.lambda3_low = e.lambda3_high = fl(16);
    }
    else if (r == 2 || r == 10)
    {
        e.lambda1 = ce(4);
        e.lambda2 = ce(8);
        e.lambda3_low = e.lambda3_high = ce(16);
    }
    else if (r == 5 || r == 11)
    {
        e.lambda1 = fl(4);
        e.lambda2 = ce(8);
        e.lambda3_low = e.lambda3_high = ce(16);
    }
    else if (r == 8)
    {
        e.lambda1 = fl(4);
        e.lambda2 = fl(8);
        e.lambda3_low = fl(16);
        e.lambda3_high = ce(16);
    }
    else
    { // r in {3,4,12,13}: theorem excludes existence
        e.invalid_by_theorem = true;
        e.lambda1 = e.lambda2 = e.lambda3_low = e.lambda3_high = -1;
    }
    return e;
}

// Check Theorem-2 modular freq conditions on (0,0) block, using BANDWIDTH tol (allow +/- tol)
bool check_modular_conditions_bandwidth1(const vector<int> &V0, const vector<int> &V1, int n, int k, const Expected &expv, int tol)
{
    if (expv.invalid_by_theorem)
        return false;
    if (k <= 1)
        return false;
    if (k > n)
        return false; // k cannot exceed n
    auto diff = generate_difference_matrix(V0, V1, n);
    auto first = first_order(diff, n);
    int maxR = min(k - 1, n - 1); // Ensure we don't exceed valid residue range (1..n-1)
    // first-order: check first[r][0][0] within +/-tol of lambda1
    for (int r = 1; r <= maxR; ++r)
    {
        int val = first[r][0][0];
        if (abs(val - (int)expv.lambda1) > tol)
            return false;
    }
    // second-order: check all pairs within tolerance
    if (maxR >= 2)
    {
        for (int i = 1; i <= maxR; ++i)
        {
            for (int j = i + 1; j <= maxR; ++j)
            {
                Mat sec = second_onepair(diff, i, j);
                int val = sec[0][0];
                if (abs(val - (int)expv.lambda2) > tol)
                    return false;
            }
        }
    }
    // third-order: allow within [lambda3_low - tol , lambda3_high + tol]
    if (maxR >= 3)
    {
        int low = (int)expv.lambda3_low - tol;
        int high = (int)expv.lambda3_high + tol;
        for (int i = 1; i <= maxR; ++i)
        {
            for (int j = i + 1; j <= maxR; ++j)
            {
                for (int l = j + 1; l <= maxR; ++l)
                {
                    Mat th = third_onetriple(diff, i, j, l);
                    int val = th[0][0];
                    if (val < low || val > high)
                        return false;
                }
            }
        }
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////
// Main: run n=8..20, brute force partitions, write results
////////////////////////////////////////////////////////////////////////////////

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int BANDWIDTH_TOL = 1; // bandwidth = 1 comparisons
    ofstream fout("results_strength4_bw1.txt");
    if (!fout)
    {
        cerr << "Cannot open results_strength4_bw1.txt\n";
        return 1;
    }

    for (int n = 8; n <= 18; ++n)
    {
        auto t0 = high_resolution_clock::now();

        int sizeV0 = (n % 2 == 0) ? (n / 2) : ((n - 1) / 2);
        int choose = sizeV0 - 1; // fixing element 1 in V0
        vector<int> elems;
        for (int i = 2; i <= n; ++i)
            elems.push_back(i);
        int m = (int)elems.size();
        if (choose < 0 || choose > m)
        {
            fout << "n=" << n << ": invalid choose\n";
            continue;
        }

        vector<int> bit(m, 0);
        for (int i = 0; i < choose; ++i)
            bit[i] = 1;
        sort(bit.begin(), bit.end(), greater<int>());

        int bestK = 1;
        vector<int> bestV0, bestV1;

        bool found_valid_by_theorem = false;
        vector<int> validV0, validV1;
        int validK = -1;

        Expected expv = expected_theorem2(n);

        // iterate partitions (1 fixed in V0)
        do
        {
            vector<int> V0;
            V0.push_back(1);
            for (int i = 0; i < m; ++i)
                if (bit[i])
                    V0.push_back(elems[i]);
            vector<int> mark(n + 1, 0);
            for (int x : V0)
                mark[x] = 1;
            vector<int> V1;
            for (int i = 1; i <= n; ++i)
                if (!mark[i])
                    V1.push_back(i);

            int k = find_k_for_partition_strength4(V0, V1, n, BANDWIDTH_TOL);
            if (k > bestK)
            {
                bestK = k;
                bestV0 = V0;
                bestV1 = V1;
            }

            if (!expv.invalid_by_theorem && k >= 4 && !found_valid_by_theorem)
            {
                if (check_modular_conditions_bandwidth1(V0, V1, n, k, expv, BANDWIDTH_TOL))
                {
                    found_valid_by_theorem = true;
                    validV0 = V0;
                    validV1 = V1;
                    validK = k;
                    // keep searching to find possibly higher k but we already have a theorem-valid partition
                }
            }
        } while (prev_permutation(bit.begin(), bit.end()));

        auto t1 = high_resolution_clock::now();
        auto elapsed = duration_cast<milliseconds>(t1 - t0).count();

        fout << "n = " << n << "\n";
        fout << "Time(ms): " << elapsed << "\n";
        fout << "Best k (by ID logic, bandwidth tol=" << BANDWIDTH_TOL << ") = " << bestK << "\n";
        fout << "Best V0: ";
        for (int x : bestV0)
            fout << x << " ";
        fout << "\nBest V1: ";
        for (int x : bestV1)
            fout << x << " ";
        fout << "\n";

        if (expv.invalid_by_theorem)
        {
            fout << "Theorem-2: n mod 16 in {3,4,12,13} -> theorem excludes existence (skipped checks)\n";
        }
        else
        {
            fout << "Theorem-2 expected (lambda1, lambda2, lambda3_range): "
                 << expv.lambda1 << ", " << expv.lambda2 << ", [" << expv.lambda3_low << "," << expv.lambda3_high << "]\n";
            fout << "(Using bandwidth tolerance = " << BANDWIDTH_TOL << " for comparisons)\n";
            if (found_valid_by_theorem)
            {
                fout << "Found partition satisfying Theorem-2 modular checks (bandwidth=1):\n";
                fout << "k = " << validK << "\nV0: ";
                for (int x : validV0)
                    fout << x << " ";
                fout << "\nV1: ";
                for (int x : validV1)
                    fout << x << " ";
                fout << "\n";
            }
            else
            {
                fout << "No partition satisfying Theorem-2 modular checks found for this n (with bandwidth tol=" << BANDWIDTH_TOL << ").\n";
            }
        }
        fout << "----------------------------------------\n";
        fout.flush();
    }

    fout.close();
    return 0;
}
