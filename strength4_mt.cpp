#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using Mat = vector<vector<int>>;

struct Expected
{
    bool invalid_by_theorem;
    long lambda1, lambda2, lambda3_low, lambda3_high;
};

Expected expected_theorem2(int n)
{
    Expected e{false, -1, -1, -1, -1};
    int r = n % 16;
    auto fl = [&](int d)
    { return n / d; };
    auto ce = [&](int d)
    { return (n + d - 1) / d; };
    if (r == 0 || r == 1 || r == 6 || r == 7 || r == 9 || r == 14 || r == 15)
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
    {
        e.invalid_by_theorem = true;
    }
    return e;
}

vector<vector<vector<vector<int>>>> generate_difference_matrix(const vector<int>& V0, const vector<int>& V1, int n)
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
                    sub[i][j] = d;
                }
            }
            diff[a][b] = move(sub);
        }
    }
    return diff;
}

vector<Mat> first_order(const vector<vector<vector<vector<int>>>>& diff, int n)
{
    int s = 2;
    int maxR = n - 1;
    vector<Mat> res(maxR + 1, Mat(s, vector<int>(s, 0)));
    for (int r = 1; r <= maxR; ++r)
    {
        for (int a = 0; a < s; ++a)
        {
            for (int b = 0; b < s; ++b)
            {
                int cnt = 0;
                const auto& sub = diff[a][b];
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

Mat second_onepair(const vector<vector<vector<vector<int>>>>& diff, int r1, int r2)
{
    int s = 2;
    Mat M(s, vector<int>(s, 0));
    for (int a = 0; a < s; ++a)
    {
        for (int b = 0; b < s; ++b)
        {
            int rows = 0;
            const auto& sub = diff[a][b];
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
                    ++rows;
            }
            M[a][b] = rows;
        }
    }
    return M;
}

Mat third_onetriple(const vector<vector<vector<vector<int>>>>& diff, int r1, int r2, int r3)
{
    int s = 2;
    Mat M(s, vector<int>(s, 0));
    for (int a = 0; a < s; ++a)
    {
        for (int b = 0; b < s; ++b)
        {
            int rows = 0;
            const auto& sub = diff[a][b];
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
                    ++rows;
            }
            M[a][b] = rows;
        }
    }
    return M;
}

bool mat_within(const Mat& A, const Mat& B, int tol)
{
    if (A.size() != B.size())
        return false;
    for (size_t i = 0; i < A.size(); ++i)
    {
        if (A[i].size() != B[i].size())
            return false;
        for (size_t j = 0; j < A[i].size(); ++j)
            if (abs(A[i][j] - B[i][j]) > tol)
                return false;
    }
    return true;
}

int prefix_equal_first(const vector<Mat>& first, int n, int tol)
{
    int maxR = n - 1;
    int p = 1;
    for (int r = 2; r <= maxR; ++r)
    {
        if (!mat_within(first[r], first[1], tol))
            break;
        p = r;
    }
    return p;
}

bool second_all_equal_up_to_ref(const vector<vector<vector<vector<int>>>>& diff, int upto, int tol)
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

bool third_all_equal_up_to_ref(const vector<vector<vector<vector<int>>>>& diff, int upto, int tol)
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

int max_t_binary(const vector<vector<vector<vector<int>>>>& diff, int p, int tol)
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

int find_k_for_partition_strength4(const vector<int>& V0, const vector<int>& V1, int n, int tol)
{
    auto diff = generate_difference_matrix(V0, V1, n);
    auto first = first_order(diff, n);
    int p = prefix_equal_first(first, n, tol);
    int t = max_t_binary(diff, p, tol);
    return min(n, t + 1);
}

bool check_modular_conditions_bandwidth1(const vector<int>& V0, const vector<int>& V1, int n, int k, const Expected& expv, int tol)
{
    if (expv.invalid_by_theorem || k <= 1 || k > n)
        return false;
    auto diff = generate_difference_matrix(V0, V1, n);
    auto first = first_order(diff, n);
    int maxR = min(k - 1, n - 1);
    for (int r = 1; r <= maxR; ++r)
    {
        int val = first[r][0][0];
        if (abs(val - (int)expv.lambda1) > tol)
            return false;
    }
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

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int BANDWIDTH_TOL = 1;
    ofstream fout("results_strength4_bw1_mt.txt");

    for (int n = 8; n <= 18; ++n)
    {
        int sizeV0 = (n % 2 == 0) ? (n / 2) : ((n - 1) / 2);
        int choose = sizeV0 - 1;
        vector<int> elems;
        for (int i = 2; i <= n; ++i)
            elems.push_back(i);
        int m = elems.size();
        if (choose < 0 || choose > m)
            continue;

        vector<vector<int>> partitions;
        vector<int> bit(m, 0);
        for (int i = 0; i < choose; ++i)
            bit[i] = 1;
        sort(bit.begin(), bit.end(), greater<int>());
        do
        {
            vector<int> V0;
            V0.push_back(1);
            for (int i = 0; i < m; ++i)
                if (bit[i])
                    V0.push_back(elems[i]);
            partitions.push_back(move(V0));
        } while (prev_permutation(bit.begin(), bit.end()));

        int bestK = 1;
        vector<int> bestV0, bestV1;
        bool found_valid = false;
        vector<int> validV0, validV1;
        int validK = -1;
        Expected expv = expected_theorem2(n);
        auto t0 = chrono::high_resolution_clock::now();

        #pragma omp parallel
        {
            int localBest = 1;
            vector<int> localBestV0, localBestV1;
            bool localFound = false;
            vector<int> localValidV0, localValidV1;
            int localValidK = -1;

            #pragma omp for schedule(dynamic, 64)
            for (int idx = 0; idx < (int)partitions.size(); ++idx)
            {
                const auto& V0 = partitions[idx];
                vector<int> mark(n + 1, 0);
                for (int x : V0)
                    mark[x] = 1;
                vector<int> V1;
                for (int i = 1; i <= n; ++i)
                    if (!mark[i])
                        V1.push_back(i);

                int k = find_k_for_partition_strength4(V0, V1, n, BANDWIDTH_TOL);
                if (k > localBest)
                {
                    localBest = k;
                    localBestV0 = V0;
                    localBestV1 = V1;
                }
                if (!localFound && !expv.invalid_by_theorem && k >= 4)
                {
                    if (check_modular_conditions_bandwidth1(V0, V1, n, k, expv, BANDWIDTH_TOL))
                    {
                        localFound = true;
                        localValidV0 = V0;
                        localValidV1 = V1;
                        localValidK = k;
                    }
                }
            }

            #pragma omp critical
            {
                if (localBest > bestK)
                {
                    bestK = localBest;
                    bestV0 = localBestV0;
                    bestV1 = localBestV1;
                }
                if (!found_valid && localFound)
                {
                    found_valid = true;
                    validV0 = localValidV0;
                    validV1 = localValidV1;
                    validK = localValidK;
                }
            }
        }

        auto t1 = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

        fout << "n = " << n << "\n";
        fout << "Time(ms): " << elapsed << "\n";
        fout << "Best k (ID, bw=1) = " << bestK << "\n";
        fout << "Best V0: ";
        for (int x : bestV0)
            fout << x << " ";
        fout << "\nBest V1: ";
        for (int x : bestV1)
            fout << x << " ";
        fout << "\n";
        if (expv.invalid_by_theorem)
        {
            fout << "Theorem-2: excluded by n mod 16\n";
        }
        else
        {
            fout << "Theorem-2 expected (lambda1, lambda2, lambda3_range): "
                 << expv.lambda1 << ", " << expv.lambda2 << ", ["
                 << expv.lambda3_low << "," << expv.lambda3_high << "]\n";
            if (found_valid)
            {
                fout << "Found theorem-valid partition (bw=1): k=" << validK << "\nV0: ";
                for (int x : validV0)
                    fout << x << " ";
                fout << "\nV1: ";
                for (int x : validV1)
                    fout << x << " ";
                fout << "\n";
            }
            else
            {
                fout << "No theorem-valid partition found (bw=1)\n";
            }
        }
        fout << "----------------------------------------\n";
        fout.flush();
    }
    return 0;
}
