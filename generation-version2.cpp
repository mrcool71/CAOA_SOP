#include <bits/stdc++.h>
using namespace std;
using Mat = vector<vector<int>>;
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
                    if (d == 0)
                        d = n;
                    sub[i][j] = d;
                }
            }
            diff[a][b] = sub;
        }
    }
    return diff;
}
vector<Mat> first_order(const vector<vector<vector<vector<int>>>> &diff, int n)
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
map<pair<int, int>, Mat> second_order(const vector<vector<vector<vector<int>>>> &diff, int n)
{
    int s = 2;
    int maxR = n - 1;
    map<pair<int, int>, Mat> res;
    for (int r1 = 1; r1 <= maxR; ++r1)
    {
        for (int r2 = r1 + 1; r2 <= maxR; ++r2)
        {
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
            res[{r1, r2}] = M;
        }
    }
    return res;
}
bool mat_eq(const Mat &A, const Mat &B)
{
    if (A.size() != B.size())
        return false;
    for (size_t i = 0; i < A.size(); ++i)
        if (A[i] != B[i])
            return false;
    return true;
}
int prefix_equal_first(const vector<Mat> &first, int n)
{
    int maxR = n - 1;
    int p = 1;
    for (int r = 2; r <= maxR; ++r)
    {
        if (!mat_eq(first[r], first[1]))
            break;
        p = r;
    }
    return p;
}
bool second_all_equal_up_to(const map<pair<int, int>, Mat> &second, int upto, const vector<Mat> &first)
{
    if (upto <= 1)
        return true;

    // For strength t >= 2, the reference matrix MUST be for the pair {1, 2}.
    // If it doesn't exist, strength 'upto' is impossible.
    auto ref_it = second.find({1, 2});
    if (ref_it == second.end())
    {
        return false; // Cannot have strength >= 2 without Phi2(1, 2).
    }
    const Mat &ref = ref_it->second;

    // Now, check that all other required pairs exist and match the reference.
    for (int i = 1; i <= upto; ++i)
    {
        for (int j = i + 1; j <= upto; ++j)
        {
            // We already have the reference for {1, 2}, so skip it.
            if (i == 1 && j == 2)
                continue;

            auto it = second.find({i, j});
            if (it == second.end())
            {
                // If a required pair {i, j} doesn't have a second-order property, fail.
                return false;
            }
            if (!mat_eq(it->second, ref))
            {
                // If any pair's matrix doesn't match the {1, 2} reference, fail.
                return false;
            }
        }
    }

    // If all checks passed, the properties are uniform up to 'upto'.
    return true;
}

int max_t_binary(const map<pair<int, int>, Mat> &second, int p, const vector<Mat> &first)
{
    if (p <= 1)
        return p;
    int lo = 1, hi = p, ans = 1;
    while (lo <= hi)
    {
        int mid = (lo + hi) / 2;
        if (second_all_equal_up_to(second, mid, first))
        {
            ans = mid;
            lo = mid + 1;
        }
        else
            hi = mid - 1;
    }
    return ans;
}

int find_k_for_partition(const vector<int> &V0, const vector<int> &V1, int n)
{
    auto diff = generate_difference_matrix(V0, V1, n);
    auto first = first_order(diff, n);
    auto second = second_order(diff, n);
    int p = prefix_equal_first(first, n);
    int t = max_t_binary(second, p, first);
    int k = t + 1;

    if (n == 20 && k == 5)
    {
        // This block can be used for targeted debugging if the issue persists.
        // For now, the logic fix should prevent this from being entered.
    }

    return k;
}

void print_generating_vector(const vector<int> &V0, int n)
{
    vector<int> gen_vec(n + 1, 1);
    for (int x : V0)
    {
        if (x >= 1 && x <= n)
        {
            gen_vec[x] = 0;
        }
    }
    cout << "Generating Vector:";
    for (int i = 1; i <= n; ++i)
    {
        cout << " " << gen_vec[i];
    }
    cout << "\n";
}

void print_caoa(const vector<int> &V0, int n, int strength_k)
{
    vector<int> base_row(n + 1, 1);
    for (int x : V0)
    {
        if (x >= 1 && x <= n)
        {
            base_row[x] = 0;
        }
    }

    int strength_t = strength_k - 1;
    // Corrected Notation: CAOA(N, k, s, t) -> CAOA(strength_k, n, 2, strength_t)
    cout << "CAOA(" << n << ", " << strength_k << ", 2, " << strength_t << "):\n";
    vector<int> current_row;
    for (int i = 1; i <= n; ++i)
        current_row.push_back(base_row[i]);

    // Print 'strength_k' rows, not 'n' rows
    for (int i = 0; i < strength_k; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << current_row[j] << (j == n - 1 ? "" : " ");
        }
        cout << "\n";
        // Cyclically shift the row for the next iteration
        rotate(current_row.begin(), current_row.begin() + 1, current_row.end());
    }
}

void write_tabular_result(ofstream &outfile, int n, int bestK, const vector<int> &bestV0, long long duration_ms)
{
    if (!outfile.is_open())
        return;

    // Column 1: n
    outfile << n << "\t";

    // Column 2: k
    outfile << bestK << "\t";

    // Column 3: vector
    vector<int> gen_vec(n + 1, 1);
    for (int x : bestV0)
    {
        if (x >= 1 && x <= n)
        {
            gen_vec[x] = 0;
        }
    }
    for (int i = 1; i <= n; ++i)
    {
        outfile << gen_vec[i] << (i == n ? "" : " ");
    }
    outfile << "\t";

    // Column 4: time
    outfile << duration_ms << "\n";
}

int main()
{
    ofstream outfile("generation_results.txt");
    if (!outfile.is_open())
    {
        cerr << "Failed to open output file." << endl;
        return 1;
    }

    cout << "=== CAOA Generation Tool ===\n";
    cout << "Generating optimal vectors for n < 20...\n";

    // Write header for tabular data
    outfile << "n\tk\tvector\ttime_ms\n";

    for (int n = 1; n <= 10; ++n)
    {
        cout << "Processing n = " << n << "..." << endl;
        auto start_time = chrono::high_resolution_clock::now();

        int sizeV0 = (n % 2 == 0) ? (n / 2) : ((n - 1) / 2);
        int choose = sizeV0 - 1;
        vector<int> elems;
        for (int i = 2; i <= n; ++i)
            elems.push_back(i);
        int m = elems.size();
        vector<int> bit(m, 0);
        for (int i = 0; i < choose; ++i)
            bit[i] = 1;
        sort(bit.begin(), bit.end(), greater<int>());

        int bestK = 0;
        set<vector<int>> best_partitions;

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

            int k = find_k_for_partition(V0, V1, n);
            if (k > bestK)
            {
                bestK = k;
                best_partitions.clear();
                best_partitions.insert(V0);
            }
            else if (k == bestK)
            {
                best_partitions.insert(V0);
            }

        } while (prev_permutation(bit.begin(), bit.end()));

        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

        // Write results to file in tabular format
        if (!best_partitions.empty())
        {
            // Pass the first found best partition to the writing function
            write_tabular_result(outfile, n, bestK, *best_partitions.begin(), duration.count());
        }
    }

    outfile.close();
    cout << "\nSearch complete!\n";
    cout << "Results written to generation_results.txt\n";

    return 0;
}