#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using Mat = vector<vector<int>>;

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
                    if (d == 0)
                        d = n;
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

map<pair<int, int>, Mat> second_order(const vector<vector<vector<vector<int>>>>& diff, int n)
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
                            ++rowsWithBoth;
                    }
                    M[a][b] = rowsWithBoth;
                }
            }
            res[{r1, r2}] = move(M);
        }
    }
    return res;
}

bool mat_eq(const Mat& A, const Mat& B)
{
    if (A.size() != B.size())
        return false;
    for (size_t i = 0; i < A.size(); ++i)
        if (A[i] != B[i])
            return false;
    return true;
}

int prefix_equal_first(const vector<Mat>& first, int n)
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

bool second_all_equal_up_to(const map<pair<int, int>, Mat>& second, int upto)
{
    if (upto <= 1)
        return true;
    auto ref_it = second.find({1, 2});
    if (ref_it == second.end())
        return false;
    const Mat& ref = ref_it->second;
    for (int i = 1; i <= upto; ++i)
    {
        for (int j = i + 1; j <= upto; ++j)
        {
            auto it = second.find({i, j});
            if (it == second.end())
                return false;
            if (!mat_eq(it->second, ref))
                return false;
        }
    }
    return true;
}

int max_t_binary(const map<pair<int, int>, Mat>& second, int p)
{
    if (p <= 1)
        return p;
    int lo = 1, hi = p, ans = 1;
    while (lo <= hi)
    {
        int mid = (lo + hi) / 2;
        if (second_all_equal_up_to(second, mid))
        {
            ans = mid;
            lo = mid + 1;
        }
        else
            hi = mid - 1;
    }
    return ans;
}

int find_k_for_partition(const vector<int>& V0, const vector<int>& V1, int n)
{
    auto diff = generate_difference_matrix(V0, V1, n);
    auto first = first_order(diff, n);
    auto second = second_order(diff, n);
    int p = prefix_equal_first(first, n);
    int t = max_t_binary(second, p);
    return t + 1;
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ofstream outfile("generation_results_mt.txt");
    outfile << "n\tk\tvector\ttime_ms\n";

    for (int n = 1; n <= 10; ++n)
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

        int bestK = 0;
        vector<int> bestV0;
        auto t0 = chrono::high_resolution_clock::now();

        #pragma omp parallel
        {
            int localBest = 0;
            vector<int> localV0;
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
                int k = find_k_for_partition(V0, V1, n);
                if (k > localBest)
                {
                    localBest = k;
                    localV0 = V0;
                }
            }
            #pragma omp critical
            {
                if (localBest > bestK)
                {
                    bestK = localBest;
                    bestV0 = localV0;
                }
            }
        }

        auto t1 = chrono::high_resolution_clock::now();
        auto ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

        vector<int> gen_vec(n + 1, 1);
        for (int x : bestV0)
            if (x >= 1 && x <= n)
                gen_vec[x] = 0;
        outfile << n << "\t" << bestK << "\t";
        for (int i = 1; i <= n; ++i)
            outfile << gen_vec[i] << (i == n ? "\t" : " ");
        outfile << ms << "\n";
    }
    return 0;
}
