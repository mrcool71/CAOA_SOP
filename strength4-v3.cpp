// caoa_search.cpp
#include <bits/stdc++.h>
using namespace std;
using hrc = chrono::high_resolution_clock;

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ofstream fout("results.txt");
    if (!fout.is_open())
    {
        cerr << "Failed to open results.txt for writing\n";
        return 1;
    }

    fout << "CAOA search results (n = 8..20). Each line shows: n, elapsed_seconds, k_found (or 0 if none), generating_vector (0/1)\n\n";
    fout.flush();
    cout << "Starting search. This may take some time for larger n...\n";
    cout.flush();

    for (int n = 1; n <= 8; ++n)
    {
        auto t0 = hrc::now();

        int zeros = n / 2; // floor(n/2)
        int ones = n - zeros;
        int max_k = min(10, n);
        bool found_for_n = false;
        int found_k = 0;
        vector<int> found_vec;
        // We'll try k from max_k down to 4
        vector<int> counts;
        for (int k = max_k; k >= 4 && !found_for_n; --k)
        {
            // Precompute all 4-subsets of rows (indices 0..k-1)
            int count = 0;
            vector<array<int, 4>> subsets;
            for (int a = 0; a < k; ++a)
                for (int b = a + 1; b < k; ++b)
                    for (int c = b + 1; c < k; ++c)
                        for (int d = c + 1; d < k; ++d)
                            subsets.push_back({a, b, c, d});

            // Prepare initial bit vector: zeros zeros ... ones ones
            vector<int> bits(n);
            for (int i = 0; i < zeros; ++i)
                bits[i] = 0;
            for (int i = zeros; i < n; ++i)
                bits[i] = 1;

            // Use next_permutation to iterate all unique balanced vectors
            bool first_perm = true;
            do
            {
                // rows[s][j] = bits[(j + s) % n]
                // Build rows as vector of vectors of uint8_t for speed
                static uint8_t rows_storage[12][24]; // k <= 10 <=12; n <= 20 <=24
                for (int s = 0; s < k; ++s)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        int idx = j + s;
                        if (idx >= n)
                            idx -= n;
                        rows_storage[s][j] = (uint8_t)bits[idx];
                    }
                }

                bool pattern_ok = true;
                // Check each 4-subset
                for (size_t si = 0; si < subsets.size(); ++si)
                {
                    const array<int, 4> &sel = subsets[si];
                    int a = sel[0], b = sel[1], c = sel[2], d = sel[3];

                    // counts of 16 possible 4-tuples
                    int cnt[16];
                    for (int t = 0; t < 16; ++t)
                        cnt[t] = 0;

                    for (int j = 0; j < n; ++j)
                    {
                        int v = (rows_storage[a][j] << 3) | (rows_storage[b][j] << 2) | (rows_storage[c][j] << 1) | (rows_storage[d][j]);
                        cnt[v]++;
                    }

                    int mn = cnt[0], mx = cnt[0];
                    for (int t = 1; t < 16; ++t)
                    {
                        if (cnt[t] < mn)
                            mn = cnt[t];
                        if (cnt[t] > mx)
                            mx = cnt[t];
                    }
                    if (mx - mn > 1)
                    {
                        pattern_ok = false;
                        break; // this pattern fails for this subset
                    }
                }

                if (pattern_ok)
                {
                    // success for this k and this generating vector
                    found_for_n = true;
                    found_k = k;
                    found_vec = bits;
                    break; // stop enumerating patterns for this k
                }

                first_perm = false;
                count++;
            } while (next_permutation(bits.begin(), bits.end()));

            counts.push_back(count);
            if (found_for_n)
                break;
            // else try next lower k
        }

        auto t1 = hrc::now();
        chrono::duration<double> elapsed = t1 - t0;
        double secs = elapsed.count();

        if (found_for_n)
        {
            fout << "n=" << n << "  elapsed=" << fixed << setprecision(4) << secs
                 << "s  k=" << found_k << "  vec=";
            for (int b : found_vec)
                fout << b;
            fout << "\n";
            fout.flush();
            cout << "n=" << n << " -> found k=" << found_k << " (elapsed " << secs << "s, count " << counts.back() << ")\n";
            cout.flush();
        }
        else
        {
            fout << "n=" << n << "  elapsed=" << fixed << setprecision(4) << secs
                 << "s  k=0  (no balanced generating vector found for k in " << min(10, n) << "..4)\n";
            fout.flush();
            cout << "n=" << n << " -> none found (elapsed " << secs << "s)\n";
            cout.flush();
        }
    }

    fout << "\nSearch finished.\n";
    fout.close();
    cout << "Done. Results written to results.txt\n";
    return 0;
}
