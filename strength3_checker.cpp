#include <bits/stdc++.h>
using namespace std;

/*
  Check if a binary generator vector g of length n generates
  a circulant CAOA(n, k, 2, 3, b) with bandwidth b <= 1.

  Construction:
    A is k x n, A[i][j] = g[(j + i) mod n].

  Test:
    For every choice of 3 rows (i0 < i1 < i2),
      - form the 3 x n submatrix
      - for each column, read the 3-tuple (b0,b1,b2) in {0,1}^3
      - count occurrences of each of the 8 patterns
      - bandwidth = max_count - min_count
    Global bandwidth = max over all 3-row choices.

    We require global bandwidth <= 1.
*/

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, k;
    if (!(cin >> n >> k)) {
        cerr << "Failed to read n and k.\n";
        return 1;
    }

    if (n <= 0 || k <= 0) {
        cerr << "n and k must be positive.\n";
        return 1;
    }

    string s;
    if (!(cin >> s)) {
        cerr << "Failed to read binary generator string.\n";
        return 1;
    }

    if ((int)s.size() != n) {
        cerr << "Error: expected binary string of length " << n
             << ", got length " << s.size() << ".\n";
        return 1;
    }

    vector<int> g(n);
    for (int i = 0; i < n; ++i) {
        if (s[i] == '0') g[i] = 0;
        else if (s[i] == '1') g[i] = 1;
        else {
            cerr << "Error: binary string must contain only '0' and '1'.\n";
            return 1;
        }
    }

    if (k > n) {
        cerr << "Warning: k > n; using first k shifts anyway.\n";
    }

    if (k < 3) {
        cout << "NO\n";
        cout << "Reason: k < 3, so strength 3 CAOA is impossible.\n";
        return 0;
    }

    // Build the k x n circulant array A
    vector<vector<int>> A(k, vector<int>(n));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            int idx = (j + i) % n; // shift by i
            A[i][j] = g[idx];
        }
    }

    int globalMaxBandwidth = 0;

    // Iterate over all 3-row subsets: i0 < i1 < i2
    for (int i0 = 0; i0 < k - 2; ++i0) {
        for (int i1 = i0 + 1; i1 < k - 1; ++i1) {
            for (int i2 = i1 + 1; i2 < k; ++i2) {

                // Count occurrences of each 3-tuple (8 possibilities)
                int cnt[8] = {0};

                for (int col = 0; col < n; ++col) {
                    int b0 = A[i0][col];
                    int b1 = A[i1][col];
                    int b2 = A[i2][col];

                    int idx = (b0 << 2) | (b1 << 1) | b2; // 0..7
                    ++cnt[idx];
                }

                int minCnt = cnt[0], maxCnt = cnt[0];
                for (int t = 1; t < 8; ++t) {
                    minCnt = min(minCnt, cnt[t]);
                    maxCnt = max(maxCnt, cnt[t]);
                }

                int bandwidth = maxCnt - minCnt;
                globalMaxBandwidth = max(globalMaxBandwidth, bandwidth);

                // Early exit if we exceed 1
                if (globalMaxBandwidth > 1) {
                    cout << "NO\n";
                    cout << "Max bandwidth observed = " << globalMaxBandwidth << "\n";
                    return 0;
                }
            }
        }
    }

    // If we get here, globalMaxBandwidth <= 1
    cout << "YES\n";
    cout << "Max bandwidth observed = " << globalMaxBandwidth << "\n";
    return 0;
}
