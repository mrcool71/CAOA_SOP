#include <bits/stdc++.h>
using namespace std;

/*
  Check if a binary generator vector g of length n generates
  a circulant CAOA(n, k, 2, 4, b) with bandwidth b <= 1.

  Definition used (Lin–Phoa–Kao type CAOA):
  - A is k x n, entries in {0,1}.
  - For every choice of 4 rows (a 4 x n submatrix),
    each ordered 4-tuple alpha in {0,1}^4 occurs lambda(alpha) times
    as a column vector.
  - Bandwidth of that submatrix = max_alpha lambda(alpha) - min_alpha lambda(alpha).
  - Global bandwidth = maximum over all 4-row choices.
  - We require global bandwidth <= 1.

  Construction of the circulant array:
    A[i][j] = g[(j + i) mod n], for i = 0..k-1, j = 0..n-1.
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
        cerr << "Warning: k > n; using k rows (first k shifts) anyway.\n";
    }

    if (k < 4) {
        cout << "NO\n";
        cout << "Reason: k < 4, so strength 4 CAOA is impossible.\n";
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

    // Iterate over all 4-row subsets: i0 < i1 < i2 < i3
    for (int i0 = 0; i0 < k - 3; ++i0) {
        for (int i1 = i0 + 1; i1 < k - 2; ++i1) {
            for (int i2 = i1 + 1; i2 < k - 1; ++i2) {
                for (int i3 = i2 + 1; i3 < k; ++i3) {

                    // Count occurrences of each 4-tuple (16 possibilities)
                    int cnt[16] = {0};

                    for (int col = 0; col < n; ++col) {
                        int b0 = A[i0][col];
                        int b1 = A[i1][col];
                        int b2 = A[i2][col];
                        int b3 = A[i3][col];

                        int idx = (b0 << 3) | (b1 << 2) | (b2 << 1) | b3;
                        // idx ranges from 0 to 15
                        ++cnt[idx];
                    }

                    int minCnt = cnt[0], maxCnt = cnt[0];
                    for (int t = 1; t < 16; ++t) {
                        minCnt = min(minCnt, cnt[t]);
                        maxCnt = max(maxCnt, cnt[t]);
                    }

                    int bandwidth = maxCnt - minCnt;
                    globalMaxBandwidth = max(globalMaxBandwidth, bandwidth);

                    // Early exit if we already exceed 1
                    if (globalMaxBandwidth > 1) {
                        cout << "NO\n";
                        cout << "Max bandwidth observed = " << globalMaxBandwidth << "\n";
                        return 0;
                    }
                }
            }
        }
    }

    // If we got here, globalMaxBandwidth <= 1
    cout << "YES\n";
    cout << "Max bandwidth observed = " << globalMaxBandwidth << "\n";
    return 0;
}
