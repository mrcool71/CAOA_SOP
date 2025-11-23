#include <bits/stdc++.h>
using namespace std;
using Mat = vector<vector<int>>;
vector<vector<vector<vector<int>>>> generate_difference_matrix(const vector<int>& V0, const vector<int>& V1, int n){
    vector<vector<vector<vector<int>>>> diff(2, vector<vector<vector<int>>>(2));
    vector<vector<int>> P = {V0, V1};
    for(int a=0;a<2;a++){
        for(int b=0;b<2;b++){
            vector<vector<int>> sub(P[a].size(), vector<int>(P[b].size()));
            for(size_t i=0;i<P[a].size();++i){
                for(size_t j=0;j<P[b].size();++j){
                    int d = (P[a][i] - P[b][j]) % n;
                    if(d<0) d += n;
                    if(d==0) d = n;
                    sub[i][j] = d;
                }
            }
            diff[a][b] = sub;
        }
    }
    return diff;
}
vector<Mat> first_order(const vector<vector<vector<vector<int>>>> &diff, int n){
    int s = 2;
    int maxR = n-1;
    vector<Mat> res(maxR+1, Mat(s, vector<int>(s,0)));
    for(int r=1;r<=maxR;++r){
        for(int a=0;a<s;++a){
            for(int b=0;b<s;++b){
                int cnt=0;
                const auto &sub = diff[a][b];
                for(size_t i=0;i<sub.size();++i) for(size_t j=0;j<sub[i].size();++j) if(sub[i][j]==r) ++cnt;
                res[r][a][b] = cnt;
            }
        }
    }
    return res;
}
map<pair<int,int>,Mat> second_order(const vector<vector<vector<vector<int>>>> &diff, int n){
    int s = 2;
    int maxR = n-1;
    map<pair<int,int>,Mat> res;
    for(int r1=1;r1<=maxR;++r1){
        for(int r2=r1+1;r2<=maxR;++r2){
            Mat M(s, vector<int>(s,0));
            for(int a=0;a<s;++a){
                for(int b=0;b<s;++b){
                    int rowsWithBoth=0;
                    const auto &sub = diff[a][b];
                    for(size_t i=0;i<sub.size();++i){
                        bool h1=false, h2=false;
                        for(size_t j=0;j<sub[i].size();++j){
                            if(sub[i][j]==r1) h1=true;
                            if(sub[i][j]==r2) h2=true;
                        }
                        if(h1 && h2) ++rowsWithBoth;
                    }
                    M[a][b] = rowsWithBoth;
                }
            }
            res[{r1,r2}] = M;
        }
    }
    return res;
}
bool mat_eq(const Mat &A, const Mat &B){
    if(A.size()!=B.size()) return false;
    for(size_t i=0;i<A.size();++i) if(A[i]!=B[i]) return false;
    return true;
}
int prefix_equal_first(const vector<Mat> &first, int n){
    int maxR = n-1;
    int p = 1;
    for(int r=2;r<=maxR;++r){
        if(!mat_eq(first[r], first[1])) break;
        p = r;
    }
    return p;
}
bool second_all_equal_up_to(const map<pair<int,int>,Mat> &second, int upto){
    if(upto<=1) return true;
    auto it = second.find({1,2});
    if(it==second.end()) return false;
    const Mat &ref = it->second;
    for(int i=1;i<=upto;++i){
        for(int j=i+1;j<=upto;++j){
            auto it2 = second.find({i,j});
            if(it2==second.end()) return false;
            if(!mat_eq(it2->second, ref)) return false;
        }
    }
    return true;
}
int max_t_binary(const map<pair<int,int>,Mat> &second, int p){
    if(p<=1) return p;
    int lo=1, hi=p, ans=1;
    while(lo<=hi){
        int mid = (lo+hi)/2;
        if(second_all_equal_up_to(second, mid)){ ans = mid; lo = mid+1; }
        else hi = mid-1;
    }
    return ans;
}
int find_k_for_partition(const vector<int>& V0, const vector<int>& V1, int n){
    auto diff = generate_difference_matrix(V0, V1, n);
    auto first = first_order(diff, n);
    auto second = second_order(diff, n);
    int p = prefix_equal_first(first, n);
    int t = max_t_binary(second, p);
    int k = t+1;
    return k;
}
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if(!(cin>>n)) return 0;
    if(n<4) return 0;
    int sizeV0 = (n % 2 == 0) ? (n/2) : ((n-1)/2);
    int choose = sizeV0 - 1;
    vector<int> elems;
    for(int i=2;i<=n;++i) elems.push_back(i);
    int m = elems.size();
    vector<int> bit(m,0);
    for(int i=0;i<choose;++i) bit[i]=1;
    sort(bit.begin(), bit.end(), greater<int>());
    int bestK = 1;
    vector<int> bestV0, bestV1;
    Mat bestPhi1, bestPhi2;
    do{
        vector<int> V0;
        V0.push_back(1);
        for(int i=0;i<m;++i) if(bit[i]) V0.push_back(elems[i]);
        vector<int> mark(n+1,0);
        for(int x: V0) mark[x]=1;
        vector<int> V1;
        for(int i=1;i<=n;++i) if(!mark[i]) V1.push_back(i);
        int k = find_k_for_partition(V0, V1, n);
        if(k > bestK){
            bestK = k;
            bestV0 = V0;
            bestV1 = V1;
            auto diff = generate_difference_matrix(V0, V1, n);
            auto first = first_order(diff, n);
            bestPhi1 = first[1];
            if(k-1 >= 2) bestPhi2 = second_order(diff, n).at({1,2});
            else bestPhi2 = Mat(2, vector<int>(2,0));
        }
    } while(prev_permutation(bit.begin(), bit.end()));
    cout << "best k = " << bestK << "\n";
    cout << "V0:";
    for(int x: bestV0) cout << " " << x;
    cout << "\nV1:";
    for(int x: bestV1) cout << " " << x;
    cout << "\nPhi1:\n";
    for(size_t i=0;i<bestPhi1.size();++i){
        for(size_t j=0;j<bestPhi1[i].size();++j) cout << bestPhi1[i][j] << " ";
        cout << "\n";
    }
    cout << "Phi2:\n";
    for(size_t i=0;i<bestPhi2.size();++i){
        for(size_t j=0;j<bestPhi2[i].size();++j) cout << bestPhi2[i][j] << " ";
        cout << "\n";
    }
    return 0;
}
