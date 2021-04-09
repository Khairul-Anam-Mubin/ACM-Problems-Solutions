/***
**      Bismillahir Rahmanir Rahim                
**              ALLAHU AKBAR
**
**     Author: Khairul Anam Mubin (__Broly__)
**     Bangladesh University of Business and Technology,
**     Dept. of CSE.
***/
#include <bits/stdc++.h>
using namespace std;
 
#define F            first
#define S            second  
#define Fin          freopen("input.txt","r",stdin)
#define Fout         freopen("output.txt","w",stdout)
#define Precision(a) cout << fixed << setprecision(a)
#define FasterIO     ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0); 
#define Test         int test;cin >> test;for(int tc = 1; tc <= test; tc++)
 
#define INF9         2147483647
#define INF18        9223372036854775806
#define eps          1e-8 
 
const double pi = 2 * acos(0.0);
typedef long long ll;
typedef unsigned long long ull;
typedef long double ld;

template <typename T> T Sqr(T x) { T n = x * x ; return n ;}
template <typename T> T Pow(T B,T P){ if(P==0) return 1; if(P&1) return B*Pow(B,P-1);  else return Sqr(Pow(B,P/2));}
template <typename T> T Abs(T a) {if(a<0)return -a;else return a;}
template <typename T> T Gcd(T a,T b){if(a<0)return Gcd(-a,b);if(b<0)return Gcd(a,-b);return (b==0)?a:Gcd(b,a%b);}
template <typename T> T Lcm(T a,T b) {if(a<0)return Lcm(-a,b);if(b<0)return Lcm(a,-b);return a*(b/Gcd(a,b));}
template <typename T> T Exgcd(T a,T b,T &x,T &y) {if(a<0){T d=Exgcd(-a,b,x,y);x=-x;return d;}   if(b<0){T d=Exgcd(a,-b,x,y);y=-y;return d;}   if(b==0){x=1;y=0;return a;}else{T d=Exgcd(b,a%b,x,y);T t=x;x=y;y=t-(a/b)*y;return d;}}
template <typename T> T BigMod (T b,T p,T m){if (p == 0) return 1;if (p%2 == 0){T s = BigMod(b,p/2,m);return ((s%m)*(s%m))%m;}return ((b%m)*(BigMod(b,p-1,m)%m))%m;}
template <typename T> T ModInvPrime (T b,T m){return BigMod(b,m-2,m);}
template <typename T> T ModInvNotPrime(T a , T m) {T x , y ;Exgcd(a , m , x , y) ;x %= m ;if(x < 0) x += m ;return x ;}
template <typename T> typename std::vector<T>::iterator Insert_sorted (std :: vector<T> & vec, T const& item) {return vec.insert (std::upper_bound( vec.begin(), vec.end(), item ), item);}
template <typename T> inline string ToBinary(T n) {string r ;while(n != 0) {r = ( n % 2 == 0 ? "0" : "1" ) + r ; n >>= 1;} return r ;} 
long long BinaryToDecimal(string s) {int len = s.size();long long n = 0, p = 1;for (int i = len - 1; i >= 0; i-- , p *= 2) n += p * (s[i] - '0');return n;}

char Uplowch(char ch){if(ch >= 'A' &&  ch <= 'Z') ch += 32; return ch;}
char Lowupch(char ch){if(ch >= 'a' &&  ch <= 'z') ch -= 32; return ch;}
bool Isalpha(char ch){if((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z')) return true; return false;}
int Strtoint(string str){stringstream ss(str);int x = 0;ss >> x ;return x ;}
string Intostr(int x){stringstream ss; ss << x; string str = ss.str(); return str;}
vector<string> Linetostr(string str){string s; vector <string> v ;istringstream is(str); while(is >> s)v.push_back(s);return v;}
template <typename T> void Print(T ar[] , int a , int b) {for (int i = a; i + 1 <= b ; i++)cout << ar[i] << " ";cout << ar[b] << "\n";}
template <typename T> void Print(T ar[] , int n) {for (int i = 0; i + 1 < n; i++)cout << ar[i] << " ";cout << ar[n - 1] << "\n";} 
template <typename T> void Print(const vector<T> &v) {for (int i = 0; i + 1 < v.size() ; i++) cout << v[i] << " ";cout << v.back() << "\n";}

//int dx[] = {0 , 0 , -1 , 1 , -1 , -1 , 1 , 1 , 0} ;  // right , left , up , down
//int dy[] = {1 ,-1 , 0 , 0 , -1 , 1 , -1 , 1 , 0} ;

/*******************String Hashing*********************/
// Double Hashing
// 1. Modular Exponentian Needed
// 2. Init must be call and set the maximum length of the string
// 3. If sub string hash required then Compute hash have to call
// 4. If prefix and suffix hash required ComputePreAndSufHash have to call

struct ModularExponentiation {
    template <typename T> T Pow(T b, T p) {
        T res = 1;
        while (p > 0) {
            if (p % 2 == 1) res = res * b;
            b = b * b;
            p /= 2;
        }
        return res;
    }
    template <typename T> T Mod(T a, T m) { 
        return (((a % m) + m) % m);
    }
    template <typename T> T BigMod(T b, T p, T m) {
        T res = 1;
        if (b > m) b %= m;
        while (p) {
            if (p % 2 == 1) res = res * b % m;
            b = b * b % m;
            p /= 2;
        }
        return res;
    }
    template <typename T> T ModInv(T b,T m) {
        return BigMod(b , m - 2 , m);
    }
}Ex;
struct DoubleHashing {
    long long base[2] = {1949313259, 1997293877};
    long long mod[2] = {2091573227, 2117566807};
    vector <long long> pow[2] , inv[2];
    vector <long long> prehash[2] , sufhash[2];
    int maxN , flag = 0 , len;
    void Init(int n) {
        maxN = n + 2;
        for (int i = 0; i < 2; i++) {
            pow[i].resize(maxN);
            inv[i].resize(maxN);
        }
        Generate();
    }
    void Generate() {
        for (int j = 0; j < 2; j++) {
            pow[j][0] = 1;
            inv[j][0] = 1;
            long long minv = Ex.ModInv(base[j] ,mod[j]);
            for (int i = 1; i < maxN; i++) {
                pow[j][i] = pow[j][i - 1] * base[j] % mod[j];
                inv[j][i] = inv[j][i - 1] * minv % mod[j];
            }
        }
    }
    long long GetHash(string &s) {
        long long hash_val[2] = {0 , 0};
        int n = s.size();
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < n; i++) {
                hash_val[j] = (hash_val[j] + s[i] * pow[j][i]) % mod[j];
            }
        }
        return (hash_val[0] << 32LL) | hash_val[1];
    }
    void ComputeHash(string &s) {
        flag = 1;
        len = s.size();
        for (int j = 0; j < 2; j++) prehash[j].resize(maxN);
        for (int j = 0; j < 2; j++) prehash[j][0] = 0;
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < len; i++) {
                prehash[j][i + 1] = (prehash[j][i] + pow[j][i] * s[i]) % mod[j];
            }
        }
    }
    long long GetSubstrHash(int l , int r) {
        if (!flag) { cout << "ComputeHash\n"; return -1;}
        long long hash_val[2];
        for (int j = 0; j < 2; j++) 
            hash_val[j] = (prehash[j][r + 1] - prehash[j][l]) * inv[j][l] % mod[j];
        for (int j = 0; j < 2; j++) if (hash_val[j] < 0) hash_val[j] += mod[j];   
        return (hash_val[0] << 32) | hash_val[1]; 
    }
    void ComputePreAndSufHash(string &s) {
        flag = 1;
        len = s.size();
        for (int j = 0; j < 2; j++) { 
            prehash[j].resize(maxN);
            sufhash[j].resize(maxN);
        }
        for (int j = 0; j < 2; j++) prehash[j][0] = sufhash[j][0] = 0;
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < len; i++) {
                prehash[j][i + 1] = (prehash[j][i] + pow[j][i] * s[i]) % mod[j];
                sufhash[j][i + 1] = (sufhash[j][i] + pow[j][len - i + 1] * s[i]) % mod[j];
            }
        }
    }
    long long GetPrefixHash(int l , int r) {
        return GetSubstrHash(l , r);
    }
    long long GetSuffixHash(int l , int r) {
        if (!flag) { cout << "ComputePreAndSufHash\n"; return -1;}
        long long hash_val[2];
        for (int j = 0; j < 2; j++) 
            hash_val[j] = (sufhash[j][r + 1] - sufhash[j][l]) * inv[j][len - r + 1] % mod[j];
        for (int j = 0; j < 2; j++) if (hash_val[j] < 0) hash_val[j] += mod[j];   
        return (hash_val[0] << 32) | hash_val[1];  
    } 
    bool IsPallindrome(int l , int r) {
        return (GetPrefixHash(l , r) == GetSuffixHash(l , r));
    }
    vector <int> RabinKarp(string &txt , string &ptrn) {
        ComputeHash(txt);
        long long ptrn_hash = GetHash(ptrn);
        vector <int> occurences;
        int txtlen = txt.size();
        int ptrnlen = ptrn.size();
        for (int i = 0; i < txtlen - ptrnlen + 1; i++) {
            long long cur_hash = GetSubstrHash(i , ((i + ptrnlen) - 1));
            // pattern match...
            if (cur_hash == ptrn_hash) 
                occurences.emplace_back(i + 1);
        }
        return occurences;
    }
} hs;
int main() {
    FasterIO
    hs.Init(100000);
    Test {
        int n; cin >> n;
        set <string> st;
        string s;
        for (int i = 0; i < n; i++) {
            cin >> s;
            st.insert(s);
        }
        n = st.size();
        vector <pair <ll , ll>> vec[26];
        for (auto it : st) {
            s = it;
            ll len = s.size();
            ll val = hs.GetHash(s);
            int ch = s[0] - 'a';
            vec[ch].push_back(make_pair(len , val));
        }
        cin >> s;
        hs.ComputeHash(s);
        n = s.size();
        ll dp[n + 2];
        for (int i = 0; i <= n; i++) dp[i] = INF9;
        dp[0] = 0;
        for (int i = 0; i < n; i++) {
            int ch = s[i] - 'a';
            for (auto it: vec[ch]) {
                ll len = it.first;
                ll val = it.second;
                if (i + len - 1 < n) {
                    ll subhash = hs.GetSubstrHash(i, i + len - 1);
                    if (subhash == val) {
                        dp[i + len] = min(dp[i + len] , dp[i] + 1);
                    }
                } else break;
            }
        }
        cout << "Case " << tc << ": ";
        if (dp[n] > n) cout << "impossible\n";
        else cout << dp[n] << "\n";
    }
    return 0;
}