#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#define ll long long

using namespace std;

int incnt = 0, outcnt = 0;

bool IS_End = false;



map<string, int> Name2Num;
map<int, string> Num2Name;

struct Command {
    bool EndJdg = false;
    bool TransposeJdg = false;
    bool InverseJdg = false;
    bool DeterminantJdg = false;
    bool RankJdg = false;
    bool Read_Finish = false;
    bool Read_Begin = false;
    bool nameJdg = false;
    bool calJdg = false;
    bool rowOPJdg = false;
};

void Conduct(Command &kmd);
Command Cmd();

struct Matrix {
    vector< vector<double> > Matx;
    vector< vector<double> > Matx_T;

    int n, m, n_T, m_T;
    bool IS_Aug = false;
    bool IS_legal = true;
    bool IS_Named = false;
    int num = 0;
    string Name;

    void ReadMatrix() {
        bool IS_End = false;
        int t = 0;
        int i = 0, j = 0;
        bool IS_Num = false;
        Command CmdRcoder;
        Matx.push_back(vector<double>{});
        Matx_T.push_back(vector<double>{});
        while(!IS_End) {
            char c = getchar();
            if(c == '\\') { CmdRcoder = Cmd(); IS_End = CmdRcoder.Read_Finish; n = i; continue;}
            else if(IS_Num && (c == '\n' || c == '\r' || c == EOF)) {
                if(Matx.size() <= i) Matx.push_back(vector<double>{});
                Matx[i].push_back((double)(t));
                if(Matx_T.size() <= j) Matx_T.push_back(vector<double>{});
                Matx_T[j].push_back((double)(t));
                ++i; j = 0;
                IS_Num = false; t = 0;
                continue;
            }
            else if(IS_Num && (c == ',' || c == ' ')) {
                if(Matx.size() <= i) Matx.push_back(vector<double>{});
                Matx[i].push_back((double)(t));
                if(Matx_T.size() <= j) Matx_T.push_back(vector<double>{});
                Matx_T[j].push_back((double)(t));
                t = 0; ++j;
                IS_Num = false;
                continue;
            }
            else if( '0' <= c && c <= '9') {
                t = (t<<3) + (t<<1) + (c^48);
                IS_Num = true;
                continue;
            }
            else if(IS_Num && c == '|') {
                if(Matx.size() <= i) Matx.push_back(vector<double>{});
                Matx[i].push_back((double)(t));
                if(Matx_T.size() <= j) Matx_T.push_back(vector<double>{});
                Matx_T[j].push_back((double)(t));
                IS_Aug = true;
                ++j;
                IS_Num = false; t = 0;
                continue;
            }
        }
        num = ++incnt + outcnt;
        m = Matx[0].size();
        n = Matx.size();
        n_T = m; m_T = n;

    }
    void PrintMatrix() {
        if(!IS_legal) { printf("Error: Illegal matrix operation.\n"); return; }
        printf("Matrix number is %d%s%s :", num, IS_Named ? ", Name is " : "", IS_Named ? Name : "");
        int width = 0;

        for (const auto& row : Matx) {
            for (double value : row) {
                int length = static_cast<int>(std::to_string(value).size());
                if (length > width)
                    width = length;
            }
        }

        for (const auto& row : Matx) {
            for (double value : row)
                printf("%*g", width + 1, value);
            printf("\n");
        }
    }
/*      for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                printf("%d%c", Matx[i][j], j == m-1 ? '\n' : ',');
            }
        }   */

    Matrix operator+(const Matrix &b) {
        Matrix Result;
        if(n != b.n || m != b.m) { printf("Error: Matrix dimensions do not match for addition.\n"); Result.IS_legal = false; return Result; }
        
        int k = 1;
        Result.n = n; Result.m = m;
        Result.Matx.resize(Result.n, vector<double>(Result.m, 0));
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                Result.Matx[i][j] = Matx[i][j] + 1.0 * k * b.Matx[i][j];
            }
        }

        return Result;
    }

    Matrix operator-(const Matrix &b) {
        Matrix Result;
        if(n != b.n || m != b.m) { printf("Error: Matrix dimensions do not match for reduction.\n"); Result.IS_legal = false; return Result; }
        
        int k = -1;
        Result.n = n; Result.m = m;
        Result.Matx.resize(Result.n, vector<double>(Result.m, 0));
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                Result.Matx[i][j] = Matx[i][j] + 1.0 * k * b.Matx[i][j];
            }
        }
        
        return Result;
    }
    Matrix operator*(const Matrix &b) {
        Matrix Result;
        if(m != b.n) { printf("Error: Matrix dimensions do not match for multiplication.\n"); Result.IS_legal = false; return Result; }
        Result.n = n; Result.m = b.m;
        Result.Matx.resize(Result.n, vector<double>(Result.m, 0));
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < b.m; ++j) {
                Result.Matx[i][j] = 0;
                for(int k = 0; k < m; ++k) {
                    Result.Matx[i][j] += Matx[i][k] * b.Matx[k][j];
                }
            }
        }
        num = incnt + ++outcnt;
        return Result;
    }
    Matrix operator*(const double c) const {
        Matrix Result = *this;

        for (auto& row : Result.Matx)
            for (double& value : row)
                value *= c;

        for (auto& row : Result.Matx_T)
            for (double& value : row)
                value *= c;

        return Result;
    }

    friend Matrix operator*(const double c, const Matrix& m) {
        return m * c;
    }
    Matrix operator^(const int &k) {
        Matrix Result;
        if(n != m) { printf("Error: Matrix is not square for exponentiation.\n"); Result.IS_legal = false; return Result; }

        Result.n = Result.m = n;
        Result.Matx.resize(Result.n, vector<double>(Result.m, 0));
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                Result.Matx[i][j] = (i == j ? 1 : 0);
            }
        }

        Matrix base = *this;
        int exp = k;
        while(exp > 0) {
            if(exp % 2 == 1) {
                Result = Result * base;
            }
            base = base * base;
            exp /= 2;
        }

        return Result;
    }
    Matrix Trans() {
        Matrix Result;

        Result.n = m; Result.m = n;
        Result.Matx = Matx_T; Result.Matx_T = Matx;

        return Result;
    }

    Matrix RowTransform(const int type, const int row1, const int row2 = 0, const double factor = 1.0) {
        
        Matrix Result = *this;

        if(!IS_legal || n <= 0 || m <= 0 || Matx.size() != static_cast<size_t>(n)) {
            printf("Error: Invalid matrix for row transformation.\n");
            Result.IS_legal = false;
            return Result;
        }
        for(const auto &row : Matx) {
            if(row.size() != static_cast<size_t>(m)) {
                printf("Error: Matrix rows have invalid dimensions.\n");
                Result.IS_legal = false;
                return Result;
            }
        }
        if(type < 1 || type > 3 || row1 < 1 || row1 > n || ((type == 1 || type == 3) && (row2 < 1 || row2 > n)) || (type == 2 && factor == 0) || (type == 3 && row1 == row2)) {
            printf("Error: Invalid elementary row transformation.\n");
            Result.IS_legal = false;
            return Result;
        }

        Matrix E = IdentityMatrix(n);
        
        int i = row1 - 1;
        int j = row2 - 1;
        if(type == 1) {
            E.Matx[i].swap(E.Matx[j]);     //swap row i, row j
        }
        else if(type == 2) {
            Matrix Row;
            Row.n = 1; Row.m = n;
            Row.n_T = n; Row.m_T = 1;
            Row.Matx.push_back(E.Matx[i]);
            Row.Matx_T.assign(n, vector<double>(1, 0));
            for(int c = 0; c < n; ++c)
                Row.Matx_T[c][0] = Row.Matx[0][c];
            Row = Row * factor;
            E.Matx[i] = Row.Matx[0];
        }
        else { E.Matx[i][j] = factor;}

        Result.Matx = (E * (*this)).Matx;
        Result.n_T = m; Result.m_T = n;
        Result.Matx_T.assign(m, vector<double>(n, 0));
        for(int r = 0; r < n; ++r)
            for(int c = 0; c < m; ++c)
                Result.Matx_T[c][r] = Result.Matx[r][c];

        return Result;
    }

    Matrix LU_Decomposition() {
        Matrix L, U;    
        //if(n != m) { printf("Error: Matrix is not square for LU decomposition.\n"); L.IS_legal = false; U.IS_legal = false; return pair<Matrix, Matrix>(L, U); }
        L.n = n; L.m = m;
        U.n = n; U.m = m;
        L.Matx.resize(n, vector<double>(m, 0));
        U = *this;
        for(int i = 0; i < n; i++ ) {
            L.Matx[i][i] = 1;
            int tot = U.Matx[i][0];
            for(int j = i + 1; j < n; j++ ) {
                L.Matx[j][i] = U.Matx[j][i] / tot;
                U.RowTransform(3, j, i, -L.Matx[j][i]);   
        }
        
        return pair<Matrix, Matrix>(L, U);
    }
    Matrix Gauss() {
        Matrix Result;
        Result.n = n; Result.m = m;
        Result.Matx = Matx; Result.Matx_T = Matx_T;
        return Gauss.LU_Decomposition().second;
    }
    Matrix Inverse() {
        Matrix Result, X, Y, L, U;
        if(n != m) { printf("Error: Matrix is not square for inversion.\n"); Result.IS_legal = false; return Result; }
        pair<Matrix, Matrix> LU = LU_Decomposition();
        L = LU.first; U = LU.second;
        Matrix I = IdentityMatrix(n);
        

        
        return Result;
    }
    Matrix Cofactor() {
        if (n == 1 && m == 1) return *this[0][0];
        Matrix cofactor = *this;
        if(n != m) { printf("Error: Matrix is not square for determinant calculation.\n"); Result.IS_legal = false; return Result; }
        

        return Result;
    }
};
vector<Matrix> M;

Matrix IdentityMatrix(const int Size) {
    Matrix Result;
    Result.n = Result.m = Size;
    Result.Matx.resize(Size, vector<double>(Size, 0));
    for(int i = 0; i < Size; ++i) { Result.Matx[i][i] = 1;}
    Result.num = incnt + ++outcnt;
    return Result;
}

Matrix cal(int num1, char c, int num2) {
    Matrix ans;
    if(c == '^') {
        Matrix &A = M[num1 - 1];
        ans = A ^ num2;
        ans.num = incnt + ++outcnt;
        ans.PrintMatrix();
        return ans;
    }
    if(c == '*') {
        Matrix &B = M[num2 - 1];
        ans = num1 * B;
        ans.num = incnt + ++outcnt;
        ans.PrintMatrix();
        return ans;
    }
    if(num1 <= 0 || num1 > M.size() || num2 <= 0 || num2 > M.size()) {
        printf("Error: Invalid matrix number for calculation.\n");
        ans.IS_legal = false;
        return ans;
    }
    Matrix &A = M[num1 - 1];
    Matrix &B = M[num2 - 1];

    switch(c) {
        case '+':
            ans = A + B;
            ans.num = incnt + ++outcnt;
            break;
        case '-':
            ans = A - B;
            ans.num = incnt + ++outcnt;
            break;
        case '@':
            ans = A * B;
            ans.num = incnt + ++outcnt;
            break;
        default:
            printf("Error: Unsupported operation '%c'.\n", c);
            ans.IS_legal = false;
            break;
    }
    ans.PrintMatrix();
    return ans;
}

Command Cmd() {
    Command kmd;
    kmd.EndJdg = false; kmd.TransposeJdg = false; kmd.InverseJdg = false; kmd.DeterminantJdg = false; kmd.RankJdg = false;
    string s, Name;
    cin >> s;
    while(1) {
        if(s == "end") { kmd.EndJdg = true; break;}
        else {
            if (s == "transpose") { kmd.TransposeJdg = true;}
            else if(s == "inverse") { kmd.InverseJdg = true;}
            else if(s == "determinant") { kmd.DeterminantJdg = true;}
            else if(s == "rank") { kmd.RankJdg = true;}
            else if(s == "read") { kmd.Read_Begin = true;}
            else if(s == "name") { kmd.nameJdg = true;}
            else if(s == "cal") { kmd.calJdg = true;}
            else if(s == "finish") { kmd.Read_Finish = true;}
            else if(s == "rowOP") {kmd.rowOPJdg = true;}
            Conduct(kmd);
            break;
        }
    }
    return kmd;
}

void Conduct(Command &kmd) {
    if(kmd.Read_Begin) {
        Matrix NewM;
        NewM.ReadMatrix();
        M.push_back(NewM);
        printf("Matrix number %d has been read.\n", NewM.num);
    }
    else if(kmd.nameJdg) {
        // Code for naming the matrix

    }
    else if(kmd.TransposeJdg) {
        int Num;
        scanf("%d", &Num);
        Matrix T = M[Num - 1].Trans();
        T.PrintMatrix();
    }
    else if(kmd.InverseJdg) {
        //scanf("%d", &num);
    }
    else if(kmd.DeterminantJdg) {
        // Code for determinant operation
    }
    else if(kmd.RankJdg) {
        // Code for rank operation
    }
    else if(kmd.calJdg) {
        char c;
        //Matrix rst;
        int num1, num2;
        scanf("%d %c %d", &num1, &c, &num2);
        M.push_back(cal(num1, c, num2));
    }
}


int main() {
    char c;
    Command CmdRcoder;
    while(!IS_End) {
        c = getchar();
        if(c == '\\') { CmdRcoder = Cmd();}
    }
    return 0;
}
