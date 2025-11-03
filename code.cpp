// code.cpp (generated for OJ submission)
#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// Do not use "using namespace std;"

namespace sjtu {
class int2048 {
private:
  // sign: false for positive or zero, true for negative
  bool neg;
  // little-endian digits in base BASE (each element 0..BASE-1)
  std::vector<int> d;
  static const int BASE = 1000; // 10^3
  static const int BASE_DIGS = 3;

  // helpers
  void trim();
  static int absCmp(const int2048 &a, const int2048 &b);
  static int2048 absAdd(const int2048 &a, const int2048 &b);
  static int2048 absSub(const int2048 &a, const int2048 &b); // assume |a|>=|b|
  static void mulFFT(const std::vector<int> &a, const std::vector<int> &b, std::vector<int> &res);
  static void divModAbs(const int2048 &a, const int2048 &b, int2048 &q, int2048 &r);

public:
  // Constructors
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &);
  // Output the stored big integer, no need for newline
  void print();

  // Add a big integer
  int2048 &add(const int2048 &);
  // Return the sum of two big integers
  friend int2048 add(int2048, const int2048 &);

  // Subtract a big integer
  int2048 &minus(const int2048 &);
  // Return the difference of two big integers
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

#endif



namespace sjtu {

// ============ helpers ============
void int2048::trim() {
  while (!d.empty() && d.back() == 0) d.pop_back();
  if (d.empty()) neg = false;
}

int int2048::absCmp(const int2048 &a, const int2048 &b) {
  if (a.d.size() != b.d.size()) return (a.d.size() < b.d.size()) ? -1 : 1;
  for (int i = (int)a.d.size() - 1; i >= 0; --i) {
    if (a.d[i] != b.d[i]) return (a.d[i] < b.d[i]) ? -1 : 1;
  }
  return 0;
}

int2048 int2048::absAdd(const int2048 &a, const int2048 &b) {
  int2048 res;
  res.neg = false;
  const int n = (int)a.d.size(), m = (int)b.d.size();
  int len = n > m ? n : m;
  res.d.assign(len, 0);
  long long carry = 0;
  for (int i = 0; i < len; ++i) {
    long long cur = carry;
    if (i < n) cur += a.d[i];
    if (i < m) cur += b.d[i];
    res.d[i] = (int)(cur % BASE);
    carry = cur / BASE;
  }
  if (carry) res.d.push_back((int)carry);
  return res;
}

int2048 int2048::absSub(const int2048 &a, const int2048 &b) { // assume |a|>=|b|
  int2048 res;
  res.neg = false;
  const int n = (int)a.d.size(), m = (int)b.d.size();
  res.d.assign(n, 0);
  long long borrow = 0;
  for (int i = 0; i < n; ++i) {
    long long cur = (long long)a.d[i] - (i < m ? b.d[i] : 0) - borrow;
    if (cur < 0) {
      cur += BASE;
      borrow = 1;
    } else borrow = 0;
    res.d[i] = (int)cur;
  }
  res.trim();
  return res;
}

static void fft(std::vector<std::complex<double>> &a, bool invert) {
  int n = (int)a.size();
  static std::vector<int> rev;
  static std::vector<std::complex<double>> roots{0, 1};

  if ((int)rev.size() != n) {
    int k = __builtin_ctz(n);
    rev.assign(n, 0);
    for (int i = 0; i < n; i++) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (k - 1));
  }
  if ((int)roots.size() < n) {
    int k = __builtin_ctz((int)roots.size());
    roots.resize(n);
    const double PI = 3.141592653589793238462643383279502884;
    while ((1 << k) < n) {
      double ang = 2 * PI / (1 << (k + 1));
      for (int i = 1 << (k - 1); i < (1 << k); ++i) {
        roots[2 * i] = roots[i];
        double ang_i = ang * (2 * i + 1 - (1 << k));
        roots[2 * i + 1] = std::complex<double>(std::cos(ang_i), std::sin(ang_i));
      }
      ++k;
    }
  }

  for (int i = 0; i < n; ++i) if (i < rev[i]) { auto tmp=a[i]; a[i]=a[rev[i]]; a[rev[i]]=tmp; }

  for (int len = 1; len < n; len <<= 1) {
    for (int i = 0; i < n; i += 2 * len) {
      for (int j = 0; j < len; ++j) {
        std::complex<double> u = a[i + j];
        std::complex<double> v = a[i + j + len] * roots[len + j];
        a[i + j] = u + v;
        a[i + j + len] = u - v;
      }
    }
  }
  if (invert) {
    // reverse a[1..n-1]
    for (int l = 1, r = n - 1; l < r; ++l, --r) { auto tmp=a[l]; a[l]=a[r]; a[r]=tmp; }
    for (int i = 0; i < n; ++i) a[i] /= n;
  }
}

void int2048::mulFFT(const std::vector<int> &A, const std::vector<int> &B, std::vector<int> &res) {
  int n1 = (int)A.size(), n2 = (int)B.size();
  if (!n1 || !n2) { res.clear(); return; }
  // small threshold
  if ((long long)n1 * n2 <= 8000) {
    std::vector<long long> tmp(n1 + n2, 0);
    for (int i = 0; i < n1; ++i) {
      long long ai = A[i];
      for (int j = 0; j < n2; ++j) tmp[i + j] += ai * B[j];
    }
    res.assign(n1 + n2, 0);
    long long carry = 0;
    for (size_t i = 0; i < tmp.size(); ++i) {
      long long cur = tmp[i] + carry;
      res[i] = (int)(cur % BASE);
      carry = cur / BASE;
    }
    while (carry) { res.push_back((int)(carry % BASE)); carry /= BASE; }
    while (!res.empty() && res.back() == 0) res.pop_back();
    return;
  }
  int n = 1;
  while (n < n1 + n2) n <<= 1;
  std::vector<std::complex<double>> fa(n), fb(n);
  for (int i = 0; i < n1; ++i) fa[i] = std::complex<double>(A[i], 0);
  for (int i = 0; i < n2; ++i) fb[i] = std::complex<double>(B[i], 0);
  for (int i = n1; i < n; ++i) fa[i] = 0;
  for (int i = n2; i < n; ++i) fb[i] = 0;
  fft(fa, false); fft(fb, false);
  for (int i = 0; i < n; ++i) fa[i] *= fb[i];
  fft(fa, true);
  res.assign(n1 + n2, 0);
  long long carry = 0;
  for (int i = 0; i < n1 + n2; ++i) {
    double rr = fa[i].real();
    long long val = (long long)((rr >= 0 ? rr + 0.5 : rr - 0.5)) + carry;
    res[i] = (int)(val % BASE);
    carry = val / BASE;
  }
  while (carry) { res.push_back((int)(carry % BASE)); carry /= BASE; }
  while (!res.empty() && res.back() == 0) res.pop_back();
}

static void mulSmallVec(std::vector<int> &v, int m, int BASE) {
  long long carry = 0;
  for (size_t i = 0; i < v.size(); ++i) {
    long long cur = 1ll * v[i] * m + carry;
    v[i] = (int)(cur % BASE);
    carry = cur / BASE;
  }
  while (carry) { v.push_back((int)(carry % BASE)); carry /= BASE; }
  while (!v.empty() && v.back() == 0) v.pop_back();
}

static int divSmallVec(std::vector<int> &v, int m, int BASE) {
  long long rem = 0;
  for (int i = (int)v.size() - 1; i >= 0; --i) {
    long long cur = v[i] + rem * BASE;
    v[i] = (int)(cur / m);
    rem = cur % m;
  }
  while (!v.empty() && v.back() == 0) v.pop_back();
  return (int)rem;
}

void int2048::divModAbs(const int2048 &A, const int2048 &B, int2048 &Q, int2048 &R) {
  // assumes |B| > 0
  if (A.d.empty()) { Q.d.clear(); Q.neg=false; R.d.clear(); R.neg=false; return; }
  if (B.d.size() == 1) {
    int b = B.d[0];
    Q.d.assign(A.d.size(), 0);
    long long rem = 0;
    for (int i = (int)A.d.size() - 1; i >= 0; --i) {
      long long cur = A.d[i] + rem * BASE;
      Q.d[i] = (int)(cur / b);
      rem = cur % b;
    }
    Q.neg=false; Q.trim();
    R.d.clear(); if (rem) R.d.push_back((int)rem); R.neg=false; R.trim();
    return;
  }
  if (absCmp(A, B) < 0) { Q.d.clear(); Q.neg=false; R = A; R.neg=false; return; }

  std::vector<int> u = A.d;
  std::vector<int> v = B.d;
  int norm = BASE / (v.back() + 1);
  if (norm > 1) { mulSmallVec(u, norm, BASE); mulSmallVec(v, norm, BASE); }
  u.push_back(0);
  int n = (int)v.size();
  int m = (int)u.size() - n - 1; // since u has an extra leading 0 now
  Q.d.assign(m + 1, 0);
  for (int i = m; i >= 0; --i) {
    unsigned long long u2 = (unsigned long long)u[i + n] * (unsigned long long)BASE + (unsigned long long)u[i + n - 1];
    unsigned long long qhat = u2 / (unsigned long long)v[n - 1];
    if (qhat >= (unsigned long long)BASE) qhat = (unsigned long long)BASE - 1ULL;
    // refine qhat using next digit when available
    if (n > 1) {
      unsigned long long v2 = (unsigned long long)v[n - 1] * (unsigned long long)BASE + (unsigned long long)v[n - 2];
      unsigned long long u3 = (unsigned long long)u[i + n] * (unsigned long long)BASE * (unsigned long long)BASE + (unsigned long long)u[i + n - 1] * (unsigned long long)BASE + (unsigned long long)((i + n - 2 >= 0) ? u[i + n - 2] : 0ULL);
      while (qhat * v2 > u3) --qhat;
    }

    // subtract qhat * v from u at position i
    unsigned long long borrow = 0;
    for (int j = 0; j < n; ++j) {
      unsigned long long prod = (unsigned long long)qhat * (unsigned long long)v[j] + borrow;
      long long cur = (long long)u[i + j] - (long long)(prod % BASE);
      borrow = prod / BASE;
      if (cur < 0) { cur += BASE; ++borrow; }
      u[i + j] = (int)cur;
    }
    long long cur2 = (long long)u[i + n] - (long long)borrow;
    if (cur2 < 0) {
      // qhat was too big, fix by adding back v and decrementing qhat
      qhat -= 1ULL;
      unsigned long long carry = 0;
      for (int j = 0; j < n; ++j) {
        unsigned long long cur = (unsigned long long)u[i + j] + (unsigned long long)v[j] + carry;
        if (cur >= (unsigned long long)BASE) { u[i + j] = (int)(cur - BASE); carry = 1; } else { u[i + j] = (int)cur; carry = 0; }
      }
      u[i + n] = (int)((unsigned long long)u[i + n] + carry);
    } else {
      u[i + n] = (int)cur2;
    }
    Q.d[i] = (int)qhat;
  }
  // remainder in u[0..n-1]; unnormalize
  if (norm > 1) divSmallVec(u, norm, BASE);
  R.d.assign(u.begin(), u.begin() + n);
  R.neg = false; R.trim();
  Q.neg = false; Q.trim();
}

// ============ ctors ============
int2048::int2048() : neg(false) {}

int2048::int2048(long long v) { *this = int2048(); if (v < 0) { neg = true; v = -v; } else neg = false; while (v) { d.push_back((int)(v % BASE)); v /= BASE; } trim(); }

int2048::int2048(const std::string &s) { read(s); }

int2048::int2048(const int2048 &o) : neg(o.neg), d(o.d) {}

// ============ basic io ============
void int2048::read(const std::string &s) {
  neg = false; d.clear();
  int i = 0; while (i < (int)s.size() && (s[i] == ' ' || s[i] == '\n' || s[i] == '\t' || s[i] == '\r')) ++i;
  bool sign = false;
  if (i < (int)s.size() && (s[i] == '+' || s[i] == '-')) { sign = s[i] == '-'; ++i; }
  while (i < (int)s.size() && s[i] == '0') ++i; // skip leading zeros
  std::vector<int> tmp;
  for (int j = (int)s.size() - 1; j >= i; j -= BASE_DIGS) {
    int x = 0;
    int l = i > (j - BASE_DIGS + 1) ? i : (j - BASE_DIGS + 1);
    for (int k = l; k <= j; ++k) x = x * 10 + (s[k] - '0');
    tmp.push_back(x);
  }
  d.swap(tmp);
  trim();
  neg = sign && !d.empty();
}

void int2048::print() {
  if (d.empty()) { std::cout << 0; return; }
  if (neg) std::cout << '-';
  int n = (int)d.size();
  std::cout << d.back();
  char buf[8];
  for (int i = n - 2; i >= 0; --i) {
    std::snprintf(buf, sizeof(buf), "%0*d", BASE_DIGS, d[i]);
    std::cout << buf;
  }
}

// ============ basic add/minus ============
int2048 &int2048::add(const int2048 &b) {
  if (b.d.empty()) return *this;
  if (this->d.empty()) { *this = b; return *this; }
  if (this->neg == b.neg) {
    int2048 sum = absAdd(*this, b);
    sum.neg = this->neg;
    *this = sum;
  } else {
    int cmp = absCmp(*this, b);
    if (cmp == 0) { d.clear(); neg = false; }
    else if (cmp > 0) { int2048 diff = absSub(*this, b); diff.neg = this->neg; *this = diff; }
    else { int2048 diff = absSub(b, *this); diff.neg = b.neg; *this = diff; }
  }
  return *this;
}

int2048 add(int2048 a, const int2048 &b) { return a.add(b); }

int2048 &int2048::minus(const int2048 &b) {
  int2048 nb = b; if (!nb.d.empty()) nb.neg = !nb.neg; return this->add(nb);
}

int2048 minus(int2048 a, const int2048 &b) { return a.minus(b); }

// ============ operators ============
int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const { int2048 t(*this); if (!t.d.empty()) t.neg = !t.neg; return t; }

int2048 &int2048::operator=(const int2048 &o) { if (this == &o) return *this; neg = o.neg; d = o.d; return *this; }

int2048 &int2048::operator+=(const int2048 &o) { return add(o); }
int2048 operator+(int2048 a, const int2048 &b) { return a += b; }

int2048 &int2048::operator-=(const int2048 &o) { return minus(o); }
int2048 operator-(int2048 a, const int2048 &b) { return a -= b; }

int2048 &int2048::operator*=(const int2048 &o) {
  if (d.empty() || o.d.empty()) { d.clear(); neg = false; return *this; }
  std::vector<int> prod;
  mulFFT(this->d, o.d, prod);
  d.swap(prod);
  neg = (this->neg != o.neg) && !d.empty();
  trim();
  return *this;
}
int2048 operator*(int2048 a, const int2048 &b) { return a *= b; }

int2048 &int2048::operator/=(const int2048 &o) {
  // floor division per README
  if (o.d.empty()) return *this; // undefined, but won't happen
  if (this->d.empty()) { neg = false; return *this; }
  int2048 A = *this; A.neg = false;
  int2048 B = o; B.neg = false;
  int2048 q, r;
  divModAbs(A, B, q, r);
  // adjust signs for floor division
  bool same = (this->neg == o.neg);
  if (same) {
    *this = q; this->neg = false; // quotient non-negative when signs same
  } else {
    if (r.d.empty()) {
      *this = q; this->neg = true; // exact division, negative quotient
    } else {
      // q = -q - 1
      q.neg = true;
      int carry = 1;
      for (size_t i = 0; i < q.d.size(); ++i) {
        int cur = q.d[i] + carry;
        if (cur >= BASE) { q.d[i] = cur - BASE; carry = 1; } else { q.d[i] = cur; carry = 0; break; }
      }
      if (carry) q.d.push_back(1);
      q.trim();
      *this = q;
    }
  }
  return *this;
}
int2048 operator/(int2048 a, const int2048 &b) { return a /= b; }

int2048 &int2048::operator%=(const int2048 &o) {
  if (o.d.empty()) return *this; // undefined
  if (this->d.empty()) return *this;
  int2048 A = *this; A.neg = false;
  int2048 B = o; B.neg = false;
  int2048 q, r;
  divModAbs(A, B, q, r);
  // adjust remainder sign to match divisor; use definition r = x - (x/y)*y
  bool same = (this->neg == o.neg);
  if (same) {
    // remainder has same sign as divisor
    *this = r; this->neg = o.neg && !r.d.empty();
  } else {
    if (r.d.empty()) { d.clear(); neg = false; }
    else {
      // r = |B| - r_abs, sign same as divisor
      int2048 t;
      t.neg = false;
      // B >= r always
      t = absSub(B, r);
      *this = t; this->neg = o.neg;
    }
  }
  return *this;
}
int2048 operator%(int2048 a, const int2048 &b) { return a %= b; }

std::istream &operator>>(std::istream &is, int2048 &x) {
  std::string s; is >> s; x.read(s); return is;
}
std::ostream &operator<<(std::ostream &os, const int2048 &x) {
  if (x.d.empty()) { os << 0; return os; }
  if (x.neg) os << '-';
  os << x.d.back();
  char buf[8];
  for (int i = (int)x.d.size() - 2; i >= 0; --i) {
    std::snprintf(buf, sizeof(buf), "%0*d", int2048::BASE_DIGS, x.d[i]);
    os << buf;
  }
  return os;
}

bool operator==(const int2048 &a, const int2048 &b) {
  return a.neg == b.neg && a.d == b.d;
}
bool operator!=(const int2048 &a, const int2048 &b) { return !(a == b); }

bool operator<(const int2048 &a, const int2048 &b) {
  if (a.d.empty() && b.d.empty()) return false;
  if (a.neg != b.neg) return a.neg;
  int cmp = int2048::absCmp(a, b);
  return a.neg ? (cmp > 0) : (cmp < 0);
}

bool operator>(const int2048 &a, const int2048 &b) { return b < a; }
bool operator<=(const int2048 &a, const int2048 &b) { return !(b < a); }
bool operator>=(const int2048 &a, const int2048 &b) { return !(a < b); }

} // namespace sjtu
