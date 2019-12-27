/*
  dp_circle_check.cpp

  Copyright (c) 2019-2020 Li Junyu <2018301050@szu.edu.cn>.

  This file is part of source code of MitoFlex.

  MitoFlex is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  MitoFlex is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MitoFlex.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include <string>

using namespace std;

struct subset {
  int start;
  int length;
};

inline subset get_start(const char* s1, const char* s2) {
  // Init, precalculate the value of two string, avoid multiple call of strlen
  int n = strlen(s1) + 1, m = strlen(s2) + 1;
  int** dp = new int*[n];
  for (int i = 0; i < n; i++) dp[i] = new int[m];

  int maxv = 0, maxi = 0;

  // DP - here use 1 as lower bound, avoid the out-of-range detection when
  // acessing 0
  for (int i = 1; i < n; i++) {
    for (int j = 1; j < m; j++) {
      dp[i][j] = s1[i - 1] == s2[j - 1] ? dp[i - 1][j - 1] + 1 : 0;
      if (dp[i][j] > maxv) {
        maxv = dp[i][j];
        maxi = i;
      }
    }
  }

  subset result;
  result.start = maxi - maxv;
  result.length = maxv;

  return result;
}

int main() {
  string s1, s2;
  cin >> s1 >> s2;

  subset res = get_start(s1.c_str(), s2.c_str());
  printf("%d %d\n", res.start, res.length);
  printf("%.*s\n", res.length, s1.c_str() + res.start);
}
