/*
  filter_pe.cpp

  Copyright (c) 2019-2020 Henry Lee <2018301050@szu.edu.cn>.

  This file is part of source code of MitoX.

  MitoX is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  MitoX is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MitoX.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <unistd.h>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <set>
#include <string>

using namespace std;

/*
  filter_pe is only a runner, which is called by filter.py to filter data.
  This is NOT a complete tool, as it only handles the filtering method, other
  things like argument processing, help printing and others are done by
  filter.py. This is for seperating the real processing part from the calling
  part, which I think could extend the toolkit's flexibility while doing no harm
  to the toolkit integrity.
*/

inline bool deny_ns(string seq, int ns) {
  return count(seq.begin(), seq.end(), 'N') >= ns;
}

inline bool deny_quality(string quaseq, int quality, double limit) {
  return count_if(quaseq.begin(), quaseq.end(), [quality](char c) {
           return c >= quality;
         }) >= quaseq.length() * limit;
}

bool has_suffix(const std::string& str, const std::string& suffix) {
  return str.size() >= suffix.size() &&
         str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline bool filter_pe(string input1, string input2, string output1,
                      string output2, string adapter1, string adapter2,
                      bool dedup, int mismatch, int align, int ns, int quality,
                      int start, int length, double limit) {
  set<string> adapters;
  if (!adapter1.empty()) {
    ifstream a1;
    string temp;
    a1.open(adapter1);
    if (!a1.is_open()) {
      printf("Unable to open adapter list file 1 %s.\n", adapter1.c_str());
      return false;
    }
    while (!a1.eof()) {
      getline(a1, temp);
      adapters.insert(temp);
    }
  }
  if (!adapter2.empty()) {
    ifstream a2;
    string temp;
    a2.open(adapter2);
    if (!a2.is_open()) {
      printf("Unable to open adapter list file 2 %s.\n", adapter2.c_str());
      return false;
    }
    while (!a2.eof()) {
      getline(a2, temp);
      adapters.insert(temp);
    }
  }
  bool filter_adpater = adapters.empty();

  FILE *ifile1 = NULL, *ifile2 = NULL;
  FILE *ofile1 = NULL, *ofile2 = NULL;

  ifile1 = has_suffix(input1, ".gz")
               ? popen(string("gzip -dc ").append(input1).c_str(), "rb")
               : fopen(input1.c_str(), "rb");

  ifile2 = has_suffix(input2, ".gz")
               ? popen(string("gzip -dc ").append(input2).c_str(), "rb")
               : fopen(input2.c_str(), "rb");

  ofile1 = has_suffix(output1, ".gz")
               ? popen(("gzip -dc " + (output1)).c_str(), "wb")
               : fopen(output1.c_str(), "rb");

  ofile2 = has_suffix(output2, ".gz")
               ? popen(("gzip -dc " + (output2)).c_str(), "wb")
               : fopen(output2.c_str(), "rb");

  if (ifile1 == NULL || ofile1 == NULL || ifile2 == NULL || ofile2 == NULL) {
    printf("Unable to open file.");
    return false;
  }

  string head1, head2, ns1, ns2, plus1, plus2, qua1, qua2;
  char temp[1024];
  while (fgets(temp, 1024, ifile1)) {
    head1 = temp;
    fgets(temp, 1024, ifile1);
    ns1 = temp;
    fgets(temp, 1024, ifile1);
    plus1 = temp;
    fgets(temp, 1024, ifile1);
    qua1 = temp;

    fgets(temp, 1024, ifile2);
    head2 = temp;
    fgets(temp, 1024, ifile2);
    ns2 = temp;
    fgets(temp, 1024, ifile2);
    plus2 = temp;
    fgets(temp, 1024, ifile2);
    qua2 = temp;

    head1.pop_back();
    head2.pop_back();
    plus1.pop_back();
    plus2.pop_back();
    qua1.pop_back();
    qua2.pop_back();
    ns1.pop_back();
    ns2.pop_back();

    if (filter_adpater && (adapters.find(ns1) != adapters.end() ||
                           adapters.find(ns2) != adapters.end()))
      continue;

    if (start != -1) {
      ns1.erase(0, start);
      ns2.erase(0, start);
      qua1.erase(0, start);
      qua2.erase(0, start);
    }
    if (length != -1) {
      ns1.erase(length);
      ns2.erase(length);
      qua1.erase(length);
      qua2.erase(length);
    }

    if (deny_ns(ns1, ns) || deny_ns(ns2, ns) ||
        deny_quality(qua1, quality, limit) ||
        deny_quality(qua2, quality, limit))
      continue;

    fputs(head1.c_str(), ofile1);
    fputs("\n", ofile1);
    fputs(head2.c_str(), ofile2);
    fputs("\n", ofile2);
    fputs(ns1.c_str(), ofile1);
    fputs("\n", ofile1);
    fputs(ns2.c_str(), ofile2);
    fputs("\n", ofile2);
    fputs(plus1.c_str(), ofile1);
    fputs("\n", ofile1);
    fputs(plus2.c_str(), ofile1);
    fputs("\n", ofile2);
    fputs(head1.c_str(), ofile1);
    fputs("\n", ofile1);
    fputs(head1.c_str(), ofile1);
    fputs("\n", ofile2);
  }

  return true;
}

int main(int argc, char** argv) {
  string i1, i2, o1, o2;
  string a1, a2;
  bool dedup;
  int mismatch = 3, align = 15;
  int start = -1, end = -1;
  int ns = 10;
  int quality = 55;
  double limit = 0.2;

  bool valid = true;
  int i;
  const char* optstring = "1:2:3:4:5:6:dm:a:s:e:n:q:l:";
  while ((i = getopt(argc, argv, optstring)) != -1) {
    switch (i) {
      case '1':
        i1 = optarg;
        break;
      case '2':
        i2 = optarg;
        break;
      case '3':
        o1 = optarg;
        break;
      case '4':
        o2 = optarg;
        break;
      case '5':
        a1 = optarg;
        break;
      case '6':
        a2 = optarg;
        break;
      case 'd':
        dedup = true;
        break;
      case 'm':
        mismatch = atoi(optarg);
        break;
      case 'a':
        align = atoi(optarg);
        break;
      case 's':
        start = atoi(optarg);
        break;
      case 'e':
        end = atoi(optarg);
        break;
      case 'n':
        ns = atoi(optarg);
        break;
      case 'q':
        quality = atoi(optarg);
        break;
      case 'l':
        limit = atof(optarg);
        break;
      case '?':
        valid = false;
        break;
    }
  }

  if (i1.empty() || i2.empty()) {
    printf("Input missing. Two input files are needed.\n");
    valid = false;
  }

  if (o1.empty()) {
    o1 = i1 + ".filtered";
    if (has_suffix(i1, ".gz")) o1.append(".gz");
  }

  if (o2.empty()) {
    o2 = i2 + ".filtered";
    if (has_suffix(i2, ".gz")) o2.append(".gz");
  }

  if (start != -1 || end != -1)
    if (start < -1 || start >= end || end < -1) {
      printf(
          "Invalid position statement with start %i and end %i. Both start and "
          "end can't be lower than -1, or start can't come latter than end.\n",
          start, end);
      valid = false;
    }

  if (quality < 0) {
    printf("Invalid quality valve %i. Can't specify a valve lower than 0.\n",
           quality);
    valid = false;
  }

  if (limit < 0 || limit > 1) {
    printf("Invalid limit %g. It should be a value between 0 and 1.\n", limit);
    valid = false;
  }

  int cal_start = start;
  int cal_length = end == -1 ? -1 : end - start + 1;

  if (!valid) {
    printf("Error(s) found in arugments, exiting.\n");
    return -1;
  } else if (!filter_pe(i1, i2, o1, o2, a1, a2, dedup, mismatch, align, ns,
                        quality, cal_start, cal_length, limit)) {
    printf("Error occured in running filter method.\n");
    return -1;
  }
  return 0;
}