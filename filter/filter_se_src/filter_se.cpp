/*
  filter_se.cpp

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

#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

using namespace std;

/*
  filter_se is only a runner, which is called by filter.py to filter data.
  This is NOT a complete tool, as it only handles the filtering method, other
  things like argument processing, help printing and others are done by
  filter.py. This is for seperating the real processing part from the calling
  part, which I think could extend the toolkit's flexibility while doing no harm
  to the toolkit integrity.
*/

bool has_suffix(const std::string& str, const std::string& suffix) {
  return str.size() >= suffix.size() &&
         str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline bool ns_check(char* bpseq, int ns_limit) {
  int ns = 0;
  for (int i = 0; i < strlen(bpseq); i++)
    if (bpseq[i] == 'N') ns++;
  return ns_limit > ns;
}

inline bool quality_check(char* quaseq, int quality, double limit) {
  int unqualified = 0;
  for (int i = 0; i < strlen(quaseq); i++)
    if (quality > quaseq[i]) unqualified++;

  return unqualified < int(strlen(quaseq) * limit);
}

inline bool filter_se(string input, string output, int start, int length,
                      int ns, int quality, double limit, int seq_c) {
  FILE *ifile = NULL, *ofile = NULL;

  char *head, *bpseq, *plus, *quaseq;
  head = (char*)malloc(1024);
  bpseq = (char*)malloc(1024);
  plus = (char*)malloc(1024);
  quaseq = (char*)malloc(1024);

  int sz = 0;

  ifile = has_suffix(input, string(".gz"))
              ? popen(string("gzip -dc ").append(input).c_str(), "r")
              : fopen(input.c_str(), "rb");

  ofile = has_suffix(output, string(".gz"))
              ? popen(string("gzip > ").append(output).c_str(), "w")
              : fopen(output.c_str(), "wb");

  if (!ifile) {
    printf("Unable to open input file %s.\n", input.c_str());
    return false;
  }

  if (!ofile) {
    printf("Unable to open output file %s.\n", output.c_str());
    return false;
  }

  while (fgets(head, 1024, ifile)) {
    fgets(bpseq, 1024, ifile);
    fgets(plus, 1024, ifile);
    fgets(quaseq, 1024, ifile);

    char *bps = bpseq, *quas = quaseq;

    strtok(bps, "\n");
    strtok(quas, "\n");

    if (start != -1) {
      bps += start;
      quas += start;
    }

    if (length != -1) {
      bps[length] = 0;
      quas[length] = 0;
    }

    if (!ns_check(bps, ns) || !quality_check(quas, quality, limit)) continue;

    fputs(head, ofile);
    fputs("\n", ofile);
    fputs(bps, ofile);
    fputs("\n", ofile);
    fputs(plus, ofile);
    fputs("\n", ofile);
    fputs(quaseq, ofile);
    fputs("\n", ofile);

    if (seq_c != -1) {
      sz++;
      if (sz >= seq_c) break;
    }
  }

  fcloseall();

  return true;
}

int main(int argc, char** argv) {
  string input, output;
  int start = -1, end = -1;
  int ns = 10;
  int quality = 55;
  double limit = 0.2;
  bool valid = true;
  int seq_c = -1;

  int i;
  const char* optstring = "i:o:s:e:n:q:l:";
  while ((i = getopt(argc, argv, optstring)) != -1) {
    switch (i) {
      case 'i':
        input = optarg;
        break;
      case 'o':
        output = optarg;
        break;
      case 's':
        start = atoi(optarg);
        break;
      case 'e':
        end = atoi(optarg);
        break;
      case 'q':
        quality = atoi(optarg);
        break;
      case 'l':
        limit = atof(optarg);
        break;
      case 'n':
        ns = atoi(optarg);
      case 'z':
        seq_c = atoi(optarg);
      case '?':
        valid = false;
        printf("Unrecognized argument.");
        break;
    }
  }

  // Test all the arguments
  if (input.empty()) {
    printf("Input file not specified.\n");
    valid = false;
  }

  if (output.empty()) {
    output = input + ".filtered";
    if (has_suffix(input, ".gz"))
      output.append(".gz");  // file.gz.filtered.gz...
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

  if (seq_c == 0) {
    printf(
        "Invalid sequence count %i. Output file should have sequences more "
        "than 1.\n",
        seq_c);
    valid = false;
  }

  int cal_start = start;
  int cal_length = end == -1 ? -1 : end - start + 1;

  if (!valid) {
    printf("Error(s) found, exiting.\n");
  } else if (!filter_se(input, output, cal_start, cal_length, ns, quality,
                        limit, seq_c)) {
    printf("Error occured when running filter function.\n");
  } else {
    return 0;
  }

  return -1;
}