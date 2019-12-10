/*
  filter_pe

  Copyright (c) 2019-2020 Henry Lee <2018301050@szu.edu.cn>.

  This file is part of MitoX.

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
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

/*
  filter_pe is only a runner, which is called by filter.py.
  This is NOT a complete tool, it only handles the filtering method,
  argument processing, help printing and others are done by filter.py.
  This is for seperating the real processing part from the calling part,
  which I think could extend the toolkit's flexibility while doing no
  harm to the toolkit integrity.
*/
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

  if (i1.empty() || i2.empty()){
    printf("Input missing. Two input files are needed.");
    valid = false;
  }

  if (o1.empty()){
    o1 = i1;
    o1.append(".filtered");
  }

  if (o2.empty()){
    o2 = i2;
    o2.append(".filtered");
  }

  if (!valid) {
    printf("Error(s) found, exiting.\n");
    return -1;
  }
  return 0;
}