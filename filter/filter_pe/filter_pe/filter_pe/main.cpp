#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>


using namespace std;

bool filter_dup() {

}

int main(int argc, char** argv)
{
	string i1, i2, o1, o2;
	string a1, a2;
	bool dedup = false;
	int ns = 10;
	int quality = 55; double limit = 0.2;
	int mismatch = 3, align = 15;
	int start = -1, end = -1;
	bool valid = true;

	int i;
	const char* optstring = "1:2:3:4:5:6:dn:q:l:m:a:s:e:";
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
		case 'q':
			quality = atoi(optarg);
			break;
		case 'l':
			quality = atof(optarg);
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
		case '?':
			valid = false;
			break;
		}

	}

	if (i1.empty() || i2.empty()) {
		printf("Input file not specified.\n");
		valid = false;
	}

	if (o1.empty()) {
		o1 = i1;
		o1.append(".filtered");
	}

	if (o2.empty()) {
		o2 = i2;
		o2.append(".filtered");
	}

	if (start != -1 || end != -1)
		if (start < -1 || start >= end || end < -1) {
			printf("Invalid position statement with start %i and end %i. Both start and end can't be lower than -1, or start can't come latter than end.\n", start, end);
			valid = false;
		}

	if (quality < 0) {
		printf("Invalid quality valve %i. Can't specify a valve lower than 0.\n", quality);
		valid = false;
	}

	if (limit < 0 || limit > 1) {
		printf("Invalid limit %g. It should be a value between 0 and 1.\n", limit);
		valid = false;
	}

	int cal_start = start;
	int cal_length = end == -1 ? -1 : end - start + 1;

	if (!valid) {
		printf("Error(s) found, exiting.");
		return -1;
	}



	return 0;
}