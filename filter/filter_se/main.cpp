#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>


using namespace std;

bool has_suffix(const std::string& str, const std::string& suffix) {
	return str.size() >= suffix.size() &&
		str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline bool ns_check(char* bpseq, int ns_limit) {
	int ns = 0;
	for (int i = 0; i < strlen(bpseq); i++)
		if (bpseq[i] == 'N')
			ns++;
	return ns_limit > ns;
}

inline bool quality_check(char* quaseq, int quality, double limit) {
	int unqualified = 0;
	for (int i = 0; i < strlen(quaseq); i++)
		if (quality > quaseq[i])
			unqualified++;

	return unqualified < int(strlen(quaseq) * limit);
}

bool filter_se(string input, string output, int start, int length, int ns, int quality, double limit) {

	FILE* ifile = NULL, * ofile = NULL;

	char* head, * bpseq, * plus, * quaseq;
	head = (char*)malloc(1024);
	bpseq = (char*)malloc(1024);
	plus = (char*)malloc(1024);
	quaseq = (char*)malloc(1024);

	ifile = has_suffix(input, string(".gz")) ?
		popen(string("gzip -dc ").append(input).c_str(), "rb") :
		fopen(input.c_str(), "rb");

	ofile = has_suffix(output, string(".gz")) ?
		popen(string("gzip > ").append(output).c_str(), "wb") :
		fopen(output.c_str(), "wb");


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

		char* bps = bpseq, * quas = quaseq;

		if (start != -1) {
			bps += start;
			quas += start;
		}

		if (length != -1) {
			bps[length] = 0;
			quas[length] = 0;
		}

		if (!ns_check(bps, ns) || !quality_check(quas, quality, limit))
			continue;

		fputs(head, ofile);
		fputs(bps, ofile);
		fputs(plus, ofile);
		fputs(quaseq, ofile);
	}
	return true;
}

int main(int argc, char** argv) {

	string input, output;
	int start = -1, end = -1;
	int ns = 10;
	int quality = 55; double limit = 0.2;
	bool valid = true;

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
		case '?':
			valid = false;
			break;
		}
	}
	//Test all the arguments

	if (input.empty()) {
		printf("Input file not specified.\n");
		valid = false;
	}

	if (output.empty()) {
		output = input;
		output.append(".filtered");
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
		printf("Error(s) found, exiting.\n");
	}
	else if (!filter_se(input, output, cal_start, cal_length, ns, quality, limit)) {
		printf("Error occured when running filter function.\n");
	}
	else {
		return 0;
	}


	return -1;
}