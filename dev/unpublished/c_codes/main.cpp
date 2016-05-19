#include <stdio.h>
#include <string.h>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <stdlib.h>
#include <time.h>

/**
 * @brief Fast string tokanizer. (thread safe too)
 */
const char *toksplit(
	       const char *src, /* Source of tokens */
	       char tokchar, /* token delimiting char */
	       char *token, /* receiver of parsed token */
	       size_t lgh) /* length token can receive */
/* not including final '\0' */
{
if (src) {
while (' ' == *src) src++;
while (*src && (tokchar != *src)) {
if (lgh) {
  *token++ = *src;
  --lgh;
}
src++;
}
if (*src && (tokchar == *src)) src++;
}
*token = '\0';
return src;
}

//Read buffer size
#define BUFFER_SIZE 30000

//Map for computing the median
//groupMap["Some_group"][column][0] will always be the sum!
//groupMap["Some_group"][column][1] will always be the count!
std::unordered_map<std::string,std::vector<std::vector<double>>> groupMap;

//Group and first data columns
int gCol = 0;
int firstCol = 0;
FILE* meanFile = NULL;
FILE* medianFile = NULL;

long entryCount = 0;
int colCount = 0;

/**
 * Naive string to double
 */
double naive(const char *p) {
    double r = 0.0;
    bool neg = false;
    if (*p == '-') {
        neg = true;
        ++p;
    }
    while (*p >= '0' && *p <= '9') {
        r = (r*10.0) + (*p - '0');
        ++p;
    }
    if (*p == '.') {
        double f = 0.0;
        int n = 0;
        ++p;
        while (*p >= '0' && *p <= '9') {
            f = (f*10.0) + (*p - '0');
            ++p;
            ++n;
        }
        r += f / std::pow(10.0, n);
    }
    if (neg) {
        r = -r;
    }
    return r;
}


/**
 * @brief Parse lines out of buffered read. Way faster than fgets()
 */
int parseBuffer(char* buf, bool final=false) {
        
	//TODO all these hardcoded values have TO GO!
	//First, get each line of the buffer
	char *line = new char[10000];
	char *tok = new char[1000];
	const char *bufferTok = buf;
	const char *prevLine = NULL;

	std::string gName;
	//Get an iterator
	std::unordered_map<std::string,std::vector<std::vector<double>>>::iterator it;
	double num = 0.0;

	while (*bufferTok) {//We have each line
		prevLine = bufferTok;
		if (!final && !strrchr(bufferTok,'\n')) {
			break;
		}
		bufferTok = toksplit(bufferTok,'\n',line,10000);
		//Break up line
		const char * pch = line;
		++entryCount;
		int k = 0;
		while (*pch)
		{
			pch = toksplit(pch,'\t',tok,1000);
			if (k==gCol) {//Group name
				gName = tok;
				//Find or add this group
				if (groupMap.find(gName)==groupMap.end()) {
					groupMap[gName] = std::vector<std::vector<double>>();
				}
				it = groupMap.find(gName);

			} else if (k >= firstCol) {
				if (it->second.size() <= (k-firstCol)) {
					num = atof(tok);
					//num = naive(tok);

					it->second.push_back(std::vector<double>(1,num));
					it->second[k-firstCol].push_back(1);
					it->second[k-firstCol].push_back(num);
				} else {
					num = atof(tok);
					//num = naive(tok);					
					it->second[k-firstCol][0] += num;
					++it->second[k-firstCol][1];	
					if (medianFile) {
						try {
							it->second[k-firstCol].push_back(num);
						} catch(const std::bad_alloc& e) {
							fprintf(stderr,"Out of space, trying a delete\n");
							return 0;
						}
					}
				}
			}
			++k;
		}
	}
	delete[] tok;
	//Check that the last line had a proper ending!
	//Otherwise must have been a truncated line, so remove the added vector and return rewind value
	if (prevLine != NULL && !strrchr(prevLine,'\n')) {
		int rev = strlen(prevLine);
		delete[] line;
		return -rev;		
	}
	delete[] line;
	return 0;
}

/**
 * @brief Dump the current map to relevant files
 */
void dumpMap() {
	//Get an iterator
	std::unordered_map<std::string,std::vector<std::vector<double>>>::iterator it;
	//Print sum
	for (it = groupMap.begin(); it != groupMap.end(); ++it) {
		std::string line = it->first;
		std::string meanLine = it->first;
		std::string medianLine = it->first;
		double median = 0;
		for (int j=0; j<colCount; ++j) {
			if (medianFile) {//Sort!
				std::vector<double>::iterator start = it->second[j].begin();
				++start;++start;
				sort(start,it->second[j].end());
				//Get median, ignoring 0's
				int move = 0;
				for (std::vector<double>::iterator vIt = start; vIt != it->second[j].end(); ++vIt, ++move) {
					if ((*vIt) != 0)
						break;				
				}
				int len = it->second[j][1]-move;
				int middle = ((it->second[j][1]-move)/2)+2/*This is the skipping due to the "special" use of [0] and [1]*/+move;
				if ((len % 2 == 1) || (len < 2))
					median = it->second[j][middle];
				else
					median = (it->second[j][middle]+it->second[j][middle+1])/2;
				
				medianLine += '\t'+std::to_string(median);
			}
			line += '\t'+std::to_string(it->second[j][0]);
			meanLine += '\t'+std::to_string(it->second[j][0]/it->second[j][1]);
		}
		fprintf(stdout,"%s\n",line.c_str());
		if (meanFile) 
			fprintf(meanFile,"%s\n",meanLine.c_str());
		if (medianFile)
			fprintf(medianFile,"%s\n",medianLine.c_str());
		
	}
}

/**
 * @brief Usage text
 */
std::string usage() {
	return("\nUsage:\n groupByColumn [-g group_column] [-f first_data_column] [-h header_rows]\n [-m median_file] [-a mean_file] <input_file>\n\
\n\
 -g                      non-unique values in this column will be grouped and rows in data rows summed up\n\
 -f                      first column with data\n\
 -m                      compute median and output to <median_file>\n\
 -a                      compute mean and output to <mean_file>\n\
 -h                      number of rows that should be ignored for grouping [1]\n");
}

/**
 * @brief Entry point
 */
int main(int argc, char* argv[]) {

	int hRows = 1;
	int setFlags = 0;
	//Parse command stuff
	opterr = 0;
    int c;
	while ((c = getopt (argc, argv, "g:f:h:m:a:")) != -1)
		switch (c)
		{
		case 'g':
			gCol = atoi(optarg);
			++setFlags;
		   	break;
		case 'm':
			medianFile = fopen(optarg,"w");
			if (medianFile == NULL) {
				fprintf(stderr,"Failed to open file: %s\n",optarg);
				return -1;
			}
			break;
		case 'a':
			meanFile = fopen(optarg,"w");
			if (meanFile == NULL) {
				fprintf(stderr,"Failed to open file: %s\n",optarg);
				return -1;
			}
			break;
		case 'f':
			firstCol = atoi(optarg);
		    ++setFlags;
		    break;
		case 'h':
			hRows = atoi(optarg);
			++setFlags;
		}

	if ((setFlags < 2) || (argc == optind)) {
		fprintf(stderr,"%s",usage().c_str());
		return -1;
	}

	//TODO usage text and all that.
	FILE *in = fopen(argv[optind],"r");
	if (in == NULL) {
		fprintf(stderr,"Could not open file: %s\n",argv[optind]);
		return -1;
	}
	//Create a decent buffer for quickly reading that file
	char* buffer = new char[BUFFER_SIZE];
	//Drop header lines
	char *tok = new char[1000];
	while(--hRows >= 0) {
		fgets(buffer,BUFFER_SIZE,in);
		std::string line;
		const char * pch = buffer;
		int k = 0;
		while (*pch)
		{
			pch = toksplit(pch,'\t',tok,1000);
			if (k==gCol) {//Group name
				line = tok;
			} else if (k >= firstCol) {
				line += '\t';
				line += tok;
			}			
			++k;			
		}
		if (colCount == 0)
			colCount = k-firstCol;

		fprintf(stdout,"%s",line.c_str());
		if (meanFile)
			fprintf(meanFile,"%s",line.c_str());
		if (medianFile)
			fprintf(medianFile,"%s",line.c_str());
	}
	delete[] tok;
	clock_t a = clock();

	while (true) {//Read entire file in chuncks
		size_t res = fread(buffer,1,BUFFER_SIZE,in);
		if (res != BUFFER_SIZE) {//Done reading. We won't have a \0 on the end of this!!! Damn you fread!
			buffer[res] = '\0';
			parseBuffer(buffer,true);
			//fprintf(stderr,"Done with reading file\nRead %ld entries\n",entryCount);
			break;
		}
		int rewind = parseBuffer(buffer);
		if (rewind > 0) {//Something went wrong, give up.	
			return -1;
		} else if (rewind<0) {
			int ret = fseek(in,rewind,SEEK_CUR);
			if (ret != 0) {
				fprintf(stderr,"Cannot seek in file\n");
			}
		}
	}	
	delete[] buffer;
	fprintf(stderr,"Time to read: %3.5f\n", ((float)(clock()-a))/CLOCKS_PER_SEC);

	clock_t t = clock();	

	dumpMap();

	fprintf(stderr,"Time to dump: %3.5f\n", ((float)(clock()-t))/CLOCKS_PER_SEC);
	if (meanFile != NULL) 
		fclose(meanFile);
	if (medianFile != NULL)
		fclose(medianFile);

	return 0;

}
