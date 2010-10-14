/*
 * parseTestConfig.C
 *
 *  Created on: Jun 10, 2010
 *      Author: anshu
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>

using namespace std;

void findTests(ifstream& configfile, string& line);
void parseTests(ifstream& configfile, string& line);
void addData(string& line);
istream& mygetline ( istream& is, string& str );
void mygetlineCheck ( istream& is, string& str );
bool reservedWordFound(string& line);
void writeConfigFiles(int testNum, int keywordNum, int val);
string trimleft(string& str);

int testCount = 0;
int lineNum = 0;
int keywordCount = 0;

//directory found in test definition
bool dirFound = false;
//END found that closes a test
bool endFound = false;

//holds directory path of each test
vector<string> dirs;

//dimension 1: tests
//dimension 2: keywords
//dimension 3: values (1st element is keyword name)
vector<vector<vector<string> > > data;

int main()
{
	ifstream configfile;
	configfile.open ("testConfig.txt");
	if(configfile.is_open())
	{
		string line;

		while(mygetline(configfile,line))
			findTests(configfile, line);

		configfile.close();
	}
	else
	{
		cerr << ("test configuration file not found or unable to open\n");
		exit(1);
	}
}

void findTests(ifstream& configfile, string& line)
{
	//if an ON test is found, begin parsing it
	if(line.find("TEST on")!=string::npos)
	{
		cout << "Found an ON test" << endl;
		testCount++;
		data.resize(testCount);

		//Until the end of the test is not found parse the test parameters
		mygetlineCheck(configfile, line);
		while(line.find("END")==string::npos)
		{
			//if beginning of another test or EOF is not found, keep parsing
			if(line.find("TEST")==string::npos)
				parseTests(configfile, line);
			else
			{
				cerr << "expected END near line "<< lineNum << endl;
				exit(1);
			}
			mygetlineCheck(configfile, line);
		}
		cout << "Found END of test on line " << lineNum << endl;

		if(dirFound == true)
		{
			//TODO: check if same number of variables for each keyword
			int keywords = data[testCount-1].size();
			int vals = data[testCount-1][0].size();
			cout << "Found " << keywords << " keywords and " << vals
				 << " values" << endl;

			for(int i=1; i<keywords; i++)
			{
				if(data[testCount-1][i].size() != vals)
				{
					cerr << "Different number of keyword values detected "
						 << "for test ending on line " << lineNum << endl;
					exit(1);
				}
			}

			//Write config files and test directories to a file
			for(int i=1; i<vals; i++)
				writeConfigFiles(testCount,keywords,i);

			dirFound = false;
			keywordCount = 0;

			cout << "Finished writing a test" << endl;
		}
		else
		{
			cerr << "Directory not defined for test ending on line "
				 << lineNum << " - test will be ignored" << endl;
			exit(1);
		}

	}
}

void parseTests(ifstream& configfile, string& line)
{
	if(line.find("dir")!=string::npos)
	{
		mygetlineCheck(configfile,line);
		cout << "Found directory location: " << line << endl;
		dirs.push_back(trimleft(line));
		dirFound = true;
	}

	if(line.find("\\")!=string::npos)
	{
		cout << "Found keyword: " << line << endl;
		keywordCount++;
		data[testCount-1].resize(keywordCount);
		addData(line);

		mygetlineCheck(configfile, line);
		while(line.find("end")==string::npos)
		{
			if(!reservedWordFound(line))
			{
				cout << "     value = " << line << endl;
				addData(line);
			}
			else
			{
				cerr<<"expected \"end\" for keyword near line "<< lineNum
						<< endl;
				//exit(1);
			}
			mygetlineCheck(configfile, line);
		}
	}
}

void addData(string& line)
{
	data[testCount-1][keywordCount-1].push_back(trimleft(line));
}

istream& mygetline ( istream& is, string& str )
{
	lineNum++;
	return getline(is,str);
}

void mygetlineCheck ( istream& is, string& str )
{
	lineNum++;
	if(!getline(is,str))
	{
		cerr<< "unexpected EOF or read error on line " << lineNum << endl;
		exit(1);
	}
}

bool reservedWordFound(string& line)
{
	if( line.find("TEST")!=line.string::npos ||
		line.find("\\")!=string::npos ||
		line.find("dir")!=string::npos ||
		line.find("END")!=string::npos)
		return true;
	else
		return false;
}

void writeConfigFiles(int testNum, int keywordNum, int val)
{
	//Read in the config file
	string buffer;
	stringstream input;
	input << dirs[testNum-1] << "cpaimd_config.p1";
	string infilename = input.str();
	ifstream testConfigfile;
	testConfigfile.open(infilename.c_str());
	if(testConfigfile.is_open())
	{
		string inputLine;
		while(getline(testConfigfile,inputLine))
		{
			buffer.append(inputLine);
			buffer.push_back('\n');
		}
		testConfigfile.close();
		cout << "Finished reading in \"" << infilename << "\"" << endl;
	}
	else
	{
		cerr << "Could not open or find \"" << infilename << "\"" << endl;
		exit(1);
	}

	//output a config file with extra parameters added
	stringstream out;
	out << "tmp/test" << testNum << "-" << val - 1 << ".config";
	string outfilename = out.str();
	ofstream configfile;
	configfile.open(outfilename.c_str());
	if(configfile.is_open())
	{
		configfile << buffer
				   << "=====================================================\n"
				   << "~charm_conf_pc_def[" << endl;


		for(int i=0; i<keywordNum; i++)
		{
				configfile << data[testNum-1][i][0] << "{"
						   << data[testNum-1][i][val] << "}" << endl;
		}

		configfile << "]" << endl;

		configfile.close();
		cout << "Finished writing \"" << outfilename << "\"" << endl;
	}
	else
	{
		cerr << "Could not make \"" << outfilename << "\"" << endl;
		exit(1);
	}

	//output the directory where the reference data lives
	stringstream dirout;
	dirout << "tmp/test" << testNum << "-" << val - 1 << ".dir";
	string dirfilename = dirout.str();
	ofstream dirfile;
	dirfile.open(dirfilename.c_str());
	if(dirfile.is_open())
	{
		dirfile << dirs[testNum-1] << endl;
		dirfile.close();
	}
	else
	{
		cerr << "Could not make \"" << dirfilename << "\"" << endl;
		exit(1);
	}

	//output the test names
	stringstream nameout;
	nameout << "tmp/test" << testNum << "-" << val - 1 << ".name";
	string namefilename = nameout.str();
	ofstream namefile;
	namefile.open(namefilename.c_str());
	if(namefile.is_open())
	{
		namefile << "testlist+= " << "test" << testNum << "-" << val - 1 << endl;
		namefile.close();
	}
	else
	{
		cerr << "Could not make \"" << namefilename << "\"" << endl;
		exit(1);
	}
}

string trimleft(string& str)
{
	size_t startpos = str.find_first_not_of(" \t");
	return str.substr( startpos );
}
