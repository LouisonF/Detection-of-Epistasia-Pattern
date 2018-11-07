/*
 * file_parsing.cpp
 *
 *  Created on: 7 nov. 2018
 *      Author: louison
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>



using namespace std;

int file_parsing(string file)
{
	 vector<string> file_content;

	ifstream input(file, ios::in);

	if(input)
	{
		string line;
		while(getline(input, line))
		{
			file_content.push_back(line);
			cout << line << endl;
		}
	}else
	{
		cout << "Can't open the file provided" << endl;
	}

	return 0;
}

int main()
{

	string file_name;
	cout << "Please enter an input file (will be deprecated when using of parameter will be effective) " << endl;
	cin >> file_name;

	file_parsing(file_name);


	return 0;
}
