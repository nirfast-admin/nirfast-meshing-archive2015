#ifndef __FileOperation
#define __FileOperation

//#include <stdio.h>
#include <fstream>
//#include <iostream>
#include <iomanip>
//#include <stdlib.h>
#include <string>
#include <vector>
//#include <cstdlib>

#include "CPoint.h"

#define adjustf std::fixed << std::setprecision(10) << std::setw(16) << std::setfill(' ') << std::right

class CFileOperation{
public:
	CFileOperation();
	CFileOperation(std::string fn) { myfilename = fn; }
	virtual ~CFileOperation();


	static void readSPRHeader(std::string inputFile, short *row, short *column, short *slice, short *dimension, std::string& dataType);
	static void readACQPHeader(std::string inputFile, short *row, short *column, short *slice, short *dimension, std::string& dataType);
	static std::string MakeFilename(std::string fileTarg, std::string OldExt, std::string NewExt);
	static void OpenOutFile(std::ofstream &fp, std::string fn, std::string mode="binary");
	static void OpenInFile(std::ifstream &fp, std::string fn, std::string mode="binary");
	static void GetLine(std::ifstream &in_file, std::string &s);
	static void CleanString(std::string &s);
	static std::string GetArg(std::string& s, char ch=',');
	static std::string GetExtension(std::string s);
	std::vector<Point *> ReadXYZ(std::string fn, std::string delimeter=" ");


protected:

private:
	std::string myfilename;
};

#endif



