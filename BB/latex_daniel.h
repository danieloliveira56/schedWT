#ifndef LATEXTABLE_H
#define LATEXTABLE_H
#include <vector>
#include <sstream>

class LatexTable
{
public:
	std::vector<std::string> header;
	std::vector<std::vector <std::string> > fields;
	std::string caption;
	std::string label;

	LatexTable();
	LatexTable(std::vector<std::string> header, std::vector<std::vector <std::string> > fields);
	~LatexTable();
	std::string print(bool, bool);
	std::string printLineSeparator();
	void appendHeader(std::string);
	void appendField(std::string);
	void appendField(int);
	void appendField(double);
	void newLine();
	void setCaption(std::string);
	void setLabel(std::string);
	std::string printExcel(bool);
private:	
	std::string latexText;
	std::string xlText;
};

#endif // LATEXTABLE_H