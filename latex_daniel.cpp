#include "./BB/latex_daniel.h"
#include <vector>
#include <sstream>

LatexTable::LatexTable()
{
	latexText = "";
	caption = "__insert_caption__";
	label = "__insert_label__";
}

LatexTable::LatexTable(std::vector<std::string> header_, std::vector<std::vector <std::string> > fields_)
{
	header = header_;
	fields = fields_;
}

LatexTable::~LatexTable()
{
}

std::string LatexTable::printLineSeparator()
{
	return "\t\\midrule\n";
}

std::string LatexTable::print(bool printHeader = true, bool printFormat = true)
{
	if (printFormat)
	{
		latexText = "\\begin{longtable} {l|";
		for (int i = 1; i < header.size(); i++)
		{
			latexText += "r";
		}
		latexText += "}\n";

		latexText += "\t\\caption{"+ caption +"}  \\label{table:" + label + "} \\\\\n";
		latexText += "\t\\toprule\n";
	}

	if (printHeader)
	{
		latexText += "\t\t";
		for (int i = 0; i < header.size(); i++)
		{
			if (i > 0)
			{
				latexText += " & ";
			}
			latexText += header[i];
		}
		latexText += "\\\\\n";
	}

	if (printFormat) {
		latexText += "\t\\midrule\n";
		latexText += "\t\\endfirsthead\n";
		latexText += "\t\\midrule\n";
		latexText += "\t\\caption*{Table \\ref{table:"+ label +"} Continued: " + caption + "} \\\\\n";
		latexText += "\t\\midrule\n";
		latexText += "\t\\endhead\n";
	}

	for (int i = 0; i < fields.size(); i++)
	{
		latexText += "\t\t";
		for (int j = 0; j < fields[i].size(); j++)
		{
			if (j > 0)
			{
				latexText += " & ";
			}
			latexText += fields[i][j];
		}
		latexText += "\\\\\n";
	}
	
	if (printFormat)
	{
		latexText += "\t\\bottomrule\n";
		latexText += "\\end{longtable}\n";
	}


	return latexText;
}

void LatexTable::appendHeader(std::string textToAppend)
{
	header.push_back(textToAppend);
}

//fazer template para appendField
void LatexTable::appendField(std::string stringTextToAppend)
{
	//std::replace( stringTextToAppend.begin(), stringTextToAppend.end(), '_', '\\_'); // replace all 'x' to 'y'
	fields[fields.size()-1].push_back(stringTextToAppend);
}

void LatexTable::appendField(int intTextToAppend)
{
	std::stringstream  strs;
	strs << intTextToAppend;
	fields[fields.size()-1].push_back(strs.str() );
}

void LatexTable::appendField(double doubleTextToAppend)
{
	std::stringstream  strs;
	strs << doubleTextToAppend;
	fields[fields.size()-1].push_back(strs.str() );
}

void LatexTable::newLine()
{
	std::vector<std::string> emptyvector;
	fields.push_back(emptyvector);
}

void LatexTable::setCaption(std::string stringCaption)
{
	caption = stringCaption;
}

void LatexTable::setLabel(std::string stringLabel)
{
	label = stringLabel;
}

std::string LatexTable::printExcel(bool printHeader = true)
{
	if (printHeader)
	{
		for (int i = 0; i < header.size(); i++)
		{
			if (i > 0)
			{
				xlText += "\t";
			}
			xlText += header[i];
		}
		xlText += "\n";
	}

	for (int i = 0; i < fields.size(); i++)
	{
		for (int j = 0; j < fields[i].size(); j++)
		{
			if (j > 0)
			{
				xlText += "\t";
			}
			xlText += fields[i][j];
		}
		xlText += "\n";
	}
	return xlText;
}