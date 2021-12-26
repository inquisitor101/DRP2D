#include "option_structure.hpp"



void AddBoolOption
(
  std::ifstream  &inputFile,
  const char     *keyword,
  bool           &value,
  std::string    &defaultValue,
  bool            ResetLine
)
 /*
  * Function that reads a boolean input, otherwise uses the defaultValue.
  */
{
  // Temporary storage for the boolean option.
  std::string NameValue;
  // Read zone conformity, if specified.
  AddScalarOption(inputFile, keyword, NameValue, defaultValue, ResetLine);

  // Make sure input is either true or false.
  if( NameValue.compare("true") && NameValue.compare("false") ){
    std::string message = keyword; message += " must be either: true of false.";
    Terminate("AddBoolOption", __FILE__, __LINE__, message);
  }
  // Convert input string to boolean.
  std::istringstream(NameValue) >> std::boolalpha >> value;
}


bool FindStringVectorFromKeyword
(
 std::ifstream            &inputFile,
 const char               *keyword,
 std::vector<std::string> &valStringVector
)
/*
 * Function that searches for and finds a string vector from input file.
 */
{
  // Length of target string.
  const int lenKeyword = std::strlen(keyword);

  // Initialize search as false.
  bool keywordFound = false;

  // Line stream.
  std::string lineBuf;

  // Read line by line, until/if value of string is matched.
  while(std::getline(inputFile, lineBuf)){
    // Remove tabs and new lines.
    ReplaceTabsAndReturns(lineBuf);

    // Ignore comments, i.e. "%"
    std::string::size_type pos = lineBuf.find("%");
    if(pos != std::string::npos ) lineBuf.erase(pos, lineBuf.length()-pos);

    // Remove any blanks in line.
    RemoveBlanks(lineBuf);

    // Evaluate current line.
    if( lineBuf.length() ){

      // Compare line with target string.
      if(lineBuf.compare(0, lenKeyword, keyword) == 0){
				// Find symbol ( position.
				pos = lineBuf.find("(");

				// Current line found begins after ( symbol.
				std::string lineFound = lineBuf.substr(pos+1);

				// Initialize string stream of current line.
				std::stringstream lineStream(lineFound);
				// Iterate while the current line is not over.
				while( lineStream.good() ) {

					// Initialize data.
					std::string data;
					// Populate data.
					std::getline(lineStream, data, ',');
					// Remove blanks in keyword.
					RemoveBlanks(data);
					// Add found value to vector.
					valStringVector.push_back(data);

					// Check whether data is corrupted or empty.
					if( data == "" )
						Terminate("FindStringVectorFromKeyword", __FILE__, __LINE__,
											"Data is corrupted and/or empty!");
				}

				// If string not empty, assert the format used is appropriate.
				if( !valStringVector.empty() ){

					// Extract final index.
					unsigned short idx = valStringVector.size()-1;

					// Extract last character of string.
					char endSymbol = valStringVector[idx].back();

					// Exit immediately, if incosistent use is found.
					if( endSymbol != ')' )
						Terminate("FindStringVectorFromKeyword", __FILE__, __LINE__,
											"Wrong format used, proper use is: ( data1, data2, etc )");

					// Remove final character: right-parenthesis symbol.
					valStringVector[idx].pop_back();
				}

        // Keyword found.
        keywordFound = true;
      }
    }
    // If target string is found, break out of loop.
    if(keywordFound) break;
  }

  // Return happily.
  return keywordFound;
}


bool FindStringFromKeyword
(
 std::ifstream &inputFile,
 const char    *keyword,
 std::string   &valString
)
/*
 * Function that searches for and finds a string from input file.
 */
{
  // Length of target string.
  const int lenKeyword = std::strlen(keyword);

  // Initialize search as false.
  bool keywordFound = false;

  // Line stream.
  std::string lineBuf;

  // Read line by line, until/if value of string is matched.
  while(std::getline(inputFile, lineBuf)){
    // Remove tabs and new lines.
    ReplaceTabsAndReturns(lineBuf);

    // Ignore comments, i.e. "%"
    std::string::size_type pos = lineBuf.find("%");
    if(pos != std::string::npos ) lineBuf.erase(pos, lineBuf.length()-pos);

    // Remove any blanks in line.
    RemoveBlanks(lineBuf);

    // Evaluate current line.
    if( lineBuf.length() ){
      // Create a lower-case version of current line.
      std::string lineBufLower = lineBuf;
      //CreateLowerCase(lineBufLower);

      // Compare line with target string.
      if(lineBufLower.compare(0, lenKeyword, keyword) == 0){
        pos = lineBuf.find("=");
        valString = lineBuf.substr(pos+1);

        // Keyword found.
        keywordFound = true;
      }
    }
    // If target string is found, break out of loop.
    if(keywordFound) break;
  }

  // Return happily.
  return keywordFound;
}


void RemoveBlanks
(
  std::string &lineBuf
)
 /*
  * Function that takes string and removes any blanks in it.
  */
{
  // Find the first non space character.
  int strLength = lineBuf.length();
  int posLeading = 0;
  while(posLeading < strLength && isspace(lineBuf[posLeading])) ++posLeading;

  // Find the last non space character.
  int posTrailing = strLength - 1;
  while(posTrailing >= 0 && isspace(lineBuf[posTrailing])) --posTrailing;

  // Determine the situation.
  if(posLeading == strLength || posTrailing < 0)
  {
    // No non-blanks in the string. Set lineBuf to an empty string.
    lineBuf = "";
  }
  else
  {
    // Non-blanks are present. Remove the blanks. First the trailing ones,
    // because otherwise the positions are screwed up.
    int eraseLenBack = strLength - posTrailing - 1;
    if( eraseLenBack ) lineBuf.erase(posTrailing+1, eraseLenBack);
    lineBuf.erase(0, posLeading);
  }
}


void ReplaceTabsAndReturns
(
  std::string &lineBuf
)
 /*
  * Function that replaces tabs and returns with blanks.
  */
{
  // Replace the tabs.
  for(;;)
  {
    std::string::size_type pos = lineBuf.find("\t");
    if(pos == std::string::npos) break;
    lineBuf.replace(pos, 1, " ");
  }

  // Replace the returns.
  for(;;)
  {
    std::string::size_type pos = lineBuf.find("\n");
    if(pos == std::string::npos) break;
    lineBuf.replace(pos, 1, " ");
  }
}


void CreateLowerCase
(
  std::string &lineBuf
)
 /*
  * Function that returns a lowered-case version of input string.
  */
{
  // Determine the length of the string and convert its elements to
  // lower case.
  std::string::size_type strLength = lineBuf.length();
  for(std::string::size_type i=0; i<strLength; ++i)
    lineBuf[i] = (char) tolower(lineBuf[i]);
}


void Terminate
(
  const char        *functionName,
  const char        *fileName,
  const int          lineNumber,
  const std::string &errorMessageIn
)
 /*
  * Function which prints an error message and exits.
  */
{
  // Insert a comment character and a space on places where a newline
  // character occurs in the string, such that the error message
  // looks nicer.
  std::string  errorMessage  = errorMessageIn;
  std::string  tmpString     = errorMessage;
  std::string::size_type off = 1;

  for(;;)
  {
    std::string::size_type loc = tmpString.find("\n");
    if(loc == std::string::npos) break;
    errorMessage.insert(loc+off, "# ");
    off += loc+3;
    tmpString.erase(0,loc+1);
  }

  // Header of the error message.
	std::cout << std::endl;
	std::cout << "#" << std::endl
            << "#=========================== !!! Error !!! "
            << "============================" << std::endl
            << "#" << std::endl;

  // Write the function name, file name and line number.
  std::cout << "#* Run-time error in " << functionName << std::endl;
  std::cout << "#* File: " << fileName
            << ", Line: "  << lineNumber << std::endl
            << "#" << std::endl;

  // Write the error message and the terminating line.
  std::cout << "# " << errorMessage << std::endl
            << "#" << std::endl;
  std::cout << "#=========================================="
            << "============================" << std::endl
            << "#" << std::endl << std::flush;

  // And exit.
  exit(1);
}


#ifdef ENABLE_NAN_CHECK
void CheckFloatingError
(
  void
)
 /*
  * Function that determines if a floating-point error is encountered.
  */
{
  std::string message;
  // Check what error has been raised.
  if( std::fetestexcept(FE_DIVBYZERO) ) message += "Error... FE_DIVBYZERO,";
  // if( std::fetestexcept(FE_INEXACT)   ) message += "Error... FE_INEXACT,";
  if( std::fetestexcept(FE_INVALID)   ) message += "Error... FE_INVALID,";
  if( std::fetestexcept(FE_OVERFLOW)  ) message += "Error... FE_OVERFLOW,";
  //if( std::fetestexcept(FE_UNDERFLOW) ) message += "Error... FE_UNDERFLOW,";

  if( !message.empty() )
    Terminate("CheckFloatingError", __FILE__, __LINE__, message);
}
#endif



