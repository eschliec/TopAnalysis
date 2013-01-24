#include "CommandLineParameters.hh"

#include <fstream>
CLParameterBase::CLParameterBase( const std::string &identifier,
                                  const std::string &helpText,
                                  const bool isRequired,
                                  const unsigned int minArguments,
                                  const unsigned int maxArguments )
  : _identifier(identifier),
    _helpText(helpText),
    _minArguments(minArguments),
    _maxArguments(maxArguments),
    _required(isRequired) {
  _wasSet = false;
  CLAnalyser::registerGlobalParameter(this);
}



void CLParameterBase::addArgument(const std::string argument){
  _arguments.push_back(argument);
  parseArgument(argument);
}

const std::string CLParameterBase::getArguments() const {
  std::stringstream example;
  for (std::vector<std::string>::const_iterator iter = _arguments.begin(); iter != _arguments.end(); ++iter)
    example << " \"" <<*iter << "\"" ;
  return example.str();
}

const std::string CLParameterBase::getCommandLineExample() const {
  std::stringstream example;
  if (! isRequired()) example << "[";
  if (_identifier.size() > 1) example << "--";
  else example << "-";
  example << _identifier << " ";

  example << getExample();

  if (! isRequired() ) example << "]";

  return example.str();
}

const std::string CLParameterBase::getSteeringFileExample() const {
  std::stringstream example;
  example << _identifier << " ";

  example << getExample();
  return example.str();
}

const std::string CLParameterBase::getExample() const {
  std::stringstream example;

  if (_minArguments > 0) {
    if (_maxArguments > 1) {
      example << "<value1> ";
      if (_minArguments > 2) example << ".. ";
      if (_minArguments > 1) example << "<value" << _minArguments << "> ";
    }
    else example << "<value> ";
  }

  if (_maxArguments > _minArguments) {
    if (_maxArguments == std::numeric_limits<unsigned int>::max()) example << "[ .. <valueN> ] ";
    else if (_maxArguments - _minArguments > 1) example << "[ .. <value" << _maxArguments << "> ] ";
    else if (_minArguments > 0) example << "[ <value" << _maxArguments << "> ] ";
    else  example << "[ <value> ]";
  }
  example << " ";
  return example.str();
}

const std::string CLParameterBase::getHelp() const {
  return _helpText;
}


CLAnalyser *CLAnalyser::_globalInstance = 0;

void CLAnalyser::registerGlobalParameter(CLParameterBase *parameter) {
  if (_globalInstance == 0) _globalInstance = new CLAnalyser();
  _globalInstance->registerParameter(parameter);
}

void CLAnalyser::interpretGlobal(const int argc, char** argv ) {
  if (_globalInstance != 0) _globalInstance->interpret(argc,argv);
}

void CLAnalyser::registerParameter(CLParameterBase *parameter) {
  _parameters.push_back(parameter);
}

CLParameterBase* CLAnalyser::findParameter(const std::string &name) const {

  for (parVec_t::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter)
    if ( (*iter)->getIdentifier() == name ) return *iter;

  std::cout << "Parameter " << name << " is unknown" << std::endl;
  exit(1);

  return 0;
}

void CLAnalyser::checkParameters() const {
  bool error = false;
  for (parVec_t::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter)
    if ( ! (*iter)->check() ) error = true;
  if (error) {
    std::cout << "Bad parameters, cannot continue." << std::endl;
    exit(1);
  }
}

int CLAnalyser::getNumberOfRequiredTokens() const {
  int noRequired = 0;
  for (parVec_t::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter)
    if ( (*iter)->isRequired() ) noRequired += (*iter)->getMinArguments();
  return noRequired;
}

void CLAnalyser::printExample(const std::string &programName) const {
  std::cout << "usage: " << programName;
  for (parVec_t::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter)
    std::cout << " " << (*iter)->getCommandLineExample();

  std::cout << std::endl;
}

void CLAnalyser::printHelp() const {
  std::cout << "CLP toolkit parameters:" << std::endl
	    << "  --help                       this help text" << std::endl
	    << "  --CLP-write-file <filename>  write steering file with current parameters" << std::endl
	    << "  --CLP-read-file <filename>   read parameters from steering file (command line settings will superseede file settings)" << std::endl
	    << std::endl
	    << "program parameters:" << std::endl;

  for (parVec_t::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter)
    std::cout << " " << (*iter)->getCommandLineExample() << " " << (*iter)->getHelp() << std::endl;
}

void CLAnalyser::readSteeringFile(const std::string &fileName) {
  std::ifstream file;

  file.open(fileName.c_str());

  while ( file.good() && ! file.eof() ) {

    std::string line;

    std::getline(file, line);

    std::stringstream stream(line);

    std::string identifier;
    if ( (stream >> identifier).fail() ) continue;
    if (identifier.size() == 0 || identifier[0] == '#') continue;

    std::string token;
    CLParameterBase *parameter = findParameter(identifier);

    if (parameter->isSet()) {
      std::cout << identifier << " already set by command line, ignoring steering file settings." << std::endl;
      continue;
    }
    else {
      std::cout << " setting parameter: " << identifier;
      parameter->set();
      while ( stream.good() && ! stream.eof() ) {
	if ( (stream >> token).fail() ) continue;
	if ( token.size() > 0 && token[0] == '"') {
	  std::stringstream quotedParameter;
	  quotedParameter << token.substr(1);
	  while ( token[token.size()-1] != '"' ) {
	    if ( (stream >> token).fail() ) {
	      std::cout << "error in steering file, missing end of quotes for parameter " << identifier << std::endl;
	      exit(1);
	    }
	    quotedParameter << " " << token;
	  }
	  token = quotedParameter.str().substr(0,quotedParameter.str().size()-1);
	}
	std::cout << " " << token ;
	parameter->addArgument(token);
      }
    }
    std::cout << std::endl;
  }

  file.close();

}

void CLAnalyser::writeSteeringFile(const std::string &fileName) const {
  std::ofstream file;

  file.open(fileName.c_str());

  for (parVec_t::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter) {
    file << "# " << (*iter)->getIdentifier() << "  " << (*iter)->getHelp() << std::endl;
    if ( (*iter)->isRequired()) file << "# " << (*iter)->getIdentifier() << " is required" << std::endl;
    file << "# " << (*iter)->getSteeringFileExample() << std::endl;
    if ( (*iter)->isSet() ) file << (*iter)->getIdentifier() << " " << (*iter)->getArguments() << std::endl;
    file << std::endl;
  }

  file.close();
}

void CLAnalyser::interpret(const int argc, char** argv) {

  if (argc-1 < getNumberOfRequiredTokens() ) {

    for (int iPar = 1; iPar < argc; ++iPar ) {
      std::string token(argv[iPar]);
      if (token == "--help" || token == "-h") {
	printHelp();
	exit(0);
      }
    }

    printExample(argv[0]);
    exit(0);
  }

  std::string currentParameter;
  std::string steeringReadFileName = "";
  std::string steeringWriteFileName = "";

  for (int iPar = 1; iPar < argc; ++iPar ) {
    std::string token(argv[iPar]);

    if (token == "--help" || token == "-h") {
      printHelp();
      exit(0);
    }

    if (token == "--CLP-read-file") {
      ++iPar;
      if (iPar < argc) {
	steeringReadFileName = argv[iPar];
	continue;
      }
      else {
	std::cout << "--CLP-read-file needs filename as argument" << std::endl;
	exit(1);
      }
    }
    if (token == "--CLP-write-file") {
      ++iPar;
      if (iPar < argc) {
	steeringWriteFileName = argv[iPar];
	continue;
      }
      else {
	std::cout << "--CLP-write-file needs filename as argument" << std::endl;
	exit(1);
      }
    }

    if (token.size() > 1 && token.substr(0,1) == "-")
      if (token.size() > 2 && token.substr(1,1) == "-") {
	currentParameter = token.substr(2);
	findParameter(currentParameter)->set();
      }
      else {
        // test if this is just a negative number
        std::istringstream stream(token);
        double dummyNumber;
        if ( ( ! (stream >> dummyNumber).fail() ) && stream.eof() )
          findParameter(currentParameter)->addArgument(token);
        else { // otherwise we have one letter options
          for (unsigned int i = 1; i < token.size(); ++i) {
            currentParameter = token[i];
            findParameter(currentParameter)->set();
          }
        }
      }
    else {
      findParameter(currentParameter)->addArgument(token);
    }
  }

  if (steeringReadFileName != "") {
    std::cout << "Will read steering file " << steeringReadFileName << " now." << std::endl;
    readSteeringFile(steeringReadFileName);
  }

  checkParameters();

  if (steeringWriteFileName != "") {
    std::cout << "Will write steering file " << steeringWriteFileName << " now, programm will not be executed." << std::endl;
    writeSteeringFile(steeringWriteFileName);
    exit(0);
  }

}

template <>
void CLParameter<std::string>::parseArgument(const std::string token) {

  _values.push_back(token);
}

