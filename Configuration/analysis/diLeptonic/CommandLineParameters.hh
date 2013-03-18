#ifndef CommandLineParameters_hh
#define CommandLineParameters_hh

#include <string>
#include <sstream>
#include <limits>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <functional>

class CLAnalyser;

/**
 * Base class of a command line parameter.
 *
 * This class contains the type independent part of the implementation.
 *
 * @author Benjamin.Lutz@desy.de
 * @version 0.1
 * @date November 2009
 */
class CLParameterBase {

  friend class CLAnalyser;

public:
  /**
   * check if this parameter was set by the command line
   *
   * @returns true if the parameter was set
   */
  bool isSet() const { return _wasSet; }
  /**
   * get the name of the parameter
   *
   * @returns identifying name
   */
  const std::string& getIdentifier() const {return _identifier;}

protected:

  /**
   * Constructor for a basic parameter. The parameter will be
   * automtically registered in the global list of parameters.
   *
   * @param[in] identifier the name of the parameter (without - or --)
   * @param[in] helpText description what this parameter does
   * @param[in] isRequired if this parameter is mandatory
   * @param[in] minArguments minimum number of required arguments
   * @param[in] maxArguments maximum number of allowed arguments
   */
  CLParameterBase( const std::string &identifier,
		   const std::string &helpText,
		   const bool isRequired = false,
		   const unsigned int minArguments = 0,
		   const unsigned int maxArguments = std::numeric_limits<unsigned int>::max() );

  /**
   * helper function to generate steering file
   */
  const std::string getSteeringFileExample() const;
  /**
   * helper function to generate steering file
   */
  const std::string getArguments() const;
  /**
   * helper function to generate usage message
   */
  const std::string getCommandLineExample() const;
  /**
   * helper function to generate usage message
   */
  const std::string getExample() const;
  /**
   * helper function to generate help message
   */
  const std::string getHelp() const;

  /**
   * function to add argument to parameter
   *
   * @param argument argument string
   */
  void addArgument(const std::string argument);
  /**
   * function to convert argument token to internal type
   *
   * this has to be implemented in the daughter class
   */
  virtual void parseArgument(const std::string token) = 0;
  /**
   * function to check if all requirements are met
   *
   * this has to be implemented in the daughter class
   */
  virtual bool check() const = 0;

  /**
   * set parameter state to set
   */
  void set() { _wasSet = true; }

  /**
   * query if this parameter is mandatory
   */
  bool isRequired() const { return _required;}
  /**
   * @returns the minimum number of arguments for this parameter
   */
  int getMinArguments() const {return _minArguments;}

  const std::string _identifier;
  const std::string _helpText;

  const unsigned int _minArguments;
  const unsigned int _maxArguments;

  const bool _required;

  bool _wasSet;

  std::vector<std::string> _arguments;
};

/**
 * Templated implementation of a command line parameter.
 *
 * This class contains the type dependent part of the implementation.
 *
 * @author Benjamin.Lutz@desy.de
 * @version 0.1
 * @date November 2009
 */
template <class T>
class CLParameter : public CLParameterBase {
public:
  /**
   * @param[in] identifier the name of the parameter (without - or --)
   * @param[in] helpText description what this parameter does
   * @param[in] isRequired if this parameter is mandatory
   * @param[in] minArguments minimum number of required arguments
   * @param[in] maxArguments maximum number of allowed arguments
   */
  CLParameter( const std::string &identifier,
	       const std::string &helpText,
	       const bool isRequired = false,
	       const unsigned int minArguments = 0,
	       const unsigned int maxArguments = std::numeric_limits<int>::max(),
               std::function<bool(T)> checkFunction = 0
             )
    : CLParameterBase( identifier, helpText, isRequired, minArguments, maxArguments ),
      _checkFunction(checkFunction)
  {}

  /**
   * get the arguments of the parameter
   *
   * @returns vector of the arguments
   */
  const std::vector<T> &getArguments() const {
    return _values;
  }

  /**
   * get the number of arguments of the parameter
   *
   * @returns number of the arguments
   */
  unsigned int size() const {
    return _values.size();
  }

  /**
   * get one argument of the parameter
   *
   * @param[in] index index of the argument to be returned
   * @returns the argument with index
   */
  T operator [] (const unsigned int index) const {
    return _values.at(index);
  }
protected:
  virtual void parseArgument(const std::string token);
  virtual bool check() const;

private:

  std::vector<T> _values;
  std::function<bool(T)> _checkFunction;
};

template <class T>
bool CLParameter<T>::check() const {
  bool result = true;
  if ( isSet() && _values.size() < _minArguments) {
    std::cout << "Parameter " << getIdentifier() << " needs at least " << _minArguments << " arguments."<< std::endl;
    result = false;
  }
  if ( _values.size() > _maxArguments) {
    if (_maxArguments==0) std::cout << "parameter " << getIdentifier() << " does not take arguments. " << std::endl;
    else std::cout << "parameter " << getIdentifier() << " does not take more than " << _maxArguments << " arguments. " << std::endl;
    result = false;
  }
  if ( _required && ! isSet() ) {
    std::cout << "Parameter " << getIdentifier() << " needs to be set. " << std::endl;
    result = false;
  }
  if (_checkFunction && isSet()) {
    for (size_t i = 0; i < _values.size(); ++i) {
      if (!_checkFunction(_values[i])) {
        std::cout << "Parameter " << getIdentifier() << " has invalid value: " << _values[i] << std::endl;
        result = false;
      }
    }
  }
  return result;
}

template <class T>
void CLParameter<T>::parseArgument(const std::string token) {

  std::stringstream conversionStream;
  conversionStream << token;

  T value;

  if ( (conversionStream >> value).fail() ) {
    std::cout << "cannot convert " << token << " to the required type for parameter " << getIdentifier() << std::endl;
    exit(1);
  }

  _values.push_back(value);
}



/**
 * Class that can search the command line for parameters.
 *
 * @author Benjamin.Lutz@desy.de
 * @version 0.1
 * @date November 2009
 */
class CLAnalyser {

public:
  /**
   * make a parameter known to the interpreter
   *
   * @param[in,out] parameter which should be checked for
   */
  void registerParameter(CLParameterBase *parameter);
  /**
   * analyse the command line
   *
   * @param[in] argc number of command line arguments
   * @param[in] argv array of comman line arguments
   */
  void interpret(const int argc, char** argv);

  /**
   * add a parameter to the global list of parameters
   *
   * @param[in,out] parameter which should be checked for
   */
  static void registerGlobalParameter(CLParameterBase *parameter);
  /**
   * analyse the command line checking for all global prameters
   *
   * @param[in] argc number of command line arguments
   * @param[in] argv array of comman line arguments
   */
  static void interpretGlobal(const int argc, char** argv);
private:

  static CLAnalyser *_globalInstance;

  CLParameterBase* findParameter(const std::string &name) const ;
  void checkParameters() const;
  int getNumberOfRequiredTokens() const;

  void printExample(const std::string &programName) const;
  void printHelp() const;

  void readSteeringFile(const std::string &fileName);
  void writeSteeringFile(const std::string &fileName) const;

  typedef std::vector<CLParameterBase*> parVec_t;
  parVec_t _parameters;
};


#endif
