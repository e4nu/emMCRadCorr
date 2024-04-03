#include "AppsUtils.h"

using namespace e4nu ;

std::string utils::GetArg(std::string op, int argc, char ** argv )
{
  const int buf_size = 2048*128;
  char *  argument   = new char[buf_size];
  strcpy(argument, "");

  while(argc>2)
    {
      if (argv[1][0] == '-' && argv[1][1] == '-') {

	char op_cur[buf_size];
	strcpy(op_cur,&argv[1][2]);

	if (strcmp(op.c_str(),op_cur)==0) {
	  if (strlen(&argv[2][0]) ) {
            strcpy(argument,&argv[2][0]);
	  }
	}
      }
      argc--;
      argv++;

    }

  std::string value = std::string(argument);
  delete [] argument;
  return value ;
}


bool utils::ExistArg(std::string op, int argc, char ** argv )
{
  const int buf_size = 2048*128;
  char *  argument   = new char[buf_size];
  strcpy(argument, "");

  while(argc>2)
    {
      if (argv[1][0] == '-' && argv[1][1] == '-') {

	char op_cur[buf_size];
	strcpy(op_cur,&argv[1][2]);

	if (strcmp(op.c_str(),op_cur)==0) {
	  return true ;
	}
      }
      argc--;
      argv++;

    }
  delete [] argument ;
  return false;
}

double utils::GetCLAS6TargetThickness( double tgt ) { 
  static const unsigned int kPdgHe3 = 1000020030; 
  static const unsigned int kPdgC12 = 1000060120 ; 
  static const unsigned int kPdgFe56 = 1000260560 ;

  const double tH = 0.027760 ;
  const double tHe = 0.005772;
  const double tC = 0.004183;
  const double tFe = 0.008532 ;
  
  if ( tgt == kPdgHe3 ) return tHe ; 
  else if ( tgt == kPdgC12 ) return tC ; 
  else if ( tgt == kPdgFe56 ) return tFe ; 
  return tH ; 
}
