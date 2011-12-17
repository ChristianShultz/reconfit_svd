// semble_file_management.cc -
//
// Saturday, October 22 2011
//

#include"semble_file_management.h"


namespace SEMBLE
{

  namespace SEMBLEIO
  {
    std::string getPath(void)
    {
      char cCurrentPath[FILENAME_MAX];

      if(!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
        {
          std::cout << "Couldnt get path in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
          exit(1);
        }

      cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; //stick in a null character
      std::string dum(cCurrentPath);
      return  dum += std::string("/");
    }


    void makeDirectoryPath(const std::string &s)
    {
      if(!!!(access(s.c_str(), 0) == 0))  //check if it already exists, theres no reason to use extra system cmds if we dont have to since they're evil
        {                                 //NB, access returns true if there was a file OR directory with that name, would be nice to check if its a directory
          std::cout << "Making path:" << s << std::endl;
          std::string cmd = "mkdir -p ";
          cmd += s;                       
	  int dum = system(cmd.c_str());  //assigning the return value of system to an integer should force it to wait till the call is finished 
	 
	  /*
	    std::cout << dum << std::endl;
	  */

	  //the return value of system is platform specific and may change with updates to the qcd machines  so we need to use access to check that the call worked not dum
	  if(!!!(access(s.c_str(), 0) == 0)) //check to see that the directory was actually made
	    {
	      std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " error " << std::endl;
	      std::cout << "System call failed to make the directory:" << s << " exiting.." << std::endl;
	      exit(1);
	    }
        }
    }

  }//SEMBLEIO
}//SEMBLE
